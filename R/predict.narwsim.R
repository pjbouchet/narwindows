#' Forecast NARW abundance
#'
#' @param obj Object returned by narw()
#' @param n Number of replicate projections
#' @param yrs Number of years for the projection
#' 
#' @export
#'
#' @author Phil J. Bouchet
#'
predict.narwsim <- function(obj,
                            n = 100,
                            yrs = 35,
                            ...) {

  # Adapted from original code by Scott Creel
  # https://www.montana.edu/screel/teaching/bioe-440r-521/course-outline/stoch_projection_new.R
  
  #'------------------------------------------------------------
  # Function checks ----
  #'------------------------------------------------------------
  
  if(!inherits(obj, "narwsim")) stop("Object must be of class <narwsim>")
  if(is.null(obj$gam$post)) warning("Posterior samples for terminal functions not available; predictions are based on the mean fitted relationships between body condition and survival/health. Run the expand() function to estimate variance through posterior simulation.")
  
  if(sum(suppressWarnings(purrr::map_lgl(.x = obj$gam$fit, .f = ~any(is.null(unlist(.x))))) > 0)){
    w1 <- "A problem occurred" 
    if(obj$param$nsim<100) w2 <- " (sample size likely insufficient)" else w2 = ""
    stop(paste0(w1, w2, " --> Cannot generate projections"))
  }
  
  if(length(obj$param$cohortID) < 6) stop("Missing cohorts --> Cannot generate projections")
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
  # Default arguments
  progress <- TRUE
  
  # Default values
  if("progress" %in% names(args)) progress <- args[["progress"]] else progress <- TRUE
  
  # User can specify the <yrs> argument either as a number of years from now
  # or as a target end year
  if(yrs > 2000){
    yrs <- yrs - lubridate::year(lubridate::now())
  }

  cohortID <- obj$param$cohortID
  cohorts.proj <- obj$param$cohorts |> dplyr::slice(1) |> 
    dplyr::mutate(name = gsub(", female", "", name), abb = gsub(",f", "", abb)) |> 
    dplyr::bind_rows(obj$param$cohorts) |> 
    dplyr::mutate(name = gsub("male, female", "female", name), abb = gsub("m,f", "f", abb))
  cohorts.proj <- cohorts.proj[c(1,3,5,2,4,6,7,8)]
  
  # Attributes to monitor during projection
  mat.attribs <- c("alive", "cohort", "female", "age", "length", "tot_mass", "lean_mass", "bc",
                   "p_surv", "min_bc", "time_rest", "birth", "reprod")
  
  # Current year
  current.yr <- lubridate::year(lubridate::now())
  
  #'------------------------------------------------------
  # TERMINAL FUNCTIONS ----
  #'------------------------------------------------------
  
  # Extract terminal functions
  # mod <- obj$gam$fit
  # mod[["gest"]] <- gam_gest
  
  # Spline interpolation of the fitted GAM relationships
  # mbc_preds <- splinefun(x = gam_gest$mod[,2], y = gam_gest$fitted.values)
  # bc_preds <- obj$gam$pred$bc
  # surv_preds <- obj$gam$pred$surv
  
  # Fitted models
  surv.model <- obj$gam$fit$surv
  bc.model <- obj$gam$fit$bc
  mbc.model <- gam_gest$gam
  
  # Inverse link functions
  ilink.surv <- family(obj$gam$fit$surv)$linkinv
  ilink.bc <- family(obj$gam$fit$bc)$linkinv
  ilink.mbc <- family(gam_gest$gam)$linkinv
  
  # Posterior samples
  if(is.null(obj$gam$post)){
    
    if(inherits(surv.model, "gam")){
      surv.mcmc <- matrix(data = coef(surv.model), nrow = 1, ncol = length(coef(surv.model)))
    } else if(inherits(surv.model, "scam")){
      surv.mcmc <- matrix(data = tsgam:::coef.scam(surv.model), nrow = 1, ncol = length(tsgam:::coef.scam(surv.model)))
    }
    
    if(inherits(bc.model, "gam")){
      bc.mcmc <- matrix(data = coef(bc.model), nrow = 1, ncol = length(coef(bc.model)))
    } else if(inherits(surv.model, "scam")){
      bc.mcmc <- matrix(data = tsgam:::coef.scam(bc.model), nrow = 1, ncol = length(tsgam:::coef.scam(bc.model)))
    }

    mbc.mcmc <- matrix(data = coef(mbc.model), nrow = 1, ncol = length(coef(mbc.model)))
    
    n.mcmc <- 1
    
  } else {
    surv.mcmc <- obj$gam$post$surv
    bc.mcmc <- obj$gam$post$bc
    mbc.mcmc <- obj$gam$post$mbc
    n.mcmc <- nrow(surv.mcmc)
  }
  
  #'------------------------------------------------------
  # Abortion rate ----
  #'------------------------------------------------------
  
  abort.rate <- sum(obj$abort[, abort])/obj$param$nsim
  
  #'------------------------------------------------------
  # INITIALIZATION ----
  #'------------------------------------------------------
  
  # Define initial population vector
  # N0 <- c(0, 7, 212, 0, 71, 0, 7, 61)
  N_0 <- N0[["yr2019"]]
  names(N_0) <- cohorts.proj[, name]
  
  console("Initializing")
  
  # Initial cohort vector
  cohort.vec <- do.call(c, sapply(X = 1:nrow(cohorts.proj), FUN = function(x) rep(cohorts.proj$id[x], each = N_0[x])))
  
  # Sterile females - 4% never reproduce (based on PCoMS estimate)
  nonrep <- 0.04
  reprod.fem <- rep(1, sum(N_0))
  reprod.fem[sample(
    x = which(cohort.vec == 6), # Assigning 0 to resting females only
    size = round(nonrep * sum(N_0[6:8]), 0) # Based on total number of mature females in the population
  )] <- 0
  
  # Number of individuals in each cohort
  narw.pop <- array(
    data = NA, c(n, yrs + 1, nrow(cohorts.proj)),
    dimnames = list(
      paste0("prj ", 1:n),
      paste0("yr ", 0:yrs),
      cohorts.proj$name
    )
  )
  narw.pop[, 1, ] <- rep(N_0, each = n)
  
  # Total population size
  tot.pop <- matrix(0, n, yrs + 1, dimnames = list(paste0("prj", 1:n), paste0("yr ", 0:yrs)))
  tot.pop[, 1] <- sum(N_0)
  
  births <- array(
    data = NA, c(yrs + 1, n),
    dimnames = list(
      paste0("yr ", 0:yrs),
      paste0("prj ", 1:n)
    )
  )
  
  #'------------------------------------------------------
  # PROJECTIONS ----
  #'------------------------------------------------------
  
  # This uses nested loops. 
  # The prj loop (outermost loop) replicates the projection <n> times.
  # The i loop is next, and steps across all years of projection from an initial population vector.
  
  # Set up progress bar
  pb <- progress::progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = n, clear = FALSE, width = 80
  )
  
  console("Initializing", TRUE)
  cat("Running projections ...\n")
  
  start.time <- Sys.time()
  
  narw.individuals <- vector(mode = "list", length = n)
  
  for(prj in 1:n){
    
    if(progress) pb$tick() # Update progress bar
    
    random.curves <- sample(n.mcmc, n.mcmc, replace = TRUE)
    
    # Create matrices and initialize them
    # rows: years <yrs>
    # columns: <attributes>
    # layers: individuals <n>
    # 4th dimension: replicate projection -> later converted to list
    narw.indiv <- array(
      data = NA, c(yrs + 1, length(mat.attribs), sum(N_0)),
      dimnames = list(
        paste0("yr ", 0:yrs), mat.attribs,
        paste0("whale ", 1:(sum(N_0)))
      )
    )
    
    # Alive and population cohort
    narw.indiv[1, "alive", ] <- 1
    narw.indiv[1, "cohort", ] <- cohort.vec
    
    # Sex
    #  -- Calves (male)
    narw.indiv[, "female", 1:N_0[1]] <- 0
    #  -- Calves (female)
    fem <- which(cohort.vec == 0)
    narw.indiv[, "female", fem[fem > N_0[1]]] <- 1
    #  -- Juveniles and adults
    narw.indiv[, "female", which(cohort.vec %in% c(2, 4, 5, 6))] <- 1
    narw.indiv[, "female", which(cohort.vec %in% c(1, 3))] <- 0
    
    # Age
    ages <- start_age_vec(cohort.vec)
    narw.indiv[1, "age", ] <- ages
    
    # Total body length
    lengths <- age2length_vec(ages)
    narw.indiv[1, "length", ] <- lengths
    
    mass <- length2mass_vec(lengths)
    narw.indiv[1, "lean_mass", ] <- mass
    
    # Body conditon
    bc <- start_bcondition_vec(cohort.vec)
    narw.indiv[1, "bc", ] <- bc
    
    # Total mass
    narw.indiv[1, "tot_mass", ] <- mass / (1-bc)
    
    # Calving interval - mean of 5.3 years, up to apparent gaps of 13 years (Kraus et al. 2001)
    # NARWC Report Card 2022: Average inter-birth interval of 7.7 yrs, median of 8, min/max (2/13)
    # This corresponds to a Normal (7.7, 1.45)
    # Mean age at first calving = 9.53 +/- 2.32 (Kraus et al. 2001)
    # 
    # Stewart et al. 2022 -- 
    # The degree to which the energetic reserves of females are depleted during lactation may govern
    # the length of the resting period between successful pregnancies (Miller et al. 2011, Marón et al. 2015).

    # Non-reproductive females
    narw.indiv[1, "reprod", ] <- reprod.fem
    
    # Reserves needed to leave resting state and initiate pregnancy
    Xp <- predict(mbc.model, data.frame(mass = narw.indiv[1, "tot_mass", ]), type = "lpmatrix")
    
    # Different function for each individual
    mbc.curves <- mbc.mcmc[random.curves[1:nrow(Xp)],]
    narw.indiv[1, "min_bc", ] <- (cohort.vec == 6) * ilink.mbc(rowSums(Xp*mbc.curves))
      # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.mbc(Xp[r,] %*% mbc.curves[r,])})
    
  
    # narw.indiv[1, "min_bc", ] <- ilink.mbc(Xp %*% mbc.mcmc[sample(n.mcmc, 1, replace = TRUE),]) * (cohort.vec == 6)
    
    narw.indiv[1, "time_rest", ] <- (cohort.vec == 6)
      
    # narw.indiv[1, "rest", ] <- ifelse(reprod.fem == 0, 1, (cohort.vec == 6))
    # * (1 - (bc >= narw.indiv[1, "min_bc", ])))

    # Calving events
    narw.indiv[1, "birth", ] <- as.numeric(cohort.vec == 4)
    
    # Survival
    Xp <- predict(surv.model, data.frame(start_bc = bc, cohort = cohort.vec), type = "lpmatrix")
    surv.curves <- surv.mcmc[random.curves[1:nrow(Xp)],]
    narw.indiv[1, "p_surv", ] <- ilink.surv(rowSums(Xp*surv.curves))
      
      # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.surv(Xp[r,] %*% surv.curves[r,])})
    # narw.indiv[1, "p_surv", ] <- ilink.surv(Xp %*% surv.mcmc[sample(n.mcmc, 1, replace = TRUE),])
    
    # narw.indiv[1, "p_surv", ] <- (surv_preds[["0"]](bc) * (cohort.vec == 0) +
    #                                 surv_preds[["1"]](bc) * (cohort.vec == 1) +
    #                                 surv_preds[["2"]](bc) * (cohort.vec == 2) +
    #                                 surv_preds[["3"]](bc) * (cohort.vec == 3) +
    #                                 surv_preds[["4"]](bc) * (cohort.vec == 4) +
    #                                 surv_preds[["5"]](bc) * (cohort.vec == 5) +
    #                                 surv_preds[["6"]](bc) * (cohort.vec == 6))
    
    # narw.indiv[1, "p_surv", ] <- 0.85

    #'------------------------------------------------------
    # Loop over years
    #'------------------------------------------------------
    
    for(i in 2:(yrs+1)){
      
      if(tot.pop[prj, i-1] > 0) {
        
        # alive <- narw.indiv[i-1, "alive", ] * (narw.indiv[i-1, "age", ] <=69)
        
        #' ----------------------------
        # ALIVE
        #' ----------------------------
        # Determine whether the animal survived
        # alive <- rbinom(n = length(ps), size = 1, prob = ps) * (narw.indiv[i-1, "age", ] <=69)
        alive <- rbinom(
          n = dim(narw.indiv)[3], 
          size = 1,
          prob = (narw.indiv[i - 1, "alive", ] * narw.indiv[i - 1, "p_surv", ])
        ) * (narw.indiv[i - 1, "age", ] <= 69)
        
        narw.indiv[i, "alive", ] <- alive
        
        #' ----------------------------
        # AGE & SEX
        #' ----------------------------
        
        # Increment age
        narw.indiv[i, "age", ] <- alive * (narw.indiv[i-1, "age", ] + 1)
        
        # Sex remains the same
        narw.indiv[i, "female", ] <- narw.indiv[i-1, "female", ]
        
        #' ----------------------------
        # COHORTS
        #' ----------------------------
        # Maturity - transitions between cohorts
        narw.indiv[i, "cohort", ] <-
          alive * increment_cohort(
            cohort = narw.indiv[i-1, "cohort", ],
            age = narw.indiv[i, "age", ],
            female = narw.indiv[i, "female", ],
            bc = narw.indiv[i-1, "bc", ],
            min_bc = narw.indiv[i-1, "min_bc", ],
            reprod = narw.indiv[i-1, "reprod", ],
            abort = abort.rate)
        
        #' ----------------------------
        # GROWTH
        #' ----------------------------
        
        # Increment length
        narw.indiv[i, "length", ] <- alive * age2length_vec(narw.indiv[i, "age", ])
        
        # Increment lean mass
        narw.indiv[i, "lean_mass", ] <- alive * length2mass_vec(narw.indiv[i, "length", ])
        
        # Predict new body condition from current body condition
        Xp <- predict(bc.model,
          data.frame(
            start_bc = narw.indiv[i - 1, "bc", ],
            cohort = narw.indiv[i - 1, "cohort", ]
          ), type = "lpmatrix")
        
        bc.curves <- bc.mcmc[random.curves[1:nrow(Xp)],]
        narw.indiv[i, "bc", ] <- alive * ilink.bc(rowSums(Xp*bc.curves))
          # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.bc(Xp[r,] %*% bc.curves[r,])})
        
        # narw.indiv[i, "bc", ] <- alive * ilink.bc(Xp %*% bc.mcmc[sample(n.mcmc, 1, replace = TRUE),])
        
        # narw.indiv[i, "bc", ] <-
        #   alive * (bc_preds[["0"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 0) +
        #              bc_preds[["1"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 1) +
        #              bc_preds[["2"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 2) +
        #              bc_preds[["3"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 3) +
        #              bc_preds[["4"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 4) +
        #              bc_preds[["5"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 5) +
        #              bc_preds[["6"]](narw.indiv[i-1, "bc", ]) * (narw.indiv[i-1, "cohort", ] == 6))
        
        # Increment total mass
        narw.indiv[i, "tot_mass", ] <- alive * narw.indiv[i, "lean_mass", ] / (1 - narw.indiv[i, "bc", ])
        
        #' ----------------------------
        # REPRODUCTION
        #' ----------------------------
        
        # Reproductive females
        narw.indiv[i, "reprod", ] <- alive * narw.indiv[i-1, "reprod", ]
        
        # Resting state
        # narw.indiv[i, "rest", ] <-
        #   (narw.indiv[i, "reprod", ] == 1) * alive * (narw.indiv[i, "cohort", ] == 6) *
        #   (1 - (narw.indiv[i, "bc", ] >= narw.indiv[i-1, "min_bc", ]))
                 
        # Which animals are juvenile females that are ready to start reproducing
        # juvenile.females.ofage <- 
        #   (narw.indiv[i-1,"cohort", ] == 2) * (narw.indiv[i-1, "age", ] >= 9) * (narw.indiv[i-1, "bc", ] >= narw.indiv[i-1, "min_bc", ])
        
        # Time spent in resting phase
        narw.indiv[i, "time_rest", ] <-
          alive * (narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i-1, "time_rest", ] + 1)
        
        # Minimum body condition needed to successfully bring fetus to term without starving
        # No evidence of reproductive senescence in right whales - Hamilton et al. (1998)
        Xp <- predict(mbc.model, data.frame(mass = narw.indiv[i, "tot_mass", ]), type = "lpmatrix")
        mbc.curves <- mbc.mcmc[random.curves[1:nrow(Xp)],]
        narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * ilink.mbc(rowSums(Xp*mbc.curves))
          # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.mbc(Xp[r,] %*% mbc.curves[r,])})
        
        # narw.indiv[i, "min_bc", ] <- alive * (narw.indiv[i, "cohort", ] == 6) * 
        #   ilink.mbc(Xp %*% mbc.mcmc[sample(n.mcmc, 1, replace = TRUE),])
        
        # narw.indiv[i, "min_bc", ] <-
        #   alive * mbc_preds(narw.indiv[i, "tot_mass", ]) * (narw.indiv[i, "cohort", ] == 6) 
        # * narw.indiv[i, "rest", ]
        
        # Birth of new calf, conditional on the mother being alive and in pregnant state
        narw.indiv[i, "birth", ] <- alive * (narw.indiv[i, "cohort", ] == 4)
        new.births <- sum(narw.indiv[i, "birth", ])
        births[i, prj] <- new.births
        
        if(new.births > 0){
          
          # Determine sex of newborn calves
          # Note: This needs to be calculated here as maintaining a % of non-reproductive females
          # in the population requires knowing the number of female calves being born
          new.females <- rbinom(new.births, 1, 0.5)
          
          # Number of non-reproductive females in the population
          Nr <- nrow(t(narw.indiv[i, ,])[t(narw.indiv[i, ,])[,"alive"] == 1 & t(narw.indiv[i, ,])[,"female"] == 1 & t(narw.indiv[i, ,])[,"reprod"] == 0, ])
          
          # Total abundance of females given new births
          Nt <- nrow(t(narw.indiv[i, ,])[t(narw.indiv[i, ,])[,"alive"] == 1 & t(narw.indiv[i, ,])[,"female"] == 1 & t(narw.indiv[i, ,])[,"cohort"] %in% c(2,4,5,6), ]) + sum(new.females)
          
          # Newborn females that are non-reproductive
          Nb <- max(0, round(nonrep * Nt) - Nr)
          newborns.nonreprod <- rep(1, new.births)
          newborns.nonreprod[sample(which(new.females == 1), size = Nb)] <- 0
          
          # Create array for newborns
          new.calves <- array(
            data = NA, c(yrs + 1, length(mat.attribs), new.births),
            dimnames = list(
              paste0("yr ", 0:yrs), mat.attribs,
              paste0("whale ", ((dim(narw.indiv)[3]+1):(dim(narw.indiv)[3] + new.births)))
            )
          )
          
          new.calves[i,,] <- add_calf(new.births, mat.attribs, new.females, newborns.nonreprod)
          narw.indiv <- abind::abind(narw.indiv, new.calves, along = 3)
          # test <- abind::abind(narw.indiv, new.calves, along = 3)
          # t(test[i, ,])[t(test[i, ,])[,"alive"] == 1 & t(test[i, ,])[,"female"] == 1 & t(test[i, ,])[,"reprod"] == 0, ]
          # nrow(t(test[i, ,])[t(test[i, ,])[,"alive"] == 1 & t(test[i, ,])[,"female"] == 1 & t(test[i, ,])[,"reprod"] == 0, ])
          # nrow(t(test[i, ,])[t(test[i, ,])[,"alive"] == 1 & t(test[i, ,])[,"female"] == 1, ])
          # t(new.calves[i, ,])[t(new.calves[i, ,])[,"alive"] == 1 & t(new.calves[i, ,])[,"female"] == 1 & t(new.calves[i, ,])[,"reprod"] == 0, ]
          # 
        }
        
        #' ----------------------------
        # SURVIVAL
        #' ----------------------------
        
        # Predict survival probability based on body condition
        Xp <- predict(surv.model,
         data.frame(
           start_bc = narw.indiv[i, "bc", ],
           cohort = narw.indiv[i, "cohort", ]
         ), type = "lpmatrix")
        
        surv.curves <- surv.mcmc[random.curves[1:nrow(Xp)],]
        narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * ilink.surv(rowSums(Xp*surv.curves))
        
        # Formerly using sapply which is much slower
        # sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.surv(Xp[r,] %*% surv.curves[r,])}) 

        # microbenchmark::microbenchmark(
        # e1 = {        ilink.surv(rowSums(Xp*surv.curves))},
        # e2 = {         sapply(seq_len(nrow(Xp)), FUN = function(r){ilink.surv(Xp[r,] %*% surv.curves[r,])}) },
        # times = 100
        # )

        # narw.indiv[i, "p_surv", ] <- narw.indiv[i, "alive", ] * ilink.surv(Xp %*% surv.mcmc[sample(n.mcmc, 1, replace = TRUE),])
        
        
        # ps <- narw.indiv[i, "alive", ] *
        #   (surv_preds[["0"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 0) +
        #      surv_preds[["1"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 1) +
        #      surv_preds[["2"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 2) +
        #      surv_preds[["3"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 3) +
        #      surv_preds[["4"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 4) +
        #      surv_preds[["5"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 5) +
        #      surv_preds[["6"]](narw.indiv[i, "bc", ]) * (narw.indiv[i, "cohort", ] == 6))
        # 
        #   # 0.95 * narw.indiv[i, "alive", ]
        # 
        # narw.indiv[i, "p_surv", ] <- ps
        # narw.indiv[i, "p_surv", ] <-  0.85 * narw.indiv[i, "alive", ]
        
        #' ----------------------------
        # TOTALS
        #' ----------------------------
        
        # Number of  in each cohort
        # Calves (male)
        narw.pop[prj, i, 1] <- sum((narw.indiv[i, "cohort", ] == 0) *
                                     (narw.indiv[i, "female", ] == 0) *
                                     (narw.indiv[i, "alive", ] == 1))
        
        # Juveniles and adults (male)
        narw.pop[prj, i, 2] <- sum((narw.indiv[i, "cohort", ] == 1) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 3] <- sum((narw.indiv[i, "cohort", ] == 3) * (narw.indiv[i, "alive", ] == 1))
        
        # Calves (female)
        narw.pop[prj, i, 4] <- sum((narw.indiv[i, "cohort", ] == 0) *
                                     (narw.indiv[i, "female", ] == 1) *
                                     (narw.indiv[i, "alive", ] == 1))
        
        # Juvenile and reproductive adults (female)
        narw.pop[prj, i, 5] <- sum((narw.indiv[i, "cohort", ] == 2) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 6] <- sum((narw.indiv[i, "cohort", ] == 4) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 7] <- sum((narw.indiv[i, "cohort", ] == 5) * (narw.indiv[i, "alive", ] == 1))
        narw.pop[prj, i, 8] <- sum((narw.indiv[i, "cohort", ] == 6) * (narw.indiv[i, "alive", ] == 1))
        
        # Total population size
        tot.pop[prj, i] <- sum(narw.indiv[i, "alive", ], na.rm = TRUE)
        
      } # End totpop >0
    } # End years
    
    narw.individuals[[prj]] <- narw.indiv
    
  } # End projections
  
  end.time <- Sys.time()
  run_time <- hms::round_hms(hms::as_hms(difftime(time1 = end.time, time2 = start.time, units = "auto")), 1)
  
  #'------------------------------------------------------
  # SUMMARY ----
  #'------------------------------------------------------

  console("Summarizing outputs")
  
  # **** Number of births per female ****
  
  births.per.female <- purrr::map(.x = narw.individuals, .f = ~ {
    # For each projection
    out <- apply(.x, 3, function(x) {
      # Filter data to retain instances of pregnant females being alive
      fem <- x[x[, "alive"] == 1 & x[, "female"] == 1 & x[, "cohort"] == 4 & x[, "reprod"] == 1, , drop = FALSE]
      # Remove NAs (from new individuals born during the simulation)
      fem <- fem[complete.cases(fem), ]
      # Calculate the number of births 
      ifelse(nrow(fem) == 0, NA, sum(fem[, "birth"], na.rm = TRUE))
    })
    if(inherits(out, "list")) out <- do.call(c, out)
  }) |> do.call(what = c)
  
  # **** Time spent in resting phase ****
  
  time.resting <- purrr::map(.x = narw.individuals, .f = ~ {
      lapply(X = seq_len(dim(.x)[3]), FUN = function(d) {
        tr <- .x[,,d]
        tr <- tr[tr[,"alive"] == 1 & tr[,"female"]==1 & tr[, "cohort"] == 6 & tr[, "reprod"] == 1,,drop =FALSE]
        tr <- tr[, "time_rest"]
        tr[tr == 0] <- NA
        unname(purrr::map_dbl(.x = split_byNA(tr), .f = ~max(.x)))
      }) |> do.call(what = c)
    }) |> do.call(what = c)

  # **** Inter-birth interval ****

  inter.birth <- purrr::map(.x = narw.individuals, .f = ~ {
    lapply(1:dim(.x)[3], FUN = function(i) {
      x <- .x[,,i]
      # Filter out NA records for individuals born during the simulation
      x <- x[complete.cases(x),,drop = FALSE]
      # Exclude juveniles
      birth.record <- x[x[, "alive"] == 1 & x[, "female"] == 1 & x[, "cohort"] > 2 & x[, "reprod"] == 1, , drop = FALSE]
      birth.record <- birth.record[, "birth"]
      if (length(birth.record) == 0 | all(birth.record == 0) | all(is.na(birth.record))) {
        NA
      } else {
        # Extract timeline of births
        birth.record <- birth.record[tapply(seq_along(birth.record), birth.record, min)[2]:tapply(seq_along(birth.record), birth.record, max)[2]]
        birth.record[birth.record == 1] <- NA
        unname(purrr::map_dbl(.x = split_byNA(birth.record), .f = ~ length(.x)))
      }
    }) |> do.call(what = c)
  }) |> do.call(what = c)
  
  # **** Non-reproductive females ****
  
  nonreprod.females <- purrr::map(.x = seq_len(n), .f = ~ {

    # Non-reproductive females (alive) by year
    nonrep <- lapply(X = seq_len(dim(narw.individuals[[.x]])[3]), FUN = function(d) {
      rep.status <- rep(0, yrs+1)
      rep.status[narw.individuals[[.x]][,,d][, "female"] == 1 & narw.individuals[[.x]][,,d][, "alive"] == 1 & narw.individuals[[.x]][,,d][, "reprod"] == 0] <- 1
      rep.status
    }) |> do.call(what = rbind) |> 
      colSums()

    # which(purrr::map_dbl(.x = nonrep, .f = ~.x[1]==1)==1)
    
    # Total females (by year)
    totfem <- rowSums(narw.pop[.x,,5:8])
    totfem[1] <- sum(N_0[6:8])
    nonrep/totfem
    
  }) |> do.call(what = rbind)
  
  # narw.out <- purrr::map(.x = 1:n, .f = ~{
  #   reshape_array(narw.indiv[[.x]], value, yr, attr, whale) |> 
  #     dplyr::mutate(attr = mat.attribs[attr]) |> 
  #     tidyr::pivot_wider(names_from = attr, values_from = value) |> 
  #     dplyr::mutate(prj = .x) |> 
  #     dplyr::relocate(prj, .before = yr)
  # }) |> do.call(what = rbind) |> 
  #   data.table::data.table()
  
  # narw.out <- narw.out[is.finite(rowSums(narw.out)),]
  
  
# microbenchmark::microbenchmark(
# e1 = { 
#   
#   m1 <- apply(narw.pop, c(2,3), FUN = function(x) mean(x, na.rm = TRUE))
#   l1 <- apply(narw.pop, c(2,3), FUN = function(x) quantile(x, 0.025, na.rm = TRUE))
#   u1 <- apply(narw.pop, c(2,3), FUN = function(x) quantile(x, 0.975, na.rm = TRUE))
#   
#   test <- data.frame(year = current.yr + as.numeric(gsub("yr ", "", 0:yrs)))
#   
#   out <- purrr::map(.x = cohorts.proj$name,
#              .f = ~{
#               cbind(test, .x, m1[, .x], l1[, .x], u1[, .x])
#              }) |> do.call(what = "rbind")
#  
  
  # 
  # test <- data.frame(year = current.yr + as.numeric(gsub("yr ", "", 0:yrs)),
  #                         cohort = "Adults (female, resting)",
  #                         mean = apply(narw.pop, c(2), FUN = function(x) mean(x, na.rm = TRUE)),
  #                         lwr = apply(narw.pop, c(2), FUN = function(x) quantile(x, 0.025, na.rm = TRUE)),
  #                         uppr = apply(narw.pop, c(2), FUN = function(x) quantile(x, 0.975, na.rm = TRUE)))


# },
# e2 = { },
# times = 10
# )
  
  narw.df <- purrr::map(.x = cohorts.proj$name, .f = ~{
    tmp <- narw.pop[,,.x] |>
      tibble::as_tibble() |> 
      tibble::rownames_to_column(var = "prj")
    if(n == 1){
      tmp |>
        dplyr::rename(year = "prj") |> 
        dplyr::mutate(prj = 1, year = as.numeric(year) - 1) |>
        dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::rename(N = "value") |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x)) |> 
        dplyr::relocate(prj, .before = year)
    } else {
      tmp |>
        tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
        dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
        dplyr::mutate(cohort = stringr::str_to_sentence(.x))
    }
  }) |> do.call(what = rbind) |> 
    data.table::data.table()
 
  tot.df <- tibble::as_tibble(tot.pop) |> 
    tibble::rownames_to_column(var = "prj") |> 
    tidyr::pivot_longer(!prj, names_to = "year", values_to = "N") |> 
    dplyr::mutate(year = current.yr + as.numeric(gsub("yr ", "", year))) |> 
    dplyr::mutate(cohort = "North Atlantic right whales") |> data.table::data.table()
  
  narw.conf <- narw.df[
    , list(
      mean = mean(N, na.rm = TRUE),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>
    dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)"
    )))
  
  tot.conf <- tot.df[
    , list(
      mean = mean(N, na.rm = TRUE),
      lwr = quantile(N, 0.025, na.rm = TRUE),
      uppr = quantile(N, 0.975, na.rm = TRUE)
    ),
    list(year, cohort)
  ] |>dplyr::mutate(cohort = factor(cohort, levels = c(
      "Calves (male)",
      "Calves (female)",
      "Juveniles (male)",
      "Juveniles (female)",
      "Adults (male)",
      "Adults (female, pregnant)",
      "Adults (female, lactating)",
      "Adults (female, resting)",
      "North Atlantic right whales"
    )))

  console("Summarizing outputs", TRUE)
  # cat("\r", "Summarizing outputs ✓  ", sep = "")
  # cat("\n")
  cat(paste0("Time elapsed: ", run_time))
  cat("\n")
  
  outproj <- list(param = list(yrs = yrs,
                               n = n,
                               current.yr = current.yr,
                               abort = abort.rate,
                               nonrep = nonrep),
                  init = list(N_0 = N_0,
                              reprod.fem = reprod.fem),
                  dat = list(ind = narw.individuals,
                             birth = list(tot = births, perfemale = births.per.female, inter = inter.birth),
                             nonrepfem = nonreprod.females,
                             rest = time.resting,
                             pop = narw.pop,
                             tot = tot.pop),
                  prj = list(proj = rbind(narw.df, tot.df) |> dplyr::mutate(prj = as.numeric(prj)), 
                             mean = rbind(narw.conf, tot.conf)))
  class(outproj) <- c("narwproj", class(outproj))
  return(outproj)
  
  ## Add % females that are not reproducing
  ## 
  
}
