#' Summary
#'
#' Summary information
#' @import data.table
#' @export
#' 
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' library(narwind)
#' 
#' animals <- run_model(10)
#' plot(animals)
#' }

summary.narwproj <- function(obj,
                             ...){
  
  options(pillar.sigfig = 7)
  
  # Function ellipsis –– optional arguments
  args <- list(...)
  
 # Preamble ---------------------------------------------------------------
  
  if(!inherits(obj, "narwproj")) stop("Object must be of class <narwproj>")
  
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n")
  cat("\n")
  cat("     NORTH ATLANTIC RIGHT WHALE (Eubalaena glacialis)\n\n")
  cat("          *** POPULATION PROJECTION SUMMARY ***\n")
  cat("\n")
  cat("-------------------------------------------------------------\n")
  cat("-------------------------------------------------------------\n\n")

  tot.df <- obj$prj$proj
  current.yr <- obj$param$current.yr
  yrs <- obj$param$yrs
  N <- obj$param$n
  births.per.female <- obj$dat$birth$perfemale
  tot.births <- obj$dat$birth$tot 
  time.resting <- obj$dat$rest
  inter.birth <- obj$dat$birth$inter
  nonreprod.females <- obj$dat$nonrepfem
  N_0 <- obj$init$N_0
  
  cat("Replicates: N =", N, "\n")
  cat("Projection horizon:", yrs, "years\n\n")
  
  cat("=============================================================\n")
  cat("ABUNDANCE\n")
  cat("=============================================================\n\n")
  
  cat("Initial population size:\n")
  cat("N =", sum(N_0))
  init.pop <- tibble::tibble(cohort = names(N_0), N = N_0) |> 
    janitor::adorn_totals()
 
  print(knitr::kable(init.pop, format = "simple"))
  cat("\n")
   
  # Find 95% confidence intervals on final population size
  final.pop <- tot.df[year == current.yr + yrs, list(N = quantile(N, probs = c(0.5, 0.025, 0.975), na.rm = TRUE)), cohort]
  final.pop[, value := rep(c("median", "lower", "upper"), length.out = nrow(final.pop))]
  final.pop <- final.pop |> tidyr::pivot_wider(names_from = value, values_from = N) |> data.table::data.table()
  final.pop[, N:= paste0(format(round(median,0), big.mark = ","), " (95% CI: ", format(round(lower,0), big.mark = ","), " – ", format(round(upper,0), big.mark = ","), ")")]
  
  cat("Final population size:\n")
  cat("N =", final.pop[cohort == "North Atlantic right whales", N])
  
  final.pop <- tibble::tibble(final.pop[, list(cohort, N)]) |> 
    dplyr::mutate(N = trimws(stringr::str_replace_all(N, "\\h+", " ")))
  
  print(knitr::kable(final.pop[!final.pop$cohort == "North Atlantic right whales", ], format = "simple"))
  
  cat("\n=============================================================\n")
  cat("FECUNDITY\n")
  cat("=============================================================")
  
  # cat("Calving events [per year]:\n")
  calving.events <- tibble::tibble(`Calving events` = c(paste0("Per year: N = ", median(tot.births, na.rm = TRUE), " ± ",
                                                      round(sd(tot.births, na.rm = TRUE),0), " [",
                                                      min(tot.births, na.rm = TRUE), "–",
                                                      max(tot.births, na.rm = TRUE), "]"),
                                              paste0("Per female: N = ", paste0(median(births.per.female, na.rm = TRUE), " ± ",
                                                                   round(sd(births.per.female, na.rm = TRUE),0), " [",
                                                                   min(births.per.female, na.rm = TRUE), "–",
                                                                   max(births.per.female, na.rm = TRUE), "]"))))
  print(knitr::kable(calving.events, format = "simple"))
  
  resting.period <- tibble::tibble(`Resting phase` = c(paste0(" t(rest): ", paste0(median(time.resting, na.rm = TRUE), " ± ",
                                                                                   round(sd(time.resting, na.rm = TRUE),0), " [",
                                                                                   min(time.resting, na.rm = TRUE), "–",
                                                                                   max(time.resting, na.rm = TRUE), "]")),
                                                       paste0(" t(inter-birth): ", paste0(median(inter.birth, na.rm = TRUE), " ± ",
                                                                                          round(sd(inter.birth, na.rm = TRUE),0), " [",
                                                                                          min(inter.birth, na.rm = TRUE), "–",
                                                                                          max(inter.birth, na.rm = TRUE), "]"))))
                                                     
  print(knitr::kable(resting.period, format = "simple"))
  
  cat("\n\n")
  cat(" Abortion rate:", 100 * obj$param$abort, "%")
  cat("\n")
  
  nrf <- nonreprod.females |> 
    tibble::as_tibble() |> 
    dplyr::mutate(proj = dplyr::row_number()) |> 
    tidyr::pivot_longer(!proj, names_to = "year", values_to = "nonrep") |> 
    dplyr::mutate(year = as.numeric(gsub("yr ", "", year))) |> 
    data.table::as.data.table()
  
  nrf.tbl <- nrf[, list(mean = round(100*mean(nonrep), 2),
                        sd = round(100*sd(nonrep), 2),
                        min = round(100*min(nonrep), 2),
                        max = round(100*max(nonrep), 2)), ] |> 
    dplyr::mutate(`Non-reproductive females` = paste0(mean, " (±", sd, ") [", min, "–", max, "]")) |> 
    dplyr::select(-mean, -sd, -min, -max) |> 
    dplyr::pull(`Non-reproductive females`)
  
  cat(" Non-reproductive females:", as.character(nrf.tbl), "%")
  cat("\n")
  

  # nrf.tbl <- dplyr::bind_rows(nrf.tbl, tibble::tibble(`% Non-reproductive females` = "Per year: See plot"))
  
  # print(knitr::kable(nrf.tbl, format = "simple"))
  
  # nrf.byyear <- nrf[, list(mean = round(100*mean(nonrep), 2),
  #                          sd = round(100*sd(nonrep), 2),
  #                          min = round(100*min(nonrep), 2),
  #                          max = round(100*max(nonrep), 2)), year] |>
  #   dplyr::mutate(` ` = paste0(format(mean, digits = 3), " (±", format(sd, digits = 3), ") [", format(min, digits = 3), " –", format(max, digits = 3), "]")) |>
  #   dplyr::select(-mean, -sd, -min, -max)
  
  # print(knitr::kable(nrf.byyear, format = "simple"))
  
  }
