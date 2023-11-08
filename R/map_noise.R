#' Noise layer
#'
#' @export
#' 
#' @author Phil J. Bouchet
#' @examples
#' \dontrun{
#' }

map_noise <- function(obj = NULL){
  
  if(!inherits(obj, "narwscenario")) stop("Input object must be of class <narwscenario>.")
  
  # Spatial support
  support.poly <- targets::tar_read(support_poly)
  
  turbines.df <- obj$locs
  ambient.db <- obj$ambient
  source.lvl <- obj$sourceLvL
  attenuation <- obj$lowerdB

  date_seq <- get_dates()
  ambient.r <- dummy_raster(value = ambient.db)
  coords <- raster::coordinates(ambient.r)
  
  if(nrow(turbines.df) > 0){
    
    # Strip years from dates
    turbines.df$date <- format(as.Date(turbines.df$date), "%d-%m")
    
    # Compute sound levels from distances to turbine locations
    turbines.xy <- as.matrix(turbines.df[, c("x", "y")])
    pb <- progress::progress_bar$new(total = nrow(turbines.xy))
    cat("Calculating noise levels ...\n")
    
    turbines.df$soundLvL <- lapply(1:nrow(turbines.xy), FUN = function(x){
      pb$tick()
      r <- raster::distanceFromPoints(ambient.r, coords[raster::cellFromXY(ambient.r, turbines.xy[x,]),]) # Assume that turbine is at centroid
      # r <- raster::distanceFromPoints(ambient.r, x) # Or calculate turbine-centroid distance
      r[which(is.na(ambient.r[]))] <- NA
      dB <- km2dB(r, SL = source.lvl, mitigation = attenuation)
      dB[dB < ambient.db] <- ambient.db
      dB
    })
    
    # Piling dates & days of the year
    piling_dates <- unique(turbines.df$date)
    
    # Split by date
    sound.layer <- split(turbines.df, f = factor(turbines.df$date, levels = unique(turbines.df$date)))
    
    # Sum sound fields
    pb <- progress::progress_bar$new(total = length(sound.layer))
    cat("Summing sound fields ...\n")
    dB <- lapply(X = 1:length(sound.layer), FUN = function(x) {
      pb$tick()
      add_SPL(sound.layer[[x]]$soundLvL, raster = TRUE)
    })
    names(dB) <- piling_dates
    
    # Create yearly timeline
    sound.out <- purrr::map(.x = date_seq, .f = ~ {
      if (.x %in% names(dB)){
        dB[[.x]]
      } else {
        ambient.r
      }
    })
    names(sound.out) <- date_seq
    
  } else {
    
    sound.out <- purrr::map(.x = date_seq, .f = ~ ambient.r)
    names(sound.out) <- date_seq
    
  }
 
  sound.out <- purrr::map(.x = sound.out, .f = ~as(.x, "SpatialGridDataFrame"))
  return(sound.out)
  
} # End function
