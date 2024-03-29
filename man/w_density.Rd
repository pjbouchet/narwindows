% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{w_density}
\alias{w_density}
\title{#' Initialize model
#'
#' @export
load <- function(){
  Rcpp::sourceCpp("src/simtools.cpp")
  # Rcpp::sourceCpp("src/bioenergetics_functions.cpp")
  source("R/run_model.R")
  assign("init.model", value = TRUE, envir = .GlobalEnv)
}}
\usage{
w_density(target = "SEUS", option = 1)
}
\description{
#' Initialize model
#'
#' @export
load <- function(){
  Rcpp::sourceCpp("src/simtools.cpp")
  # Rcpp::sourceCpp("src/bioenergetics_functions.cpp")
  source("R/run_model.R")
  assign("init.model", value = TRUE, envir = .GlobalEnv)
}
}
