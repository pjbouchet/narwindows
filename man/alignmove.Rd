% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{alignmove}
\alias{alignmove}
\title{Decrease in milk feeding efficiency as a function of calf age}
\usage{
alignmove(
  mobj,
  res = "min",
  digit = "min",
  unit = "days",
  spaceMethod = "greatcircle"
)
}
\arguments{
\item{a}{Calf age in days.}

\item{Tstop}{Time at which milk stops nursing and becomes nutritionally independent (in days).}

\item{Tdecrease}{Time at which milk comsuption starts to decrease (in days).}

\item{E}{Steepness of the non-linear relationship.}

\item{starvation}{Starvation threshold
Title}

\item{eta}{Steepness of the non-linear relationship}

\item{target_condition}{}

\item{maintenance_mass}{}

\item{blubber_mass}{Feeding effort as function of copepod density}

\item{gamma}{Feeding threshold}

\item{D}{Coepepod density}
}
\description{
Decrease in milk feeding efficiency as a function of calf age
}
