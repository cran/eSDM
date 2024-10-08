% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eSDM-package.R
\docType{package}
\name{eSDM-package}
\alias{eSDM-package}
\alias{eSDM}
\title{Ensemble tool for predictions from Species Distribution Models}
\description{
eSDM: A tool for creating and exploring ensembles of predictions from Species
Distribution Models
}
\details{
eSDM provides functionality for overlaying SDM predictions onto a single base
geometry and creating and evaluating ensemble predictions. This can be done
manually in R, or using the eSDM GUI (an R Shiny app) opened through
\link{eSDM_GUI}

eSDM allows users to overlay SDM predictions onto a single base
  geometry, create ensembles of these predictions via weighted or unweighted
  averages, calculate performance metrics for each set of predictions and for
  resulting ensembles, and visually compare ensemble predictions with
  original predictions. The information provided by this tool can assist
  users in understanding spatial uncertainties and making informed
  conservation decisions.

  The GUI ensures that the tool is accessible to non-R users, while also
  providing a user-friendly environment for functionality such as loading
  other polygons to use and visualizing predictions. However, user choices
  are restricted to the workflow provided by the GUI.
}
\seealso{
\url{https://swfsc.github.io/eSDM/}
}
\author{
Sam Woodman \email{sam.woodman@noaa.gov}
}
\keyword{package}
