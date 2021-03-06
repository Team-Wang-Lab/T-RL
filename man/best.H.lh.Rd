% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/best.H.lh.R
\name{best.H.lh}
\alias{best.H.lh}
\title{Search for best covariate for splitting based on lookahead.}
\usage{
best.H.lh(H, A, mus.hat, minsplit = 20)
}
\arguments{
\item{H}{A matrix of covariates before assigning final treatment, excluding previous treatment variables.}

\item{A}{A vector of observed treatment options.}

\item{mus.hat}{Estimated conditional mean outcome.}

\item{minsplit}{Minimal node size.}
}
\description{
Search for best covariate for splitting based on lookahead.
}
