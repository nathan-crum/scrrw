% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helper_functions.R
\name{meanMaxMove}
\alias{meanMaxMove}
\title{Calculate the mean of the maximum pariwise distance between each individuals' observations}
\usage{
meanMaxMove(y.pix, distMat)
}
\arguments{
\item{y.pix}{A matrix; NA when individual -row- not observed on -occasion-; cell ID where observed otherwise}

\item{distMat}{A matrix; holds pairwise distances between all grid cells in the habitat matrix}
}
\description{
Calculate the mean of the maximum pariwise distance between each individuals' observations
}
