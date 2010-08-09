\name{print.treemox}
\alias{print.treemox}
\title{Printing treemox objects}
\description{
  Printing method for objects of class \code{"treemox"}
}
\usage{
  \method{print}{treemox}(x, \dots)
}
\arguments{
  \item{x}{Object of class \code{"treemox"} returned by either \code{\link{pathmox}} or \code{\link{techmox}} .}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{print.treemox} displays the results of PATHMOX and TECHMOX trees. 
}
\seealso{\code{\link{pathmox}}, \code{\link{techmox}}, \code{\link{plot.treemox}}}