\name{print.bootnodes}
\alias{print.bootnodes}
\title{Printing bootnodes objects}
\description{
  Printing method for bootstrap nodes of PATHMOX and TECHMOX trees, i.e. objects of class \code{"bootnodes"}
}
\usage{
  \method{print}{bootnodes}(x, \dots)
}
\arguments{
  \item{x}{Object of class \code{"bootnodes"} returned by \code{\link{treemox.boot}}.}
  \item{\dots}{Further arguments are ignored}
}
\details{
  The function \code{print.bootnodes} displays the results of the bootstrap validation on PATHMOX and TECHMOX trees. 
}
\seealso{\code{\link{pathmox}}, \code{\link{techmox}},\code{\link{treemox.pls}},\code{link{treemox.boot}}}