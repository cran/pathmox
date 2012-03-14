\name{treemox.pls}
\alias{treemox.pls}
\title{PLS-PM results of terminal nodes from a PATHMOX or TECHMOX tree}
\description{
  Calculates basic PLS-PM results for the terminal nodes of PATHMOX and TECHMOX trees
}
\usage{
  treemox.pls(pls, treemox, X = NULL)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"} returned by \code{\link{plspm}}.}
  \item{treemox}{An object of class \code{"treemox"} returned by \code{\link{pathmox}} or \code{\link{techmox}}.}
  \item{X}{Optional dataset (matrix or data frame) used when argument \code{dataset=NULL} inside \code{pls}.}
}
\details{
  The argument \code{pls} must be the same used for calculating the \code{treemox} object.
  When the object \code{pls} does not contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or data frame in \code{X}.
}
\value{
  An object of class \code{"treemox.pls"}. Basically a list with the following results:
  \item{weights}{Matrix of outer weights for each terminal node}
  \item{loadings}{Matrix of loadings for each terminal node}
  \item{paths}{Matrix of path coefficients for each terminal node}
  \item{r2}{Matrix of r-squared coefficients for each terminal node}
}
\author{
  Gaston Sanchez
}
\seealso{\code{\link{pathmox}}, \code{\link{techmox}}, \code{\link{plot.treemox}}.}
\examples{
  \dontrun{
  ## example of PLS-PM in customer satisfaction analysis
  ## model with seven LVs and reflective indicators
  data(csimobile)
  data.mobile <- csimobile[,8:33]  ## select the manifest variables
  IMAG <- c(0,0,0,0,0,0,0)
  EXPE <- c(1,0,0,0,0,0,0)
  QUAL <- c(0,1,0,0,0,0,0)
  VAL <- c(0,1,1,0,0,0,0)
  SAT <- c(1,1,1,1,0,0,0)
  COM <- c(0,0,0,0,1,0,0)
  LOY <- c(1,0,0,0,1,1,0)
  inner.mob <- rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)
  outer.mob <- list(1:5,6:9,10:15,16:18,19:21,22:24,25:26)
  modes.mob <- rep("A",7)
  pls.mob <- plspm(data.mobile, inner.mob, outer.mob, modes.mob, scheme="factor", scaled=FALSE)

  # re-ordering those segmentation variables with ordinal scale (Age and Education)
  csimobile$Education <- factor(csimobile$Education, 
    levels=c("basic","highschool","university"),
    ordered=TRUE)
  # select the segmentation variables
  EXEV.mob = csimobile[,1:7]
  # Pathmox Analysis
  mob.pathmox <- pathmox(pls.mob, EXEV.mob, signif=.10, size=.10, deep=2)
  # applying function treemox.pls
  mob.nodes <- treemox.pls(pls.mob, mob.pathmox)
  # comparative barplots
  plot(mob.nodes)
  }
}