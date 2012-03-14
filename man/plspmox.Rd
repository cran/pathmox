\name{plspmox}
\alias{plspmox}
\title{PLSPM of a node from a PATHMOX or TECHMOX tree}
\description{
  Calculates a plspm on a selected node from a pathmox or techmox tree
}
\usage{
  plspmox(pls, treemox, X = NULL, node, boot.val = FALSE, br = NULL, dataset = FALSE)
}
\arguments{
  \item{pls}{An object of class \code{"plspm"}}
  \item{treemox}{An object of class \code{"treemox"}}
  \item{X}{Optional argument for data table}
  \item{node}{An integer value indicating the number of the node. Must be an integer larger than 1.}
  \item{boot.val}{A logical value indicating whether bootstrap validation is performed (\code{FALSE} by default).}
  \item{br}{An integer indicating the number bootstrap resamples. Used only when \code{boot.val=TRUE}.}
  \item{dataset}{A logical value indicating whether the data matrix should be included in the list of results (\code{FALSE} by default).}
}
\details{
  Performs a PLS-PM analysis with the elements contained in \code{node}, by calling the function \code{\link{plspm}}.
  The rest of the parameters to perform the PLS-PM analysis (i.e. \code{inner}, \code{outer}, \code{modes}, \code{scheme}, \code{scaled}, \code{tol}, \code{iter}) are inherited from the object in argument \code{pls}. \cr
  When the object \code{pls} does not contain a data matrix (i.e. \code{pls$data=NULL}), the user must provide the data matrix or data frame in \code{X}.
}
\value{
  An object of class \code{"plspm"}
}
\author{
  Gaston Sanchez
}
\seealso{
  See Also as \code{\link{plspm}}
}
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

  node2 <- plspmox(pls.mob, mob.pathmox, node=2)
  node3 <- plspmox(pls.mob, mob.pathmox, node=3)
  }
}
