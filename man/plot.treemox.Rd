\name{plot.treemox}
\alias{plot.treemox}
\title{Plot PATHMOX and TECHMOX trees}
\description{
  Plot method for objects of class \code{"treemox"}  
}
\usage{
  \method{plot}{treemox}(x, root.col = "grey", node.col = "orange", leaf.col = "green2", shadow.size = 0.003, 
      node.shadow = "red", leaf.shadow = "darkgreen", cex = 0.7, seg.col = "blue3", lwd = 1,
      show.pval = TRUE, pval.col = "blue", main = NULL, cex.main = 1, \dots)
}
\arguments{
  \item{x}{An object of class \code{"treemox"} returned by \code{\link{pathmox}} or \code{\link{techmox}}.}
  \item{root.col}{Fill color of root node.}
  \item{node.col}{Fill color of child nodes.}
  \item{leaf.col}{Fill color of leaf nodes.}
  \item{shadow.size}{Relative size of shadows.}
  \item{node.shadow}{Color of shadow of child nodes.}
  \item{leaf.shadow}{Color of shadow of leaf nodes.}
  \item{cex}{A numerical value indicating the magnification to be used for plotting text.}
  \item{seg.col}{The color to be used for the labels of the segmentation variables.}
  \item{lwd}{Line width of the tree branches.}
  \item{show.pval}{Logical value indicating whether the p-values should be plotted.}
  \item{pval.col}{The color to be used for the labels of the p-values.}
  \item{main}{A main title for the plot.}
  \item{cex.main}{The magnification to be used for the main title.}
  \item{\dots}{Further arguments are ignored.}
}
\details{
  The function \code{plot.treemox} allows to display binary trees of PATHMOX and TECHMOX analyses. \cr
  If \code{shadow.size=0}, no shadows are drawn. \cr
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
  mob.pathmox <- pathmox(pls.mob, EXEV.mob, signif=.10, size=.10, deep=2, tree=FALSE)
  
  plot(mob.pathmox, root.col="lightblue", node.col="turquoise", leaf.col="skyblue3", 
       shadow.size=0, seg.col="blue2", pval.col="magenta")
  }
}
