\name{plot.treemox.pls}
\alias{plot.treemox.pls}
\title{Comparative plot between nodes from a PATHMOX or TECHMOX tree}
\description{
  Plot method for objects of class \code{"treemox.pls"}. Barplots of path coefficients of terminal nodes with respect to those of the global (root) model
}
\usage{
  \method{plot}{treemox.pls}(x, comp.by = "nodes", nodes.names = NULL, ordered = TRUE, decreasing = FALSE, color = NULL, 
     show.box = TRUE, border = NA, cex.names = 0.75, cex.axis = 0.75, short.labs = TRUE, short.min = NULL, \dots)
}
\arguments{
  \item{x}{An object of class \code{"treemox.pls"} returned by \code{\link{treemox.pls}}.}
  \item{comp.by}{One of "nodes" or "latents". This argument indicates the type of barplots comparison.}
  \item{nodes.names}{Optional vector of names for the terminal nodes (must be a vector of length equal to the number of terminal nodes).}
  \item{ordered}{A logical value indicating whether the barplots are shown in increasing ordered.} 
  \item{decreasing}{A logical value indicating if the sort order should be increasing or decreasing.}
  \item{color}{Optional vector of colors for the bars. When \code{color=NULL} rainbow colors are used.}
  \item{show.box}{A logical value indicating whether a box is drawn around each barplot.}
  \item{border}{The color to be used for the border of the bars. Use \code{border=NA} to omit borders.}
  \item{cex.names}{Expansion factor for axis names (bar labels).}
  \item{cex.axis}{Expansion factor for numeric axis labels.}
  \item{short.labs}{Logical value indicating if the labels of the barplots should be abbreviated (\code{TRUE} by default).}
  \item{short.min}{Integer number indicating the minimum length of the abbreviations for the labels. Only used when \code{short.labs=TRUE}.}
  \item{\dots}{Arguments to be passed to/from other methods.}
}
\details{
  This function aims to visualize the comparison between path coefficients of the terminal nodes aginst the path coefficients of the global model in the root node. \cr
  When \code{comp.by="nodes"} a graphic window is displayed for each endogenous latent variable of the PLS model, and barplots of nodes are shown. \cr
  When \code{comp.by="latents"} a graphic window is displayed for each endogenous relationship of the PLS model, and barplots of independent latent variables are shown.
}
\seealso{
  \code{\link{treemox.pls}}, \code{\link{pathmox}}, \code{\link{techmox}}
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
  # applying function treemox.pls
  mob.nodes <- treemox.pls(pls.mob, mob.pathmox)

  # comparative barplots of endogenous latent variables between nodes
  plot(mob.nodes, comp.by="nodes")  ## default plot

  # comparative barplots of nodes between latent regressors
  plot(mob.nodes, comp.by="latents", decreasing=TRUE)
  }
}