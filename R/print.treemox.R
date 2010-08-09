print.treemox <-
function(x, ...)
{
    cat("SEGMENTATION TREES IN PLS-PM", "\n\n")
    cat("----------------------------------------------------------", "\n")    
    cat("TREE SPECIFICATION", "\n")
    cat("1   Mox Algorithm", "\t", x$model$mox, "\n")
    cat("2   Threshold signif", "\t", x$model$signif, "\n")
    if (x$model$size<1) {
        cat("3   Node size limit(%)",  "\t", x$model$size, "\n")
    } else {
        cat("3   Node size limit(#)",  "\t", x$model$size, "\n")
    }
    cat("4   Tree depth level",  "\t", x$model$deep, "\n")
    cat("\n")
    cat("Segmentation Variables", "\n")
    print(x$model$df.exev, print.gap=2)
    cat("----------------------------------------------------------", "\n")    
    cat("$MOX", "\n")
    print(x$MOX, digits=3, print.gap=2)
    invisible(x)
}

