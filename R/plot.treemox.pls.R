plot.treemox.pls <-
function(x, comp.by="nodes", nodes.names=NULL, ordered=TRUE, decreasing=FALSE, 
       color=NULL, show.box=TRUE, border=NA, cex.names=.75, cex.axis=.75, short.labs=TRUE, short.min=NULL, ...)
{
    # ARGUMENTS
    # x: an object of class treemox.pls
    # comp.by: one of "nodes" or "latents"
    # nodes.names: optional vector of names for the nodes (length of terminal nodes)
    # ordered: logical value indicating whether bars are ordered (increasing order)
    # decreasing: logical value indicating if the sort order should be increasing or decreasing
    # color: optional vector of colors for the bars
    # show.box: logical value indicating whether a box is drawn around each barplot
    # border: color of the borders
    # cex.names: expansion factor of labels (x-axis)
    # cex.axis: expansion factor of y-axis
    # short.labs: logical value indicating if labels should be abbreviated
    # short.min: integer number indicating the minimum length of the abbreviations

    if (length(union(comp.by,c("nodes","latents")))>2)
        stop("Invalid argument 'comp.by'")
    if (!is.null(nodes.names)) {
        if (length(nodes.names)!=(ncol(x$paths)-1)) {
            cat("NOTICE: Invalid argument 'nodes.names'. Default names were used", "\n")   
        } else {
            colnames(x$paths) <- c("Root_Node",nodes.names)
        }        
    }
    if (short.labs)
        if (mode(short.min)!="numeric" || length(short.min)!=1 || (short.min%%1)!=0)
            short.min <- 7

    IDM <- x$IDM
    lvs <- nrow(IDM)
    idm <- as.vector(IDM)
    idm[idm==1] <- 1:sum(IDM)
    SDM <- matrix(idm,lvs,lvs)
    B <- x$paths[,-1] - x$paths[,1]
    rs <- c(1,1,1,2,2,2,2,2,3,2,3,3,4,4)
    cs <- c(1,2,3,2,3,3,4,4,3,5,4,4,4,4)
    index.mat <- cbind(1:14,rs,cs)
    endo <- rowSums(IDM)
    endo[endo!=0]<-1
    who.endo <- which(endo!=0)

    if (comp.by=="nodes")
    {    
        nn <- ncol(B)
        if (is.null(color))  colors=rainbow(lvs,s=.5,v=.7)  else  colors=rep(color,length.out=lvs)
        for (i in who.endo)
        {
            n.ind <- which(IDM[i,]==1)
            indep <- SDM[i,which(SDM[i,]!=0)]
            arg.names <- colnames(IDM)[n.ind]
            if (short.labs) 
                arg.names <- abbreviate(arg.names, minlength=short.min)
            cols <- colors[n.ind]            
            dev.new()
            par(mfrow=index.mat[nn+1,2:3])
            par(mar=c(3,3,3,3)) 
            barplot(x$paths[indep,1], names.arg=arg.names, main=paste("Global:",rownames(IDM)[i]), col=cols,
                    border=border)
            abline(h=0)
            if (show.box) box("figure")
            ylim <- 1.15 * c(min(B[indep,]),max(B[indep,]))  
            for (n in 1:nn)
            {   
                if (ordered) {
                    barplot(sort(B[indep,n], decreasing), names.arg=arg.names[order(B[indep,n], decreasing=decreasing)], 
                            main=colnames(B)[n], border=border, col=cols[order(B[indep,n], decreasing=decreasing)], 
                            cex.names=cex.names, cex.main=1, cex.axis=cex.axis, ylim=ylim, ...)
                } else {
                     barplot(B[indep,n], names.arg=arg.names, main=colnames(B)[n], border=border, 
                            col=cols, cex.names=cex.names, cex.main=1, cex.axis=cex.axis, ylim=ylim, ...)
                }
                abline(h=0)
                if (show.box) box("figure")
            }   
        }
    } else {
        # comp.by="latents"
        ncols <- ncol(B)
        if (is.null(color))  colors=rainbow(ncols,s=.5,v=.7)  else  colors=rep(color,length.out=ncols)
        arg.names <- colnames(B)
        if (short.labs) 
            arg.names <- abbreviate(arg.names, minlength=short.min)
        for (i in who.endo)
        {
            indep <- SDM[i,which(SDM[i,]!=0)]    
            dev.new()
            par(mfrow=index.mat[length(indep),2:3])
            par(mar=c(3,3,3,3)) 
            ylim <- 1.15 * c(min(B[indep,]),max(B[indep,]))
            for (k in indep)
            {
                if (ordered) {
                    barplot(sort(B[k,], decreasing), names.arg=arg.names[order(B[k,], decreasing=decreasing)], 
                            main=rownames(B)[k], border=NA, cex.main=1,cex.names=cex.names, cex.axis=cex.axis, 
                            ylim=ylim, col=colors[order(B[k,], decreasing=decreasing)], ...)
                } else {
                    barplot(B[k,], names.arg=arg.names, main=rownames(B)[k], border=NA, cex.main=1,
                            cex.names=cex.names, cex.axis=cex.axis, ylim=ylim, col=colors, ...)
                }                  
                abline(h=0)
                if (show.box) box("figure")
            }   
        }
    } 
}

