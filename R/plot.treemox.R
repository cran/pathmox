plot.treemox <-
function(x, root.col="grey", node.col="orange", leaf.col="green2",
       shadow.size=0.003, node.shadow="red", leaf.shadow="darkgreen", cex=.7, 
       seg.col="blue3", lwd=1, show.pval=TRUE, pval.col="blue", main=NULL, cex.main=1, ...)
{
    MOX <- x$MOX
    last <- nrow(MOX)
    last.level <- MOX$Depth[last]
    num.levels <- rep(1,last.level+1)
    for (i in 1:last.level)
         num.levels[i+1] <- 2^i
    ## plot parameters
    dev.new()
    par(mar=c(.4,.4,1,1.5))
    openplotmat()
    elpos <- coordinates(num.levels)# estimates coordinates of elements, arranged on a plot
    fromto <- cbind(MOX[-1,2], MOX[-1,1])
    nr <- nrow(fromto)
    arrpos <- matrix(ncol=2, nrow=nr)
    ## drawing the binary tree structure
    for (i in 1:nr)
        arrpos[i,]<- straightarrow(to=elpos[fromto[i,2],], from=elpos[fromto[i,1],],
             lwd=lwd, arr.pos=0.6, arr.length=0)
    ## drawing the root node (i.e. parent node)
    textellipse(elpos[1,], 0.045, 0.045, lab=c("Root",MOX[1,6]), box.col=root.col,
                shadow.size=shadow.size, cex=cex)
    ## drawing the child nodes
    for (i in 2:last)
    {
         posi <- MOX$Node[i]
         nodlab <- c(paste("Node",posi), MOX$Size[i])
         if (MOX$Type[i]=="node") {
             textellipse(elpos[posi,], 0.05, 0.03, lab=nodlab, box.col=node.col,
                         shadow.col=node.shadow, shadow.size=shadow.size, cex=cex)
         } else { # "leaf"
             textrect(elpos[posi,], 0.045, 0.025, lab=nodlab, box.col=leaf.col,
                      shadow.col=leaf.shadow, shadow.size=shadow.size, cex=cex)
         }
    }
    ## adding segmentation variables and p.values
    aux <- 1
    for (i in seq(1,nr,by=2))
    {   
        if (i==1) k<-1 else k<-1.15     
        x1 <- (arrpos[i,1] + arrpos[i+1,1]) / 2
        text(x1, k*arrpos[i,2], MOX$Variable[i+1], cex=cex, col=seg.col)
        if (show.pval)
        {
            if (x$model$mox=="pathmox")
                text(x1, k*arrpos[i,2], paste("p.val=",round(x$FT$p.val[aux],4),sep=""), cex=.9*cex, col=pval.col, pos=1) 
            if (x$model$mox=="techmox")
            {
                pvals <- unlist(lapply(x$FT, function(x) exp(mean(log(x[,4]))) ))
                text(x1, k*arrpos[i,2], paste("pv.gm=",round(pvals[aux],4),sep=""), cex=.9*cex, col=pval.col, pos=1) 
            }
        }
        aux <- aux + 1
    }     
    ## adding segmentation categories
    for (i in 1:nr)
    {
        posi <- MOX$Node[i+1]        
        seg.cat <- as.character(MOX$Category[i+1])
        seg.cat <- unlist(strsplit(seg.cat,"/"))
        if (posi%%2==0) {
            for (h in 1:length(seg.cat))
                text(arrpos[i,1]-0.03, arrpos[i,2]+h/55, seg.cat[h], cex=cex)
        } else {
            for (h in 1:length(seg.cat))
                text(arrpos[i,1]+0.03, arrpos[i,2]+h/55, seg.cat[h], cex=cex)
        }
    }
    ## adding main title
    if (is.null(main)) {
        if (x$model$mox=="pathmox")
            text(.5, .95, c("PATHMOX Tree"), cex=cex.main)
        if (x$model$mox=="techmox")
            text(.5, .95, c("TECHMOX Tree"), cex=cex.main)
    } else {
        text(.5, .95, main, cex=cex.main)
    }
}

