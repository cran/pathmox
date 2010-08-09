treemox.pls <-
function(pls, treemox, X=NULL)
{
    # ====================== ARGUMENTS =====================
    # pls: an object of class plspm
    # treemox: an object of class treemox 
    # X: optional data matrix or dataframe

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm")
        stop("Argument 'pls' must be an object of class 'plspm'")
    if (class(treemox)!="treemox")
        stop("Argument 'treemox' must be an object of class 'treemox'")
    if (nrow(pls$latents)!=treemox$MOX$Size[1]) 
        stop("Arguments 'pls' and 'treemox' are incompatible. Different number of observations")
    if (!is.null(X)) # if X available
    {
        if (is.null(pls$data))
        {
            if (!is.matrix(X) && !is.data.frame(X))
                stop("Invalid object 'X'. Must be a numeric matrix or data frame.")
            if (nrow(X)!=nrow(pls$latents))
                stop("Argument 'pls' and 'X' are incompatible. Different number of rows.")
        }
    } else { # if no X
        if (is.null(pls$data)) 
            stop("Argument 'X' is missing. No dataset available.")
    }

    # ========================== INPUTS SETTING ==========================
    IDM <- pls$model$IDM
    blocks <- pls$model$blocks
    modes <- pls$model$modes
    scheme <- pls$model$scheme
    scaled <- pls$model$scaled
    tol <- pls$model$tol
    iter <- pls$model$iter
    outer <- pls$model$outer
    blocklist <- outer
    for (k in 1:length(blocks))
         blocklist[[k]] <- rep(k,blocks[k])
    blocklist <- unlist(blocklist)
    # data matrix DT
    if (!is.null(pls$data)) {
        DT <- pls$data
    } else {         
        # building data matrix 'DT'
        DT <- matrix(NA, nrow(pls$latents), sum(blocks))
        for (k in 1:nrow(IDM))
            DT[,which(blocklist==k)] <- as.matrix(X[,outer[[k]]])
        dimnames(DT) <- list(rownames(pls$latents), names(pls$out.weights))
    }
    MOX <- treemox$MOX
    list.nodes <- treemox$list.nodes
    term.nodes <- which(MOX$Terminal=="yes") - 1# identifying terminal nodes    
    tn <- length(term.nodes)# number of terminal nodes
    tn.labs <- paste("Node", MOX$Node[term.nodes+1], sep="_")# labels for terminal nodes
    lvs <- ncol(IDM)
    obs <- nrow(DT)
    mvs <- sum(blocks)
    mvs.names <- names(pls$out.weights)
    lvs.names <- row.names(IDM)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1

    TOW <- matrix(NA, mvs, 1+tn)# TOW: table of outer weights
    TL <- matrix(NA, mvs, 1+tn)# TL: table of loadings
    TR2 <- matrix(NA, lvs, 1+tn)# TR2: table of R2
    # obtaining labels for the relations among latent variables
    path.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))    
    TPC <- matrix(NA, length(path.labs), 1+tn)# TPC: table of path coefficients   

    ## model for root node
    TOW[,1] <- pls$out.weights
    TL[,1] <- round(pls$loadings, 3)
    betas <- pls$path.coefs
    lvs.paths <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1)
                 lvs.paths <- c(lvs.paths, round(betas[i,j],3))
    TPC[,1] <- lvs.paths
    TR2[,1] <- round(pls$r.sqr,3)

    # calculating plspm results for terminal nodes
    for (k in 1:tn)
    {        
        elems.tn <- list.nodes[[term.nodes[k]]]# elements of k-th terminal node
        DTN <- DT[elems.tn,]
        pls.tn <- .pls.basic(DTN, IDM, blocks, modes, scheme, scaled, tol, iter)
        TOW[,k+1] <- pls.tn$out.weights# adding outer weights
        TL[,k+1] <- round(pls.tn$loadings, 3)# adding loadings
        # obtaining relations among latent variables
        lvs.efs <- NULL
        betas <- pls.tn$path.coefs
        for (j in 1:lvs)
            for (i in j:lvs)
                 if (IDM[i,j]==1)
                     lvs.efs <- c(lvs.efs, round(betas[i,j],3))
        TPC[,k+1] <- lvs.efs
        TR2[,k+1] <- round(pls.tn$R2,3)
    }

    dimnames(TOW) <- list(mvs.names, c("Root_Node",tn.labs))
    dimnames(TL) <- list(mvs.names, c("Root_Node",tn.labs))
    dimnames(TPC) <- list(path.labs, c("Root_Node",tn.labs))
    dimnames(TR2) <- list(lvs.names, c("Root_Node",tn.labs))
    resul <- list(weights=TOW, loadings=TL, paths=TPC, r2=TR2, IDM=IDM)
    class(resul) <- "treemox.pls"
    resul
}

