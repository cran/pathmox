treemox.boot <-
function(pls, treemox, X=NULL, br=100)
{
    # ============================= ARGUMENTS ============================
    # pls: an object of class plspm
    # treemox: an object of class treemox 
    # X: optional matrix or dataframe used when pls$data is NULL
    # br: number of bootstraps resamples, (50 <= br <= 500)

    # ==================== Checking function arguments ===================
    if (class(pls)!="plspm")
        stop("Argument 'pls' must be an object of class 'plspm'")
    if (class(treemox)!="treemox")
        stop("Argument 'treemox' must be an object of class 'treemox'")
    if (nrow(pls$latents)!=treemox$MOX$Size[1]) 
        stop("Arguments 'pls' and 'pathmox' are incompatible. Different number of observations")
    if (!is.null(X)) # if X available
    {
        if (is.null(pls$data))
        {
            if (!is.matrix(X) && !is.data.frame(X))
                stop("Invalid object 'X'. Must be a numeric matrix or data frame.")
            if (nrow(X)!=nrow(pls$latents))
                stop("Argument 'pls' and 'X' are incompatible. Different number of rows.")
            if (nrow(X)!=pathmox$MOX$Size[1])
                stop("Arguments 'X' and 'pathmox' are incompatible. Different number of observations")
        }
    } else { # if no X
        if (is.null(pls$data))
            stop("Argument 'X' is missing. No dataset available.")
    }
    if (mode(br)!="numeric" || length(br)!=1 || (br%%1)!=0 ||
        br<50 || br>500) {
        warning("argument 'br' must be an integer between 50 and 500. Default 'br=100' is used.")   
        br <- 100
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
    lvs <- ncol(IDM)
    obs <- nrow(DT)
    endo <- rowSums(IDM)
    endo[endo!=0] <- 1
    mvs <- sum(blocks) 
    obs.names = rownames(DT)
    mvs.names = colnames(DT)
    lvs.names = colnames(IDM)

    ### Parameters from pathmox
    MOX <- treemox$MOX
    list.nodes <- treemox$list.nodes
    term.nodes = which(MOX$Terminal=="yes") - 1# identifying terminal nodes    
    tn = length(term.nodes)# number of terminal nodes
    tn.labs = paste("Node", MOX$Node[term.nodes+1], sep="_")# labels for terminal 

    path.labs <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1) 
                 path.labs <- c(path.labs, paste(lvs.names[j],"->",lvs.names[i],sep=""))   
    TPC <- matrix(NA, length(path.labs), 1+tn)# TPC: table of path coefficients   
 
    ### Path coefficients for Root Node
    betas <- pls$path.coefs
    lvs.paths <- NULL
    for (j in 1:lvs)
        for (i in j:lvs)
             if (IDM[i,j]==1)
                 lvs.paths <- c(lvs.paths, round(betas[i,j],3))
    TPC[,1] <- lvs.paths

    # calculating plspm results for terminal nodes
    for (k in 1:tn)
    {        
        elems.tn <- list.nodes[[term.nodes[k]]]# elements of k-th terminal node
        DTN <- DT[elems.tn,]
        pls.tn <- .pls.basic(DTN, IDM, blocks, modes, scheme, scaled, tol, iter)
        # obtaining relations among latent variables
        lvs.efs <- NULL
        betas <- pls.tn$path.coefs
        for (j in 1:lvs)
            for (i in j:lvs)
                 if (IDM[i,j]==1)
                     lvs.efs <- c(lvs.efs, round(betas[i,j],3))
        TPC[,k+1] <- lvs.efs
    }
    dimnames(TPC) <- list(path.labs, c("Root_Node",tn.labs))

    # ======================================================= #
    # ===================== BOOTSTRAPING ==================== #
    # default number of samples br=100
    bootnum <- br
    Path.meanboot = matrix(NA, sum(IDM), tn+1)
    Path.steboot = matrix(NA, sum(IDM), tn+1)
    Path.perc05 = matrix(NA, sum(IDM), tn+1)
    Path.perc95 = matrix(NA, sum(IDM), tn+1)

    # =============== Bootstrapping for Root Node ===========
    PATHS <- matrix(NA, bootnum, sum(IDM))
    i <- 1
    while (i <= bootnum)
    {
        boot.obs <- sample.int(nrow(DT), size=nrow(DT), replace=TRUE)
        DM.boot <- DT[boot.obs,]
        # scaling boot sample
        if (scaled) {
            sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
            X.boot <- scale(DM.boot, scale=sd.XB)
        } else {
            X.boot <- scale(DM.boot, scale=FALSE)
        }
        colnames(X.boot) <- mvs.names
        # calculating boot model parameters 
        w.boot <- .pls.weights(X.boot, IDM, blocks, modes, scheme, tol, iter)
        if (is.null(w.boot)) {
            i <- i - 1
            next
        }
        Y.boot <- X.boot %*% w.boot[[2]]
        pathmod <- .pls.paths(IDM, Y.boot, plsr=FALSE)
        P.boot <- pathmod[[2]]
        PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
        i <- i + 1
    }
    Path.meanboot[,1] <- apply(PATHS, 2, mean)
    Path.steboot[,1] <- apply(PATHS, 2, sd)
    Path.perc05[,1] <- apply(PATHS, 2, function(x) quantile(x,p=.05))
    Path.perc95[,1] <- apply(PATHS, 2, function(x) quantile(x,p=.95))

    # ============= Bootstrapping for terminal nodes ===============
    for (k in 1:tn)
    {        
        elems.tn = list.nodes[[term.nodes[k]]]# elements of k-th terminal node
        DT.new = DT[elems.tn,]
        PATHS = matrix(NA, bootnum, sum(IDM))
        i <- 1
        while (i <= bootnum)
        {            
            boot.obs <- sample.int(nrow(DT.new), size=nrow(DT.new), replace=TRUE)
            DM.boot <- DT.new[boot.obs,]
            # scaling boot sample
            if (scaled) {
                sd.XB <- sqrt((nrow(DM.boot)-1)/nrow(DM.boot)) * apply(DM.boot, 2, sd)
                X.boot <- scale(DM.boot, scale=sd.XB)
            } else {
                X.boot <- scale(DM.boot, scale=FALSE)
            }
            colnames(X.boot) <- mvs.names
            # calculating boot model parameters 
            w.boot <- .pls.weights(X.boot, IDM, blocks, modes, scheme, tol, iter)
            if (is.null(w.boot)) {
                i <- i - 1
                next
            }
            Y.boot <- X.boot %*% w.boot[[2]]
            pathmod <- .pls.paths(IDM, Y.boot, plsr=FALSE)
            P.boot <- pathmod[[2]]
            PATHS[i,] <- as.vector(P.boot[which(IDM==1)])
            i <- i + 1
        }
        Path.meanboot[,(k+1)] = apply(PATHS, 2, mean)
        Path.steboot[,(k+1)] = apply(PATHS, 2, sd)
        Path.perc05[,(k+1)] = apply(PATHS, 2, function(x) quantile(x,p=.05))
        Path.perc95[,(k+1)] = apply(PATHS, 2, function(x) quantile(x,p=.95))
    } # end for 'tn' terminal nodes

    dimnames(Path.meanboot) <- list(path.labs, c("Root_Node",tn.labs))
    dimnames(Path.steboot) <- list(path.labs, c("Root_Node",tn.labs))
    dimnames(Path.perc05) <- list(path.labs, c("Root_Node",tn.labs))
    dimnames(Path.perc95) <- list(path.labs, c("Root_Node",tn.labs))

    bootnodes <- list(PC=TPC, PMB=Path.meanboot, PSB=Path.steboot, 
                     PP05=Path.perc05, PP95=Path.perc95)
    class(bootnodes) <- "bootnodes"
    bootnodes
}

