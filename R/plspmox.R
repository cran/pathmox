plspmox <-
function(pls, treemox, X=NULL, node, boot.val=FALSE, br=NULL, dataset=FALSE)
{
    # ====================== ARGUMENTS =====================
    # pls: an object of class plspm
    # treemox: an object of class treemox (pathmox or techmox)
    # X: optional data matrix or dataframe
    # node: number of node (integer larger than 1)
    # boot.val:a logical value indicating whether bootstrap validation is done 
    # br: an integer indicating the number of bootstraps resamples, used 
    #     only when boot.val=TRUE, (100 <= br <= 1000)
    # dataset: a logical value indicating whether the data should be retrieved

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
    if (missing(node))
        stop("argument 'node' (number of node) is missing")
    if (mode(node)!="numeric" || length(node)!=1 || node<=1 || (node%%1)!=0)
        stop("Invalid number of 'node'. Must be an integer larger than 1")
    if (length(which(treemox$MOX$Node==node))==0)
        stop("Invalid number of 'node'")

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
    list.nodes <- treemox$list.nodes
    node <- which(treemox$MOX$Node==node) - 1
    node.obs <- list.nodes[[node]]
    DT.node <- DT[node.obs,]
    res <- plspm(DT.node, IDM, outer, modes, scheme, scaled, boot.val=boot.val, 
                 br=br, tol=tol, iter=iter, dataset=FALSE)
    return(res)
}

