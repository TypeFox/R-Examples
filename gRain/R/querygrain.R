querygrain <- function(object, nodes=nodeNames(object), type="marginal",
                       evidence=NULL, exclude=TRUE, normalize=TRUE,
                       result="array", details=0)
{
  UseMethod("querygrain")
}

qgrain <- querygrain

querygrain.grain <- function(object, nodes = nodeNames(object), type = "marginal",
                             evidence=NULL, exclude=TRUE, normalize=TRUE,
                             result="array", details=0){

    if (!is.null(evidence)){
        if (details>=1) cat(" Inserting (additional) evidence\n")
        object <- setEvidence(object, evidence=evidence)
    }

    type <- match.arg(type, c("marginal","joint","conditional"))
    result <- match.arg(result, c("array","data.frame"))
    t0 <- proc.time()

    if (is.null(nodes))
        return(invisible(NULL))

    if (!object$isCompiled){
        if (details>=1) cat("  Compiling (and propagating) model ...\n")
        object <- compile(object, propagate=TRUE)
    } else {
        if (!object$isPropagated){
            if (details>=1) cat("  Propagating model...\n")
            object <- propagate(object)
        }
    }

    type = match.arg(type, choices=c("marginal","joint","conditional"))
    switch(type,
           "marginal"={
               ans <- .nodeMarginal(object, nodes=nodes, exclude=exclude,
                                    details=details)
               if (result=="data.frame")
                   ans <- lapply(ans, as.data.frame.table)
           },
           "joint"={
               ans <- .nodeJoint(object, nodes=nodes, exclude=exclude,
                                 normalize=normalize, details=details)
               if (result=="data.frame")
                   ans <- as.data.frame.table(ans)
         },
           "conditional"={
               #' qobject <- querygrain.grain(object, nodes=nodes,
               #'                             type="joint", exclude=exclude,
               #'                             result="data.frame")
               #' nst     <- nodeStates(object)[nodes]
               #' ans     <- parray(nodes, nst, values=qobject$Freq,
               #'                   normalize="first")

               qobject <- querygrain.grain(object, nodes=nodes,
                                           type="joint", exclude=exclude)
               #' qq<<- qobject

               ans <- tabDiv( qobject, tabMarg(qobject, nodes[-1]) )


               if (result=="data.frame")
                   ans <- as.data.frame.table(ans)
           })
    if (object$control$timing)
        cat("Time: query", proc.time()-t0, "\n")
    ans
}



.nodeJoint <- function(object, nodes=NULL, exclude=TRUE,
                       normalize=TRUE, details=1){

    if (is.null(nodes)){
        nodes  <- object$rip$nodes
    } else {
        nodes <- intersect(object$rip$nodes, nodes)
    }

    if (exclude)
        nodes <- setdiff(nodes, getEvidence(object)$nodes)


    cliq  <- object$rip$cliques
    ## FIXME: This is potentially slow:
    ## FIXME: gRbase (1.7-1) introduces is_subsetof_ to be used instead.
    ## FIXME: Has been fixed but not checked!!!
    ##idxb <- sapply(cliq, function(cq) subsetof(nodes, cq))
    idxb <- sapply(cliq, function(cq) is_subsetof_(nodes, cq))

    if (any(idxb)){
        ## cat(".Calculating directly from clique\n")
        tab   <- object$equipot[[ which(idxb)[1] ]]
        value <- tableMargin(tab, nodes)
        if (!normalize){
            value$values <- value$values * pEvidence(object)
        }
    } else {
        ## cat(".Calculating brute force\n")
        nnodes <- length( nodes )
        dn    <- object$universe$levels[nodes]
        value <- parray(names(dn), dn)

        nodes2  <- nodes[2:nnodes]
        dn2   <- object$universe$levels[nodes2]
        gr2   <- as.matrix( expand.grid( dn2 ) )

        object  <- absorbEvidence( object )
        zz <- lapply(1:nrow(gr2), function(i){
            tmp <- setFinding(object, nodes=nodes2, states=gr2[i,])
            r <- .nodeMarginal( tmp, nodes[1] )
            v <- r[[1]] * pEvidence(tmp)
            v
        })
        zz <- unlist(zz, use.names=FALSE)
        zz[ is.na( zz )] <- 0
        if (normalize)
            zz <- zz / sum( zz )
        value[] <- zz
    }
    return(value)
}


.nodeMarginal <- function(object, nodes=NULL, exclude=TRUE, details=1){

    if (is.null(nodes)){
        nodes  <- object$rip$nodes
    } else {
        nodes <- intersect(object$rip$nodes, nodes)
    }

    if (exclude)
        nodes <- setdiff(nodes, getEvidence(object)$nodes)

    .rip      <- object$rip
    idx <- match(nodes, .rip$nodes)
    host.cq <- .rip$host[idx]

    if (length(nodes)>0){

        out <- vector("list", length(nodes))
        names(out) <- nodes

        for (i in 1:length(nodes)){
            cvert  <- nodes[i]
            idx    <- host.cq[i]
            ## querygrain - .nodeMarginal: Calculations based on equipot
            cpot   <- object$equipot[[ idx ]]
            mtab   <- tableMargin( cpot, cvert )
            mtab   <- mtab / sum( mtab )
            out[[ i ]] <- mtab
        }
        return( out )
    }
}
