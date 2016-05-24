## #################################################################
##
## Everything is about setEvidence etc here. The old setFinding type
## functions (in the Finding.R file) are just remapped.
##
## Distinguish between (marginal) evidence on single nodes
## (written evidence) and (joint) evidence on several nodels
## (j.evidence).
##
#################################################################

pEvidence <- function(object)
  attr(object$equipot,"pEvidence")

getEvidence <- function(object){
    object$evidence
}

## #################################################################
##
## setEvidence: With this syntax (nodes/states) setEvidence has been
## used in the bnlearn book
##
## #################################################################


setEvidence <- function(object, nodes=NULL, states=NULL, evidence=NULL, nslist=NULL,
                        propagate=TRUE, details=0){
    if (!is.null(nslist))
        stop("Argument 'nslist' has been deprecated; please use 'evidence' instead\n")

    .setEvidence_wrap( object, nodes=nodes,
                 states=states, evidence=evidence, propagate=propagate, details=details)
}


.setEvidence_wrap <- function(object, nodes=NULL, states=NULL, evidence=NULL, propagate=TRUE, details=0){
    if (!is.null( evidence )){
        setEvidence_(object, evidence, propagate=propagate, details=details)
    } else {
        if ( !is.null( nodes ) ){
            if (!is.null( states ) && length( nodes )==length( states )){
                evidence <- as.vector(states, "list")
                names(evidence) <- nodes
                setEvidence_(object, evidence, propagate=propagate, details=details)
            } else {
                stop( "Either states are not given or nodes and states do not have same length" )
            }
        }  else {
            stop( "Evidence is not given; nothing to do...")
        }
    }
}


## setEvidence_: The workhorse
setEvidence_ <- function(object, evidence=NULL, propagate=TRUE, details=0){

    #' details=1
    #' cat("++++ setEvidence_\n"); print(evidence)

    old.ev <- getEvidence(object)

    if (.is.JointEvidence( old.ev ))
        stop("JointEvidence (multivariate evidence) has been given;\n   to update evidence use setJointEvidence\n")

    if ( !.is.Evidence(evidence) ){

        ## Relevant if someone sticks in a row of a dataframe; not sure if this is a good solution.
        if (class(evidence)=="data.frame"){
            if (nrow(evidence)>1)
                stop("evidence is a data.frame with more than one row; only one row is allowed\n")
            evidence <- lapply(evidence, as.character)
        }

        ## Strip any NA's in the evidence; hmmm - maybe I do this twice ?? FIXME ??
        idx <- !unlist(lapply(evidence, function(e) any(is.na(e))), use.names = FALSE)
        if (length(idx)>0)
            evidence <- evidence[ idx ]

        new.ev <- .evidence2parray(evidence, object$universe$levels)
    } else {
        new.ev <- evidence
    }

    #' cat("Evidence - later\n"); print(new.ev)

    tot.ev <- new.ev # can be changed below

    if (details>0){
        cat("new.ev (nodes):"); print(new.ev$summary$nodes)
        cat("old.ev (nodes):"); print(old.ev$summary$nodes)
    }

    if (!object$isCompiled){
        object <- compile(object)
    }
    object$isInitialized  <- FALSE
    object$isPropagated   <- FALSE

    if ( length(old.ev) > 0 ){
        ## nodes on which there is already evidence will not be given
        ## new evidence
        idx <- match( intersect(new.ev$summary$nodes, old.ev$summary$nodes), new.ev$summary$nodes)
        if ( length( idx ) > 0 ){
            new.ev <- .delete.evidence( new.ev, idx)
            if (details>0){
                cat("new.ev (nodes - updated):"); print(new.ev$summary$nodes)
            }
        }
        tot.ev <- .append.evidence( old.ev, new.ev )
    }

    if ( length( new.ev ) > 0 ){
        if (details>0){
            cat("tot.ev (nodes):"); print(tot.ev$summary$nodes)
            cat("new.ev (to be inserted):"); print(new.ev$summmary$nodes)
        }
        host  <- .get.host.clique( new.ev$evidence, object$rip )
        object$temppot <- .insert.evidence.in.potential( object$temppot, new.ev, host )
        object$evidence <- tot.ev
    }

    if(details>0){
        cat("after insertion:\n"); print( object$evidence )
    }

    if (propagate){
        propagate(object)
    } else {
        object
    }
}


.insert.evidence.in.potential <- function( pot, evi.list, hostclique ){
    #' if (any(is.na(hostclique))) stop("NAs in hostclique...")
    for (i in seq_along( evi.list$evidence ) ){
        j <- hostclique[ i ]
        p <- evi.list$evidence[[ i ]]
        pot[[j]] <- tabMult__( pot[[ j ]], p )
    }
    pot
}

retractEvidence <- function(object, nodes=NULL, propagate=TRUE){
    .retractEvidence_internal(object, nodes=nodes, propagate=propagate)
}

.retractEvidence_internal <- function(object, nodes=NULL, propagate=TRUE){
    #cat("++++ retractEvidence_\n")
    .resetgrain <- function(x){
        x$temppot       <- x$origpot
        x$evidence       <- NULL
        ## FIXME Do we isInitialized? I doubt! Status: Now removed!
        ## x$isInitialized <- TRUE
        x$isPropagated  <- FALSE
        x
    }

    if ( is.null( nodes ) ){
        object <- .resetgrain( object )
    } else {
        old.ev <- getEvidence(object)
        if ( .is.JointEvidence(old.ev) )
            stop("JointEvidence (Multivariate evidence) has been set;\n   to retract retractJointEvidence instead...")

        #cat("old.ev:\n"); print(old.ev)
        idx <- match(intersect(old.ev$summary$nodes, nodes), old.ev$summary$nodes)
        #print(idx)
        if (length(idx)>0){
            new.ev <- .delete.evidence( old.ev, idx )
            #cat("new.ev:\n"); print(new.ev)
            object <- .resetgrain(object)
            if ( length(new.ev$summary$nodes) > 0 ){
                object <- setEvidence_(object, evidence=new.ev, propagate=FALSE)
            }
        }
    }

    if (propagate){
        propagate(object)
    } else {
        object
    }
}

## #################################################################
##
## setJointEvidence EXPERIMENTAL - MULTIDIMENSIONAL EVIDENCE
##
## #################################################################

setJointEvidence <- function(object, evidence=NULL, propagate=TRUE, details=0){

    #print(evidence)
    if (!object$isCompiled){
        object <- compile(object)
    }
    object$isInitialized <- FALSE
    object$isPropagated  <- FALSE

    new.ev <- if (!.is.JointEvidence(evidence))
        .j.evidence2parray(evidence, object$universe$levels)
    else
        evidence

    #print(new.ev)
    if( details>0 ){
        cat("new.ev:\n"); print( new.ev )
    }

    n.nodes <- length( new.ev$evidence )
    if (n.nodes>0){
        curr.evidence    <- getEvidence( object )
        ##  Insert new.ev to temppot
        if (length(nodes)>0){
            object$temppot <- .insert.jevi.in.potential( new.ev, #new.ev$nodes, new.ev$evidence,
                                                   object$temppot, object$rip )
            total.evidence <- .append.j.evidence( curr.evidence, new.ev )
            object$evidence <- total.evidence
        }
    }

    if (propagate){
        propagate(object)
    } else {
        object
    }
}

.insert.jevi.in.potential <- function(j.evi.list, pot, rip, details=0){
    host.idx <- .get.host.clique(j.evi.list$evidence, rip)
    #' cat("host.idx: \n"); print(host.idx)
    #' nodes <- lapply(j.evi.list$evidence, .namesDimnames)
    #' cat("nodes  : \n"); print( nodes )

    for (i in seq_along(j.evi.list$evidence)){
        j <- host.idx[ i ]
        #' cat("j.evi.list[i]: "); print( j.evi.list$evidence[[i]] );
        #' cat("pot: "); print( ftable( pot[[h.idx]] ))
        pot[[ j ]]  <- tabMult__( pot[[ j ]], j.evi.list$evidence[[ i ]] )
    }
    pot
}

retractJointEvidence <- function(object, items=NULL, propagate=TRUE){
    ##cat("++++ retractJointEvidence_\n")
    .resetgrain <- function(x){
        x$temppot       <- x$origpot
        x$evidence       <- NULL
        ## FIXME Do we isInitialized? I doubt! Status: Now removed!
        ## x$isInitialized <- TRUE
        x$isPropagated  <- FALSE
        x
    }

    if (is.null(items)){
        object <- .resetgrain( object )
    } else {
        old.ev <- getEvidence(object)
        if ( .is.Evidence(old.ev) )
            stop("Evidence (univariate evidence) has been set;\n   to retract retractEvidence instead...")
        if ( !is.numeric(items) ){
            cat(sprintf("Invalid items: %s\n", toString(items)))
            stop( "items must be numeric" )
        }
        idx <- items
        new.ev <- .delete.j.evidence( old.ev, idx )

        object <- .resetgrain(object)
        if ( length(new.ev$summary$nodes) > 0 ){
            object <- setJointEvidence(object, evidence=new.ev, propagate=FALSE)
        }
    }

    if (propagate){
        propagate(object)
    } else {
        object
    }
}


.evidence2parray <- function(evi.list, levels){

    ## First remove all evidence specified as NA
    not.na <- !unlist(lapply(lapply(evi.list,is.na), any), use.names=FALSE)
    if (length( not.na ) > 0)
        evi.list <- evi.list[ not.na ]

    evidence           <- vector("list", length(evi.list))
    is.hard.evidence   <- rep.int(TRUE,  length(evi.list))
    hard.state         <- rep.int(NA,    length(evi.list))

    #' cat(".evidence2parray\n"); print(evi.list)

    for (i in seq_along(evi.list)){
        ev <- evi.list[i]
        v <- ev[[1]]

        if( is.array(v)){
            n <- names(dimnames(v))
            is.hard.evidence[i]    <- FALSE
            evidence[[i]] <- v
            next
        }

        if (is.character(v)){
            n <- names(evi.list)[i]
            hard.state[i]  <- v
            evidence[[i]]  <- .hard.state2parray(n, v, levels[[n]])
            next
        }

        if (is.numeric(v)){
            n <- names(evi.list)[i]
            is.hard.evidence[i] <- FALSE
            evidence[[i]] <- .soft.state2parray(n, v, levels[[n]])
        }
    }

    ## print(evidence)
    ## If evidence is zero on all states or negative on some (or all) states then it is invalid
    keep <- unlist(lapply(evidence, function(e){ sum(e) !=0 && all(e>=0) }), use.names=FALSE)
    ## print(keep)

    nodes <- unique.default( unlist(lapply(evidence, .namesDimnames)),
                            use.names=FALSE )

    out <- list(summary=list(
                    nodes=nodes[keep],
                    is.hard.evidence=is.hard.evidence[keep],
                    hard.state=hard.state[keep]),
                evidence=evidence[keep])

    class( out ) <- "grainEvidence_"
    ## print.default( out )
    out
}

.j.evidence2parray <- function(evi.list, levels){

    ## First we remove all evidence specified as NA
    not.na <- !unlist(lapply(lapply(evi.list,is.na), any), use.names=FALSE)
    if (length( not.na ) > 0)
        evi.list <- evi.list[not.na]

    evidence         <- vector("list", length(evi.list))
    is.hard.evidence <- rep.int(TRUE, length(evi.list))
    hard.state       <- rep.int(NA_character_, length(evi.list))

    for (i in seq_along(evi.list)){
        n <- names(evi.list)[i]
        v <- evi.list[[i]]

        if (is.array( v )){
            is.hard.evidence[i]  <- FALSE
            evidence[[i]]  <- v
            next
        }

        if (is.character(v)){
            hard.state[i]  <- v
            evidence[[i]]  <- .hard.state2parray(n, v, levels[[n]])
            next
        }

        if (is.numeric(v)){
            is.hard.evidence[i]  <- FALSE
            ## Maximum value of soft evidence is 1, minimum is 0
            v[ v>1 ] <- 1
            v[ v<0 ] <- 0
            evidence[[i]]  <- .soft.state2parray(n, v, levels[[n]])
        }
    }

    nodes <- unique.default( unlist(lapply(evidence, .namesDimnames)), use.names=FALSE )
    #out <- list(nodes=nodes, evidence=evidence)
    out <- list(summary=list(nodes=nodes), evidence=evidence)

    class( out ) <- "grainJEvidence_"
    out
}


absorbEvidence <- function(object, propagate=TRUE ){
    # Update 'object' as
    # 1) set origpot <- temppot
    # 2) ignore any finding
    object$origpot <- object$temppot
    object$evidence <- NULL
    object$isPropagated <- FALSE

    if (propagate){
        propagate(object)
    } else {
        object
    }
}


print.grainEvidence_ <- function(x, ...){
    cat("Univariate evidence\n")
    cat(sprintf("pEvidence=%f\n", x$pEvidence))
    print(as.data.frame(x$summary))
    print(x$evidence)
    #print(list(overview=as.data.frame(x[c(1,3,4)]), evidence=x[[2]]) )
    invisible( x )
}

print.grainJEvidence_ <- function(x, ...){
    cat("General (multinode) evidence\n")
    print(x$evidence)
}


## ####################################################################
##
## Evidence related utility functions
##
## ####################################################################


.is.Evidence <- function(x){
    !is.null(x) && class(x)=="grainEvidence_"
}

.is.JointEvidence <- function(x){
    !is.null(x) && class(x)=="grainJointEvidence_"
}

.hard.state2parray <- function(n, v, lev){
    #str(list(n,v,lev))
    tab <- .fast.parray(n, list(lev), rep.int(0, length(lev)))
    tab[ match( v, lev )] <- 1
    tab
}

.soft.state2parray <- function(n, v, lev){
    #str(list(n,v,lev))
    .fast.parray(n, list(lev), v)
}

.fast.parray <- function(varNames, levels, values=1){
    #str(list(varNames, levels, values))
    dn <- if(is.list(levels)) levels else list(levels)
    names(dn) <- varNames
    array(values, dimnames=dn)
}

.append.evidence <- function(ev1, ev2){
    if (is.null(ev1))
        return (ev2)
    if (is.null(ev2))
        return (ev1)

    #cat(".append.evidence:\n"); cat("ev1:\n"); print(ev1); cat("ev2:\n"); print(ev2)
    summary <- lapply(seq_along(ev1$summary),
                      function(i) c(ev1$summary[[i]], ev2$summary[[i]]))
    names(summary) <- names(ev1$summary)
    evidence <- c(ev1$evidence, ev2$evidence)
    out <- list(summary=summary, evidence=evidence)
    class(out) <- "grainEvidence_"
    out
}

.append.j.evidence <- function(ev1, ev2){
    if (is.null(ev1))
        return (ev2)
    if (is.null(ev2))
        return (ev1)

    #cat(".append.j.evidence\n"); print(ev1); print(ev2)

    evidence <- c(ev1$evidence, ev2$evidence)
    nodes <- unique.default( unlist(lapply(evidence, .namesDimnames)) )
    out   <- list(summary=list(nodes=nodes), evidence=evidence)
    class(out) <- "grainJointEvidence_"
    out
}

.delete.evidence <- function(ev, idx){
    if (length(idx)>0){
        summary <- lapply(ev$summary, function(x) x[-idx])
        names(summary) <- names(ev$summary)
        ev <- list(summary=summary, evidence=ev$evidence[-idx])
        #ev <- lapply(ev, function(x) x[-idx])
        class(ev) <- "grainEvidence_"
    }
    ev
}

.delete.j.evidence <- function(ev, idx){
    if (length(idx)>0){
        summary <- lapply(ev$summary, function(x) x[-idx])
        names(summary) <- names(ev$summary)
        ev <- list(summary=summary, evidence=ev$evidence[-idx])
        #ev <- lapply(ev, function(x) x[-idx])
        class(ev) <- "grainJointEvidence_"
    }
    ev
}




## ####################################################################
##
## Ting der skal omdøbes / et andet sted hen ved lejlighed
##
## ####################################################################

## FIXME: .get.host.clique: ineffektiv implementation;
## FIXME: gRbase (1.7-1) introduces get_superset_ to be used instead
## FIXME: Has been fixed but not checked!!!

#' .get.host.clique <- function(evidence, rip){
#'     unlist(lapply(evidence,
#'                   function(x){
#'                       n <- names(dimnames(x))
#'                       which( isin( rip$cliques, n, index=TRUE) !=0 )[1]} ),
#'            use.names = FALSE)
#' }


.get.host.clique <- function(evidence, rip){
    unlist(lapply(evidence,
                  function(x){
                      n <- names(dimnames(x))
                      gRbase::get_superset_(n, rip$cliques, all=FALSE)
                  }),
                  use.names = FALSE)
}




