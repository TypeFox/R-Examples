#Null models for community data
#Internal use only; mostly based on picante's code
#' @importFrom picante randomizeMatrix
#' @export
.eco.null <- function(comm, method=c("taxa.labels", "richness", "frequency", "independentswap", "trialswap"), swap.iter=1000){
  #Checks
  method <- match.arg(method)
  
  if(method == "taxa.labels")
    curr.rnd <- comm[,sample(seq(ncol(comm)))] else{
        curr.rnd <- randomizeMatrix(comm, method)
    }
  return(curr.rnd)
}

#Prearing regressions for output
#' @importFrom stats coef
.prepare.regression.output <- function(observed, randomisations=NULL, method=c("lm", "quantile", "mantel"), permute=0, class=NULL){
  #Checks
  method <- match.arg(method)
  
  if(method == "lm"){
    if(permute > 0)
      rnd.slopes <- sapply(randomisations, function(x) coef(x)[2]) else rnd.slopes <- NULL
    obs.slope <- coef(observed)[2]
  }
  
  if(method == "quantile"){
    #Beware quantile regressions with only one tau!
    if(is.null(dim(coef(observed))))
      obs.slope <- coef(observed)[2] else obs.slope <- coef(observed)[2,]
    if(permute > 0){
      if(is.null(dim(coef(observed))))
        rnd.slopes <- sapply(randomisations, function(x) coef(x)[2]) else rnd.slopes <- sapply(randomisations, function(x) coef(x)[2,])
    } else rnd.slopes <- NULL
  }
  
  if(method == "mantel"){
    obs.slope <- observed$statistic
    if(permute > 0)
      rnd.slopes <- sapply(randomisations, function(x) x$statistic) else rnd.slopes <- NULL
  }
  
  output <- list(observed=observed, randomisations=randomisations, obs.slope=obs.slope, rnd.slopes=rnd.slopes,
                 method=method, permute=permute, randomisation=randomisations)
  if(!is.null(class)) output$type <- class
  class(output) <- "eco.xxx.regression"
  return(output)
}

#Printing summaries of regressions
#ADD RETURN STATEMENTS FOR FALL-THROUGH LIKE OTHERS
#' @method summary eco.xxx.regression
#' @param object \code{eco.xxx.regression} object
#' @export
#' @rdname eco.xxx.regression
#' @importFrom quantreg summary.rq print.rq
#' @importFrom stats sd
summary.eco.xxx.regression <- function(object, ...){
    x <- object  ## FIXME:  lazy hack
    cat("\n", x$type, "\n", sep="")
    cat("Method: ", x$method, "\n")
    if(x$permute > 0)
        cat("Randomisation: ", x$method, "; Permutations: ", x$permute, "\n") else cat("Randomisation: NONE\n")
    if(x$method == "quantile" && length(x$obs.slope)>1)
        cat("Observed slopes (at specified taus): ", paste(round(x$obs.slope,2), collapse=","), "\n") else cat("Observed slope: ", round(x$obs.slope,2), "\n")
    if(x$permute > 0 ){
        if(x$method == "quantile" && length(x$obs.slope)>1){
            cat("Random slope means (at specified taus) +/- SD:\n")
            for(i in seq(nrow(x$rnd.slopes)))
                cat(round(mean(x$rnd.slopes[i,]),2), " +/- ", round(sd(x$rnd.slopes[i,]),4), "\n")
        } else cat("Random slope mean +/-SD: ", round(mean(x$rnd.slopes),2), " +/- ", round(sd(x$rnd.slopes),4), "\n")
    }
    cat("Observed model summary:\n")
    if(x$method == "mantel")
        print(x$observed)
    if(x$method == "quantile")
        print.rq(x$observed)
    if(x$method == "lm")
        print(summary(x$observed))
    cat("\n")
}

#' @method print eco.xxx.regression
#' @param x \code{eco.xxx.regression} object
#' @export
#' @rdname eco.xxx.regression
print.eco.xxx.regression <- function(x, ...){
    summary(x, ...)
}

#' @importFrom graphics plot abline
#' @importFrom stats coef
.plot.regression <- function(x, y, observed, randomisations,
                             method=c("quantile", "lm", "mantel"),
                             permute=0, ...){
    method <- match.arg(method)
    plot(y ~ x, ...)
    if(method == "lm"){
        abline(observed, lwd=3)
                                        # Easiest way to silence
                                        # lapply...
        if(permute>0 && method=="lm")
            silent<-lapply(randomisations, abline, col="red")
    }
    if(method == "quantile"){
                                        # Check to see if we've got
                                        # more than one tau value and
                                        # plot accordingly
        if(is.null(dim(coef(observed)))){
            abline(coef(observed), lwd=3)
            for(j in seq(from=1,length.out=permute))
                abline(coef(randomisations[[j]]), col="red")
        } else {
            for(j in seq(ncol(coef(observed)))){
                abline(coef(observed)[,j], lwd=3)
                for(k in seq_along(randomisations))
                    abline(coef(randomisations[[k]])[,j], col="red")
            }
        }
    }
}

#' @method plot eco.xxx.regression
#' @export
#' @rdname eco.xxx.regression
#' @importFrom stats as.dist
#' @importFrom ape cophenetic.phylo
plot.eco.xxx.regression <- function(x, ...){
    if(x$type == "eco.env.regression"){
        #Beware calling without knowing which environmental distance matrix we should be using
        if(!x$altogether) stop("Cannot call 'plot' on an element of an eco.env.regression.list - uses plot(list, no.trait) instead")
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        env.matrix <- as.numeric(pianka.dist(x$data, TRUE))
        .plot.regression(env.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Pianka's Distance", ylab="Ecological Co-existnce", ...)
    }
    
    if(x$type == "eco.trait.regression"){
        #Beware calling without knowing which environmental distance matrix we should be using
        if(!x$altogether) stop("Cannot call 'plot' on an element of an eco.trait.regression.list - uses plot(list, no.trait) instead")
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        trait.matrix <- as.numeric(traits.dist(x$data, TRUE))
        .plot.regression(trait.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
    }

    if(x$type == "eco.phy.regression"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        phy.matrix <- as.numeric(as.dist(cophenetic.phylo(x$data$phy)))
        .plot.regression(phy.matrix, eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                         xlab="Phylogenetic Distance", ylab="Ecological Co-existnce", ...)
    }
}

#' List of eco.xxx.regressions
#' @method summary eco.xxx.regression.list
#' @param object \code{eco.xxx.regression.list} object
#' @param ... additional arguments to plotting functions
#' @export
#' @rdname eco.xxx.regression.list
#' @name eco.xxx.regression.list
summary.eco.xxx.regression.list <- function(object, ...){
    x <- object ## FIXME:  lazy hack
    if(x$type == "eco.env.regression.list"){
        cat("\neco.env.regression list:\n")
        for(i in seq(ncol(object$data$env))){
            cat("\n\n**", names(object$data$env)[i], "**\n")
            summary(object[[i]], ...)
        }
        return()
    }

    if(x$type == "eco.trait.regression.list"){
        cat("\neco.trait.regression list:\n")
        for(i in seq(ncol(object$data$data))){
            cat("\n\n**", names(object$data$data)[i], "**\n")
            summary(x[[i]], ...)
        }
        return()
    }

}

#' @method print eco.xxx.regression.list
#' @param x \code{eco.xxx.regression.list} object
#' @export
#' @rdname eco.xxx.regression.list
print.eco.xxx.regression.list <- function(x, ...){
    if(x$type == "eco.env.regression.list"){
        cat("\neco.env.regression list:\n")
        cat("Traits:\n")
        cat(colnames(x$data$env), sep=", ")
        cat("\nTo examine each regression of each trait, use something like 'x[[1]]', or 'print(x[[2]])'\n")
        cat("To display all at once, call something like 'summary(regression.list)'\n")
        return()
    }
    if(x$type == "eco.trait.regression.list"){
        cat("\neco.trait.regression list:\n")
        cat("Traits:\n")
        cat(colnames(x$data$data), sep=", ")
        cat("\nTo examine each regression of each trait, use something like 'x[[1]]', or 'print(x[[2]])'\n")
        cat("To display all at once, call something like 'summary(regression.list)'\n")
        return()
    }
}

#' @method plot eco.xxx.regression.list
#' @export
#' @rdname eco.xxx.regression.list
plot.eco.xxx.regression.list <- function(x, ...){
    if(x$type == "eco.env.regression.list"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        env.matrix <- pianka.dist(x$data, FALSE)
        if(is.null(which)){
            for(i in seq(ncol(x$data$env)))
                .plot.regression(as.numeric(as.dist(env.matrix[,,i])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                 xlab="Pianka's Distance", ylab="Ecological Co-existnce", main=names(x$data$env)[i], ...)
        } else .plot.regression(as.numeric(as.dist(env.matrix[,,which])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                xlab="Pianka's Distance", ylab="Ecological Co-existnce", ...)
        return()
    }

    if(x$type == "eco.trait.regression.list"){
        eco.matrix <- as.numeric(comm.dist(x$data$comm))
        trait.matrix <- traits.dist(x$data, FALSE)
        if(is.null(which)){
            for(i in seq(ncol(x$data$data)))
                .plot.regression(as.numeric(as.dist(trait.matrix[,,i])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                 xlab="Trait Distance", ylab="Ecological Co-existnce", main=names(x$data$data)[i], ...)
        } else .plot.regression(as.numeric(as.dist(trait.matrix[,,which])), eco.matrix, x$observed, x$randomisations, x$method, x$permute,
                                xlab="Trait Distance", ylab="Ecological Co-existnce", ...)
        return()
    }
}

##' Trim a phylogeny
##'
##' This is a weak wrapper around \code{ape}'s
##' \code{\link{drop.tip}}. Importantly, if asked to drop no species
##' from a phylogeny, it will just return the phylogeny (not an empty
##' phylogeny, as \code{\link{drop.tip}}) will.
##'
##' @param tree An \code{\link[ape:phylo]{phylo}} object
##' @param spp A vector of species (one, many, or none) to be removed
##' from \code{tree}
##' @return \code{\link[ape:phylo]{phylo}} object
##' @seealso \code{\link[ape:drop.tip]{drop.tip}} \code{\link[ape:extract.clade]{extract.clade}}
##' @export
drop_tip <- function(tree, spp)
  if(length(setdiff(tree$tip.label, spp)) >0) return(drop.tip(tree, spp)) else return(tree)

.removeErrors <- function(object) {
    if(inherits(object, "try-error")) return(NULL)
    object
}

#' Manipulate internals of \code{comparative.comm} object
#' @param expr expression to be evaluated within the scope of
#' \code{data}
#' @param ... ignored
#' @export
#' @rdname cc.manip
within.comparative.comm <- function(data, expr, ...){
	
	parent <- parent.frame()
	e <- evalq(environment(), data, parent)
	eval(substitute(expr), e)
	
	cc <- as.list(e)
	class(cc) <- class(data)
        return(cc)
}

.check.ext.dist <- function(ext.dist, species, n.sites){
    if(!inherits(ext.dist, "dist"))
          stop("'ext.dist' must be a distance matrix")
      if(attr(ext.dist, "Size") != n.sites)
          stop("'ext.dist' must have dimensions matching comparative.comm object's species'")
      if(!identical(attr(ext.dist, "Labels"), species))
          warning("'ext.dist' names do not match species data; continuing regardless")
    return(as.matrix(ext.dist))
}
