# shared reshaping portion of MFClus and MFClusBoot
reshapeCluster <- function(data, formula, compare, envir){
	assign("test", 5, envir = envir)
	
## remove assignment to global; this is bad practice!!
    #assign('cluster', function(x){return(x)}, envir = .GlobalEnv)
	cluster <- function(x){return(x)}
    this.call <- match.call()
    Terms <- terms(formula, specials = 'cluster', data = data)
	environment(Terms) <- environment()
    A <- model.frame(formula = Terms, data = data)
## subset data for only the comparison groups - mcv 09/03/13
	A <- A[A[,attr(Terms, "term.labels")[1]] %in% compare, ]
    dat <- A[, 1]
    group <- A[, 2]
## remove any group levels that aren't present; don't want to evaluate for empty groups - mcv 08/27/13
	group <- factor(group)
    clusters <- A[, 3]
    strat <- unique(as.character(clusters))
	
	assign('dat', dat, envir = envir)
	assign('group', group, envir = envir)
	assign('clusters', clusters, envir = envir)
	assign('strat', strat, envir = envir)
}

# used in the bootstrapping functions MFClusBoot  MFBoot HLBoot
emp.hpd <- function (X, alpha){
    # empirical hpd by shortest length interval
    X <- sort(X)
    probs <- cbind(low=seq(0,alpha,.001),high=seq(1-alpha,1,.001))
    int.len <- quantile(X,prob=probs[,'high'])-quantile(X,prob=probs[,'low'])
    shortest <- min(int.len)
    first <- which(int.len==shortest)[1]
    hpd <- quantile(X,prob=probs[first,],type=7)
    # see documentation for quantile() for type
    return(hpd)
}
