.getSingleTransition=function(clickstream, i, states) {
    template=matrix(rep(0, length(states)^2), nrow=length(states))
    if (length(clickstream) > 1) {
        dimnames(template)=list(states, states)
        app=rep("[[-]]", i)
        clickstream=c(clickstream, app)
        clicks=unlist(clickstream) 
        clicks2=c(clicks[-(1:i)], rep(NA, i))
        dat=data.table(clicks, clicks2) 
        transition=as.data.frame(dcast.data.table(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2"))
        row.names(transition)=transition[,1]        
        pos1=which(!(names(transition) %in% states))
        pos2=which(!(row.names(transition) %in% states)) 
        if (length(pos1)>0 && length(pos2)>0) {
            transition=transition[-pos2,-pos1, drop=F]
            sums=rowSums(transition)
            sums=ifelse(sums==0, 1, sums)
            transition=t(transition/sums)
            template[dimnames(transition)[[1]], dimnames(transition)[[2]]]=transition
        } 
    }
    singleTransition = as.vector(template)
    nms = expand.grid(states, states)
    names(singleTransition) = paste(nms$Var1, nms$Var2, sep=",")
    return(singleTransition)
}

.getSingleFrequencies=function(clickstream, states) {
    transition=rep(0, length(states))
    names(transition)=states
    freq=table(clickstream)
    transition[names(freq)]=freq
    transition=transition/sum(transition)
    return(transition)
}





#' Performs K-Means Clustering on a List of Clickstreams
#' 
#' Performs k-means clustering on a list of clickstreams. For each clickstream a
#' transition matrix of a given order is computed. These transition matrices
#' are used as input for performing k-means clustering.
#' 
#' 
#' @param clickstreamList A list of clickstreams for which the cluster analysis
#' is performed.
#' @param order The order of the transition matrices used as input for
#' clustering (default is 0; 0 and 1 are possible).
#' @param centers The number of clusters.
#' @param ...  Additional parameters for k-means clustering (see
#' \code{\link{kmeans}}).
#' @return This method returns a \code{ClickstreamClusters} object (S3-class).
#' It is a list with the following components: \item{clusters}{ The resulting list of
#' \code{Clickstreams} objects.}
#' \item{centers}{ A matrix of cluster centres.  } \item{states}{ Vector of states} 
#' \item{totss}{ The total sum of squares.  } \item{withinss}{ Vector of within-cluster 
#' sum of squares, one component per cluster.  } \item{tot.withinss}{ Total within-cluster sum of
#' squares, i.e., \code{sum(withinss)}.  } \item{betweenss}{ The
#' between-cluster sum of squares, i.e., \code{totss - tot.withinss}.  }
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{print.ClickstreamClusters}},
#' \code{\link{summary.ClickstreamClusters}}
#' @examples
#' 
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf <- tempfile()
#' writeLines(clickstreams, csf)
#' cls <- readClickstreams(csf, header = TRUE)
#' clusters <- clusterClickstreams(cls, order = 0, centers = 2)
#' print(clusters)
#' 
#' @export clusterClickstreams
clusterClickstreams=function(clickstreamList, order=0, centers, ...) {
    states=unique(as.character(unlist(clickstreamList, use.names = FALSE)))
    if (order==0) {
        transitionData=laply(.data=clickstreamList, .fun=.getSingleFrequencies, states)
    } else {        
        transitionData=ldply(.data=clickstreamList, .fun=.getSingleTransition, order, states)
        transitionData=transitionData[,-1]
    }
    fit=kmeans(transitionData, centers, ...)
    clusterList=list()
    for (i in 1:centers) {
        clusterList[[i]]=clickstreamList[which(fit$cluster==i)]
        class(clusterList[[i]])="Clickstreams"        
    }
    clusters=list(clusters=clusterList, centers=fit$centers, states=states, totss=fit$totss,
                  withinss=fit$withinss, tot.withinss=fit$tot.withinss, betweenss=fit$betweenss,
                  order=order)
    class(clusters)="ClickstreamClusters"
    return(clusters)
} 

.rotate=function(x, n) {
    l=length(x)
    n=n %% l
    if (n == 0) {
        return(x)
    }
    tmp=x[(l-n+1):l]
    x[(n+1):l]=x[1:(l-n)]
    x[1:n]=tmp
    return(x)
}

.getQ=function(i, clickstreamList) { 
    clicks=clickstreamList
    for (j in 1:i) {
        clicks=rbind(clicks, "[[-]]")
    }
    clicks=unlist(clicks, use.names=F)
    clicks2=.rotate(clicks, -i)
    dat=data.table(clicks, clicks2) 
    transition=as.data.frame(dcast.data.table(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2"))
    transition=transition[,-1]  
    pos=which(names(transition)=="[[-]]")
    rnames=names(transition)[-pos]
    transition=transition[,-pos]
    transition=transition[-pos,]
    sums=colSums(t(transition))
    sums[sums==0]=1
    ll=sum(transition*log(transition/sums), na.rm=T)
    transition=as.data.frame(t(transition/sums))
    names(transition)=rnames
    rownames(transition)=rnames
    return(list(ll=ll, transition=transition))
}

.foo=function(params) {
    QX=get("QX")
    X=get("X")    
    error=0
    for (i in 1:length(QX)) {
        error=error+(params[i]*QX[[i]]-X)
    }
    return(sum(error^2))
}

.constr=function(params) {
    return(sum(params))
}

.sollp=function(X, QX, lps = FALSE) {
    if (lps) {
        sumDir = "=="
    } else {
        sumDir = "<="
    }
    f.obj = c(rep(1, length(X)), rep(0, length(QX)))
    f.con = t(sapply(X=seq(1, 2*length(X), 1), FUN=function(x) {
        sign=(x %% 2)*2 -1; 
        i=ceiling(x/2);
        pre=rep(0, length(X))
        pre[i]=1
        c(pre, sign * sapply(X=QX, FUN=function(x) x[i]))
    }))
    f.con = rbind(f.con, c(rep(0, length(X)), rep(1, length(QX))))
    f.con = rbind(f.con, diag(length(X) + length(QX)))
    f.dir = c(rep(">=", 2*length(X)), sumDir, rep(">=", length(X) + length(QX)))
    f.rhs = sapply(X=seq(1, 2*length(X), 1), FUN=function(x) {
        sign=(x %% 2)*2 -1; 
        i=ceiling(x/2);
        sign * X[i]
    })
    f.rhs = c(f.rhs, 1, rep(0, length(X) + length(QX)))
    result = solveLP(cvec = as.numeric(f.obj), 
            bvec = as.numeric(f.rhs),
            Amat = as.matrix(f.con),
            maximum = FALSE,
            const.dir = as.character(f.dir),
            lpSolve = lps)
    return(as.numeric(result$solution))
}



#' Fits a List of Clickstreams to a Markov Chain
#' 
#' This function fits a list of clickstreams to a Markov chain. Zero-order,
#' first-order as well as higher-order Markov chains are supported. For
#' estimating higher-order Markov chains this function solves the following
#' linear or quadratic programming problem:\cr \deqn{\min ||\sum_{i=1}^k X-\lambda_i
#' Q_iX||}{min ||\sum X-\lambda_i Q_iX||} \deqn{\mathrm{s.t.}}{s.t.}
#' \deqn{\sum_{i=1}^k \lambda_i = 1}{sum \lambda_i = 1} \deqn{\lambda_i \ge
#' 0}{\lambda_i \ge 0} The distribution of states is given as \eqn{X}.
#' \eqn{\lambda_i} is the lag parameter for lag \eqn{i} and \eqn{Q_i} the
#' transition matrix.
#' 
#' For solving the quadratic programming problem of higher-order Markov chains,
#' an augmented Lagrange multiplier method from the package
#' \code{\link{Rsolnp}} is used.
#' 
#' @param clickstreamList A list of clickstreams for which a Markov chain is
#' fitted.
#' @param order (Optional) The order of the Markov chain that is fitted from
#' the clickstreams. Per default, Markov chains with \code{order=1} are fitted.
#' It is also possible to fit zero-order Markov chains (\code{order=0}) and
#' higher-order Markov chains.
#' @param verbose (Optional) An optimal logical variable to indicate whether warnings
#' and infos should be printed.
#' @param control (Optional) The control list of optimization parameters. Parameter
#' \code{optimizer} specifies the type of solver used to solve the given
#' optimization problem. Possible values are "linear" (default) and "quadratic".
#' Parameter \code{use.lpSolve} determines whether lpSolve or linprog is used as
#' linear solver.
#' @return Returns a \code{MarkovChain} object.
#' @note At least half of the clickstreams need to consist of as many clicks as
#' the order of the Markov chain that should be fitted.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link[=MarkovChain-class]{MarkovChain}},
#' \code{\link[=Rsolnp]{Rsolnp}}
#' @references This method implements the parameter estimation method presented
#' in Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and
#' Applications}, 2nd edition, Springer, 2013.
#' @examples
#' 
#' # fitting a simple Markov chain
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf <- tempfile()
#' writeLines(clickstreams, csf)
#' cls <- readClickstreams(csf, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' show(mc)
#' 
#' @export fitMarkovChain
fitMarkovChain=function(clickstreamList, order=1, verbose=TRUE, control=list()) { 
    ans=list()
    params=unlist(control)
    if (is.null(params)) {
        ans$optimizer="linear"
        ans$use.lpSolve=FALSE
    } else {
        npar = tolower(names(unlist(control)))
        names(params) = npar
        if(any(substr(npar, 1, 9) == "optimizer")) ans$optimizer = params["optimizer"] else ans$optimizer = "linear"
        if(any(substr(npar, 1, 11) == "use.lpsolve")) ans$use.lpSolve = params["use.lpSolve"] else ans$use.lpSolve = FALSE
    }
    clickVector=unlist(clickstreamList, use.names = FALSE)
    states=unique(as.character(clickVector))
    #### DECIDE ####
    if (as.numeric(R.Version()$major) > 3 || 
        (as.numeric(R.Version()$major) == 3 && as.numeric(R.Version()$minor) >= 2)) {
        lens=as.numeric(lengths(clickstreamList))
    } else {
        lens=as.numeric(laply(.data=clickstreamList, .fun=function(x) c(length(x))))
    }
    n=sum(lens)
    ratio=length(lens[lens>order])/length(lens)
    if (ratio<0.5) {
        stop("The order is to high for the specified click streams.")
    } else if (ratio<1 && verbose) {        
        warning(paste("Some click streams are shorter than ", order, ".", sep=""))
    } 
    start=table(clickVector[c(1, cumsum(lens[-length(lens)])+1)])
    start=start/sum(start)
    end=table(clickVector[cumsum(lens)])
    end=end/sum(end)
    
    q1=.getQ(1, clickstreamList)
    transitions=q1$transition
    diag(transitions)=0
    cs=colSums(transitions)
    rs=rowSums(transitions)
    absorbingPositions=which(cs==0)
    nonAbsorbingPositions=which(cs>0)
    absorbingStates=names(absorbingPositions)
    Q=transitions[nonAbsorbingPositions, nonAbsorbingPositions]
    R=transitions[absorbingPositions, nonAbsorbingPositions]
    It=diag(length(nonAbsorbingPositions))
    N=solve(It-Q)
    B=N %*% t(R)
    if (length(absorbingStates) > 0) {
        absorbingProbabilities=data.frame(state=names(nonAbsorbingPositions), B/rowSums(B))
        row.names(absorbingProbabilities)=nonAbsorbingPositions
    } else {
        absorbingProbabilities=data.frame()
    }
    infTransitions=transitions * t(transitions)
    csInf=colSums(infTransitions)
    transientStates=names(which(csInf>0))
    x=as.data.frame(table(clickVector))
    names(x)=c("states", "frequency") 
    probability=x$frequency/sum(x$frequency)
    x=cbind(x, probability)
    if (order==0) {
        transitions=list(x)
        lambda=0
        logLikelihood=sum(x$frequency*log(x$probability))
    } else {
        X=as.numeric(x$probability)
        ll=vector()
        Q=alply(.data=seq(1,order,1), .margins=1, .fun=.getQ, clickstreamList)
        transitions=llply(.data=Q, .fun=function(q) q$transition)
        ll=laply(.data=Q, .fun=function(q) q$ll)
        QX=llply(.data=transitions, .fun=function(tr) as.matrix(tr)%*%X) 
        environment(.foo)=environment()
        if (ans$optimizer == "quadratic") {
            params=rep(1/order, order)
            model=solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
            lambda=model$pars
        } else {
            solution=.sollp(X, QX, ans$use.lpSolve)
            lambda=solution[seq(length(solution) - length(QX) + 1, length(solution), 1)]
        }
        logLikelihood=as.numeric(lambda %*% ll)
    }
    markovChain=new("MarkovChain", states=states, order=order, start=start, end=end, transitions=transitions,
                    lambda=lambda, logLikelihood=logLikelihood, observations=n,
                    transientStates=transientStates, absorbingStates=absorbingStates,
                    absorbingProbabilities=absorbingProbabilities)
    return(markovChain)
}
