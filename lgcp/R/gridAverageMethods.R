###
# Generic functions for the class gridAverage    
###

##' GAinitialise function
##'
##' Generic function defining the the initialisation step for the \code{gridAverage} class of functions. 
##' The function is called invisibly within \code{MALAlgcp} and facilitates the computation of
##' Monte Carlo Averages online.
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GAinitialise
##' @seealso \link{setoutput}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAinitialise <- function(F,...){
    UseMethod("GAinitialise")
}



##' GAupdate function
##'
##' Generic function defining the the update step for the \code{gridAverage} class of functions.
##' The function is called invisibly within \code{MALAlgcp} and facilitates the computation of
##' Monte Carlo Averages online.
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GAupdate
##' @seealso \link{setoutput}, \link{GAinitialise}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAupdate <- function(F,...){
    UseMethod("GAupdate")
}



##' GAfinalise function
##'
##' Generic function defining the the finalisation step for the \code{gridAverage} class of functions.
##' The function is called invisibly within \code{MALAlgcp} and facilitates the computation of
##' Monte Carlo Averages online.
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GAfinalise
##' @seealso \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAreturnvalue}
##' @export

GAfinalise <- function(F,...){
    UseMethod("GAfinalise")
}



##' GAreturnvalue function
##'
##' Generic function defining the the returned value for the \code{gridAverage} class of functions.
##' The function is called invisibly within \code{MALAlgcp} and facilitates the computation of
##' Monte Carlo Averages online.
##'
##' @param F an object    
##' @param ... additional arguments  
##' @return method GAreturnvalue
##' @seealso \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}
##' @export

GAreturnvalue <- function(F,...){
    UseMethod("GAreturnvalue")
}



###
# Functions to facilitate online computation of Monte Carlo averages
###



##' nullAverage function
##'
##' A null scheme, that does not perform any computation in the running of \code{lgcpPredict}, it is the default
##' value of \code{gridmeans} in the argument \code{output.control}.
##'
##' @return object of class nullAverage
##' @seealso \link{setoutput}, \link{lgcpPredict}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

nullAverage <- function(){
    obj <- "NULL"
    class(obj) <- c("nullAverage","gridAverage")
    return(obj)
}


##' MonteCarloAverage function
##'
##' This function creates an object of class \code{MonteCarloAverage}. The purpose of the function is to compute 
##' Monte Carlo expectations online in the function \code{lgcpPredict}, it is set in the argument \code{gridmeans}
##' of the argument \code{output.control}.
##'
##' 
##' A Monte Carlo Average is computed as:
##' \deqn{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_1:t_2})] \approx \frac1n\sum_{i=1}^n g(Y_{t_1:t_2}^{(i)})}{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_1:t_2})] \approx \frac1n\sum_{i=1}^n g(Y_{t_1:t_2}^{(i)})}
##' where \eqn{g}{g} is a function of interest, \eqn{Y_{t_1:t_2}^{(i)}}{Y_{t_1:t_2}^{(i)}} is the \eqn{i}{i}th retained sample from the target  
##' and \eqn{n}{n} is the total number of retained iterations. For example, to compute the mean of \eqn{Y_{t_1:t_2}}{Y_{t_1:t_2}} set,
##' \deqn{g(Y_{t_1:t_2}) = Y_{t_1:t_2},}{g(Y_{t_1:t_2}) = Y_{t_1:t_2},}
##' the output from such a Monte Carlo average would be a set of \eqn{t_2-t_1}{t_2-t_1} grids, each cell of which 
##' being equal to the mean over all retained iterations of the algorithm (NOTE: this is just an example computation, in
##' practice, there is no need to compute the mean on line explicitly, as this is already done by defaul in \code{lgcpPredict}).
##' For further examples, see below. The option \code{last=TRUE} computes,
##' \deqn{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_2})],}{E_{\pi(Y_{t_1:t_2}|X_{t_1:t_2})}[g(Y_{t_2})],}
##' so in this case the expectation over the last time point only is computed. This can save computation time.
##'
##' @param funlist a character vector of names of functions, each accepting single argument Y
##' @param lastonly compute average using only time T? (see ?lgcpPredict for definition of T)
##' @return object of class MonteCarloAverage
##' @seealso \link{setoutput}, \link{lgcpPredict}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}, \link{exceedProbs}
##' @examples
##' fun1 <- function(x){return(x)}   # gives the mean
##' fun2 <- function(x){return(x^2)} # computes E(X^2). Can be used with the 
##'                                  # mean to compute variances, since 
##'                                  # Var(X) = E(X^2) - E(X)^2
##' fun3 <- exceedProbs(c(1.5,2,3))  # exceedance probabilities, 
##'                                  #see ?exceedProbs
##' mca <- MonteCarloAverage(c("fun1","fun2","fun3"))
##' mca2 <- MonteCarloAverage(c("fun1","fun2","fun3"),lastonly=TRUE)
##' @export

MonteCarloAverage <- function(funlist,lastonly=TRUE){
    obj <- list()
    obj$funlist <- as.list(funlist)
    fl <- unlist(funlist)
    for(i in 1:length(fl)){
        if(!inherits(get(fl[i]),"function")){
            stop(paste("Function",fl[i],"not defined."))
        }
    }
    obj$lastonly <- lastonly
    result <- list()
    iter <- 0
    itinc <- function(){
        iter <<- iter + 1
    }
    retit <- function(){
        return(iter)
    }
    ini <- function(Y){ # initialise result object
        if (lastonly){
            for(i in 1:length(funlist)){
                result[[i]] <<- eval(call(funlist[[i]],Y[[length(Y)]]))
            }            
        }
        else{
            for(i in 1:length(funlist)){
                result[[i]] <<- list()
                for (j in 1:length(Y)){
                    result[[i]][[j]] <<- eval(call(funlist[[i]],Y[[j]]))
                }
            }
        }
    } 
    upd <- function(Y){ # update result object
        if (lastonly){
            for(i in 1:length(funlist)){
                result[[i]] <<- result[[i]] + eval(call(funlist[[i]],Y[[length(Y)]]))
            } 
        }
        else{
            for(i in 1:length(funlist)){
                for (j in 1:length(Y)){
                    result[[i]][[j]] <<- result[[i]][[j]] + eval(call(funlist[[i]],Y[[j]]))
                }
            }
        }
    }      
    fin <- function(){ # take mean of result object
        if (lastonly){
            for(i in 1:length(funlist)){
                result[[i]] <<- result[[i]] / iter 
            } 
        }
        else{
            for(i in 1:length(result)){
                for (j in 1:length(result[[i]])){
                    result[[i]][[j]] <<- result[[i]][[j]] / iter
                }
            }
        }
    }    
    retres <- function(){
        return(result)
    }
    obj$iterinc <- itinc
    obj$returniter <- retit
    obj$initialise <- ini
    obj$update <- upd
    obj$finalise <- fin
    obj$returnresult <- retres
    class(obj) <- c("MonteCarloAverage","gridAverage")
    return(obj)
}



##' GAinitialise.nullAverage function
##'
##' This is a null function and performs no action. 
##'
##' @method GAinitialise nullAverage
##' @param F an object of class nullAverage    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAinitialise.nullAverage <- function(F,...){
    return(NULL)
}



##' GAinitialise.MonteCarloAverage function
##'
##' Initialise a Monte Carlo averaging scheme.
##'
##' @method GAinitialise MonteCarloAverage
##' @param F an object of class MonteCarloAverage    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{MonteCarloAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAinitialise.MonteCarloAverage <- function(F,...){
    return(NULL)
}



##' GAupdate.nullAverage function
##'
##' This is a null function and performs no action.
##'
##' @method GAupdate nullAverage
##' @param F an object of class nullAverage    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAupdate.nullAverage <- function(F,...){
    return(NULL)
}



##' GAupdate.MonteCarloAverage function
##'
##' Update a Monte Carlo averaging scheme. This function performs the Monte Carlo sum online.
##'
##' @method GAupdate MonteCarloAverage
##' @param F an object of class MonteCarloAverage   
##' @param ... additional arguments 
##' @return updates Monte Carlo sums
##' @seealso \link{MonteCarloAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAupdate.MonteCarloAverage <- function(F,...){
    F$iterinc()
    M <- get("M",envir=parent.frame()) 
    N <- get("N",envir=parent.frame()) 
    if(get("SpatialOnlyMode",envir=parent.frame())){
        if (F$returniter()==1){
            F$initialise(Y=list(get("oldtags",envir=parent.frame())$Y[1:M,1:N])) # note Y converted to a list to make the other functions work   
        }
        else{
            F$update(Y=list(get("oldtags",envir=parent.frame())$Y[1:M,1:N])) 
        }
    }
    else{ # otherwise in space-time mode 
        if (F$returniter()==1){
            F$initialise(Y=lapply(get("oldtags",envir=parent.frame())$Y,function(x){x[1:M,1:N]}))   
        }
        else{
            F$update(Y=lapply(get("oldtags",envir=parent.frame())$Y,function(x){x[1:M,1:N]})) 
        }
    }    
}



##' GAfinalise.nullAverage function
##'
##' This is a null function and performs no action.
##'
##' @method GAfinalise nullAverage
##' @param F an object of class nullAverage    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAfinalise.nullAverage <- function(F,...){
    return(NULL)
}



##' GAfinalise.MonteCarloAverage function
##'
##' Finalise a Monte Carlo averaging scheme. Divide the sum by the number of iterations.
##'
##' @method GAfinalise MonteCarloAverage
##' @param F an object of class MonteCarloAverage    
##' @param ... additional arguments 
##' @return computes Monte Carlo averages
##' @seealso \link{MonteCarloAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAfinalise.MonteCarloAverage <- function(F,...){
    F$finalise()
}



##' GAreturnvalue.nullAverage function##'
##'
##' This is a null function and performs no action.
##'
##' @method GAreturnvalue nullAverage
##' @param F an object of class nullAverage    
##' @param ... additional arguments 
##' @return nothing
##' @seealso \link{nullAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAreturnvalue.nullAverage <- function(F,...){
    return(NULL)
}



##' GAreturnvalue.MonteCarloAverage function
##'
##' Returns the required Monte Carlo average.
##'
##' @method GAreturnvalue MonteCarloAverage
##' @param F an object of class MonteCarloAverage   
##' @param ... additional arguments 
##' @return results from MonteCarloAverage
##' @seealso \link{MonteCarloAverage}, \link{setoutput}, \link{GAinitialise}, \link{GAupdate}, \link{GAfinalise}, \link{GAreturnvalue}
##' @export

GAreturnvalue.MonteCarloAverage <- function(F,...){
    return(list(funlist=F$funlist,return=F$returnresult()))
}



##' exceedProbs function
##'
##' This function can be called using \code{MonteCarloAverage} (see \code{fun3} the examples in the help file for
##' \link{MonteCarloAverage}). It computes exceedance probabilities,
##' \deqn{P[\exp(Y_{t_1:t_2})>k],}{P[\exp(Y_{t_1:t_2})>k],}
##' that is the probability that the relative reisk exceeds threshold \eqn{k}{k}. Note that it is possible
##' to pass vectors of tresholds to the function, and the exceedance probabilities will be computed for each
##' of these.
##'  
##' @param threshold vector of threshold levels for the indicator function
##' @param direction default 'upper' giving exceedance probabilities, alternative is 'lower', which gives 'subordinate probabilities'
##' @return a function of Y that computes the indicator function I(exp(Y)>threshold) evaluated for each cell of a matrix Y
##' If several tresholds are specified an array is returned with the [,,i]th slice equal to I(exp(Y)>threshold[i])
##' @seealso \link{MonteCarloAverage}, \link{setoutput}
##' @export

exceedProbs <- function(threshold,direction="upper"){
    fun <- function(Y){
        EY <- exp(Y)
        d <- dim(Y)
        len <- length(threshold)
        if(len==1){
            if(direction=="upper"){
                return(matrix(as.numeric(EY>threshold),d[1],d[2]))
            }
            else{
                return(matrix(as.numeric(EY<threshold),d[1],d[2]))
            }            
        }
        else{
            A <- array(dim=c(d[1],d[2],len))
            
            for(i in 1:len){
                if(direction=="upper"){
                    A[,,i] <- matrix(as.numeric(EY>threshold[i]),d[1],d[2])
                }
                else{
                    A[,,i] <- matrix(as.numeric(EY<threshold[i]),d[1],d[2])
                }
            }
            return(A)
        }
    }
    attr(fun,"threshold") <- threshold
    attr(fun,"direction") <- direction
    return(fun)
}

##' exceedProbsAggregated function
##'
##' NOTE THIS FUNCTION IS IN TESTING AT PRESENT
##'
##' This function computes regional exceedance probabilities after MCMC has finished, it requires the information to have been dumped to disk, and
##' to have been computed using the function lgcpPredictAggregated
##' \deqn{P[\exp(Y_{t_1:t_2})>k],}{P[\exp(Y_{t_1:t_2})>k],}
##' that is the probability that the relative risk exceeds threshold \eqn{k}{k}. Note that it is possible
##' to pass vectors of tresholds to the function, and the exceedance probabilities will be computed for each
##' of these.
##'  
##' @param threshold vector of threshold levels for the indicator function
##' @param lg an object of class aggregatedPredict
##' @param lastonly logical, whether to only compute the exceedances for the last time point. default is TRUE  
##' @return a function of Y that computes the indicator function I(exp(Y)>threshold) evaluated for each cell of a matrix Y, but with values aggregated to regions
##' If several tresholds are specified an array is returned with the [,,i]th slice equal to I(exp(Y)>threshold[i])
##' @seealso \link{lgcpPredictAggregated}
##' @export

exceedProbsAggregated <- function(threshold,lg=NULL,lastonly=TRUE){

    if(!is.null(lg)){
        M <- lg$M
        N <- lg$N
        verifyclass(lg,"aggregatedPredict")
        regpop <- lg$app$spdf$population
        nreg <- length(regpop)
        if(is.null(regpop)){
            stop("Reguire regional population denominators to compute exceedances. Add data $population to the SpatialPolygonsDataFrame.")
        }
        olay <- lg$overlay
        cellarea <- (lg$mcens[2]-lg$mcens[1])*(lg$ncens[2]-lg$ncens[1])
        lambda <- lg$grid
        mu <- lg$temporal
        nt <- 1
        if(!lastonly){
            nt <- length(mu)
        }        
        cts <- lg$RegCounts
    }
    else{
        stop("Currently, this is only implemented for post processing.")
    }
    len <- length(threshold)
    tempfunc <- function(i,EY,d,cnt){
        if(cts[i]==0){
            return(rep(NA,len))
        }
        else{
            A <- rep(NA,len)            
            for(j in 1:len){
                A[j] <- as.numeric((sum((cellarea*mu[cnt]*lambda[[j]][1:M,1:N][olay==i]*EY[olay==i]))/regpop[i])>threshold[j])
            }
            return(A)
        }
    }
    fun <- function(Y){
        EY <- exp(Y)
        d <- dim(Y)
        cnt <- 1
        if(lastonly){
           cnt <- nt 
        }
        return(t(sapply(1:nreg,function(i){ans<-tempfunc(i,EY=EY,d=d,cnt=cnt);cnt<<-cnt+1;if(cnt>nt){cnt<<-1};return(ans)})))
    }
    attr(fun,"threshold") <- threshold
    class(fun) <- "exceedProbs"
    return(fun)
}

