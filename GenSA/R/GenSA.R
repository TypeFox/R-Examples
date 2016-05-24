# File: GenSA.R
# 
# Author: Sylvain Gubian
# Aim: Function for General Simulated Annealing

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################


GenSA <- function(par=NULL, fn, lower, upper, control=list(), ...)
{
# Do some checks
    jc = NULL
        if (!is.function(fn) || is.null(fn)) {
            stop("'fn' has to be a R function.")
        }

# Create an environment for sharing data between R and C
    genSA.env <- new.env(hash=TRUE, parent=emptyenv())

        fn1 <- function(par) {
            ret <- fn(par, ...)
                return(ret)
        }


    if (!is.null(jc)) {
        if (!is.function(jc)) {
            stop("'jc' has to be a R function.")
        }
        jc1 <-  function(par, ...) {
            return(jc(par,...))
        }
    }
    else
    {
        jc1 <- NULL
    }

    LSE <- function(theta, ui, ci, mu,xlow,xhigh,count) {
        assign("xlow", xlow, envir=genSA.env)
            assign("xhigh", xhigh, envir=genSA.env)
            assign("count", count, envir=genSA.env)

            res <- constrOptim(theta=theta, f=fn2, ui=ui, ci=ci,mu=mu, grad=NULL,outer.eps = 1e-06,)
            counts <- get("count", genSA.env)
            ret <- list(value=res$value, convergence=res$convergence, par=res$par,counts=as.integer(counts))
            return(ret)
    }

    fn2 <- function(par, ...) {
        if (!is.null(jc1)){
            in.constraint <- jc1(par, ...) ;
            if (!in.constraint) {
                return(1e10)
            }
        }
        else {
            penalty <- 0
                xlow <- get("xlow", genSA.env)
                xhigh <- get("xhigh", genSA.env)
                count <- get("count", genSA.env)
                counts <- count
                for(i in (length(lower))) {
                    if (par[i] >= xlow[i] && par[i] <= xhigh[i]) {
                        delta.energy <- 0
                    } else {
                        if (par[i] < xlow[i]) {
                            delta.energy <-  abs(par[i] - xlow[i]) * 1e11
                        }
                        if (par[i] > xhigh[i]) {
                            delta.energy <- abs(par[i] - xhigh[i]) * 1e11
                        }
                    }
                    penalty <- penalty + delta.energy			
                }
            if (penalty > 1.e-10) {
                to.return <- penalty+1.e10
                    return(to.return)
            } else {
                to.return <- fn1(par, ...)
                    counts <- counts + 1
                    to.return <- to.return + penalty
                    if (is.nan(to.return)) {	
                        to.return = count * 1e5 + 1e10;
                        assign("count", counts, envir=genSA.env)
                            return(to.return)
                    }  else {
                        assign("count", counts, envir=genSA.env)
                            return(to.return)
                    }
            }
        }
    }


    assign("LSE", LSE, envir=genSA.env)


        con <- list(
                maxit = 5000,
                threshold.stop = NULL,
                temperature = 5230,
                visiting.param = 2.62,
                acceptance.param = -5.0,
                max.time = NULL,
                nb.stop.improvement = 1e6,
                smooth = TRUE,
                max.call = 10000000,
                simple.function = FALSE,
                trace.fn = NULL,
                verbose = FALSE,
                trace.mat=TRUE
                )
        con$high.dim = TRUE
        con$markov.length = 2 * length(lower)
        con$tem.restart = .1

# Perform some checks before callinc C code
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if(length(noNms <- namc[!namc %in% nmsC]))
            warning("unknown names in control: ", paste(noNms,collapse=", "))


                if (!exists("par") && (length(lower)==0 || length(upper)==0)) {
                    stop("There is no par or no lower/upper bounds defined")
                }

    if (length(lower) != length(upper)) {
        stop("Lower and upper bounds vector do not have the same length")
    }

    cmp <- unique(lower < upper) 
        if (length(cmp)!=1 || !cmp) {
            stop("Lower and upper bounds are not consistent (lower >= upper)")
        }

    if (!is.null(par) ) {

        if (length(lower)==0 || length(lower) != length(par)) {
            stop("Lower bounds vector size does not match with par size, using -Inf")
        }
        if (length(upper)==0 || length(upper) != length(par)) {
            stop("Upper bounds vector size does not match with par size, using -Inf")
        }
        if (any(is.na(par)) || any(is.nan(par)) || any(is.infinite(par))) {
            stop("par contains NA, NAN or Inf")
        }
    }
    else {
        if (con$verbose) {
            cat("Initializing par with random data inside bounds\n")
        }
        par <- vector()
#initialize par with random values in the bounds
            par <- lower + runif(length(lower))*(upper-lower)
    }

    ret <- list()
# Create instance of the GenSACaller
        instance <- .Call("createInstance", PACKAGE="GenSA")
        if (is.null(instance)) {
            stop("Can not create GenSACaller instance!")
        }

# Call execute on the instance
    res <- .Call("execute", par, lower, upper, fn1, jc1, con, genSA.env, instance, PACKAGE="GenSA")
        if (is.null(res)) {
            stop("Can not call execute function on instance")
        }

# Get the results in a list
    res <- .Call("getREnergy", instance, PACKAGE="GenSA")
        if (is.null(res)) {
            message("Can not get minimum function value")
        }
        else
        {
            ret$value <- res
        }

    res <- .Call("getRXMiniVector", instance, PACKAGE="GenSA")
        if (is.null(res)) {
            message("Can not get calculated par values")
        }
        else
        {
            ret$par <- res
        }

    nr <- .Call("getRTraceMatSize", instance, PACKAGE="GenSA")
        if (nr >0) {
            ret$trace.mat <- matrix(NA, nr, 4)
                ret$trace.mat[,1] <- as.integer(.Call("getRTraceMat", instance, "nSteps", PACKAGE="GenSA"))
                ret$trace.mat[,2] <- as.numeric(.Call("getRTraceMat", instance, "temperature", PACKAGE="GenSA"))
                ret$trace.mat[,3] <- as.numeric(.Call("getRTraceMat", instance, "currentEnergy", PACKAGE="GenSA"))
                ret$trace.mat[,4] <- as.numeric(.Call("getRTraceMat", instance, "minEnergy", PACKAGE="GenSA"))
                colnames(ret$trace.mat) <- c("nb.steps", "temperature", "function.value", "current.minimum")
        }
    if (!is.null(con$trace.fn)) {
        if (con$verbose) {
            cat(paste('Writing trace.mat data into file:', con$trace.fn,'\n'))
        }
        write.table(ret$trace.mat, file=con$trace.fn, col.names=TRUE, row.names=FALSE)
            ret$trace.mat <- paste('trace.mat is written in file:', con$trace.fn)
    }


    res <- .Call("getRNbFuncCall", instance, PACKAGE="GenSA")
        if (is.null(res)) {
            message("Can not get number of function calls")
        }
        else
        {
            ret$counts <- res
        }
    .Call("releaseInstance", instance, PACKAGE="GenSA")
        ret
}
