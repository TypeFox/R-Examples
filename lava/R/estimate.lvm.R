###{{{ estimate.lvm

##' Estimation of parameters in a Latent Variable Model (lvm)
##'
##' Estimate parameters. MLE, IV or user-defined estimator.
##'
##' A list of parameters controlling the estimation and optimization procedures
##' is parsed via the \code{control} argument. By default Maximum Likelihood is
##' used assuming multivariate normal distributed measurement errors. A list
##' with one or more of the following elements is expected:
##'
##' \describe{
##' \item{start:}{Starting value. The order of the parameters can be shown by
##' calling \code{coef} (with \code{mean=TRUE}) on the \code{lvm}-object or with
##' \code{plot(..., labels=TRUE)}. Note that this requires a check that it is
##' actual the model being estimated, as \code{estimate} might add additional
##' restriction to the model, e.g. through the \code{fix} and \code{exo.fix}
##' arguments. The \code{lvm}-object of a fitted model can be extracted with the
##' \code{Model}-function.}
##'
##' \item{starterfun:}{Starter-function with syntax
##' \code{function(lvm, S, mu)}.  Three builtin functions are available:
##' \code{startvalues}, \code{startvalues0}, \code{startvalues1}, ...}
##'
##' \item{estimator:}{ String defining which estimator to use (Defaults to
##' ``\code{gaussian}'')}
##'
##' \item{meanstructure}{Logical variable indicating
##' whether to fit model with meanstructure.}
##'
##' \item{method:}{ String pointing to
##' alternative optimizer (e.g. \code{optim} to use simulated annealing).}
##'
##' \item{control:}{ Parameters passed to the optimizer (default
##' \code{stats::nlminb}).}
##'
##' \item{tol:}{ Tolerance of optimization constraints on lower limit of
##' variance parameters.  } }
##'
##' @aliases estimate
##' @param x \code{lvm}-object
##' @param data \code{data.frame}
##' @param estimator String defining the estimator (see details below)
##' @param control control/optimization parameters (see details below)
##' @param missing Logical variable indiciating how to treat missing data.
##' Setting to FALSE leads to complete case analysis. In the other case
##' likelihood based inference is obtained by integrating out the missing data
##' under assumption the assumption that data is missing at random (MAR).
##' @param weight Optional weights to used by the chosen estimator.
##' @param weightname Weight names (variable names of the model) in case
##' \code{weight} was given as a vector of column names of \code{data}
##' @param weight2 Optional additional dataset used by the chosen
##' estimator.
##' @param id Vector (or name of column in \code{data}) that identifies
##' correlated groups of observations in the data leading to variance estimates
##' based on a sandwich estimator
##' @param fix Logical variable indicating whether parameter restriction
##' automatically should be imposed (e.g. intercepts of latent variables set to
##' 0 and at least one regression parameter of each measurement model fixed to
##' ensure identifiability.)
##' @param index For internal use only
##' @param graph For internal use only
##' @param silent Logical argument indicating whether information should be
##' printed during estimation
##' @param quick If TRUE the parameter estimates are calculated but all
##' additional information such as standard errors are skipped
##' @param method Optimization method
##' @param param set parametrization (see \code{help(lava.options)})
##' @param cluster Obsolete. Alias for 'id'.
##' @param p Evaluate model in parameter 'p' (no optimization)
##' @param ... Additional arguments to be passed to the low level functions
##' @return A \code{lvmfit}-object.
##' @author Klaus K. Holst
##' @seealso \code{score}, \code{information}, ...
##' @keywords models regression
##' @export
##' @method estimate lvm
##' @examples
##' dd <- read.table(header=TRUE,
##' text="x1 x2 x3
##' 1.0 2.0 1.4
##' 2.1 4.1 4.0
##' 3.1 3.4 7.0
##' 4.2 6.1 3.5
##' 5.3 5.2 2.3
##' 1.1 1.6 2.9")
##' e <- estimate(lvm(c(x1,x2,x3)~u),dd)
##'
##' ## Simulation example
##' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
##' covariance(m) <- v1~v2+v3+v4
##' dd <- sim(m,10000) ## Simulate 10000 observations from model
##' e <- estimate(m, dd) ## Estimate parameters
##' e
##'
##' ## Using just sufficient statistics
##' n <- nrow(dd)
##' e0 <- estimate(m,data=list(S=cov(dd)*(n-1)/n,mu=colMeans(dd),n=n))
##'
##' ## Multiple group analysis
##' m <- lvm()
##' regression(m) <- c(y1,y2,y3)~u
##' regression(m) <- u~x
##' d1 <- sim(m,100,p=c("u,u"=1,"u~x"=1))
##' d2 <- sim(m,100,p=c("u,u"=2,"u~x"=-1))
##'
##' mm <- baptize(m)
##' regression(mm,u~x) <- NA
##' covariance(mm,~u) <- NA
##' intercept(mm,~u) <- NA
##' ee <- estimate(list(mm,mm),list(d1,d2))
##'
##' ## Missing data
##' d0 <- makemissing(d1,cols=1:2)
##' e0 <- estimate(m,d0,missing=TRUE)
##' e0
`estimate.lvm` <-
    function(x, data=parent.frame(),
             estimator="gaussian",
             control=list(),
             missing=FALSE,
             weight, weightname,
             weight2,
             id,
             fix,
             index=TRUE,
             graph=FALSE,
             silent=lava.options()$silent,
             quick=FALSE,
             method,
             param,
             cluster,
             p,
             ...) {

        if (length(exogenous(x)>0)) {
            catx <- categorical2dummy(x,data)
            x <- catx$x; data <- catx$data
        }
        cl <- match.call()
        if (!base::missing(param)) {
            oldparam <- lava.options()$param
            lava.options(param=param)
            on.exit(lava.options(param=oldparam))
        }
        if (!base::missing(method)) {
            control["method"] <- list(method)
        }

        optim <- list(
            iter.max=lava.options()$iter.max,
            trace=ifelse(lava.options()$debug,3,0),
            gamma=lava.options()$gamma,
            gamma2=1,
            ngamma=lava.options()$ngamma,
            lambda=0.05,
            abs.tol=1e-9,
            epsilon=1e-10,
            delta=1e-10,
            rel.tol=1e-10,
            S.tol=1e-5,
            stabil=FALSE,
            start=NULL,
            constrain=lava.options()$constrain,
            method=NULL,
            starterfun="startvalues0",
            information="E",
            meanstructure=TRUE,
            sparse=FALSE,
            tol=lava.options()$tol)

        defopt <- lava.options()[]
        defopt <- defopt[intersect(names(defopt),names(optim))]
        optim[names(defopt)] <- defopt
        if (length(control)>0) {
            optim[names(control)] <- control
        }

        if (is.environment(data)) {
            innames <- intersect(ls(envir=data),vars(x))
            data <- as.data.frame(lapply(innames,function(x) get(x,envir=data)))
            names(data) <- innames
        }

        if (!lava.options()$exogenous) exogenous(x) <- NULL

        redvar <- intersect(intersect(parlabels(x),latent(x)),colnames(data))
        if (length(redvar)>0)
            warning(paste("Latent variable exists in dataset",redvar))
        ## Random-slopes:
        xfix <- setdiff(colnames(data)[(colnames(data)%in%parlabels(x,exo=TRUE))],latent(x))
        if (base::missing(fix)) {
            fix <- ifelse(length(xfix)>0,FALSE,TRUE)
        }
        Debug(list("start=",optim$start))

        if (!missing(cluster)) id <- cluster
        if (!missing & (is.matrix(data) | is.data.frame(data))) {
            includelist <- c(manifest(x),xfix)
            if (!base::missing(weight) && is.character(weight)) includelist <- c(includelist,weight)
            if (!base::missing(weight2) && is.character(weight2)) includelist <- c(includelist,weight2)
            if (!base::missing(id) && is.character(id)) includelist <- c(includelist,id)
            data <- na.omit(data[,intersect(colnames(data),includelist),drop=FALSE])
        }

        ## Weights...
        if (!base::missing(weight)) {
            if (is.character(weight)) {
                weight <- data[,weight,drop=FALSE]
                if (!base::missing(weightname)) {
                    colnames(weight) <- weightname
                } else {
                    yvar <- index(x)$endogenous
                    nw <- seq_len(min(length(yvar),ncol(weight)))
                    colnames(weight)[nw] <- yvar[nw]
                }
            }
            weight <- cbind(weight)
        } else {
            weight <- NULL
        }
        if (!base::missing(weight2)) {
            if (is.character(weight2)) {
                weight2 <- data[,weight2]
            }
        } else {
            weight2 <- NULL
        }
        ## Correlated clusters...
        if (!base::missing(id)) {
            if (is.character(id)) {
                id <- data[,id]
            }
        } else {
            id <- NULL
        }

        Debug("procdata")

        dd <- procdata.lvm(x,data=data)
        S <- dd$S; mu <- dd$mu; n <- dd$n
        ## Debug(list("n=",n))
        ## Debug(list("S=",S))
        ## Debug(list("mu=",mu))

        ##  if (fix)
        {
            var.missing <- setdiff(vars(x),colnames(S))
            if (length(var.missing)>0) {## Convert to latent:
                new.lat <- setdiff(var.missing,latent(x))
                if (length(new.lat)>0)
                    x <- latent(x, new.lat)
            }
        }

        ## Run hooks (additional lava plugins)
        myhooks <- gethook()
        for (f in myhooks) {
            res <- do.call(f, list(x=x,data=data,weight=weight,weight2=weight2,estimator=estimator,optim=optim))
            if (!is.null(res$x)) x <- res$x
            if (!is.null(res$data)) data <- res$data
            if (!is.null(res$weight)) weight <- res$weight
            if (!is.null(res$weight2)) weight2 <- res$weight2
            if (!is.null(res$optim)) optim <- res$optim
            if (!is.null(res$estimator)) estimator <- res$estimator
            rm(res)
        }

        checkestimator <- function(x,...) {
            ffname <- paste0(x,c("_objective","_gradient"),".lvm")
            exists(ffname[1])||exists(ffname[2])
        }
        if (!checkestimator(estimator)) { ## Try down/up-case version
            estimator <- tolower(estimator)
            if (!checkestimator(estimator)) {
                estimator <- toupper(estimator)
            }
        }
        ObjectiveFun  <- paste0(estimator, "_objective", ".lvm")
        GradFun  <- paste0(estimator, "_gradient", ".lvm")
        if (!exists(ObjectiveFun) & !exists(GradFun))
            stop("Unknown estimator.")

        Method <-  paste0(estimator, "_method", ".lvm")
        if (!exists(Method)) {
            Method <- "nlminb1"
        } else {
            Method <- get(Method)
        }
        NoOptim <- "method"%in%names(control) && is.null(control$method)
        if (is.null(optim$method) && !(NoOptim)) {
            optim$method <- if (missing) "nlminb1" else Method
        }

        if (!quick & index) {
            ## Proces data and setup some matrices
            x <- fixsome(x, measurement.fix=fix, S=S, mu=mu, n=n,debug=!silent)
            if (!silent)
                message("Reindexing model...\n")
            if (length(xfix)>0) {
                index(x) <- reindex(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE)
            } else {
                x <- updatelvm(x,sparse=optim$sparse,zeroones=TRUE,deriv=TRUE,mean=TRUE)
            }
        }
        if (is.null(estimator) || estimator==FALSE) {
            return(x)
        }

        if (length(index(x)$endogenous)==0) stop("No observed outcome variables. Check variable names in model and data.")
        if (!optim$meanstructure) {
            mu <- NULL
        }

        nparall <- index(x)$npar + ifelse(optim$meanstructure, index(x)$npar.mean+index(x)$npar.ex,0)
        ## Get starting values
        if (!missing(p)) {
            start <- p
            optim$start <- p
        } else {
            myparnames <- coef(x,mean=TRUE)
            paragree <- FALSE
            paragree.2 <- c()
            if (!is.null(optim$start)) {
                paragree <- myparnames%in%names(optim$start)
                paragree.2 <- names(optim$start)%in%myparnames
            }
            if (sum(paragree)>=length(myparnames))
                optim$start <- optim$start[which(paragree.2)]

            if (! (length(optim$start)==length(myparnames) & sum(paragree)==0))
                if (is.null(optim$start) || sum(paragree)<length(myparnames)) {
                    if (is.null(optim$starterfun) && lava.options()$param!="relative")
                        optim$starterfun <- startvalues0
                    start <- suppressWarnings(do.call(optim$starterfun, list(x=x,S=S,mu=mu,debug=lava.options()$debug,silent=silent,data=data,...)))
                    if (!is.null(x$expar) && length(start)<nparall) {
                        ii <- which(index(x)$e1==1)
                        start <- c(start, structure(unlist(x$expar[ii]),names=names(x$expar)[ii]))
                    }
                    ## Debug(list("start=",start))
                    if (length(paragree.2)>0) {
                        start[which(paragree)] <- optim$start[which(paragree.2)]
                    }
                    optim$start <- start
                }
        }

        coefname <- coef(x,mean=optim$meanstructure,fix=FALSE);
        names(optim$start) <- coefname

        ## Missing data
        if (missing) {
            return(estimate.MAR(x=x,data=data,fix=fix,control=optim,debug=lava.options()$debug,silent=silent,estimator=estimator,weight=weight,weight2=weight2,cluster=id,...))
        }

        ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
        constr <- lapply(constrain(x), function(z)(attributes(z)$args))
        xconstrain <- intersect(unlist(constr), manifest(x))
        xconstrainM <- TRUE
        XconstrStdOpt <- TRUE
        if (length(xconstrain)>0) {
            constrainM <- names(constr)%in%unlist(x$mean)
            for (i in seq_len(length(constr))) {
                if (!constrainM[i]) {
                    if (constr[[i]]%in%xconstrain) {
                        xconstrainM <- FALSE
                        break;
                    }
                }
            }
            if (xconstrainM & ( (is.null(control$method) || optim$method=="nlminb0") & (lava.options()$test & estimator=="gaussian") ) ) {
                XconstrStdOpt <- FALSE
                optim$method <- "nlminb0"
                if (is.null(control$constrain)) control$constrain <- TRUE
            }
        }

        ## Setup optimization constraints
        lowmin <- -Inf
        lower <- rep(lowmin,length(optim$start))
        if (length(optim$constrain)==1 & optim$constrain)
            lower[variances(x)+index(x)$npar.mean] <- optim$tol
        if (any(optim$constrain)) {
            if (length(optim$constrain)!=length(lower))
                constrained <- is.finite(lower)
            else
                constrained <- optim$constrain
            lower[] <- -Inf
            optim$constrain <- TRUE
            constrained <- which(constrained)
            nn <- names(optim$start)
            CS <- optim$start[constrained]
            CS[CS<0] <- 0.01
            optim$start[constrained] <- log(CS)
            names(optim$start) <- nn
        }
        ## Fix problems with starting values?
        optim$start[is.nan(optim$start)] <- 0
        ## Debug(list("lower=",lower))

        ObjectiveFun  <- paste0(estimator, "_objective", ".lvm")
        GradFun  <- paste0(estimator, "_gradient", ".lvm")
        if (!exists(ObjectiveFun) & !exists(GradFun)) stop("Unknown estimator.")

        InformationFun <- paste0(estimator, "_hessian", ".lvm")

        mymodel <- x
        myclass <- "lvmfit"

        ## Random slopes?
        if (length(xfix)>0 | (length(xconstrain)>0 & XconstrStdOpt | !lava.options()$test)) { ## Yes
            x0 <- x

            if (length(xfix)>0) {
                myclass <- c("lvmfit.randomslope",myclass)
                nrow <- length(vars(x))
                xpos <- lapply(xfix,function(y) which(regfix(x)$labels==y))
                colpos <- lapply(xpos, function(y) ceiling(y/nrow))
                rowpos <- lapply(xpos, function(y) (y-1)%%nrow+1)
                myfix <- list(var=xfix, col=colpos, row=rowpos)
                x0 <- x
                for (i in seq_along(myfix$var))
                    for (j in seq_len(length(myfix$col[[i]])))
                        regfix(x0, from=vars(x0)[myfix$row[[i]][j]],
                               to=vars(x0)[myfix$col[[i]][j]]) <-
                                   colMeans(data[,myfix$var[[i]],drop=FALSE])
                x0 <- updatelvm(x0,zeroones=TRUE,deriv=TRUE)
                x <- x0
                yvars <- endogenous(x0)
                ## Alter start-values/constraints:
                new.par.idx <- which(coef(mymodel,mean=TRUE,fix=FALSE)%in%coef(x0,mean=TRUE,fix=FALSE))
                if (length(optim$start)>length(new.par.idx))
                    optim$start <- optim$start[new.par.idx]
                lower <- lower[new.par.idx]
                if (optim$constrain) {
                    constrained <- match(constrained,new.par.idx)
                }
            }
            mydata <- as.matrix(data[,manifest(x0)])

            myObj <- function(pp) {
                if (optim$constrain) {
                    pp[constrained] <- exp(pp[constrained])
                }
                myfun <- function(ii) {
                    if (length(xfix)>0)
                        for (i in seq_along(myfix$var)) {
                            x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
                        }
                    if (is.list(weight2)) {
                        res <- do.call(ObjectiveFun, list(x=x0, p=pp, data=mydata[ii,], n=1, weight=weight[ii,], weight2=weight2[ii,]))
                    } else {
                        res <- do.call(ObjectiveFun, list(x=x0, p=pp, data=mydata[ii,], n=1, weight=weight[ii,], weight2=weight2))
                    }
                    return(res)
                }
                sum(sapply(seq_len(nrow(data)),myfun))
            }

            myGrad <- function(pp) {
                if (optim$constrain) {
                    pp[constrained] <- exp(pp[constrained])
                }
                myfun <- function(ii) {
                    if (length(xfix)>0)
                        for (i in seq_along(myfix$var)) {
                            x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
                        }
                    if (is.list(weight2)) {
                        rr <- do.call(GradFun, list(x=x0, p=pp, data=mydata[ii,,drop=FALSE], n=1, weight=weight[ii,], weight2=weight2))
                    } else
                        {
                            rr <- do.call(GradFun, list(x=x0, p=pp, data=mydata[ii,,drop=FALSE], n=1, weight=weight[ii,], weight2=weight2[ii,]))
                        }
                    return(rr)
                }
                ss <- rowSums(rbind(sapply(seq_len(nrow(data)),myfun)))
                if (optim$constrain) {
                    ss[constrained] <- ss[constrained]*pp[constrained]
                }
                return(ss)
            }


            myInfo <- function(pp) {
                myfun <- function(ii) {
                    if (length(xfix)>0)
                        for (i in seq_along(myfix$var)) {
                            x0$fix[cbind(rowpos[[i]],colpos[[i]])] <- index(x0)$A[cbind(rowpos[[i]],colpos[[i]])] <- data[ii,xfix[i]]
                        }
                    if (is.list(weight2)) {
                        res <- do.call(InformationFun, list(p=pp, obj=myObj, x=x0, data=data[ii,],
                                                            n=1, weight=weight[ii,], weight2=weight2))
                    } else {
                        res <- do.call(InformationFun, list(p=pp, obj=myObj, x=x0, data=data[ii,],
                                                            n=1, weight=weight[ii,], weight2=weight2[ii,]))
                    }
                    return(res)
                }
                L <- lapply(seq_len(nrow(data)),function(x) myfun(x))
                val <- apply(array(unlist(L),dim=c(length(pp),length(pp),nrow(data))),c(1,2),sum)
                if (!is.null(attributes(L[[1]])$grad)) {
                    attributes(val)$grad <- colSums (
                        matrix( unlist(lapply(L,function(i) attributes(i)$grad)) , ncol=length(pp), byrow=TRUE)
                        )
                }
                return(val)
            }

##################################################
        } else { ## No, standard model

            ## Non-linear parameter constraints involving observed variables? (e.g. nonlinear regression)
            xconstrain <- c()
            for (i in seq_len(length(constrain(x)))) {
                z <- constrain(x)[[i]]
                xx <- intersect(attributes(z)$args,manifest(x))
                if (length(xx)>0) {
                    warg <- setdiff(attributes(z)$args,xx)
                    wargidx <- which(attributes(z)$args%in%warg)
                    exoidx <- which(attributes(z)$args%in%xx)
                    parname <- names(constrain(x))[i]
                    y <- names(which(unlist(lapply(intercept(x),function(x) x==parname))))
                    el <- list(i,y,parname,xx,exoidx,warg,wargidx,z)
                    names(el) <- c("idx","endo","parname","exo","exoidx","warg","wargidx","func")
                    xconstrain <- c(xconstrain,list(el))
                }
            }
            yconstrain <- unlist(lapply(xconstrain,function(x) x$endo))
            iconstrain <- unlist(lapply(xconstrain,function(x) x$idx))

            MkOffset <- function(pp,grad=FALSE) {
                if (length(xconstrain)>0) {
                    Mu <- matrix(0,nrow(data),length(vars(x))); colnames(Mu) <- vars(x)
                    M <- modelVar(x,p=pp,data=data)
                    M$parval <- c(M$parval,  x$mean[unlist(lapply(x$mean,is.numeric))])
                    for (i in seq_len(length(xconstrain))) {
                        pp <- unlist(M$parval[xconstrain[[i]]$warg]);
                        myidx <- with(xconstrain[[i]],order(c(wargidx,exoidx)))
                        mu <- with(xconstrain[[i]],
                                   apply(data[,exo,drop=FALSE],1,
                                         function(x) func(
                                             unlist(c(pp,x))[myidx])))
                        Mu[,xconstrain[[i]]$endo] <- mu
                    }
                    offsets <- Mu%*%t(M$IAi)[,endogenous(x)]
                    return(offsets)
                }
                return(NULL)
            }

            myObj <- function(pp) {
                if (optim$constrain) {
                    pp[constrained] <- exp(pp[constrained])
                }
                offset <- MkOffset(pp)
                mu0 <- mu; S0 <- S; x0 <- x
                if (!is.null(offset)) {
                    x0$constrain[iconstrain] <- NULL
                    data0 <- data[,manifest(x0)]
                    data0[,endogenous(x)] <- data0[,endogenous(x)]-offset
                    pd <- procdata.lvm(x0,data=data0)
                    S0 <- pd$S; mu0 <- pd$mu
                    x0$mean[yconstrain] <- 0
                }
                do.call(ObjectiveFun, list(x=x0, p=pp, data=data, S=S0, mu=mu0, n=n, weight=weight
                                          ,weight2=weight2, offset=offset
                                           ))
            }

            myGrad <- function(pp) {
                if (optim$constrain)
                    pp[constrained] <- exp(pp[constrained])
                ##  offset <- MkOffset(pp)
                ##  mu0 <- mu; S0 <- S; x0 <- x
                ## if (!is.null(offset)) {
                ##   x0$constrain[iconstrain] <- NULL
                ##   data0 <- data[,manifest(x0)]
                ##   data0[,endogenous(x)] <- data0[,endogenous(x)]-offset
                ##   pd <- procdata.lvm(x0,data=data0)
                ##   S0 <- pd$S; mu0 <- pd$mu
                ## }
                S <- do.call(GradFun, list(x=x, p=pp, data=data, S=S, mu=mu, n=n, weight=weight
                                         , weight2=weight2##, offset=offset
                                           ))
                if (optim$constrain) {
                    S[constrained] <- S[constrained]*pp[constrained]
                }
                if (is.null(mu) & index(x)$npar.mean>0) {
                    return(S[-c(seq_len(index(x)$npar.mean))])
                }
                if (length(S)<length(pp))  S <- c(S,rep(0,length(pp)-length(S)))

                return(S)
            }
            myInfo <- function(pp) {
                I <- do.call(InformationFun, list(p=pp,
                                                  obj=myObj,
                                                  x=x, data=data,
                                                  S=S, mu=mu,
                                                  n=n,
                                                  weight=weight, weight2=weight2,
                                                  type=optim$information
                                                  ))
                if (is.null(mu) && index(x)$npar.mean>0) {
                    return(I[-seq_len(index(x)$npar.mean),-seq_len(index(x)$npar.mean)])
                }
                return(I)
            }

        }

        myHess <- function(pp) {
            p0 <- pp
            if (optim$constrain)
                pp[constrained] <- exp(pp[constrained])
            I0 <- myInfo(pp)
            attributes(I0)$grad <- NULL
            D <- attributes(I0)$grad
            if (is.null(D)) {
                D <- myGrad(p0)
                attributes(I0)$grad <- D
            }
            if (optim$constrain) {
                I0[constrained,-constrained] <- apply(I0[constrained,-constrained,drop=FALSE],2,function(x) x*pp[constrained]);
                I0[-constrained,constrained] <- t(I0[constrained,-constrained])
                if (sum(constrained)==1) {
                    I0[constrained,constrained] <- I0[constrained,constrained]*outer(pp[constrained],pp[constrained]) - D[constrained]
                } else {
                    I0[constrained,constrained] <- I0[constrained,constrained]*outer(pp[constrained],pp[constrained]) - diag(D[constrained],ncol=length(constrained))
                }
            }
            return(I0)
        }
        if (is.null(tryCatch(get(InformationFun),error = function (x) NULL)))
            myInfo <- myHess <- NULL
        if (is.null(tryCatch(get(GradFun),error = function (x) NULL)))
            myGrad <- NULL

        if (!silent) message("Optimizing objective function...")
        if (optim$trace>0 & !silent) message("\n")
        ## Optimize with lower constraints on the variance-parameters
        if ((is.data.frame(data) | is.matrix(data)) && nrow(data)==0) stop("No observations")

        if (!missing(p)) {
            opt <- list(estimate=p)
            ## if (!is.null(myGrad))
            ##     opt <- c(opt,list(gradient=myGrad(p)))
            ## if (!is.null(myObj))
            ##     opt <- c(opt,list(objective=myObj(p)))

        } else {
            if (!is.null(optim$method)) {
                opt <- do.call(optim$method,
                               list(start=optim$start, objective=myObj, gradient=myGrad, hessian=myHess, lower=lower, control=optim, debug=debug))
                if (is.null(opt$estimate))
                    opt$estimate <- opt$par
                if (optim$constrain) {
                opt$estimate[constrained] <- exp(opt$estimate[constrained])
                }

                if (XconstrStdOpt & !is.null(myGrad))
                    opt$gradient <- as.vector(myGrad(opt$par))
                else {
                    opt$gradient <- numDeriv::grad(myObj,opt$par)
                }
            } else {
                if (!NoOptim) {
                    opt <- do.call(ObjectiveFun, list(x=x,data=data,control=control,...))
                    opt$gradient <- rep(0,length(opt$estimate))
                } else {
                    opt <- list(estimate=optim$start,
                                gradient=rep(0,length(optim$start)))
                }
            }
        }
        if (!is.null(opt$convergence)) {
            if (opt$convergence!=0) warning("Lack of convergence. Increase number of iteration or change starting values.")
        } else if (!is.null(opt$gradient) && mean(opt$gradient)^2>1e-3) warning("Lack of convergence. Increase number of iteration or change starting values.")
        if (quick) {
            return(opt$estimate)
        }
        ## Calculate std.err:
        pp <- rep(NA,length(coefname)); names(pp) <- coefname
        if (!is.null(names(opt$estimate))) {
            pp[names(opt$estimate)] <- opt$estimate
            pp.idx <- na.omit(match(coefname,names(opt$estimate)))
        } else {
            pp[] <- opt$estimate
            pp.idx <- seq(length(pp))
        }

        suppressWarnings(mom <- tryCatch(modelVar(x, pp, data=data),error=function(x)NULL))
        if (NoOptim) {
            asVar <- matrix(NA,ncol=length(pp),nrow=length(pp))
        } else {

            if (!silent) message("\nCalculating asymptotic variance...\n")
            asVarFun  <- paste0(estimator, "_variance", ".lvm")

            if (!exists(asVarFun)) {
                if (is.null(myInfo)) {
                if (!is.null(myGrad))
                    myInfo <- function(pp)
                        numDeriv::jacobian(myGrad,pp,method=lava.options()$Dmethod)
                else
                    myInfo <- function(pp)
                        numDeriv::hessian(myObj,pp)
            }
                I <- myInfo(opt$estimate)
                asVar <- tryCatch(Inverse(I),
                                  error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
            } else {
                asVar <- tryCatch(do.call(asVarFun,
                                          list(x=x,p=opt$estimate,data=data,opt=opt)),
                                  error=function(e) matrix(NA, length(opt$estimate), length(opt$estimate)))
            }

            if (any(is.na(asVar))) {
                warning("Problems with asymptotic variance matrix. Possibly non-singular information matrix!")
            }
            if (!is.null(attributes(asVar)$pseudo) && attributes(asVar)$pseudo) {
                warning("Near-singular covariance matrix, using pseudo-inverse!")
            }
            diag(asVar)[diag(asVar)==0] <- NA
        }

        mycoef <- matrix(NA,nrow=nparall,ncol=4)
        mycoef[pp.idx,1] <- opt$estimate ## Will be finished during post.hooks

### OBS: v = t(A)%*%v + e
        res <- list(model=x, call=cl, coef=mycoef,
                    vcov=asVar, mu=mu, S=S, ##A=A, P=P,
                    model0=mymodel, ## Random slope hack
                    estimator=estimator, opt=opt,expar=x$expar,
                    data=list(model.frame=data, S=S, mu=mu,
                        C=mom$C, v=mom$v, n=n,
                        m=length(latent(x)), k=length(index(x)$manifest), weight2=weight2),
                    weight=weight, weight2=weight2,
                    cluster=id,
                    pp.idx=pp.idx,
                    graph=NULL, control=optim)

        class(res) <- myclass

        myhooks <- gethook("post.hooks")
        for (f in myhooks) {
            res0 <- do.call(f,list(x=res))
            if (!is.null(res0))
                res <- res0
        }

        if(graph) {
            res <- edgelabels(res,type="est")
        }

        return(res)
    }

###}}} estimate.lvm

###{{{ estimate.formula

##' @export
estimate.formula <- function(x,data=parent.frame(),pred.norm=c(),unstruct=FALSE,silent=TRUE,id=NULL,distribution=NULL,estimator="gaussian",...) {
    cl <- match.call()
    formulaId <- union(Specials(x,c("cluster")),Specials(x,c("id")))    
    formulaSt <- paste0("~.-cluster(",formulaId,")-id(",formulaId,")")
    if (!is.null(formulaId)) {
        id <- formulaId
        x <- update(x,as.formula(formulaSt))
    }
    if (!is.null(id))
        x <- update(x,as.formula(paste(".~.+",id)))
    varnames <- all.vars(x)
    mf <- model.frame(x,data)
    mt <- attr(mf, "terms")
    yvar <- names(mf)[1]
    y <- mf[,yvar]
    opt <- options(na.action="na.pass")
    mm <- model.matrix(x,data)
    options(opt)
    covars <- colnames(mm)
    covars <- unlist(lapply(covars, function(x) gsub("[^a-zA-Z0-9._]","",x)))
    colnames(mm) <- covars

    if (attr(terms(x),"intercept")==1) {
        covars <- covars[-1]
        it <- c()
    } else {
        it <- "0"
    }
    
    if (!is.null(id)) covars <- setdiff(covars,id)
    model <- lvm(toformula(yvar,c(it,covars)),silent=TRUE)
    if (!is.null(distribution)) {
        lava::distribution(model,yvar) <- distribution
        estimator <- "glm"
    }
    mydata <- na.omit(as.data.frame(cbind(y,mm))); names(mydata)[1] <- yvar
    exogenous(model) <- setdiff(covars,pred.norm)
    if (unstruct) {
        model <- covariance(model,pred.norm,pairwise=TRUE)
    }
    estimate(model,mydata,silent=silent,id=id,estimator=estimator,...)
}

###}}} estimate.formula
