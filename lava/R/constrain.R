##' Define range constraints of parameters
##'
##' @aliases Range.lvm
##' @title Define range constraints of parameters
##' @param a Lower bound
##' @param b Upper bound
##' @return function
##' @author Klaus K. Holst
##' @export
Range.lvm <- function(a=0,b=1) {
  if (b==Inf) {
    f <- function(x) {
      res <- a+exp(x)
      attributes(res)$grad <- exp
      res
    }
    return(f)
  }
  if (a==-Inf) {
    f <- function(x) {
      res <- -exp(x)+b
      attributes(res)$grad <- function(x) -exp(x)
      res
    }
    return(f)
  }
  f <- function(x) {
    res <- (a+b*exp(x))/(1+exp(x))
    attributes(res)$grad <- function(x) exp(x)*(b-a-a*b*exp(x))/(1+exp(x))^2
    res
  }
  return(f)
}

##' Add non-linear constraints to latent variable model
##'
##' Add non-linear constraints to latent variable model
##'
##' Add non-linear parameter constraints as well as non-linear associations
##' between covariates and latent or observed variables in the model (non-linear
##' regression).
##'
##' As an example we will specify the follow multiple regression model:
##'
##' \deqn{E(Y|X_1,X_2) = \alpha + \beta_1 X_1 + \beta_2 X_2} \deqn{V(Y|X_1,X_2)
##' = v}
##'
##' which is defined (with the appropiate parameter labels) as
##'
##' \code{m <- lvm(y ~ f(x,beta1) + f(x,beta2))}
##'
##' \code{intercept(m) <- y ~ f(alpha)}
##'
##' \code{covariance(m) <- y ~ f(v)}
##'
##' The somewhat strained parameter constraint \deqn{ v =
##' \frac{(beta1-beta2)^2}{alpha}}
##'
##' can then specified as
##'
##' \code{constrain(m,v ~ beta1 + beta2 + alpha) <- function(x)
##' (x[1]-x[2])^2/x[3] }
##'
##' A subset of the arguments \code{args} can be covariates in the model,
##' allowing the specification of non-linear regression models.  As an example
##' the non-linear regression model \deqn{ E(Y\mid X) = \nu + \Phi(\alpha +
##' \beta X)} where \eqn{\Phi} denotes the standard normal cumulative
##' distribution function, can be defined as
##'
##' \code{m <- lvm(y ~ f(x,0)) # No linear effect of x}
##'
##' Next we add three new parameters using the \code{parameter} assigment
##' function:
##'
##' \code{parameter(m) <- ~nu+alpha+beta}
##'
##' The intercept of \eqn{Y} is defined as \code{mu}
##'
##' \code{intercept(m) <- y ~ f(mu)}
##'
##' And finally the newly added intercept parameter \code{mu} is defined as the
##' appropiate non-linear function of \eqn{\alpha}, \eqn{\nu} and \eqn{\beta}:
##'
##' \code{constrain(m, mu ~ x + alpha + nu) <- function(x)
##' pnorm(x[1]*x[2])+x[3]}
##'
##' The \code{constraints} function can be used to show the estimated non-linear
##' parameter constraints of an estimated model object (\code{lvmfit} or
##' \code{multigroupfit}). Calling \code{constrain} with no additional arguments
##' beyound \code{x} will return a list of the functions and parameter names
##' defining the non-linear restrictions.
##'
##' The gradient function can optionally be added as an attribute \code{grad} to
##' the return value of the function defined by \code{value}. In this case the
##' analytical derivatives will be calculated via the chain rule when evaluating
##' the corresponding score function of the log-likelihood. If the gradient
##' attribute is omitted the chain rule will be applied on a numeric
##' approximation of the gradient.
##' @aliases constrain constrain<- constrain.default constrain<-.multigroup
##' constrain<-.default constraints parameter<-
##' @return A \code{lvm} object.
##' @author Klaus K. Holst
##' @seealso \code{\link{regression}}, \code{\link{intercept}},
##' \code{\link{covariance}}
##' @keywords models regression
##' @examples
##' ##############################
##' ### Non-linear parameter constraints 1
##' ##############################
##' m <- lvm(y ~ f(x1,gamma)+f(x2,beta))
##' covariance(m) <- y ~ f(v)
##' d <- sim(m,100)
##' m1 <- m; constrain(m1,beta ~ v) <- function(x) x^2
##' ## Define slope of x2 to be the square of the residual variance of y
##' ## Estimate both restricted and unrestricted model
##' e <- estimate(m,d,control=list(method="NR"))
##' e1 <- estimate(m1,d)
##' p1 <- coef(e1)
##' p1 <- c(p1[1:2],p1[3]^2,p1[3])
##' ## Likelihood of unrestricted model evaluated in MLE of restricted model
##' logLik(e,p1)
##' ## Likelihood of restricted model (MLE)
##' logLik(e1)
##'
##' ##############################
##' ### Non-linear regression
##' ##############################
##'
##' ## Simulate data
##' m <- lvm(c(y1,y2)~f(x,0)+f(eta,1))
##' latent(m) <- ~eta
##' covariance(m,~y1+y2) <- "v"
##' intercept(m,~y1+y2) <- "mu"
##' covariance(m,~eta) <- "zeta"
##' intercept(m,~eta) <- 0
##' set.seed(1)
##' d <- sim(m,100,p=c(v=0.01,zeta=0.01))[,manifest(m)]
##' d <- transform(d,
##'                y1=y1+2*pnorm(2*x),
##'                y2=y2+2*pnorm(2*x))
##'
##' ## Specify model and estimate parameters
##' constrain(m, mu ~ x + alpha + nu + gamma) <- function(x) x[4]*pnorm(x[3]+x[1]*x[2])
##' \donttest{ ## Reduce Ex.Timings
##' e <- estimate(m,d,control=list(trace=1,constrain=TRUE))
##' constraints(e,data=d)
##' ## Plot model-fit
##' plot(y1~x,d,pch=16); points(y2~x,d,pch=16,col="gray")
##' x0 <- seq(-4,4,length.out=100)
##' lines(x0,coef(e)["nu"] + coef(e)["gamma"]*pnorm(coef(e)["alpha"]*x0))
##' }
##'
##' ##############################
##' ### Multigroup model
##' ##############################
##' ### Define two models
##' m1 <- lvm(y ~ f(x,beta)+f(z,beta2))
##' m2 <- lvm(y ~ f(x,psi) + z)
##' ### And simulate data from them
##' d1 <- sim(m1,500)
##' d2 <- sim(m2,500)
##' ### Add 'non'-linear parameter constraint
##' constrain(m2,psi ~ beta2) <- function(x) x
##' ## Add parameter beta2 to model 2, now beta2 exists in both models
##' parameter(m2) <- ~ beta2
##' ee <- estimate(list(m1,m2),list(d1,d2),control=list(method="NR"))
##' summary(ee)
##'
##' m3 <- lvm(y ~ f(x,beta)+f(z,beta2))
##' m4 <- lvm(y ~ f(x,beta2) + z)
##' e2 <- estimate(list(m3,m4),list(d1,d2),control=list(method="NR"))
##' e2
##' @export
##' @usage
##'
##' \method{constrain}{default}(x,par,args,...) <- value
##'
##' \method{constrain}{multigroup}(x,par,k=1,...) <- value
##'
##' constraints(object,data=model.frame(object),vcov=object$vcov,level=0.95,
##'                         p=pars.default(object),k,idx,...)
##'
##' @param x \code{lvm}-object
##' @param par Name of new parameter. Alternatively a formula with lhs
##' specifying the new parameter and the rhs defining the names of the
##' parameters or variable names defining the new parameter (overruling the
##' \code{args} argument).
##' @param args Vector of variables names or parameter names that are used in
##' defining \code{par}
##' @param k For multigroup models this argument specifies which group to
##' add/extract the constraint
##' @param value Real function taking args as a vector argument
##' @param object \code{lvm}-object
##' @param data Data-row from which possible non-linear constraints should be
##' calculated
##' @param vcov Variance matrix of parameter estimates
##' @param level Level of confidence limits
##' @param p Parameter vector
##' @param idx Index indicating which constraints to extract
##' @param \dots Additional arguments to be passed to the low level functions
"constrain<-" <- function(x,...,value) UseMethod("constrain<-")
##' @export
"constrain" <- function(x,...) UseMethod("constrain")

##' @export
constrain.default <- function(x,fun, idx, level=0.95, vcov, estimate=FALSE, ...) {
  if (estimate) {
    return(constraints(x,...))
  }
  if (missing(fun)) {
    if (inherits(Model(x),"multigroup")) {
      res <- list()
      for (m in Model(x)$lvm) {
        if (length(constrain(m))>0)
          res <- c(res, constrain(m))
      }
      return(res)
    }
    return(Model(x)$constrain)
  }
  if (is.numeric(x)) {
     b <- x
   } else {
     b <- pars(x)
   }
  if (missing(vcov)) {
    S <- stats::vcov(x)
  } else {
    S <- vcov
  }
  if (!missing(idx)) {
    b <- b[idx]; S <- S[idx,idx,drop=FALSE]
  }
  fb <- fun(b)
  pl <- 1-(1-level)/2
  D <- rbind(numDeriv::grad(fun,b))
  se <- (D%*%S%*%t(D))^0.5
  res <- c(fb,se,fb+c(-1,1)*qnorm(pl)*se)
  pstr <- paste0(format(c(round(1000-1000*pl),round(pl*1000))/10),"%")
  names(res) <- c("Estimate","Std.Err",pstr)
  res
}

##' @export
"constrain<-.multigroupfit" <-
  "constrain<-.multigroup" <- function(x,par,k=1,...,value) {
    constrain(Model(x)$lvm[[k]],par=par,...) <- value
    return(x)
}

##' @export
"constrain<-.default" <- function(x,par,args,...,value) {
    if (inherits(par,"formula")) {
        lhs <- getoutcome(par)
        xf <- attributes(terms(par))$term.labels
        par <- lhs
        if (par%in%vars(x)) {
            if (is.na(x$mean[[par]])) {
                intercept(x,par) <- par
            } else {
                par <- x$mean[[par]]
            }
        }
        args <- xf
    }
    if (is.null(value) || suppressWarnings(is.na(value))) {
        if (!is.null(par)) {
            Model(x)$constrain[[par]] <- NULL
            Model(x)$constrainY[[par]] <- NULL
        } else {
            Model(x)$constrain[[args]] <- NULL
        }
        return(x)
    }
    for (i in args) {
        if (!(i%in%c(parlabels(Model(x)),vars(Model(x)),
                     names(constrain(x))))) {
            if (!lava.options()$silent)
                message("\tAdding parameter '", i,"'\n",sep="")
            parameter(x,silent=TRUE) <- i
        }
    }

    if (par%in%vars(x)) {
        if (!"..."%in%names(formals(value))) {
            formals(value) <- c(formals(value),alist(...=))
        }
        Model(x)$constrainY[[par]] <- list(fun=value,args=args)
    } else {
        ## Wrap around do.call, since functions are not really
        ## parsed as call-by-value in R, and hence setting
        ## attributes to e.g. value=cos, will be overwritten
        ## if value=cos is used again later with new args.
        Model(x)$constrain[[par]] <- function(x) do.call(value,list(x))
        attributes(Model(x)$constrain[[par]])$args <- args
        index(Model(x)) <- reindex(Model(x))
    }
  return(x)
}

##' @export
constraints <- function(object,data=model.frame(object),vcov=object$vcov,level=0.95,
                        p=pars.default(object),k,idx,...) {
  if (class(object)[1]=="multigroupfit") {
    if (!missing(k)) {
      if (class(data)[1]=="list") data <- data[[k]]
      parpos <- modelPar(object, seq_len(length(p)))$p[[k]]
      if (nrow(data)>1 & !missing(idx)) {
        res <- t(apply(data,1,function(x) constraints(Model(object)$lvm[[k]],data=x,p=p[parpos],vcov=vcov[parpos,parpos],level=level)[idx,]))
        return(res)
      }
      return(constraints(Model(object)$lvm[[k]],data=data,p=p[parpos],vcov=vcov[parpos,parpos],level=level))
    }
    return(attributes(CoefMat.multigroupfit(object,data=data,vcov=vcov,...))$nlincon)
  }

  if (NROW(data)>1 & !missing(idx)) {
    res <- t(apply(data,1,function(x) constraints(object,data=x,p=p,vcov=vcov,level=level)[idx,],...))
    return(res)
  }

  if (length(index(object)$constrain.par)<1) return(NULL)
  parpos <- Model(object)$parpos
  if (is.null(parpos)) {
      parpos <- with(index(object),matrices2(Model(object),seq_len(npar+npar.mean+npar.ex)))
      parpos$A[index(object)$M0==0] <- 0
      parpos$P[index(object)$P0==0] <- 0
      parpos$v[index(object)$v1==0] <- 0
      parpos$e[index(object)$e1==0] <- 0
  }
  myidx <- unlist(lapply(parpos$parval, function(x) {
    if (!is.null(attributes(x)$reg.idx)) {
      return(parpos$A[attributes(x)$reg.idx[1]])
    } else if (!is.null(attributes(x)$cov.idx)) {
      return(parpos$P[attributes(x)$cov.idx[1]])
    } else if (!is.null(attributes(x)$m.idx)) {
      return(parpos$v[attributes(x)$m.idx[1]])
    } else if (!is.null(attributes(x)$e.idx))
        return(parpos$e[attributes(x)$e.idx[1]])
    else NA
  }))
  names(myidx) <- names(parpos$parval)
  mynames <- c()
  N <- length(index(object)$constrain.par)
  if (N>0)
  res <- c()
  count <- 0
  mydata <- rbind(numeric(length(manifest(object))))
  colnames(mydata) <- manifest(object)
  data <- rbind(data)
  iname <- intersect(colnames(mydata),colnames(data))
  mydata[1,iname] <- unlist(data[1,iname])
  for (pp in index(object)$constrain.par) {
    count <- count+1
    myc <- constrain(Model(object))[[pp]]
    mycoef <- numeric(6)
    val.idx <- myidx[attributes(myc)$args]
    val.idx0 <- na.omit(val.idx)
    M <- modelVar(Model(object),p=p,data=as.data.frame(mydata))
    vals <- with(M,c(parval,constrainpar))[attributes(myc)$args]
    fval <- mycoef[1] <- myc(unlist(vals))
    myc0 <- function(theta) {
      theta0 <- unlist(vals);
##    theta0[val.idx0] <- theta[val.idx0];
      theta0[!is.na(val.idx)] <- theta
      res <- myc(theta0)
      return(res)
    }
    vals0 <- unlist(vals)[!is.na(val.idx)]
##  vals0 <- unlist(vals)[na.omit(val.idx)]
    if (length(vals0)==0)
      mycoef[2] <- NA
    else {
      if (!is.null(attributes(fval)$grad)) {
        Gr <- cbind(attributes(fval)$grad(unlist(vals0)))
      } else {
        Gr <- cbind(as.numeric(numDeriv::jacobian(myc0, unlist(vals0))))
      }
      V <- vcov[val.idx0,val.idx0]
      mycoef[2] <- (t(Gr)%*%V%*%Gr)^0.5
    }
    ## if (second) {
    ##   if (!is.null(attributes(fval)$hessian)) {
    ##     H <- attributes(fval)$hessian(unlist(vals))
    ##   } else {
    ##     H <- hessian(myc, unlist(vals))
    ##   }
    ##   HV <- H%*%vcov[val.idx,val.idx]
    ##   mycoef[1] <- mycoef[1] + 0.5*sum(diag(HV))
    ##   mycoef[2] <- mycoef[2] + 0.5*sum(diag(HV%*%HV))
    ## }
    mycoef[3] <- mycoef[1]/mycoef[2]
    mycoef[4] <- 2*(pnorm(abs(mycoef[3]),lower.tail=FALSE))
    mycoef[5:6] <- mycoef[1] + c(1,-1)*qnorm((1-level)/2)*mycoef[2]
    res <- rbind(res,mycoef)
    mynames <- c(mynames,pp)
    if (!is.null(attributes(fval)$inv)){
      res2 <- attributes(fval)$inv(mycoef[c(1,5,6)])
      res <- rbind(res, c(res2[1],NA,NA,NA,res2[2],res2[3]))
      mynames <- c(mynames,paste0("inv(",pp,")"))
    }
  }
  rownames(res) <- mynames
  colnames(res) <- c("Estimate","Std. Error", "Z value", "Pr(>|z|)", "2.5%", "97.5%")
  return(res)
}
