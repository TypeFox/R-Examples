"Survr" <-
function (id,time,event) 
{
 
    if(length(unique(id))!=length(event[event==0]))
      {
        stop("Data doesn't match. Every subject must have a censored time")
      }

    if(length(unique(event))>2 | max(event)!=1 | min(event)!=0)
      {
        stop("event must be 0-1")
      }

    ans<-cbind(id,time,event)

    oldClass(ans) <- "Survr"
    invisible(ans)

}

"is.Survr" <-
function(x)
inherits(x, "Survr")


"psh.fit" <-
function (x, tvals) 
{
    if (!is.Survr(x)) {
        stop("\n x must be a Survr object")
    }
    n <- length(unique(x[, 1]))
    failed <- c(x[, 2][x[, 3] == 1])
    censored <- c(x[, 2][x[, 3] == 0])
    m <- table(x[, 1]) - 1
    sfailed <- sort(failed)
    nfailed <- length(failed)
    summ <- .Fortran("distinctfailed", as.integer(n), as.integer(m), 
        as.double(failed), as.double(sfailed), as.integer(nfailed), 
        as.double(censored), as.integer(0), as.double(rep(0, 
            nfailed)), as.integer(rep(0, nfailed)), as.integer(rep(0, 
            n * nfailed)), PACKAGE = "survrec")
    numdistinct <- summ[[7]]
    distinct <- summ[[8]][1:numdistinct]
    numdeaths <- summ[[9]][1:numdistinct]
    vAtRisk <- summ[[10]][1:(n * numdistinct)]
    AtRisk <- matrix(vAtRisk, nrow = n, ncol = numdistinct)
    survfuncPSHple <- vector("numeric", numdistinct)
    AtRiskTotals <- t(AtRisk) %*% c(rep(1, n))
    survfuncPSHple.i <- 1 - (numdeaths/AtRiskTotals)
    survfuncPSHple <- cumprod(survfuncPSHple.i)

    se.NA<-cumsum(numdeaths/(AtRiskTotals)^2)
    se.PLE<-sqrt(se.NA*survfuncPSHple^2) 

    if (!missing(tvals)) {
        tvalslen <- length(tvals)
        tvals.o <- sort(tvals)
        PSHpleAttvals <- surv.search(tvals.o, distinct, survfuncPSHple)
    }
    else {
        tvals <- NA
        PSHpleAttvals <- NA
    }
    ans <- list(n = n, m = m, failed = failed, censored = censored, 
        time = distinct, n.event = numdeaths, AtRisk = AtRisk, 
        survfunc = survfuncPSHple, std.error=se.PLE,tvals = tvals, 
        PSHpleAttvals = PSHpleAttvals)
    oldClass(ans) <- "survfitr"
    ans
}




"wc.fit" <-
function (x, tvals) 
{
    if (!is.Survr(x)) {
        stop("\n x must be a Survr object")
    }
    n <- length(unique(x[, 1]))
    failed<-c(x[,2][x[,3]==1]) 
    gap <- c(x[, 2])
    cen.gap <- c(x[, 2][x[, 3] == 0])
    event <- c(x[,3])
    tot <- length(gap) 
    distinct <- sort(unique(c(x[, 2][x[, 3] == 1])))    
    ndistinct <- length(distinct)
    uid <- unique(x[, 1])
    m <- as.integer(table(x[, 1]))
    mMax<-max(m)
    wt<-rep(1,n)
        
    summ <- .Fortran("wc2",
		  as.integer(n),
		  as.double(matrix(0,nrow=n,ncol=mMax)),
		  as.double(wt),
		  as.double(m),
		  as.integer(mMax),
		  as.integer(m),
		  as.double(matrix(0,nrow=n,ncol=mMax)),
		  as.double(m - 1),
		  as.integer(ndistinct),
		  as.double(distinct),
		  as.integer(tot),
		  as.double(gap),		  
                  as.double(event),		  
                  r=as.double(rep(0,ndistinct)),
		  d=as.double(rep(0,ndistinct)),
		  surv=as.double(rep(0,ndistinct)),
		  var=as.double(rep(0,ndistinct)),  PACKAGE = "survrec")
 
    survfuncWCple <- summ$surv

    if (!missing(tvals)) {
        tvalslen <- length(tvals)
        tvals.o <- sort(tvals)
        WCpleAttvals <- surv.search(tvals.o, distinct, survfuncWCple)
    }
    else {
        tvals <- NA
        WCpleAttvals <- NA
    }
    ans <- list(n = n, m = m, failed = failed, censored = cen.gap, 
        time = distinct, n.event = summ$d, AtRisk = summ$r, 
        survfunc = survfuncWCple, std.error=summ$var, tvals = tvals, 
        WCpleAttvals = WCpleAttvals)
    oldClass(ans) <- "survfitr"
    ans
}


"mlefrailty.fit"<-
function (x, tvals, lambda = NULL, alpha = NULL, alpha.min, alpha.max, 
    tol = 1e-07, maxiter = 500, alpha.console=TRUE) 
{
    if (!is.Survr(x)) {
        stop("\n x must be a Survr object")
    }
    n <- length(unique(x[, 1]))
    failed <- c(x[, 2][x[, 3] == 1])
    censored <- c(x[, 2][x[, 3] == 0])
    m <- table(x[, 1])-1
    sfailed <- sort(failed)
    nfailed <- length(failed)
    summ <- .Fortran("distinctfailed", as.integer(n), as.integer(m), 
        as.double(failed), as.double(sfailed), as.integer(nfailed), 
        as.double(censored), as.integer(0), as.double(rep(0, 
            nfailed)), as.integer(rep(0, nfailed)), as.integer(rep(0, 
            n * nfailed)),PACKAGE="survrec")
    numdistinct <- summ[[7]]
    distinct <- summ[[8]][1:numdistinct]
    numdeaths <- summ[[9]][1:numdistinct]
    vAtRisk <- summ[[10]][1:(n * numdistinct)]
    AtRisk <- matrix(vAtRisk, nrow = n, ncol = numdistinct)
    if (is.null(lambda[1])) {
        lambda <- as.vector(numdeaths/apply(AtRisk,2,sum))
    }
    if (is.null(alpha)) {
 
        if(alpha.console)
              cat("\nNeeds to Determine a Seed Value for Alpha")
        if (missing(alpha.min)) {
            alpha.min <- 0.5
        }
        if (missing(alpha.max)) {
            alpha.max <- max(distinct)
        }
        tol.max <- (alpha.max - alpha.min)/50
        Seed <- .Fortran("searchforseed", as.integer(n), as.integer(m), 
            as.integer(numdistinct), as.double(distinct), as.integer(numdeaths), 
            as.integer(AtRisk), as.double(lambda), as.double(alpha.min), 
            as.double(alpha.max), as.double(tol.max), as.double(0), 
            as.integer(0),PACKAGE="survrec")
        alpha <- Seed[[11]]
        if(alpha.console)
               cat("\n Seed Alpha: ", alpha)
        IER <- Seed[[12]]
        if (IER != 0 && alpha.console) {
            warning("Problem with the seed value for alpha")
            if (IER == 129) {
                warning("bL is greater than or equal to bR. Minimum as bL")
            }
            if (IER == 130) {
                warning("tol is greater than the interval bL to bR")
            }
            if (IER == 131) {
                warning("the function is not unimodal. Check your results")
            }
        }
    }
    alphadel <- alpha/4
    alphaseeds <- c(alpha, alpha - alphadel, alpha - 2 * alphadel, 
        alpha - 3 * alphadel, alpha + alphadel, alpha + 2 * alphadel, 
        alpha + 3 * alphadel)
    status <- 0
    ind <- 0
    while ((status == 0) && (ind < 7)) {
        ind <- ind + 1
        alpha <- alphaseeds[ind]
        Estimates <- .Fortran("emalgo", as.integer(n), as.integer(m), 
            as.integer(numdistinct), as.double(distinct), as.integer(numdeaths), 
            as.integer(AtRisk), as.double(lambda), as.double(alpha), 
            as.double(tol), as.integer(maxiter), as.integer(status),PACKAGE="survrec")
        status <- Estimates[[11]]
    }
    alpha <- Estimates[[8]]
    if (alpha.console)
      { 
        cat("\n ")
        cat("\n Alpha estimate=", alpha)
        cat("\n ")
      }
    lambda <- Estimates[[7]]
    if (!missing(tvals)) 
        tvalslen <- length(tvals)
    if (!(status == 1)) {
        cat("\n\n WARNING: No estimates will be provided!")
        cat("\n Value of (status,alpha) from iteration is ", 
            c(status, alpha), "\n\n")
        alpha <- NA
        survfuncMLE <- c(rep(NA, numdistinct))
        if (!missing(tvals)) 
            MLEAttvals <- c(rep(NA, tvalslen))
        else {
            tvals <- NA
            MLEAttvals <- NA
        }
    }
    else {
        if (alpha >= 1e+05) 
            alpha <- 1e+05
        temp <- .Fortran("mlevalue", as.integer(numdistinct), 
            as.double(alpha), as.double(lambda), as.double(rep(0, 
                numdistinct)),PACKAGE="survrec")
        survfuncMLE <- temp[[4]]
        if (!missing(tvals)) {
            tvalslen <- length(tvals)
            tvals.o <- sort(tvals)
            MLEAttvals <- surv.search(tvals.o, distinct, survfuncMLE)
        }
        else {
            tvals <- NA
            MLEAttvals <- NA
        }
    }
    ans <- list(n = n, m = m, failed = failed, censored = censored, 
        time = distinct, n.event = numdeaths, AtRisk = AtRisk, 
        status = status, alpha = alpha, lambda=lambda, survfunc = survfuncMLE, 
        tvals = tvals, MLEAttvals = MLEAttvals)
    oldClass(ans) <- "survfitr"
    ans
}




"surv.search" <-
function (tvals,time,surv)
{
  time.c<-c(0,time)
  surv.c<-c(1,surv)
  mm<-outer(tvals,time.c,">=")
  pos<-apply(mm,1,sum)
  return(surv.c[pos])
}


"q.search" <-
function (f, q=0.5) 
{
     tt <- c(0, f$time)
     ss <- c(1, f$surv)
     if(ss[length(ss)] > q)
        stop(paste("\noverall survival estimate does not fall below ",q))
     ans<-min(tt[ss <= q])
     return(ans)
}



"survfitr" <-
function (formula, data, type="MLEfrailty",...) 
{
   method <- charmatch(type, c("pena-strawderman-hollander", 
             "wang-chang", "MLEfrailty"), nomatch= 0)
   if(method == 0)
	{
	  stop("estimator must be pena-strawderman-hollander wang-chang or MLEfrailty")
        }

    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) || 
        inherits(formula, "Survr")) {

    stop("formula.default(object): invalid formula")
     }

    m <- match.call(expand.dots = FALSE)
    m$type<- m$... <- NULL
    Terms <- terms(formula, "strata")
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Survr(Y)) 
        stop("Response must be a survival recurrent object")
    ll <- attr(Terms, "term.labels")
    
    if(method==1) FUN<-psh.fit
    if(method==2) FUN<-wc.fit
    if(method==3) FUN<-mlefrailty.fit

    if (ncol(m)>1) {
        group <- m[ll][, 1]  
        k <- levels(group)
        ans <- NULL
        for (i in 1:length(k)) {
            temp <- Y[group == k[i], ]
            temp1 <- Survr(temp[, 1], temp[, 2], temp[, 3])
            ans[[i]] <- FUN(temp1,...)
        }
        names(ans) <- k
        oldClass(ans) <- "survfitr"
        attr(ans, "strata") <- length(k)
        attr(ans, "group") <- ll
    }
    else {
        temp<-Survr(Y[,1],Y[,2],Y[,3]) 
        ans <- FUN(temp,...)
    }
    ans
}

"survdiffr" <-
function (formula, data, q, B=500, boot.F="WC",boot.G="none",...) 
{
   aux <- runif(1)

   method.F <- charmatch(boot.F, c("PSH","WC", "semiparametric"), nomatch= 0)
   if(method.F == 0)
	{
	  stop("bootstrap froom F must be PSH WC or semiparametric")
        }
   if(method.F == 3)
        {
          stop("Problems with Fortran code. Please contact with mantainer") 
        }

   method.G <- charmatch(boot.G, c("none","empirical"), nomatch= 0)
   if(method.G == 0)
	{
	  stop("bootstrap from G must be none or empirical")
        }

   if (method.F==1) {
        type.boot <- 2 + method.G - 1
        type <- "p"
     }  
   if (method.F==2) {
        type.boot <- 4 + method.G - 1 
        type <- "w"
     }  
   if (method.F==3) {
        type.boot <- 6 + method.G - 1
        type <- "M"
     }  


    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Survr")) || 
        inherits(formula, "Survr")) {

    stop("formula.default(object): invalid formula")
     }

    m <- match.call(expand.dots = FALSE)
    m$q<-  m$B<-  m$boot.F<- m$... <- NULL
    Terms <- terms(formula, "strata")
    ord <- attr(Terms, "order")
    if (length(ord) & any(ord != 1)) 
        stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)
    Y <- model.extract(m, "response")
    if (!is.Survr(Y)) 
        stop("Response must be a survival recurrent object")
    ll <- attr(Terms, "term.labels")
    group <- m[ll][, 1]


    if (!is.null(group)) {
        k <- levels(group)
        ans <-list(NULL)
        for (i in 1:length(k)) {
            temp <- Y[group == k[i], ]
            Sr <- Survr(temp[, 1], temp[, 2], temp[, 3])
            n <- length(unique(Sr[,1]))
            failed<-c(Sr[,2][Sr[,3]==1])
            censored<-c(Sr[,2][Sr[,3]==0])
            m <- table(Sr[, 1])-1

           sfailed<-sort(failed)
           nfailed<-length(failed)

           tau<-rep(NA,n)
           id.unique<-unique(Sr[,1])
           for (j in 1:n)
            {
             tau[j]<-sum(Sr[,2][Sr[,1]==id.unique[j]])
            }
           summ <- .Fortran("bootmedian",
                         as.integer(n),
                         as.integer(m),
                         as.double(failed),
                         as.double(sfailed),
                         as.integer(nfailed),
                         as.double(censored),
                         as.double(tau),
                         as.integer(B),
                         as.integer(type.boot), 
                         as.double(q),
                         as.double(rep(0,B)), PACKAGE = "survrec")

           ans[[i]]<-list(NULL)
           ans[[i]]$t0<-q.search(survfitr(Sr~1,type=type),q=q)
           ans[[i]]$t<-cbind(summ[[11]])
           ans[[i]]$R<-B
           ans[[i]]$data<-unclass(Sr)
           ans[[i]]$seed<-.Random.seed
           ans[[i]]$statistic<-NULL
           ans[[i]]$sim<-c("ordinary")
           ans[[i]]$call<-call
           ans[[i]]$stype<-c("i")
           ans[[i]]$strata<-rep(1,nrow(Sr))
           ans[[i]]$weights<-rep(1/nrow(Sr),nrow(Sr))
           oldClass(ans[[i]])<-"boot"
         }
         names(ans)<-k
    }
    else {
            Sr<-Survr(Y[,1],Y[,2],Y[,3]) 
            n <- length(unique(Sr[,1]))
            failed<-c(Sr[,2][Sr[,3]==1])
            censored<-c(Sr[,2][Sr[,3]==0])
            m <- table(Sr[, 1])-1

           sfailed<-sort(failed)
           nfailed<-length(failed)

           tau<-rep(NA,n)
           id.unique<-unique(Sr[,1])
           for (i in 1:n)
            {
             tau[i]<-sum(Sr[,2][Sr[,1]==id.unique[i]])
            }
           summ <- .Fortran("bootmedian",
                         as.integer(n),
                         as.integer(m),
                         as.double(failed),
                         as.double(sfailed),
                         as.integer(nfailed),
                         as.double(censored),
                         as.double(tau),
                         as.integer(B),
                         as.integer(type.boot), 
                         as.double(q),
                         as.double(rep(0,B)), PACKAGE = "survrec")
                     
           ans<-NULL
           ans$t0<-q.search(survfitr(Sr~1,type=type),q=q)
           ans$t<-cbind(summ[[11]])
           ans$R<-B
           ans$data<-unclass(Sr)
           ans$seed<-.Random.seed
           ans$statistic<-NULL
           ans$sim<-c("ordinary")
           ans$call<-call
           ans$stype<-c("i")
           ans$strata<-rep(1,nrow(Sr))
           ans$weights<-rep(1/nrow(Sr),nrow(Sr))
           oldClass(ans)<-"boot"
    
    }

ans

}



"plot.survfitr" <-
function (x, conf.int=TRUE, prob = FALSE, ...) 
{
    dostep <- function(x, y) {
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }
    y.lab <- ifelse(prob, "Probability Estimates", "Survivor Probability Estimates")
    if (!prob) {
        if (!is.null(attr(x, "strata"))) {
            plot(dostep(x[[1]]$time, x[[1]]$surv), type = "n", 
                xlab = "Time", ylab = y.lab, ...)
            for (i in 1:attr(x, "strata")) {
                y <- x[[i]]$surv
                if((is.null(x[[i]]$std.error)) || (!conf.int))
                   e<-1
                else 
                   e <- x[[i]]$std.error
                lines(dostep(x[[i]]$time, y), col = i)
                lines(dostep(x[[i]]$time, y+qnorm(0.975)*e), col=i, lty=2)
                lines(dostep(x[[i]]$time, y-qnorm(0.975)*e), col=i, lty=2)
            }
        }
        else {
            y <- x$surv
            if((is.null(x$std.error)) || (!conf.int))
              e<-1
            else 
              e <- x$std.error
            plot(dostep(x$time, y), type = "l", ylim = c(0, max(y)), 
                xlab = "Time", ylab = y.lab)
            lines(dostep(x$time, y+qnorm(0.975)*e),lty=2)
            lines(dostep(x$time, y-qnorm(0.975)*e),lty=2)

        }
    }
    else {
        if (!is.null(attr(x, "strata"))) {
            plot(dostep(x[[1]]$time, 1 - x[[1]]$surv), type = "n", 
                xlab = "Time", ylab = y.lab, ...)
            for (i in 1:attr(x, "strata")) {
                y <- x[[i]]$surv
                if((is.null(x[[i]]$std.error)) || (!conf.int))
                   e<-1
                else 
                   e <- x[[i]]$std.error
                lines(dostep(x[[i]]$time, 1-y), col = i)
                lines(dostep(x[[i]]$time, 1-y+qnorm(0.975)*e), col=i, lty=2)
                lines(dostep(x[[i]]$time, 1-y-qnorm(0.975)*e), col=i, lty=2)

            }
        }
        else {
            y <- x$surv
            if((is.null(x$std.error)) || (!conf.int))
                e<-1
            else 
                e <- x$std.error
            plot(dostep(x$time, 1 - y), type = "l", ylim = c(0, 
                max(y)), xlab = "Time", ylab = y.lab)
            lines(dostep(x$time, 1 - y+qnorm(0.975)*e),lty=2)
            lines(dostep(x$time, 1 - y-qnorm(0.975)*e),lty=2)
        }
    }
    return(invisible())
}



"lines.survfitr"<-
function (x, prob=FALSE, ...) 
{
    dostep <- function(x, y) {
        if (is.na(x[1] + y[1])) {
            x <- x[-1]
            y <- y[-1]
        }
        n <- length(x)
        if (n > 2) {
            dupy <- c(TRUE, diff(y[-n]) != 0, TRUE)
            n2 <- sum(dupy)
            xrep <- rep(x[dupy], c(1, rep(2, n2 - 1)))
            yrep <- rep(y[dupy], c(rep(2, n2 - 1), 1))
            list(x = xrep, y = yrep)
        }
        else if (n == 1) 
            list(x = x, y = y)
        else list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
    }
if(!prob)
    lines(dostep(x$time, x$survfunc), ...)
else
    lines(dostep(x$time, 1-x$survfunc), ...)
    return(invisible())
}



"print.survfitr" <-
function (x,scale=1,digits = max(options()$digits - 4, 3), ...) 
{

  savedig <- options(digits = digits)
  on.exit(options(savedig))
  plab<-c("n","events","mean","se(mean)","median","recurrences: min","max","median")

  pfun <- function(x)
    #compute the mean, se(mean) and median survival
	{
          minmin <- function(y, xx)
            {
	     if(any(!is.na(y) & y == 0.5)) {
		if(any(!is.na(y) & y < 0.5))
		  0.5 * (min(xx[!is.na(y) & y == 0.5]) + min(xx[!is.na(y) & y < 0.5]))
		else 0.5 * (min(xx[!is.na(y) & y == 0.5]) + max(xx[!is.na(y) & y == 0.5]))
			}
			else min(xx[!is.na(y) & y <= 0.5])
		}

                stime<-x$time/scale
           	    n <- length(stime)
          	    if (is.matrix(x$AtRisk))
                      {n.risk<-apply(x$AtRisk,2,sum)}
                    else
                      {n.risk<-x$AtRisk}
                hh <- c(x$n.event[ - n]/(n.risk[ - n] * (n.risk[ - n] - x$n.event[ - n])), 0)
	       
                med <- minmin(x$survfunc, x$time)
               

                dif.time <- c(diff(c(0, stime)), 0)
                mean <- dif.time * c(1, x$survfunc)

                varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)                 

                ans<-c(x$n,sum(x$m),sum(mean),sqrt(varmean),med,min(x$m),max(x$m),median(x$m))
                ans
}


if(is.null(attr(x,"strata")))
  {
   # no strata 
    x1<-rbind(pfun(x))
    cat("Survival for recurrent event data")
    cat("\n")
    dimnames(x1)<-list(" ",plab)
    print(x1)
    cat("\n")
  }

else
  {
      cat("Survival for recurrent event data. Group=" ,attr(x,"group"))
      cat("\n")  
      x1<-NULL
      for (i in 1:attr(x,"strata"))
        {
         temp<-pfun(x[[i]])
         x1<-rbind(x1,temp)
        }
      temp<-names(x)
      dimnames(x1)<-list(temp,plab)
      print(x1)
      cat("\n")
   }

invisible(x)

}


"summary.survfitr" <-
function (object, ...) 
{
    x <- object
    if (!inherits(x, "survfitr")) 
        stop("Invalid data")
    if (is.null(x[[1]]$std.error)) {
        plab <- c("time", "n.event", "n.risk", "surv")
        if (!is.null(attr(x, "strata"))) {
            ans <- list(NA)
            for (i in 1:attr(x, "strata")) {
                if (is.matrix(x[[i]]$AtRisk))
                    {n.risk<-apply(x[[i]]$AtRisk,2,sum)}
                else
                    {n.risk<-x[[i]]$AtRisk}
                temp <- cbind(x[[i]]$time, x[[i]]$n.event, n.risk, 
                  x[[i]]$surv)
                dimnames(temp) <- list(rep("", nrow(temp)), plab)
                ans[[i]] <- temp
            }
            names(ans) <- names(x)
            oldClass(ans) <- "summary.survfitr"
            attr(ans, "strata") <- attr(x, "strata")
        }
        else {
            if (is.matrix(x$AtRisk))
              n.risk<-apply(x$AtRisk,2,sum)
            else
              n.risk<-x$AtRisk
            temp <- cbind(x$time, x$n.event, n.risk, x$surv)
            dimnames(temp) <- list(rep("", nrow(temp)), plab)
            ans <- temp
        }
    }
    else {
        plab <- c("time", "n.event", "n.risk", "surv", "std.error")
        if (!is.null(attr(x, "strata"))) {
            ans <- list(NA)
            for (i in 1:attr(x, "strata")) {
                if (is.matrix(x[[i]]$AtRisk))
                    {n.risk<-apply(x[[i]]$AtRisk,2,sum)}
                else
                    {n.risk<-x[[i]]$AtRisk}
                temp <- cbind(x[[i]]$time, x[[i]]$n.event, n.risk, 
                  x[[i]]$surv, x[[i]]$std.error)
                dimnames(temp) <- list(rep("", nrow(temp)), plab)
                ans[[i]] <- temp
            }
            names(ans) <- names(x)
            oldClass(ans) <- "summary.survfitr"
            attr(ans, "strata") <- attr(x, "strata")
        }
        else {
            if (is.matrix(x$AtRisk))
              n.risk<-apply(x$AtRisk,2,sum)
            else
              n.risk<-x$AtRisk
            temp <- cbind(x$time, x$n.event, n.risk, x$surv, 
                x$std.error)
            dimnames(temp) <- list(rep("", nrow(temp)), plab)
            ans <- temp
        }
    }
    ans
}



"print.summary.survfitr" <-
function (x,scale=1,digits = max(options()$digits - 4, 3), ...) 
{
  savedig <- options(digits = digits)
  on.exit(options(savedig))

if(is.null(attr(x,"strata")))
  {
   # no strata 
    cat("\n")
    temp<-x
    oldClass(temp)<-NULL
    print(temp)
    cat("\n")
  }

else
 { 
  for (i in 1:attr(x,"strata"))
    {
     cat("\n      Group=",names(x)[i])
     cat("\n")
     print(x[[i]])
     cat("\n")
    }
 }
 invisible(x)
}



############ First.lib ###############

.onLoad <- function(lib, pkg){
   library.dynam("survrec", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("survrec", libpath)


############ End of .First.lib ###############
