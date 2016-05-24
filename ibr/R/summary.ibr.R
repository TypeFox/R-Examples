summary.ibr <- function(object, criteria="call", ...) {
  r <- object$residuals
  tt <- terms(object)
  colnamesx <- delete.response(tt)
  n <- length(r)
  sigma2 <- sum(r^2)/(n)
  stderr <- sqrt(n*sum(r^2)/(n-object$finaldf))
  if (any(criteria=="call")) {
      criteria <- object$parcall$criterion
      if (criteria=="user") {
          anscrit <- NULL
      } else  {
          anscrit <- object$allcriteria
          names(anscrit) <- criteria
      }
    itercrit <- object$alliter
  } else {
    crit <-c("aic","aicc","gcv","bic","gmdl")
    if (all(!(criteria%in%crit))) stop(paste("criteria are:",crit,"\n"))
    criteria <- criteria[criteria%in%crit]
    anscrit <- NULL
    if (any(criteria=="gcv"))  anscrit <- c(anscrit,log(sigma2)-2*log(1-object$finaldf/n))
    if (any(criteria=="aic"))  anscrit <- c(anscrit,log(sigma2)+2*object$finaldf/n)
    if (any(criteria=="aicc"))  anscrit <- c(anscrit,log(sigma2)+1+(2*(object$finaldf+1))/(n-object$finaldf-2))
    if (any(criteria=="bic"))  anscrit <- c(anscrit,log(sigma2) + log(n)*(object$finaldf)/n)
    if (any(criteria=="gmdl")) {
        Sbul <-   n*sigma2/(n-object$finaldf)
        mf <- object$call
        m <- match(c("formula", "data", "subset"), names(mf), 0L)
        mf <- mf[c(1L, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1L]] <- quote(stats::model.frame)
        mf <- eval(mf, environment(tt))
        y <- model.response(mf, "numeric")
        anscrit <- c(anscrit,log(Sbul)+object$finaldf/n*log((sum(y^2)-n*sigma2)/(object$finaldf*Sbul)))
    }
    criteria <- criteria[criteria!="user"]
    names(anscrit) <- criteria
    itercrit <- rep(NA,length(anscrit))
    names(itercrit) <- criteria
    if ((all(criteria!="user"))&&(any(criteria %in%object$parcall$criterion))) {
      itercrit[criteria%in%object$parcall$criterion] <- object$alliter[object$parcall$criterion%in%criteria]
      names(itercrit)[criteria%in%object$parcall$criterion] <- object$parcall$criterion[object$parcall$criterion%in%criteria]
    }
  }
  if (object$parcall$critmethod=="aggregation") crit4iter <- paste("aggregation of:",paste(object$parcall$criterion,collapse=", ")) else crit4iter <- object$parcall$criterion[1]
  ans <- list(residuals=r,Std.Error=stderr,Initial.Df=object$initialdf,
              Final.Df=object$finaldf,Resid.Df=n-object$finaldf,criteria=anscrit,
              iterations=itercrit,kernel=object$parcall$kernel, iter=object$iter,
              crit4iter=crit4iter,bandwidth=object$bandwidth,
              smoother=object$parcall$smoother,m=object$parcall$m,s=object$parcall$s)
  if (object$parcall$smoother=="k") {
    if (is.null(colnamesx)) {
      names(ans$bandwidth) <- paste("X",1:object$parcall$p,sep="")
    } else {
      names(ans$bandwidth) <- paste(colnamesx,1:object$parcall$p,sep="")
    }
  }
  class(ans) <- "summary.ibr"
  ans
}
