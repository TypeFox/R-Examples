print.DirichletRegConfint <- function(x, digits=3, ...){

  e <- x$e
  repar <- x$repar
  cat("\n", paste(100*x$level, collapse="%, ", sep=""), "% Confidence Intervals ", sep="")
  if(e) cat("(exponentiated)\n\n") else cat("(original form)\n\n")

  ctab <- x$coefficients

  if(repar){

    if(e){
      for(v in 1:length(ctab)){
        ctab[[v]] <- lapply(ctab[[v]], function(l) if(is.null(l)) NULL else exp(l))
      }
    }

    iter <- 1
    for(cc in rev(seq_along(x$ci[[1]]))[-1]){
        iter <- iter + 1
      for(ll in seq_along(x$ci)){
        if(is.null(x$ci[[1]][[cc]])) next
        ctab[[1]][[cc]] <- cbind(ctab[[1]][[cc]], x$ci[[ll]][[cc]])
      }
    }
    for(ll in seq_along(x$ci)){
      ctab[[2]][[1]] <- cbind(ctab[[2]][[1]], x$ci[[ll]][[iter]])
    }

    for(cc in seq_along(ctab[[1]])){
      if(is.null(ctab[[1]][[cc]])) next
      ll <- ncol(ctab[[1]][[cc]])
      ind <- c(rev(seq(2, ll, by=2)), seq(1, ll, by=2))

      ctab[[1]][[cc]] <- ctab[[1]][[cc]][,ind,drop=FALSE]
    }
    ll <- ncol(ctab[[2]][[1]])
    ind <- c(rev(seq(2, ll, by=2)), seq(1, ll, by=2))
    ctab[[2]][[1]] <- ctab[[2]][[1]][,ind,drop=FALSE]

  } else {

  if(e){
    ctab <- lapply(ctab, exp)
  }

    for(cc in rev(seq_along(x$ci[[1]]))){
      for(ll in seq_along(x$ci)){
        ctab[[cc]] <- cbind(ctab[[cc]], x$ci[[ll]][[cc]])
      }
    }
    for(cc in seq_along(x$ci[[1]])){
      ll <- ncol(ctab[[cc]])
      ind <- c(rev(seq(2, ll, by=2)), seq(1, ll, by=2))

      ctab[[cc]] <- ctab[[cc]][ , ind, drop=FALSE ]
    }

  }


  lo_lab <- paste(format(100*(1 - rev(x$level))/2, format="f"), "%", sep="")
  hi_lab <- paste(format(100*(x$level + (1 - x$level)/2), format="f"), "%", sep="")

  if(repar){
    for(tt in 1:2){
    cat("- ",ifelse(tt == 1, "Beta-Parameters:\n", "Gamma-Parameters\n"),sep="")
      for(i in seq_along(ctab[[tt]])){
        if(tt == 1) cat("Variable: ",names(ctab[[1]])[i],"\n",sep="")
        if(is.null(ctab[[tt]][[i]])){cat("  variable omitted\n\n"); next}
        ttab <- round(ctab[[tt]][[i]], digits)
        colnames(ttab) <- c(lo_lab, ifelse(e, "exp(Est.)", "Est."), hi_lab)
        print(ttab, digits=digits, print.gap=2)
        cat("\n")
      }
    }
  } else {
    for(i in seq_along(ctab)){
      cat("Variable: ",names(x$coefficients)[i],"\n",sep="")
      ttab <- round(ctab[[i]], digits)
      colnames(ttab) <- c(lo_lab, ifelse(e, "exp(Est.)", "Est."), hi_lab)
      print(ttab, digits=digits, print.gap=2)
      cat("\n")
    }
  }

}
