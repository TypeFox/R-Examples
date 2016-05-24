"summary.longiPenal"<-
  function(object,level=.95, len=6, d=2, lab=c("coef","hr"), ...)
  {
    x <- object
    if (!inherits(x, "longiPenal"))
      stop("Object must be of class 'longiPenal'")

    z<-abs(qnorm((1-level)/2))
    co <- x$coef
    se <- sqrt(diag(x$varH))
    or <- exp(co)
    li <- exp(co-z * se)
    ls <- exp(co+z * se)
    p <-  signif(1 - pchisq((co/se)^2, 1), 5)

    rl <- cbind(co, se, p)

    dimnames(rl) <- list(names(co), c(lab[1], "SE","p"))

    ddl <- dim(rl)

    al <- formatC(rl, d, len,format="f")

    dim(al) <- ddl
    if(length(ddl) == 1){
      ddl<-c(1,ddl)
      dim(al)<-ddl
      labl<-" "
    }
    else
      labl <- dimnames(rl)[[1]]

    mxl <- max(nchar(labl)) + 1

    al[which(al[,3]==formatC(0, d, len,format="f")),3]<-formatC("<1e-16", d, len,format="f")

    cat("Longitudinal outcome:\n")
    cat("------------- \n")
    cat(paste(rep(" ",mxl),collapse=""),paste("  ",dimnames(rl)[[2]]),"\n")
    for(i in (x$nvar[1]+1):ddl[1])
    {
      labl[i] <- paste(c(rep(" ", mxl - nchar(labl[i])), labl[i]),collapse = "")
      cat(labl[i], al[i, 1:3]," \n")
    }



    r <- cbind(or, li, ls)

    dimnames(r) <- list(names(co), c(lab[2], paste(level*100,"%",sep=""), "C.I."))

    n<-r

    dd <- dim(n)
    n[n > 999.99] <- Inf
    a <- formatC(n, d, len,format="f")

    dim(a) <- dd
    if(length(dd) == 1){
      dd<-c(1,dd)
      dim(a)<-dd
      lab<-" "
    }
    else
      lab <- dimnames(n)[[1]]

    mx <- max(nchar(lab)) + 1

    cat("\n")
    cat("Terminal event:\n")
    cat("--------------- \n")
    cat(paste(rep(" ",mx),collapse=""),paste("   ",dimnames(n)[[2]]),"\n")
    for(i in 1:x$nvar[1])
    {
      lab[i] <- paste(c(rep(" ", mx - nchar(lab[i])), lab[i]),collapse = "")
      cat(lab[i], a[i, 1], "(", a[i, 2], "-", a[i, 3], ") \n")
    }


  }
