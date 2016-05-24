summary.mcprofile <-
function(object, margin=0, adjust="single-step", alternative=c("two.sided","less","greater"), ...){
  CM <- object$CM
  est <- coefficients(object$object)
  df <- object$df
  if (!(alternative[1] %in% c("two.sided", "less", "greater"))) stop("alternative has to be one of 'two.sided', 'less', or 'greater'")
  pam <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")  
  if (!(adjust[1] %in% c(pam, "single-step"))) stop(paste("adjust has to be one of:", paste(c(pam, "single-step"), collapse=", ")))
  ##### Find stats
  sdlist <- object$srdp
  spl <- lapply(sdlist, function(x){
    x <- na.omit(x)
    try(interpSpline(x[,1], x[,2]))
  })
  ptest <- function(spl, delta){
    pst <- try(predict(spl, delta)$y, silent=TRUE)
    if (class(pst) == "try-error") NA else pst   
  }
  stat <- sapply(spl, function(x) -1*ptest(x, delta=margin))
  naid <- is.na(stat)
  stat[is.na(stat)] <- 0
  if (is.null(df)){
    switch(alternative[1], less = {
      praw <- pnorm(stat, lower.tail=TRUE)
    }, greater = {
      praw <- pnorm(stat, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pnorm(abs(stat), lower.tail=FALSE)*2)
    })
  } else {
    switch(alternative[1], less = {
      praw <- pt(stat, df=df, lower.tail=TRUE)
    }, greater = {
      praw <- pt(stat, df=df, lower.tail=FALSE)
    }, two.sided = {
      praw <- pmin(1, pt(abs(stat), df=df, lower.tail=FALSE)*2)
    })
  }
  if (length(praw) > 1){
    if (adjust[1] %in% pam) padj <- p.adjust(praw, method=adjust[1])
    pfct <- function(q) {
      switch(alternative[1], two.sided = {
       low <- rep(-abs(q), dim)
       upp <- rep(abs(q), dim)
      }, less = {
       low <- rep(q, dim)
       upp <- rep(Inf, dim)
      }, greater = {
       low <- rep(-Inf, dim)
       upp <- rep(q, dim)
      })
      if (is.null(df)){
       pmvnorm(lower = low, upper = upp, corr = cr)
      } else {
       pmvt(lower = low, upper = upp, df=df, corr = cr)
      }
    }
    if (adjust[1] == "single-step"){
      vc <- vcov(object$object)
      VC <- CM %*% vc %*% t(CM)
      d <- 1/sqrt(diag(VC))
      dd <- diag(d)
      cr <- dd %*% VC %*% dd
      dim <- ncol(cr)
      padj <- numeric(length(stat))
      error <- 0
      for (i in 1:length(stat)) {
       tmp <- pfct(stat[i])
       if (error < attr(tmp, "error")) error <- attr(tmp, "error")
       padj[i] <- tmp
      }
      padj <- 1 - padj
      attr(padj, "error") <- error
    }
  } else { padj <- praw }
  stat[naid] <- NA
  padj[naid] <- NA

  out <- list()
  out$CM <- CM
  out$statistic <- stat
  out$p.values <- padj
  out$margin <- margin
  out$estimate = est
  out$alternative <- alternative[1]
  out$adjust <- adjust[1]
  #out$cr <- cr
  out$df <- df
  class(out) <- "mcpSummary"
  return(out)
}

