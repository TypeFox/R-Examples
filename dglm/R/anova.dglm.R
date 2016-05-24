anova.dglm <- function(object, ..., dispersion = NULL, test = NULL)
{
  #  ANOVA for double glm (likelihood ratio tests for mean and dispersion models)
  #  GKS  16 Jan 98
  #
  response <- as.character(stats::formula(object)[2])
  mterms <- object$terms
  dterms <- object$dispersion.fit$terms
  df <- c(length(mterms),length(dterms))
  #
  #  Extract null mean model formula
  #
  if(df[1]>0) {
    o <- attr(mterms, "offset")
    i <- attr(mterms, "intercept")
    factor.labels <- dimnames(attr(mterms,"factors"))[[1]]
    fm <- paste(response,"~")
    if(!is.null(o)) {
      fm <- paste(fm,factor.labels[o])
      if(i == 0)
        fm <- paste(fm,"-1")
    }
    else {
      if(i == 0)
        fm <- paste(fm,"-1")
      else
        fm <- paste(fm,"1")
    }
    fm <- parse(text=fm)
    mode(fm) <- "call"
    df[1] <- object$rank - i
  }
  #
  #  Extract null dispersion model formula
  #
  if(df[2]>0) {
    o <- attr(dterms, "offset")
    i <- attr(dterms, "intercept")
    factor.labels <- dimnames(attr(dterms,"factors"))[[1]]
    if(!is.null(o)) {
      fd <- paste("~",factor.labels[o])
      if(i == 0)
        fd <- paste(fd,"-1")
    }
    else {
      if(i == 0)
        fd <- "~ -1"
      else
        fd <- "~ 1"
    }
    fd <- parse(text=fd)
    mode(fd) <- "call"
    df[2] <- object$dispersion.fit$rank - i
  }
  #
  #  Fit null models and store likelihoods
  # 
  lik <- rep(object$m2loglik,4)
  names(lik) <- c("Null","Mean","Full","Disp")
  ncall <- object$call
  if(df[2]>0) {   
    ncall["dformula"] <- fd
    lik["Mean"] <- eval(ncall)$m2loglik
    if(df[1]>0) {
      ncall["formula"] <- fm
      lik["Null"] <- eval(ncall)$m2loglik
      ncall["dformula"] <- object$call["dformula"]
      ncall["formula"] <- fm
      lik["Disp"] <- eval(ncall)$m2loglik
    }
    else {
      lik["Null"] <- lik["Mean"]
      lik["Disp"] <- lik["Full"]
    }
  }
  else {
    lik["Mean"] <- lik["Full"]
    if(df[1]>0) {
      ncall["formula"] <- fm
      lik["Null"] <- eval(ncall)$m2loglik
      lik["Disp"] <- lik["Null"]
    }
    else
      lik["Disp"] <- lik["Null"] <- lik["Full"]
  }
  seqdev <- c(lik["Null"] - lik["Mean"], lik["Mean"] - lik["Full"])
  adjdev <- c(lik["Disp"] - lik["Full"], lik["Mean"] - lik["Full"])
  #
  #  Output table
  #
  heading <- c("Analysis of Deviance Table",
               paste("\n",object$family$family," double generalized linear model",sep = ""),
               paste("\nResponse: ", response,"\n", sep = "") )
  aod <- data.frame(row.names = c("Mean model","Dispersion model"),
                    DF = df,
                    Seq.Chisq = seqdev,
                    Seq.P = 1 - stats::pchisq(seqdev,ifelse(df > 0,df,NA)),
                    Adj.Chisq = adjdev,
                    Adj.P = 1 - stats::pchisq(adjdev,ifelse(df > 0,df,NA)))
  structure(aod, heading = heading, class = c("anova","data.frame") )
  
}
