##' Performs Likelihood-ratio, Wald and score tests
##' @title Statistical tests
##' @aliases compare
##' @export
##' @param object \code{lvmfit}-object
##' @param \dots Additional arguments to low-level functions
##' @return Matrix of test-statistics and p-values
##' @author Klaus K. Holst
##' @seealso \code{\link{modelsearch}}, \code{\link{equivalence}}
##' @keywords htest
##' @examples
##' m <- lvm();
##' regression(m) <- c(y1,y2,y3) ~ eta; latent(m) <- ~eta
##' regression(m) <- eta ~ x
##' m2 <- regression(m, c(y3,eta) ~ x)
##' set.seed(1)
##' d <- sim(m,1000)
##' e <- estimate(m,d)
##' e2 <- estimate(m2,d)
##'
##' compare(e)
##'
##' compare(e,e2) ## LRT, H0: y3<-x=0
##' compare(e,scoretest=y3~x) ## Score-test, H0: y3~x=0
##' compare(e2,par=c("y3~x")) ## Wald-test, H0: y3~x=0
##'
##' B <- diag(2); colnames(B) <- c("y2~eta","y3~eta")
##' compare(e2,contrast=B,null=c(1,1))
##'
##' B <- rep(0,length(coef(e2))); B[1:3] <- 1
##' compare(e2,contrast=B)
##'
##' compare(e,scoretest=list(y3~x,y2~x))
compare <- function(object,...) UseMethod("compare")

##' @export
compare.default <- function(object,...,par,contrast,null,scoretest,Sigma,level=.95,df=NULL) {
  if (!missing(par) || (!missing(contrast) && is.character(contrast))) {
      if (!missing(contrast) && is.character(contrast)) par <- contrast
      contrast <- rep(0,length(coef(object)))
      myidx <- parpos(Model(object),p=par)
      contrast[myidx] <- 1
      contrast <- diag(contrast,nrow=length(contrast))[which(contrast!=0),,drop=FALSE]
      if (!missing(null) && length(null)>1) null <- null[attributes(myidx)$ord]
  }
  ### Wald test
  if (!missing(contrast)) {
    B <- contrast
    p <- coef(object)
    pname <- names(p)
    B <- rbind(B);
    colnames(B) <- if (is.vector(contrast)) names(contrast) else colnames(contrast)
    if (missing(Sigma)) {
      Sigma <- vcov(object)
    }
    if (ncol(B)<length(p)) {
      nn <- colnames(B)
      myidx <- parpos(Model(object),p=nn)
      B0 <- matrix(0,nrow=nrow(B),ncol=length(coef(object)))
      B0[,myidx] <- B[,attributes(myidx)$ord]
      B <- B0
    }
    if (missing(null)) null <- rep(0,nrow(B))
    if (length(null)==1) null <- rep(null,nrow(B))
    Bp <- B%*%p
    V <- B%*%Sigma%*%t(B)
    ct <- cbind(Bp,diag(V)^.5)
    p <- 1-(1-level)/2
    qp <- if(!is.null(df)) qt(p,df=df) else qnorm(p)
    ct <- cbind(ct,ct[,1] + qp*cbind(-1,1)%x%ct[,2])
    colnames(ct) <- c("Estimate","Std.Err",paste0(c(1-p,p)*100,"%"))
    rownames(ct) <- rep("",nrow(ct))
    Q <- t(Bp-null)%*%Inverse(V)%*%(Bp-null)
    df <- qr(B)$rank; names(df) <- "df"
    attributes(Q) <- NULL; names(Q) <- "chisq";
    pQ <- ifelse(df==0,NA,pchisq(Q,df,lower.tail=FALSE))

    method = "- Wald test -";
    cnames <- c()
    if (!is.null(pname)) {
      msg <- c()
      for (i in seq_len(nrow(B))) {
        Bidx <- which(B[i,]!=0)
        Bval <- abs(B[i,Bidx]); Bval[Bval==1] <- ""
        sgn  <- rep(" + ",length(Bval)); sgn[sign(B[i,Bidx])==-1] <- " - ";
        if (sgn[1]==" + ") sgn[1] <- "" else sgn[1] <- "-"
        cnames <- c(cnames,paste0(sgn,Bval,paste0("[",pname[Bidx],"]"),collapse=""))
        msg <- c(msg,paste0(cnames[i]," = ",null[i]))
      }
      method <- c(method,"","Null Hypothesis:",msg)
##      method <- c(method,"","Observed:",paste(formatC(as.vector(Bp)),collapse=" "))
    }

    res <- list(##data.name=hypothesis,
                statistic = Q, parameter = df,
                p.value=pQ, method = method, estimate=ct, vcov=V, coef=ct[,1],
                null=null, cnames=cnames
                )
    class(res) <- "htest"
    attributes(res)$B <- B
    return(res)
  }

  ### Score test
  if (!missing(scoretest)) {
    altmodel <- Model(object)
    if (inherits(scoretest,"formula")) scoretest <- list(scoretest)
    for (i in scoretest) {
      regression(altmodel) <- i
    }
    p0 <- numeric(length(coef(altmodel)))
    idx <-  match(coef(Model(object)),coef(altmodel))
    p0[idx] <- coef(object)
    Sc2 <- score(altmodel,p=p0,data=model.frame(object),weigth=Weight(altmodel),
                 estimator=object$estimator,...)
    I <- information(altmodel,p=p0,n=object$data$n,
                     data=model.frame(object),weigth=Weight(object),
                     estimator=object$estimator,...
                     )
    iI <- try(solve(I), silent=TRUE)
    Q <- ifelse (inherits(iI, "try-error"), NA, ## Score test
                 ## rbind(Sc)%*%iI%*%cbind(Sc)
                 (Sc2)%*%iI%*%t(Sc2)
                 )
    attributes(Q) <- NULL; names(Q) <- "chisq"
    df <- length(p0)-length(coef(object)); names(df) <- "df"
    pQ <- ifelse(df==0,NA,pchisq(Q,df,lower.tail=FALSE))
    res <- list(data.name=as.character(scoretest),
                statistic = Q, parameter = df,
                p.value=pQ, method = "- Score test -")
    class(res) <- "htest"
    return(res)
  }

  ### Likelihood ratio test
  objects <- list(object,...)
  if (length(objects)<2) {
    if (!(inherits(object,"lvmfit"))) {
      cc <- rbind(logLik(object),AIC(object))
      rownames(cc) <- c("logLik","AIC")
      colnames(cc) <-  ""
      return(cc)
    }
    L0 <- logLik(object)
    L1 <- satmodel(object,logLik=TRUE)
    df <- attributes(L1)$df-attributes(L0)$df; names(df) <- "df"
    Q <- abs(2*(L0-L1));
    attributes(Q) <- NULL; names(Q) <- "chisq";
    pQ <- ifelse(df==0,NA,pchisq(Q,df,lower.tail=FALSE))

    values <- c(L0,L1); names(values) <- c("log likelihood (model)", "log likelihood (saturated model)")
    res <- list(statistic = Q, parameter = df,
                p.value=pQ, method = "- Likelihood ratio test -",
                estimate = values)
    class(res) <- "htest"
    return(res)
  }
  if (length(objects)==2)
    return(comparepair(objects[[1]],objects[[2]]))
  res <- list()
  for (i in seq_len(length(objects)-1)) {
    res <- c(res, list(comparepair(objects[[i]],objects[[i+1]])))
  }
    return(res)
}


comparepair <- function(x1,x2) {
    l1 <- do.call("logLik",list(x1),envir=parent.frame(2))
    l2 <- do.call("logLik",list(x2),envir=parent.frame(2))
    df1 <- attributes(l1)$df;  df2 <- attributes(l2)$df;
    if (is.null(df1)) {
        df1 <- length(do.call("coef",list(x1),envir=parent.frame(2)))
        df2 <- length(do.call("coef",list(x2),envir=parent.frame(2)))
    }
    Q <- abs(2*(l1-l2))
    names(Q) <- "chisq"
    df <- abs(df1-df2); names(df) <- "df"
    p <- pchisq(Q,df=df,lower.tail=FALSE)
    values <- c(l1,l2); names(values) <- c("log likelihood (model 1)", "log likelihood (model 2)")

    res <- list(statistic = Q, parameter = df,
                p.value= p, method = "- Likelihood ratio test -",
                estimate = values)
    class(res) <- "htest"
    return(res)
}
