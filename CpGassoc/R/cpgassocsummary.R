cpgassocsummary <-
function (object)    {
  if(!("aov" %in% class(object)))  {
    z <- object
    rdf <- object$df.residual
    p <- object$rank
    if (is.null(object$terms))
        stop("invalid 'lm' object:  no 'terms' component")
    Qr <- object$qr
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1L:p
    r <- z$residuals
    cpgsum<-function(x) sum(x**2)
    w <- z$weights
    if(is.null(ncol(r))) {
        r<-as.matrix(r)
        z$coefficients<-as.matrix(z$coefficients)
        }      
    if (is.null(w)) {
        rss<-apply(r,2,cpgsum)
    }
    else {
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])

    se <- sqrt(diag(R)[2] * resvar)
    if(ncol(r)==1)  est<-z$coefficients[Qr$pivot[p1],][2]
    else est <- z$coefficients[Qr$pivot[p1],][2,]
    tval <- est/se
    pval<-2 * pt(abs(tval),rdf, lower.tail = FALSE)
    cbind(tval,se,pval)
}
else {
  asgn <- object$assign[object$qr$pivot[1L:object$rank]]
    uasgn <- unique(asgn)
    nterms <- length(uasgn)
    effects <- object$effects
    if (!is.null(effects))
        effects <- as.matrix(effects)[seq_along(asgn), , drop = FALSE]
    rdf <- object$df.residual
    nmeffect <- c("(Intercept)", attr(object$terms, "term.labels"))
    coef <- as.matrix(object$coefficients)
    resid <- as.matrix(object$residuals)
    wt <- object$weights
    df <- ss <- numeric()
    ai<-(asgn==uasgn[nterms])
    df<-sum(ai)
    if(df>1 ) {
      ss<-c(ss,colSums(effects[ai,]**2))
             }
    else{
         lasteffects<-as.matrix(effects[ai,])
         ss<-c(ss,rowSums(lasteffects**2))
       }

    cpgsum<-function(x) sum(x**2)
    mse<-apply(resid,2,cpgsum)/rdf
    ms <- ss/df
    f.val<-matrix(nrow=length(ms))
    p.val<-matrix(nrow=length(ms))
       if (rdf > 0) {
           TT <- ms/mse
           TP <- pf(TT, df, rdf, lower.tail = FALSE)
           f.val <- TT
           p.val <- TP
              }
    cbind(f.val,p.val)
    
  }
  
}
