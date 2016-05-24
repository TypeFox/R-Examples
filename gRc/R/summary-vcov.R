
summary.rcox <- function(object, type="coef",...){
  type  <- match.arg(type,c("coef","K","KC","ACA"))
  m     <- object
  vn <- unlist(lapply(getcc(object),names))
  cctype  <- c(rep("vcc", length(getSlot(m,"vcc"))), rep("ecc", length(getSlot(m,"ecc"))))

  ans   <- list(type=type,
                self          = object,
                vcc           = getSlot(m,"vcc"),
                ecc           = getSlot(m,"ecc"), 
                logL          = fitInfo(m,"logL"),
                dimension     = dimension(m),
                method        = getSlot(m,"method"),
                time          = fitInfo(m,"time"),
                short         = object$control$short
                )
  switch(type,
         "coef"={
           cm    <- as.numeric(coef(m))
           if (!is.null(fitInfo(m)$J)){
             V     <- vcov(m)
             vv    <- (diag(V))             
             X2    <- cm^2/vv
             p     <- 1-pchisq(X2, 1)           
             v     <- data.frame(cctype=cctype, cc=vn, estimate=cm,stderr=sqrt(vv),
                                 X2=X2,p=p)
           } else {
             v     <- data.frame(cctype=cctype, cc=vn, estimate=cm)
           }
           ans$coefficients <- v
           rownames(v) <- NULL           
         },
         "K"={
           ans$K <- fitInfo(m,"K")
         },
         "KC"={
           KC <- fitInfo(m,"K")
           dKC <- diag(KC)
           CorMatrix <- -cov2cor(KC)
           KC[lower.tri(KC)] <- CorMatrix[lower.tri(CorMatrix)]
           diag(KC) <- sqrt(dKC)
           ans$KC <- KC
         },
         "ACA"={
           K <- fitInfo(m,"K")
           C <- -cov2cor(K)
           diag(C) <- 1
           ans$A <- sqrt(diag(K))
           ans$C <- C

         }
         
         )
  class(ans) <- "summary.rcox"
  ans
}


print.summary.rcox <- function(x, ...){
#   cat("vcc: ", cc2str(x$vcc),"\n")
#   cat("ecc: ", cc2str(x$ecc),"\n")

#   cat("logL: ", x$logL, "dimension:", x$dimension, "\n")
#   cat("method:", x$method, "time taken:", x$time, "\n")

  print(x$self)
  switch(x$type,
         "coef"={
           cat("\n")
           print(x$coefficients)
           if (!x$short){
             cat("\n");print(getvcc(x));print(getecc(x))
           }
         },
         "K"={
           print(x$K)
         },
         "KC"={
           print(x$KC)
         },
         "ACA"={
           print(x$A)
           print(x$C)
         }
         )
  return(invisible(x))
}

vcov.rcox <- function(object, ...){
  f1     <- fitInfo(object)
  solve(f1$J)
}
