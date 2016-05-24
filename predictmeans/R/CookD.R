CookD <- function (model, group=NULL, plot=TRUE, idn=3, newwd=TRUE) {

  if (class(model)[1] %in% c("gls", "lme")) model <- update(model, method="ML") 
  if (class(model)[1]=="lmerMod") model <- update(model, REML=FALSE)
  if (class(model)[1] =="gls") {
    mdf <- nlme::getData(model)
  }else{
    mdf <- model.frame(model)
  }
  mp <- mymodelparm(model)
  beta0  <- mp$coef
  vcovb <- mp$vcov
  vb.inv <- solve(vcovb)

  if (is.null(group) || group%in%c("NULL", "")) {
    rn <- rownames(mdf)
    LOOmp <- lapply(rn, function(x)  mymodelparm(update(model, data=mdf[rn!=x, ])))
  }else{
    rn <- unique(mdf[, group])
    LOOmp <- lapply(rn, function(x)  {
  	rind <- mdf[, group]!=x
  	mymodelparm(update(model, data=mdf[rind, ]))
    })
  }

  LOObeta <- sapply(LOOmp, function(x) x$coef)          # Matrix with beta(-i)  
  rK <- t(LOObeta-beta0)  
  CookD <- diag(rK %*% tcrossprod(vb.inv, rK)/length(beta0))
  names(CookD) <- rn
  if (plot) {
    if (newwd) dev.new()
    outD <- CookD >= sort(CookD, decreasing =TRUE)[idn]                        # Outlying Di's
    labid <- names(CookD)
    plot(CookD,  xlab="Obs. number", col="blue", ylim=c(0, max(CookD)+0.005),
         main="Cook's Distance", ylab = "Cook's distance", type = "h")
    text((1:length(CookD))[outD], CookD[outD], labid[outD], pos=3)          # Annotation
    points((1:length(CookD))[outD], CookD[outD], pch=16, col="blue")
  }
  return(invisible(CookD))
}
