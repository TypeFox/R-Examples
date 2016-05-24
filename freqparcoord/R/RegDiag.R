
# regression fit diagnostics based on freqparcoord()

# one axis is the "divergences," the differences beween the parametric
# and nonparametric estimates of the regression function, while the
# other axes are the predictor variables; the goal is to use the graph
# to determine regions of under- and overfit in predictor space

# the divergences are grouped into 3 quantile segments; for example, if
# the user sets 'tail' to 0.10, then 3 groups will be formed, as
# follows:  the lowest 10% of the divergences, the highest 10%, and the
# remainder, i.e. the middle 80%; 3 parallel coordinates plots of the
# predictors will be displayed

# the parallel coordinates plot then can be used to identify regions in
# predictor space in which the parametric model tends to either under-
# or overpredict the response

# if the output of lm() or glm() is available and the model formula
# consisted only of '+' operations, not ':' etc., then call regdiag();
# otherwise, call regdiagbas()

# arguments:
#    regout:  output of lm() or glm()
#    tail:  what % of the divergences to plot at each tail
#    k, checkna etc.:  see freqparcoord()

regdiag <- function(regout,tail=0.10,k=NULL,m=5,
      checkna = TRUE,cls = NULL,nchunks = length(cls)){
   clregout <- class(regout)
   if(!"lm" %in% clregout)
      stop("The regout argument not of class 'lm'; us regdiagbas()")
   preds <- regout$model[,-1]
   resp <- regout$model[,1]
   parest <- regout$fitted.values
   p <- regdiagbas(preds,resp,parest,tail=tail,k=k,m=m,
      checkna=checkna,cls=cls,nchunks=nchunks)
   if (length(clregout) == 1) {
      p$paramr2 <- summary(regout)$adj.r.squared
      p$nonparamr2 <- cor(p$nonparest,resp)^2
   }
   p
}

regdiagbas <- function(preds,resp,parest, tail=0.10, k=NULL, m=5,
      checkna = TRUE, cls = NULL, nchunks = length(cls)){
   if (is.null(k)) k <- min(200,floor(sqrt(nrow(preds))))
   nonparest <- smoothz(cbind(preds,resp),k,sf=knnreg,
      checkna=checkna,cls=cls,nchunks=nchunks,scalefirst=T)
   divs <- parest - nonparest
   qs <- quantile(divs,probs=c(tail,1-tail))
   x <- cbind(divs,preds)
   grp <- ifelse(x[,1] < qs[1],1,2)
   grp <- ifelse(x[,1] > qs[2],3,grp)
   x <- cbind(x,grp)
   x <- x[grp > 0,]
   nc <- ncol(x)
   p <- freqparcoord(x,m,dispcols=1:(nc-1),grpvar=nc,cls=cls)
   p$nonparest <- nonparest
   p
}

