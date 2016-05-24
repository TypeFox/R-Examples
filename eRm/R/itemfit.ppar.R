`itemfit.ppar` <-
function(object)
# computes Chi-square based itemfit statistics
# for object of class "ppar" (from person.parameter)
{
  if (length(object$pers.ex)==0) {
    X <- object$X
  } else {
    X <- object$X[-object$pers.ex,]
  }

  VE <- pifit.internal(object)                  #compute expectation and variance term
  Emat <- VE$Emat
  Vmat <- VE$Vmat
  Cmat <- VE$Cmat

  st.res <- (X-Emat)/sqrt(Vmat)
  sq.res <- st.res^2                            #squared standardized residuals
  ifit <- colSums(sq.res,na.rm=TRUE)

  idf <- apply(X,2,function(x) {length(na.exclude(x))})

  i.outfitMSQ <- ifit/idf

  qsq.outfitMSQ <- (colSums(Cmat/Vmat^2, na.rm=TRUE)/idf^2) - 1/idf
  q.outfitMSQ <- sqrt(qsq.outfitMSQ)

  isumVmat<-colSums(Vmat)
  i.infitMSQ <- colSums(sq.res*Vmat, na.rm = TRUE)/isumVmat

  qsq.infitMSQ <- colSums(Cmat-Vmat^2, na.rm=TRUE)/isumVmat^2
  q.infitMSQ <- sqrt(qsq.infitMSQ)

  i.outfitZ <- (i.outfitMSQ^(1/3) - 1)*(3/q.outfitMSQ)+(q.outfitMSQ/3) # corr. rh 2011-06-15
  i.infitZ  <- (i.infitMSQ^(1/3)  - 1)*(3/q.infitMSQ) +(q.infitMSQ/3)  # hint from rainer alexandrowicz

  result <- list(i.fit=ifit,i.df=idf,st.res=st.res,i.outfitMSQ=i.outfitMSQ,i.infitMSQ=i.infitMSQ,i.outfitZ=i.outfitZ,i.infitZ=i.infitZ)

  class(result) <- "ifit"
  result
}
