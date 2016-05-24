personfit.ppar <- function(object) {
# computes Chi-square based itemfit statistics (Smith, p.77ff)
# for object of class "ppar" (from person.parameter)

  excl_obs_num <- object$pers.ex                     #mjm 2014-09-07
  excl_obs_chr <- rownames(object$X)[excl_obs_num]   #
                                                     #
  if(length(excl_obs_num) > 0L){                 # remove obs. to be excluded, but
    X <- object$X[-excl_obs_num,]                # store information to use in
  } else {                                           # subsequent functions
    X <- object$X                                    #
  }                                                  #

  VE <- pifit.internal(object)   # compute expectation and variance term
  Emat <- VE$Emat
  Vmat <- VE$Vmat
  Cmat <- VE$Cmat

  st.res <- (X-Emat)/sqrt(Vmat)
  #st.res <- (X[!TFrow,]-Emat)/sqrt(Vmat)

  sq.res <- st.res^2                            #squared standardized residuals
  pfit <- rowSums(sq.res,na.rm=TRUE)

  pdf <- apply(X, 1L, function(x){ length(na.exclude(x)) })

  #pdf <- apply(X[!TFrow,],1,function(x) {length(na.exclude(x))})   #degress of freedom (#of persons per item)

  p.outfitMSQ <- pfit/pdf

  qsq.outfitMSQ <- (rowSums(Cmat/Vmat^2, na.rm=TRUE)/pdf^2) - 1/pdf
  q.outfitMSQ <- sqrt(qsq.outfitMSQ)

  psumVmat<-rowSums(Vmat)
  p.infitMSQ <- rowSums(sq.res*Vmat, na.rm = TRUE)/psumVmat

  qsq.infitMSQ <- rowSums(Cmat-Vmat^2, na.rm=TRUE)/psumVmat^2
  q.infitMSQ <- sqrt(qsq.infitMSQ)

  p.outfitZ <- ((p.outfitMSQ)^(1/3)-1)*(3/q.outfitMSQ)+(q.outfitMSQ/3)
  p.infitZ <- ((p.infitMSQ)^(1/3)-1)*(3/q.infitMSQ)+(q.infitMSQ/3)

  result <- structure(
    list("p.fit"        = pfit,
         "p.df"         = pdf,
         "st.res"       = st.res,
         "p.outfitMSQ"  = p.outfitMSQ,
         "p.infitMSQ"   = p.infitMSQ,
         "p.outfitZ"    = p.outfitZ,
         "p.infitZ"     = p.infitZ,
         "excl_obs_num" = excl_obs_num,
         "excl_obs_chr" = excl_obs_chr),
    class = "pfit")
  return(result)

}
