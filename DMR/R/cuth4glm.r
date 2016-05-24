cuth4glm <- function(heig, ind, models, y, Data, cont, K, fam){
 sp <- list()
 form <- c('y ~ ')
 if (heig[ind] == max(heig)){
  form <- 'y ~ 1'
 } else {
  if (length(models) != 0){
    for (i in 1:length(models)){
      w <- Data[, i]
      CUT <- cutree(models[[i]], h = heig[ind])
      sp[[i]] <- CUT
      levels(w) <- CUT
      if (sd(w) == 0){
        w <- as.numeric(w)
      }
      names(sp)[i] <- colnames(Data)[i]
      if (sd(w) != 0) form <- paste(form, '+', names(sp)[i])
        Data[, i] <- w
    }
  }
  if (cont > 0){
   form <- paste(form, '+', paste(names(which((heig > heig[ind]) & (names(heig) != 'fac') )), collapse = '+'))      # tu sie na koncu formuly niepotrzebnie dokleja +
  }
 }
 m <- glm(as.formula(form), data = Data, x = TRUE, y = TRUE, family = fam)
 p <- length(m$coeff)
 llik <- logLik(m)
 bic <- -2*llik + p*K
 return(list(model = m, Crit = bic, LogL = llik, SPart = sp))
}
