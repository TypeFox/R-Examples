# library(COUNT)
# data(medpar)
# mdpar <- glm.nb(los ~ hmo+white+type2+type3, data=medpar, y=TRUE, model=TRUE)

nb2.obs.pred <- function(len, model)  {
  mu <- fitted(model)
  trun.y <- model$y[model$y < len+1]
  obs <- as.data.frame(table(trun.y) / length(model$y) * 100)
  names(obs) <- c("Count","propObsv")
  alpha <- 1/model$theta
  amu <- alpha*mu
  pred <- data.frame(Count = 0:len,
                     propPred = sapply(0:len, function(i)
                       mean(exp(i*log(amu/(1+amu)) - (1/alpha) * log(1+amu)  + 
                                log(gamma(i + 1/alpha) ) - log(gamma(i + 1) ) - 
                                log(gamma(1 / alpha))))) * 100)
  out <- merge(pred, obs, all=TRUE)
  out$propObsv[is.na(out$propObsv)] <- 0
  out$Diff <- with(out, propObsv - propPred)
  return(out[,c(1,3,2,4)])
}

# nb2.obs.pred(len=25, model=mdpar)








