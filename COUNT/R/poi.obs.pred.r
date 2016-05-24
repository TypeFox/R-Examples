# Create a table of observed vs predicted Poisson counts, and difference
#  following glm()  :  see poi.obs.pred.Rd for usage
# See Hilbe, J.M (2011), Negative Binomial Regression, 2nd ed, Cambridge Univ Press


poi.obs.pred <- function(len, model)  {
  mu <- fitted(model)
  trun.y <- model$y[model$y < len+1]
  obs <- as.data.frame(table(trun.y) / length(model$y) * 100)
  names(obs) <- c("Count","propObsv")
  pred <- data.frame(Count = 0:len,
                     propPred = sapply(0:len, function(i)
                       mean(exp(-mu)*(mu^i)/factorial(i))) * 100)
  out <- merge(pred, obs, all=TRUE)
  out$propObsv[is.na(out$propObsv)] <- 0
  out$Diff <- with(out, propObsv - propPred[1:(len+1)])
  return(out[,c(1,3,2,4)])
}

#library(COUNT)
#data(medpar)
#
#mdpar <- glm(los ~ hmo+white+type2+type3, family=poisson, data=medpar,
#             y=TRUE, model=TRUE)
#poi.obs.pred(len=25, model=mdpar)
#
#data(affairs)
#myglm <- glm(naffairs ~ kids, family=poisson, data=affairs, y=TRUE, model=TRUE)
#poi.obs.pred(len=8, model=myglm)
