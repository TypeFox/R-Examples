onecor.wtd <- function(x, y, weight=NULL){
  if(sum((!is.na(x))*(!is.na(y)))>0){
    if(is.null(weight)){
      weight <- rep(1, length(x))
    }
    use <- !is.na(y) & !is.na(x)
    x <- x[use]
    y <- y[use]
    weight <- weight[use]
    #r1 <- lm(stdz(y, weight=weight)~stdz(x, weight=weight), weight=weight)
    #corcoef <- coef(summary(r1))[2,]
    corcoef <- coef(summary(lm(stdz(y, weight=weight)~stdz(x, weight=weight), weights=weight)))[2,]
  }
  else
    corcoef <- rep(NA, 4)
  names(corcoef) <- rep("", 4)
  corcoef
}
