ag.test <- function(y, group, na.rm = TRUE) {

 if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"
  METHOD <- "Alexander-Govern Test"


  n <- length(y)
  x.levels <- levels(factor(group))
  y.se <- w <- y.means <- y.n <- NULL
  for (i in x.levels) {
   
    y.se[i] <- sqrt(var(y[group==i])/length(y[group==i]))
 
    y.means[i] <- mean(y[group==i])
 
    y.n[i] <- length(y[group==i])
  }

  w <- (1/y.se^2)/sum(1/y.se^2)


  y.plus <- sum(w * y.means)


  t <- (y.means - y.plus)/y.se


  a <- y.n - 1.5
  b <- 48*a^2
  c <- (a * log(1 + t^2/(y.n-1)))^(.5)
  z <- c + (c^3+3*c)/b - (4*c^7+33*c^5+240*c^3+855*c)/(10*b^2+8*b*c^4+1000*b)


  approx <- sum(z^2)
  df <- length(x.levels)-1
  p.value <- pchisq(approx, df=df, lower.tail=FALSE)

  names(approx) <- "X-squared"
  PARAMETER <- c(df)
  names(PARAMETER) <- c("df")
  
  if(sum(!completeObs) > 0){
if(sum(!completeObs)==1){cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations was removed due to missingness."), "\n")
}else  {cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations were removed due to missingness."), "\n")}
}
  
  structure(list(statistic = approx, parameter = PARAMETER, 
                 p.value = p.value, method = METHOD, data.name = DNAME), class = "htest")



}
