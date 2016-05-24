kw.test <- function(y, group, na.rm = TRUE){

 if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"
  METHOD <- "Kruskal-Wallis Rank Sum Test"

  n.i<- Rmean.i<- KW.noties<- ties<-correction<- KW.stat<- p.value<-NULL
  ranks = rank(y)
  N<-length(y)
  x.levels <- levels(factor(group))
  t<-as.data.frame(table(rank(y)))
  ties<-t$Freq
  for (i in x.levels) {
    n.i[i] <- length(y[group==i])
    Rmean.i[i]<-mean(ranks[group==i])
  }
  

  KW.noties = 12/(N*(N+1)) * sum( n.i*(Rmean.i - (N+1)/2)^2 )

  correction <- (1 - sum( ties^3 - ties )/(N^3 - N) )
  approx<- KW.noties/correction
  p.value<-pchisq(approx, df=(length(x.levels)-1),lower.tail = FALSE)
  df=(length(x.levels)-1)

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