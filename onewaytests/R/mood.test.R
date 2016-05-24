mood.test <- function(y, group, na.rm = TRUE) {
  
if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"
  METHOD <- "Mood's Median Test"

  y.median <-Fisher<-Chisquare<-Chisquare_pvalue<-Chisquare_df<-NULL
  
  y.median=median(y)
 
    approx<-chisq.test(y < y.median, group)$statistic
    df<-chisq.test(y < y.median, group)$parameter
    p.value<-chisq.test(y < y.median, group)$p.value

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
