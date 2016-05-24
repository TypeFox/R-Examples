aov.test <- function(y, group, na.rm = TRUE) {

 if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"
  METHOD <- "One-Way Analysis of Variance"


n <- length(y)
x.levels <- levels(factor(group))
y.sums <- y.n <- NULL

sst=sum(y^2)-(sum(y)^2)/n


for (i in x.levels) {

y.sums[i] <- sum(y[group==i])
  
y.n[i] <- length(y[group==i])

}


ssb<- sum(y.sums^2/y.n)-sum(y)^2/n

ssw=sst-ssb 


df1=length(x.levels)-1
df2=n-length(x.levels)


Ftest=ssb/df1*df2/ssw



p.value=pf(Ftest,df1,df2,lower.tail = F)

 names(Ftest) <- "F"
  PARAMETER <- c(df1, df2)
  names(PARAMETER) <- c("num df", "denom df")
  
 if(sum(!completeObs) > 0){
if(sum(!completeObs)==1){cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations was removed due to missingness."), "\n")
}else  {cat("\n", paste("NOTE: ", sum(!completeObs), " of ", sum(!completeObs)+length(y), "observations were removed due to missingness."), "\n")}
}
  
  structure(list(statistic = Ftest, parameter = PARAMETER, 
                 p.value = p.value, method = METHOD, data.name = DNAME), class = "htest")
}


