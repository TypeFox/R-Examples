
welch.test <- function(y, group, rate = 0, na.rm = TRUE) {


 if (na.rm){
    completeObs <- complete.cases(y, group)
    y <- y[completeObs]
    group <- group[completeObs]
  }
  df <- data.frame(Response = y, Group = group)
  DNAME <- "y vs group"

if (rate==0){METHOD <- "Welch's Heteroscedastic F Test"
}else{METHOD <- "Welch's Heteroscedastic F Test with Trimmed Mean and Winsorized Variance"}
  


trim=function(x,rate){
n=length(x)
xx=sort(x)
lambda=round(n*rate)
xx[(lambda+1):(n-lambda)]
}

wins=function(x,rate){
n=length(x)
xx=sort(x)
lambda=round(n*rate)
xxx=c(rep(xx[lambda+1],lambda),xx[(lambda+1):(n-lambda)],rep(xx[n-lambda],lambda))
xxx
}

n <- length(y)
x.levels <- levels(factor(group))
y.vars <- y.means <-lambda<- m <- y.n <- w <- b <- q <-NULL

for (i in x.levels) {

y.n[i] <- length(y[group==i])

lambda[i]=round(y.n[i]*rate)

b[i]=y.n[i]-2*lambda[i]

y.vars[i] <- var(wins(y[group==i],rate))

q[i]=(y.n[i]-1)*y.vars[i]/b[i]/(b[i]-1)

w[i] <- 1/q[i]

y.means[i] <- mean(trim(y[group==i],rate))
  
}

U=sum(w)

w_y=sum(w*y.means)/U

J=length(x.levels)

A=sum(w*(y.means-w_y)^2)/(J-1)

B=2*(J-2)/(J^2-1)*sum((1-w/U)^2/(b-1))



Ftest=A/(B+1)

df1=J-1
df2=(3/(J^2-1)*sum((1-w/U)^2/(b-1)))^(-1)


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

