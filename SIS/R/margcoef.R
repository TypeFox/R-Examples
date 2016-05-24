margcoef <- function(x, y, condind = NULL, family, null.model = FALSE, iterind){
n = dim(x)[1]; p = dim(x)[2]; ones = rep(1, n)
candind = setdiff(1:p, condind)
if(iterind == 0){
   if(family == "cox") 
      margcoef = abs(cor(x,y[,1]))
   else
      margcoef = abs(cor(x,y))
}else{
    if(null.model == TRUE){
       if(is.null(condind) == TRUE)  {x = x[sample(1:n),]}
       if(is.null(condind) == FALSE) {x[,candind] = x[sample(1:n),candind]}
    }
    margcoef = abs(sapply(candind, mg, x, y, ones, family, condind)) 
}
return(margcoef)
}
  
mg <- function(index, x=x, y=y, ones=ones, family=family, condind=condind){
margfit = switch(family, gaussian = coef(glm.fit(cbind(ones, x[,index], x[,condind]), y, family=gaussian()))[2],
                         binomial = coef(glm.fit(cbind(ones, x[,index], x[,condind]), y, family=binomial()))[2],
                         poisson = coef(glm.fit(cbind(ones, x[,index], x[,condind]), y, family=poisson()))[2],
                         cox = coef(coxph(y ~ cbind(x[,index], x[,condind])))[1]
                         )                        
}
