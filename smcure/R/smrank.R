smrank <-
function(beta,Time,X,n,w,Status){
     error <- drop(log(Time)-beta%*%t(X))
     tp <- numeric()     
     for(i in 1:n){
        tp[i] <- sum(as.numeric((error[i]-error)<0)*abs(error[i]-error)*w*Status[i])
        }
    sum(tp)/n
}
