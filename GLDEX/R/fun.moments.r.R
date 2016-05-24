"fun.moments.r" <-
function(x,normalise="N")
{
x<-na.omit(x)
a1 <- mean(x)
a2 <- sum((x - a1)^2)/length(x)
a3 <- mean((x - mean(x))^3)/(a2^(3/2))
a4 <- mean((x - mean(x))^4)/(a2^(4/2))

  if (normalise == "Y") {
        a4 <- a4 - 3
    }

return(c("mean"=a1, "variance"=a2, "skewness"=a3, 
"kurtosis"=a4))
}

