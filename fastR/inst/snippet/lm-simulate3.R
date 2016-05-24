sim <- function(b0=3,b1=5,lambda=1,
                x=rep(1:5,each=4)    # 4 observations at each of 5 values
    ){
    # shift to give a mean of 0.
    e <- rexp(length(x),rate=1/lambda) - lambda
    y <- b0 + b1*x + e
    model <- lm(y~x)  
    ci <- confint(model,2)
    return(b1 > ci[1] && b1 < ci[2])
}
t <- table(replicate(10000,sim())); t / 10000
prop.test(t[2],sum(t),p=0.95)
