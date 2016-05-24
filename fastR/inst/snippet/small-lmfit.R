y <- c(0,2,0,2,-1,6)
x1 <- c(1,1,2,2,3,3); x2 <- c(0,1,1,2,1,3)
v0 <- rep(1,length(y))
v1 <- x1 - mean(x1); v2=x2 - mean(x2)
w1 <- v1 - project(v1,v2)
w2 <- v2 - project(v2,v1)
#
# obtaining model fits by projection
#
p0 <- project(y,v0,type='vec'); p0
p1 <- project(y,v1,type='vec'); p1
p2 <- project(y,v2,type='vec'); p2
q1 <- project(y,w1,type='vec'); q1
q2 <- project(y,w2,type='vec'); q2
#
# this won't be a correct fit because dot(v1,v2) != 0
#
p0 + p1 + p2  
#
# here is the correct fit 
#
p0 + q1 + p2
p0 + p1 + q2
#
# we can compare the results with those from lm()
#
model <- lm(y~x1+x2); fitted(model)        
#
# this won't work to get the coefficients:
#
b1.wrong <- (p1/v1); b1.wrong   
b2.wrong <- (p2/v2); b2.wrong   
#
# now let's get the coefficients correctly:
#
b1 <- (q1/w1); b1   
b2 <- (q2/w2); b2  
a0 <- (p0/v0); a0
b0 <- a0 - b1*mean(x1) - b2*mean(x2); b0
coef(model)
