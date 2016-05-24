numgh <- function(funcD, para){

eps <- 1e-06
eps2 <- eps*eps 

para1 <- para - eps/2
para2 <- para + eps/2 
f1 <- funcD(para1)
f2 <- funcD(para2)

fi <- (f2 - f1)/ eps   

t01 <- t10 <- t11 <- para + eps
t11 <- t11 + eps 

f00 <- funcD(para) 
f01 <- funcD(t01) 
f10 <- funcD(t10) 
f11 <- funcD(t11) 

se <- (f11 - f01 - f10 + f00)/eps2 


list(fi = fi, se = se)


}
    
  