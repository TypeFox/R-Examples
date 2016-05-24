numch <- function(funcD, para1, para2){

eps <- 1e-06
eps2 <- eps*eps 

f00 <- funcD(para1, para2)  

t01 <- t10 <- t11 <- cbind(para1,para2)

t01[,1] <- t01[,1] + eps 
t10[,2] <- t10[,2] + eps 
t11[,1] <- t11[,1] + eps 
t11[,2] <- t11[,2] + eps 

f01 <- funcD(t01[,1],t01[,2]) 
f10 <- funcD(t10[,1],t10[,2]) 
f11 <- funcD(t11[,1],t11[,2]) 

(f11 - f01 - f10 + f00)/eps2  
 

}
    
  