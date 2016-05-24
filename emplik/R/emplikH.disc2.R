#####################################################
#### 2 sample, discrete hazard type constraint. #####
#### Only handles one constraint.               #####
#####################################################
emplikH.disc2 <- function(x1, d1, y1 = -Inf, x2, d2, y2 = -Inf,
                     theta, fun1, fun2, tola = 1e-6, maxi, mini) {
########Sample One########
    n1 <- length(x1)
    if (n1 <= 2)
        stop("Need more observations")
    if (length(d1) != n1)
        stop("length of x1 and d1 must agree")
    if (any((d1 != 0) & (d1 != 1)))
        stop("d1 must be 0/1's for censor/not-censor")
    if (!is.numeric(x1))
        stop("x1 must be numeric values --- observed times")
        
    newdata1 <- Wdataclean2(z=x1, d=d1)
    temp1 <- DnR(newdata1$value, newdata1$dd, newdata1$weight, y=y1)
    jump1 <- (temp1$n.event)/temp1$n.risk
    k1 <- temp1$n.event - temp1$n.risk
    index1 <- (jump1 < 1)
    k1 <- k1[index1]
    eve1 <- temp1$n.event[index1]
    tm1 <- temp1$times[index1]
    rsk1 <- temp1$n.risk[index1]
    jmp1 <- jump1[index1]
    funtime1 <- fun1(tm1)
    funh1 <- funtime1/rsk1
    
########Sample two########
    n2 <- length(x2)
    if (n2 <= 2)
        stop("Need more observations for sample 2")
    if (length(d2) != n2)
        stop("length of x2 and d2 must agree")
    if (any((d2 != 0) & (d2 != 1)))
        stop("d2 must be 0/1's for censor/not-censor")
    if (!is.numeric(x2))
        stop("x2 must be numeric values -- observed times")

    newdata2 <- Wdataclean2(z=x2, d=d2)
    temp2 <- DnR(newdata2$value, newdata2$dd, newdata2$weight, y=y2)
    jump2 <- (temp2$n.event)/temp2$n.risk
    k2 <- temp2$n.event - temp2$n.risk
    index2 <- (jump2 < 1)
    k2 <- k2[index2]  
    eve2 <- temp2$n.event[index2]
    tm2 <- temp2$times[index2]
    rsk2 <- temp2$n.risk[index2]
    jmp2 <- jump2[index2]
    funtime2 <- fun2(tm2)
    funh2 <- funtime2/rsk2
    
######## constrain function ########
inthaz <- function(x, funt1, evt1, rsk1, funt2, evt2, rsk2, tht) {
 sum(funt1*log(1-(evt1/(rsk1+x*funt1)))) -
 sum(funt2*log(1-(evt2/(rsk2-x*funt2)))) - tht }

m0 <- inthaz(0, funtime1, eve1, rsk1, funtime2, eve2, rsk2, theta)
if (m0 >0) maxi <- 0 
if (m0 <0) mini <- 0
  m1 <- inthaz(mini+tola, funtime1, eve1, rsk1, funtime2, eve2, rsk2, theta)
  m2 <- inthaz(maxi-tola, funtime1, eve1, rsk1, funtime2, eve2, rsk2, theta)
 print(c(m0, m1, m2))
 #  print(m1)
 #  print(m2)
   
 temp <- uniroot(inthaz, c(mini, maxi), tol=tola, funt1=funtime1,
   evt1=eve1, rsk1=rsk1, funt2=funtime2, evt2=eve2, rsk2=rsk2, tht=theta)

 lam <- temp$root
 onePlam <- 1+lam*funh1
 weights1 <- jmp1/onePlam
 oneMlam <- 1-lam*funh2
 weights2 <- jmp2/oneMlam
 
 loglik1 <- sum(eve1*log(onePlam))+sum((-k1)*log((1-jmp1)/(1-weights1)))
 loglik2 <- sum(eve2*log(oneMlam))+sum((-k2)*log((1-jmp2)/(1-weights2)))
 loglikR <- 2*(loglik1+loglik2)
#MZ <- inthaz(lam, funtime1, eve1, rsk1, funtime2, eve2, rsk2, theta)
#print(MZ)

list("-2LLR" = loglikR, lambda = lam, times1 = tm1, times2=tm2,
     wts1 = weights1, wts2 = weights2)
}

