# power for simple logistic regression based on
# http://personal.health.usf.edu/ywu/logistic.pdf
# Hsieh et al. STATISTICS IN MEDICINE
# Statist. Med. 17, 1623-1634 (1998)
# 

# created on Sept. 24, 2012

# OR = exp(beta.star) odds ratio
SSizeLogisticCon<-function(p1, OR, alpha=0.05, power=0.8)
{
   beta.star<-log(OR)
   za<-qnorm(1-alpha/2)
   zb<-qnorm(power)
   n<-(za+zb)^2/(p1*(1-p1)*beta.star^2)
   n.int<-ceiling(n)

   return(n.int)
}


# OR = exp(beta.star) odds ratio
powerLogisticCon<-function(n, p1, OR, alpha=0.05)
{
   beta.star=log(OR)
   za<-qnorm(1-alpha/2)
   power<-1-pnorm(za-sqrt(n*beta.star^2*p1*(1-p1)))
   return(power)
  
}

###########################
# logistic regression logit(p) = beta0+ beta1*X
# B=pr(X=1)
# p1 = Pr(D|X=0) # event rate at X=0
# p2 = Pr(D|X=1) # event rate at X=1
# alpha - type I error rate
SSizeLogisticBin<-function(p1, p2, B, alpha=0.05, power=0.8)
{
   za<-qnorm(1-alpha/2)
   zb<-qnorm(power)

   p=(1-B)*p1+B*p2
   part1 = za*sqrt(p*(1-p)/B) 
   part2 = zb*sqrt( p1*(1-p1)+p2*(1-p2)*(1-B)/B)
   part3 = (p1-p2)^2*(1-B)
   n<-(part1+part2)^2/part3
   n.int<-ceiling(n)

   return(n.int)
}

###########################
# logistic regression logit(p) = beta0+ beta1*X
# B=pr(X=1)
# p1 = Pr(D|X=0) # event rate at X=0
# p2 = Pr(D|X=1) # event rate at X=1
# alpha - type I error rate
powerLogisticBin<-function(n, p1, p2, B, alpha=0.05)
{
   za<-qnorm(1-alpha/2)

   p=(1-B)*p1+B*p2
   a = za*sqrt(p*(1-p)/B) 
   b = sqrt( p1*(1-p1)+p2*(1-p2)*(1-B)/B)
   myc = (p1-p2)^2*(1-B)
   power = pnorm( (sqrt(n*myc) - a)/b )

   return(power)
}




###########################


### Example in Table II Design (Balanced design (1)) of Hsieh et al. (1998 )
### the sample size is 317
#cat("Sample size>>\n")
#print(SSizeLogisticCon(p1=0.5, OR=exp(0.405), alpha=0.05, power=0.95))
#cat("\n")
#
#cat("Power>>\n")
#print(powerLogisticCon(n=317, p1=0.5, OR=exp(0.405), alpha=0.05))
#cat("\n")
#
# n=1281
# SSizeLogisticBin(p1=0.4, p2=0.5, B=0.5, alpha=0.05, power=0.95)
# power=0.95
# powerLogisticBin(n=1281, p1=0.4, p2=0.5, B=0.5, alpha=0.05)

##
