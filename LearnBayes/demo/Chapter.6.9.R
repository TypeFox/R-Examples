###################################################
# Section 6.9 Modeling Data with Cauchy Errors
###################################################

library(LearnBayes)

 data(darwin)
 attach(darwin)
 mean(difference)

 log(sd(difference))

 laplace(cauchyerrorpost,c(21.6,3.6),difference)

 laplace(cauchyerrorpost,.1*c(21.6,3.6),difference)$mode
 
 c(24.7-4*sqrt(34.96),24.7+4*sqrt(34.96))
 c(2.77-4*sqrt(.138),2.77+4*sqrt(.138))

 mycontour(cauchyerrorpost,c(-10,60,1,4.5),difference,
   xlab="mu",ylab="log sigma")

S=readline(prompt="Type  <Return>   to continue : ")

 fitlaplace=laplace(cauchyerrorpost,c(21.6,3.6), difference)
 windows()
 mycontour(lbinorm,c(-10,60,1,4.5),list(m=fitlaplace$mode,
   v=fitlaplace$var), xlab="mu",ylab="log sigma")

 proposal=list(var=fitlaplace$var,scale=2.5)
 start=c(20,3)
 m=1000
 s=rwmetrop(cauchyerrorpost,proposal,start,m,difference)

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 mycontour(cauchyerrorpost,c(-10,60,1,4.5),difference,
   xlab="mu",ylab="log sigma")
 points(s$par[,1],s$par[,2])

 fitgrid=simcontour(cauchyerrorpost,c(-10,60,1,4.5),difference,
  50000)
 proposal=list(var=fitlaplace$var,scale=2.5)
 start=c(20,3)
 fitrw=rwmetrop(cauchyerrorpost,proposal,start,50000,
   difference)
 proposal2=list(var=fitlaplace$var,mu=t(fitlaplace$mode))
 fitindep=indepmetrop(cauchyerrorpost,proposal2,start,50000,
  difference)
 fitgibbs=gibbs(cauchyerrorpost,start,50000,c(12,.75),
  difference)

 apply(fitrw$par,2,mean)

 apply(fitrw$par,2,sd)


