`power.equivalence.md.plot` <-
function (alpha,logscale,theta1,theta2,sigma,n,nu,title2) 
{
# plot power curves for TOST
# a vector n of sample sizes and corresponding degrees of freedom, nu, is input
# power is calculated at 201 values between theta1 and theta2

# parameters on the ordinary scale are:
# alpha:       alpha level for the 2 t-tests (usually 0.05).  Confidence interval for full test is at level 1 - 2 alpha 
# logscale:    FALSE=use regular scale; TRUE=use logarithmic scale
# theta1:      lower limit of equivalence interval
# theta2:      upper limit of equivalence interval
# sigma:       root MSE from ANOVA 
# n:           vector of number of observations
# nu:          vector of degrees of freedom for sigma
# title2:      subtitle appearing at bottom of graph

nn<-length(n)
 
power_array<-matrix(c(1:(201*nn)),nrow=nn)  #  output matrix of power values with differences appended
xdiff<-c(1:201)
delta<-(theta2-theta1)/200

for (inn in 1:nn)
{
 diff<-theta1 - delta

  for (idiff in 1:201)
   {
   diff<-diff+delta
   xdiff[idiff] <-diff
   power_array[inn,idiff]<-power.equivalence.md(alpha,logscale,theta1,theta2,diff,sigma,n[inn],nu[inn])[1]
   }
}

if (logscale){
  plot(xdiff,power_array[nn,],type='l',
  xlab='Ratio of Test to Reference',
  ylab='Power',
  main='Power of TOST',
  sub=title2,
  xaxp=c(theta1,theta2,9),
  ylim=c(0,1),
  xlim=c(theta1,theta2))}   

else{
  plot(xdiff,power_array[nn,],type='l',
  xlab='Difference of Test and Reference',
  ylab='Power',
  main='Power of TOST',
  sub=title2,
  xaxp=c(theta1,theta2,10),
  ylim=c(0,1),
  xlim=c(theta1,theta2))}   

for (inn in 1:nn)
{ lines(xdiff,power_array[inn,]) }

abline(h=.9)
abline(h=.8)
abline(h=.7)
abline(h=.6)
abline(h=.5)
abline(h=.4)
abline(h=.3)
abline(h=.2)
abline(h=.1)
abline(h=.05)
par(ps=7)

if (logscale){
  abline(v=sqrt(theta2*theta1))
  xdel <-sqrt(theta2*theta1)*.01
  xleft<-theta1+10*xdel
  text(xleft,.02,paste('Equivalence:( ',theta1,',',theta2,' )','    Sigma=',sigma))

  y<-c(1:nn)
  x<-0*c(1:nn) +  (theta1/theta2)^.25
  for (inn in 1:nn){
    y[inn]<-power.equivalence.md(alpha,logscale,theta1,theta2,x[inn],sigma,n[inn],nu[inn])[1]  }
  y<-y 
  text(x,y,n)
}

else{
  xdel <- (theta2-theta1)*.01
  xleft<-theta1+20*xdel
  text(xleft,.02,paste('Equivalence:( ',theta1,',',theta2,' )','    Sigma=',sigma))

  y<-c(1:nn)
  x<-0*c(1:nn) + theta1 + .25*(theta2-theta1)
  for (inn in 1:nn){
    y[inn]<-power.equivalence.md(alpha,logscale,theta1,theta2,x[inn],sigma,n[inn],nu[inn])[1]  }
  y<-y 
  text(x,y,n)
}
p_return<- cbind(xdiff,t(power_array)) 
return(p_return)}
