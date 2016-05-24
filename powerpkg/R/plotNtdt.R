"plotNtdt" <-
function(g,m,alpha=0.00000005,power=0.80,...)
{
# Duplicate Figures 1 and 2 from Abel and Muller-Myhsok (1998)
# Am J Hum Genet 63:664-667
# plot.ntdt(g=2,m=0.10) gives Figure 1A
# plot.ntdt(g=2,m=0.50) gives Figure 1B
tb <- ntdt.q(g=g,m=m,alpha=alpha,power=power)
plot(tb$q,tb$log.dmax.50,type="l",xlab="q",ylab="log10(N)",lty=4,
ylim=c(min(tb$log.dmax),max(tb$log.dmax.50)),
main=paste("g=",g,"m=",m,"alpha=",alpha,"power=",power),...)
lines(tb$q,tb$log.dmax.75,lty=2)
lines(tb$q,tb$log.dmax,lty=1)
legend(mean(tb$q),(max(tb$log.dmax.50)-0.05),c("0.50 dmax","0.75 dmax","dmax"),lty=c(4,2,1))
}

