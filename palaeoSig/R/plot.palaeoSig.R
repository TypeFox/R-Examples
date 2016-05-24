#plot.palaeoSig <-
#function(x, names, pos=2,...){
#  if(missing(names))names=names(x$EX)
#  with(x,{
#    hist(sim.ex, breaks=seq(min(sim.ex),max(sim.ex),length=20),xlim=c(0,MAX[1]*1.1), main="", xlab="Proportion variance explained", col="grey80", border=NA,...)
#    abline(v=MAX, col=1, lwd=2, lty=3)
#    abline(v=EX, col=1)     
#    text(EX,par()$usr[4]*.9,label=names, srt=90, pos=pos)
#  })
#}

plot.palaeoSig<-function(x, vnames, top=0.7, adj=c(0,0.5),p.val=0.95,...){
  if(missing(vnames))vnames=names(x$EX)
  with(x,{
    hist(sim.ex, breaks=seq(min(sim.ex),max(sim.ex),length=20),xlim=c(0,MAX[1]*1.1), main="", xlab="Proportion variance explained", col="grey80", border=NA,...)
    tops<-par()$usr[4]*top
    sapply(EX,function(z)lines(rep(z,2),c(0,tops)))
    #abline(v=MAX, col=1, lwd=2, lty=3)
    lines(rep(MAX,2),c(0,tops), col=1, lwd=2, lty=3)
    lines(rep(quantile(sim.ex,p.val),2),c(0,tops), col=2, lwd=1, lty=3)
    #abline(v=EX, col=1)     
    putEX<-spread.labs(EX,1.2*strwidth('A', cex=.8))

    text(putEX,par()$usr[4]*.71,label=vnames, srt=90, adj=adj, cex=.8)
  })
}

