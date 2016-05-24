#-------------------------- normal distribution--------------------------------#

shade.norm<-function(x=NULL,from=NULL,to=NULL,sigma=1,mu=0,tail="lower",show.p=TRUE,show.d=FALSE,show.dist=TRUE,digits=5,legend.cex=.9,shade.col="gray",...){
xv<-seq(mu-4*sigma,mu+4*sigma,sigma/1000)
yv<-dnorm(xv,mean=mu,sd=sigma)
curve(dnorm(x,mu,sigma),from=mu-4*sigma,to=mu+4*sigma,ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")),xlab = expression(italic(x)),...)

if(tail=="lower"){
polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==mu-4*sigma]),col=shade.col)
p<-round(pnorm(x,mu,sigma,lower.tail=TRUE),digits)
d<-round(dnorm(x,mean=mu,sd=sigma),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="upper"){
polygon(c(x,xv[xv>=x]),c(yv[xv==mu+4*sigma],yv[xv>=x]),col=shade.col)
p<-round(pnorm(x,mu,sigma,lower.tail=FALSE),digits)
d<-round(dnorm(x,mean=mu,sd=sigma),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="two"){
polygon(c(xv[xv<=-abs(x)],-abs(x)),c(yv[xv<=-abs(x)],yv[xv==mu-4*sigma]),col=shade.col)
polygon(c(abs(x),xv[xv>=abs(x)]),c(yv[xv==mu+4*sigma],yv[xv>=abs(x)]),col=shade.col)
p<-round(2*pnorm(abs(x),mu,sigma,lower.tail=FALSE),digits)
if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(2%*%italic(P),"(", italic(X)>="|",.(x),"|) = ",.(p), sep = "")))
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="two.custom"){
polygon(c(xv[xv<=x[1]],x[1]),c(yv[xv<=x[1]],yv[xv==mu-4*sigma]),col=shade.col)
polygon(c(x[2],xv[xv>=x[2]]),c(yv[xv==mu+4*sigma],yv[xv>=x[2]]),col=shade.col)
p<-round(pnorm(x[1],mu,sigma),digits) + round(pnorm(x[2],mu,sigma, lower.tail = FALSE),digits)
if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(x[1]) <= "", italic(X)>=.(x[2]),") = ",.(p), sep = "")))
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="middle"){
polygon(c(xv[xv<=mu+4*sigma],mu+4*sigma),c(yv[xv<=mu+4*sigma],yv[xv==mu-4*sigma]),col=shade.col)
polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==mu-4*sigma]),col="white")
polygon(c(to,xv[xv>=to]),c(yv[xv==mu+4*sigma],yv[xv>=to]),col="white")
p<-round(pnorm(to,mu,sigma)-pnorm(from,mu,sigma),digits)
if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(N),"(",.(mu)," , ",.(sigma^2),")", sep = "")),bty="n",cex=legend.cex)
}
}

#-------------------------- t-distribution --------------------------------#

shade.t<-function(x=NULL,from=NULL,to=NULL,nu=3,tail="lower",show.p=TRUE,show.d=FALSE,show.dist=TRUE,digits=5,legend.cex=.9,shade.col="gray",...){
sigma<-qt(.975,nu)
xv<-seq(-4*sigma,4*sigma,sigma/1000)
yv<-dt(xv,df=nu)
curve(dt(x,nu),from=-4*sigma,to=4*sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")),...)

if(tail=="lower"){
polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==-4*sigma]),col=shade.col)
p<-round(pt(x,nu,lower.tail=TRUE),digits)
d<-round(dt(x,nu),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="upper"){
polygon(c(x,xv[xv>=x]),c(yv[xv==4*sigma],yv[xv>=x]),col=shade.col)
p<-round(pt(x,nu,lower.tail=FALSE),digits)
d<-round(dt(x,nu),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="two"){
polygon(c(xv[xv<=-abs(x)],-abs(x)),c(yv[xv<=-abs(x)],yv[xv==4*sigma]),col=shade.col)
polygon(c(abs(x),xv[xv>=abs(x)]),c(yv[xv==4*sigma],yv[xv>=abs(x)]),col=shade.col)
p<-round(2*pt(abs(x),nu,lower.tail=FALSE),digits)
if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(2%*%italic(P),"(", italic(X)>="|",.(x),"|) = ",.(p), sep = "")))
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
}

if(tail=="middle"){
polygon(c(xv[xv<=4*sigma],4*sigma),c(yv[xv<=4*sigma],yv[xv==-4*sigma]),col=shade.col)
polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==-4*sigma]),col="white")
polygon(c(to,xv[xv>=to]),c(yv[xv==4*sigma],yv[xv>=to]),col="white")
p<-round(pt(to,nu)-pt(from,nu),digits)
if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X)," ~ ",italic(t),"(",.(nu),")", sep = "")),bty="n",cex=legend.cex)
}
}

#------------------------F distribution---------------------------#

shade.F<-function(x=NULL,from=NULL,to=NULL,nu1=1,nu2=5,tail="lower",show.p=TRUE,show.d=FALSE,show.dist=TRUE,prob.to.each.tail=0.025,digits=5,legend.cex=.9,shade.col="gray",...){
sigma<-qf(.999,nu1,nu1)
xv<-seq(0,sigma,sigma/1000)
yv<-df(xv,nu1,nu2)
curve(df(x,nu1,nu2),from=0,to=sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")), ...)

if(tail=="lower"){
polygon(c(xv[xv<=x],x),c(yv[xv<=x],yv[xv==0]),col=shade.col)
p<-round(pf(x,nu1,nu2,lower.tail=TRUE),digits)
d<-round(df(x,nu1,nu2),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
}
if(tail=="upper"){
polygon(c(x,xv[xv>=x]),c(yv[xv==sigma],yv[xv>=x]),col=shade.col)
p<-round(pf(x,nu1,nu2,lower.tail=FALSE),digits)
d<-round(df(x,nu1,nu2),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
}
if(tail=="two"){
qs<-round(qf(c(prob.to.each.tail,1-prob.to.each.tail),nu1,nu2),1) 
polygon(c(xv[xv<=qs[1]],qs[1]),c(yv[xv<=qs[1]],yv[xv==sigma]),col=shade.col)
polygon(c(qs[2],xv[xv>=qs[2]]),c(yv[xv==sigma],yv[xv>=qs[2]]),col=shade.col)
p<-round(prob.to.each.tail*2,digits)
if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(qs[1])<=italic(X), " and ", italic(X)>=.(qs[2]),") = ",.(p))))
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
}
if(tail=="middle"){
polygon(c(xv[xv<=sigma],sigma),c(yv[xv<=sigma],yv[xv==0]),col=shade.col)
polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==0]),col="white")
polygon(c(to,xv[xv>=to]),c(yv[xv==sigma],yv[xv>=to]),col="white")
p<-round(pf(to,nu1,nu2)-pf(from,nu1,nu2),digits)
if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="", italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",italic(F),"(",.(nu1),",",.(nu2),")", sep = "")),bty="n",cex=legend.cex,adj=0.5)
}
}

#-------------------------Chi-square distribution------------------------------#

shade.chi<-function(x=NULL,from=NULL,to=NULL,nu=1,tail="lower",show.p=TRUE,show.d=FALSE,show.dist=TRUE,prob.to.each.tail=0.025,digits=5,legend.cex=.9,shade.col="gray",...){
sigma<-qchisq(.9999,nu)
xv<-seq(0,sigma,sigma/10000)
yv<-dchisq(xv,nu)
curve(dchisq(x,nu),from=0,to=sigma,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep ="")), ...)

if(tail=="lower"){
    if(nu<3){
        ax <- abs(xv-x)
        xn <- ax==min(ax)
       polygon(c(0, x, x, sort(xv[xv <= x], decreasing = T)), c(0, 0, yv[xn], sort(yv[xv <= x])), col = shade.col) 
     }
    if(nu>=3)polygon(c(xv[xv <= x], x), c(yv[xv <= x], yv[xv == 0]), col = shade.col)

d<-round(dchisq(x,nu),digits)
p<-round(pchisq(x,nu,lower.tail=TRUE), digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
}

if(tail=="upper"){
polygon(c(x,xv[xv>=x]),c(yv[xv==sigma],yv[xv>=x]),col=shade.col)
p<-round(pchisq(x,nu,lower.tail=FALSE),digits)
d<-round(dchisq(x,nu),digits)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==FALSE)legend("topright",legend=bquote(paste("",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.d==TRUE&show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(p),",  ",italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
}
if(tail=="two"){
qs<-round(qchisq(c(prob.to.each.tail,1-prob.to.each.tail),nu),1) 
if(nu >=3){
polygon(c(xv[xv<=qs[1]],qs[1]),c(yv[xv<=qs[1]],yv[xv==sigma]),col=shade.col)
polygon(c(qs[2],xv[xv>=qs[2]]),c(yv[xv==sigma],yv[xv>=qs[2]]),col=shade.col)
}
p<-round(prob.to.each.tail*2,digits)
if(show.p==TRUE)legend("topright",bty="n",cex=legend.cex,legend=bquote(paste(italic(P),"(",.(qs[1])<=italic(X), " and ", italic(X)>=.(qs[2]),") = ",.(p))))
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
}
if(tail=="middle"){
if(nu >= 3){
polygon(c(xv[xv<=sigma],sigma),c(yv[xv<=sigma],yv[xv==0]),col=shade.col)
polygon(c(xv[xv<=from],from),c(yv[xv<=from],yv[xv==0]),col="white")
polygon(c(to,xv[xv>=to]),c(yv[xv==sigma],yv[xv>=to]),col="white")
}
p<-round(pchisq(to,nu)-pchisq(from,nu),digits)
if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="",italic(X)<=.(to),") = ",.(p), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("top",legend=bquote(paste(italic(X)," ~ ",chi^2,"(",.(nu),")")),bty="n",cex=legend.cex,adj=0.5)
}
}

#------------------------------ Weibull distribution ----------------------------#

shade.wei <- function (x = NULL, from = NULL, to = NULL, theta = 1, beta = 1, tail = "lower", 
    show.p = TRUE, show.d = FALSE, show.dist = TRUE, prob.to.each.tail = 0.025, 
    digits = 5, legend.cex = 0.9, shade.col = "gray", ...) 
{
    sigma <- qweibull(0.9999, theta, beta)
    xv <- seq(0, sigma, sigma/10000)
    yv <- dweibull(xv, theta, beta)
    curve(dweibull(x, theta, beta), from = 0, to = sigma, xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")), ...)
    if (tail == "lower") {
         if(theta<=1){
            ax <- abs(xv-x)
            xn <- ax==min(ax)
           polygon(c(0, x, x, sort(xv[xv <= x], decreasing = T)), c(0, 0, yv[xn], sort(yv[xv <= x])), col = shade.col) 
         }
        if(theta > 1)polygon(c(xv[xv <= x], x), c(yv[xv <= x], yv[xv == 0]), col = shade.col)
        d <- round(dweibull(x, theta, beta), digits)
        p <- round(pweibull(x, theta, beta, lower.tail = TRUE), digits)
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(" ", italic(f),"(", .(x),") = ", .(d), sep = "")), bty = "n", cex = legend.cex)
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), "(", italic(X)<=.(x), ") = ", .(p), sep = "")), bty = "n", cex = legend.cex)
        if(show.p == TRUE & show.d == TRUE)
            legend("topright", legend = bquote(paste(italic(P), "(", italic(X)<=.(x), ") = ", .(p), ",  ", italic(f),"(", .(x), 
                ") = ", .(d), sep = "")), bty = "n", cex = legend.cex)  
            if (show.dist == TRUE) 
            legend("top", legend = bquote(paste(italic(X)," ~ ", italic(WEI),"(", .(theta),",", .(beta), ")", sep = "")), bty = "n", cex = legend.cex, adj = 0.5)
    }
    if (tail == "upper"){
        polygon(c(x, xv[xv >= x]), c(yv[xv == sigma], yv[xv >= 
            x]), col = shade.col)
        p <- round(pweibull(x, theta, beta, lower.tail = FALSE), digits)
        d <- round(dweibull(x, theta, beta), digits)
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(" ", italic(f),"(", .(x), 
                ") = ", .(d), sep = "")), bty = "n", cex = legend.cex)
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), "(", italic(X)>=.(x), ") = ", .(p), sep = "")), bty = "n", cex = legend.cex)
        if(show.p == TRUE & show.d == TRUE)
            legend("topright", legend = bquote(paste(italic(P), "(", italic(X)>=.(x), ") = ", .(p), ",  ",italic(f),"(", .(x), 
                ") = ", .(d), sep = "")), bty = "n", cex = legend.cex)  
                
        if (show.dist == TRUE) 
            if (show.dist == TRUE) 
           legend("top", legend = bquote(paste(italic(X)," ~ ", italic(WEI),"(", .(theta),",", .(beta), ")", sep = "")), bty = "n", cex = legend.cex, adj = 0.5)
    }
    if (tail == "two"){
        qs <- round(qweibull(c(prob.to.each.tail, 1 - prob.to.each.tail), 
            theta, beta), 1)
        polygon(c(xv[xv <= qs[1]], qs[1]), c(yv[xv <= qs[1]], 
            yv[xv == sigma]), col = shade.col)
        polygon(c(qs[2], xv[xv >= qs[2]]), c(yv[xv == sigma], 
            yv[xv >= qs[2]]), col = shade.col)
        p <- round(prob.to.each.tail * 2, digits)
        if (show.p == TRUE) 
            legend("topright", bty = "n", cex = legend.cex, legend = bquote(paste(italic(P),"(",.(qs[1])<=italic(X), " and ", italic(X)>=.(qs[2]), ") = ", .(p), sep = "")))
        if (show.dist == TRUE) 
            if (show.dist == TRUE) 
            legend("top", legend = bquote(paste(italic(X)," ~ ", italic(WEI),"(", .(theta),",", .(beta), ")", sep = "")), bty = "n", cex = legend.cex, adj = 0.5)
    }
    if (tail == "middle") {
        if(theta > 1){
        polygon(c(xv[xv <= sigma], sigma), c(yv[xv <= sigma], 
            yv[xv == 0]), col = shade.col)
        polygon(c(xv[xv <= from], from), c(yv[xv <= from], yv[xv == 
            0]), col = "white")
        polygon(c(to, xv[xv >= to]), c(yv[xv == sigma], yv[xv >= 
            to]), col = "white")}
        p <- round(pweibull(to, theta) - pweibull(from, beta), digits)
        if (show.p == TRUE) 
            legend("topright", legend = bquote(paste(italic(P),"(", .(from)<="",italic(X)<=.(to), ") = ", .(p), sep = "")), bty = "n", cex = legend.cex)
        if (show.dist == TRUE) 
            legend("top", legend = bquote(paste(italic(X)," ~ ", italic(WEI),"(", .(theta),",", .(beta), ")", sep = "")), bty = "n", cex = legend.cex, adj = 0.5)
    }
}



#-------------------------Binomial distribution------------------------------#

shade.bin <- function (x = NULL, from = NULL, to = NULL, n = 1, p = 0.5, tail = "X=x", 
    show.p = TRUE, show.dist = TRUE, show.d = FALSE, legend.cex = 0.9, 
    digits = 5, ...) 
{
    xv <- seq(0, n, length = n + 1)
    yv <- dbinom(xv, n, p)
    plot(xv, yv, type = "h", xlab = expression(italic(x)), ylab = expression(paste(italic(f), 
        "(", italic(x), ")", sep = "")), col = gray(0.4), lty = 2, 
        ...)
    points(xv, yv, col = gray(0.4), pch = 19)
    if (tail == "X=x") {
        d <- round(dbinom(x, n, p), digits)
        segments(x, 0, x, d, lwd = 2)
        points(x, d, pch = 19)
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(italic(f), 
                "(", .(x), ") = ", .(d), sep = "")), bty = "n", 
                cex = legend.cex)
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X), " = ", .(x), ") = ", .(d), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.p == TRUE & show.d == TRUE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X), " = ", .(x), ") = ", .(d), ",  ", 
                italic(f), "(", .(x), ") = ", .(d), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.dist == TRUE) 
            legend("topleft", legend = bquote(paste(italic(X), 
                "  ~ ", italic(BIN), "(", .(n), ", ", .(p), ")", sep = "")), 
                bty = "n", cex = legend.cex)
    }
    if (tail == "lower") {
        xvnew <- xv[xv <= x]
        yvnew <- dbinom(xvnew, n, p)
        segments(xvnew, rep(0, length(xvnew)), xvnew, yvnew, 
            lwd = 2)
        points(xvnew, yvnew, pch = 19)
        d <- round(dbinom(x, n, p), digits)
        pv <- round(pbinom(x, n, p, lower.tail = TRUE), digits)
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(italic(f), 
                "(", .(x), ") = ", .(d), sep = "")), bty = "n", 
                cex = legend.cex)
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X) <= .(x), ") = ", .(pv), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.p == TRUE & show.d == TRUE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X) <= .(x), ") = ", .(pv), ",  ", 
                italic(f), "(", .(x), ") = ", .(d), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.dist == TRUE) 
            legend("topleft", legend = bquote(paste(italic(X), 
                "  ~ ", italic(BIN), "(", .(n), ", ", .(p), ")", sep = "")), 
                bty = "n", cex = legend.cex)
    }
    if (tail == "upper") {
        xvnew <- xv[xv >= x]
        yvnew <- dbinom(xvnew, n, p)
        segments(xvnew, rep(0, length(xvnew)), xvnew, yvnew, 
            lwd = 2)
        points(xvnew, yvnew, pch = 19)
        d <- round(dbinom(x, n, p), digits)
        pv <- round(pbinom(x, n, p, lower.tail = FALSE), digits)
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(italic(f), 
                "(", .(x), ") = ", .(d), sep = "")), bty = "n", 
                cex = legend.cex)
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X) >= .(x), ") = ", .(pv), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.p == TRUE & show.d == TRUE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X) >= .(x), ") = ", .(pv), ",  ", 
                italic(f), "(", .(x), ") = ", .(d), sep = "")), 
                bty = "n", cex = legend.cex)
        if (show.dist == TRUE) 
            legend("topleft", legend = bquote(paste(italic(X), 
                "  ~ ", italic(BIN), "(", .(n), ", ", .(p), ")", sep = "")), 
                bty = "n", cex = legend.cex)
    }
    if (tail == "two") {
        qs <- sort(c(x, n - x))
        xvnew <- c(xv[xv <= qs[1]], xv[xv >= qs[2]])
        yvnew <- dbinom(xvnew, n, p)
        segments(xvnew, rep(0, length(xvnew)), xvnew, yvnew, 
            lwd = 2)
        points(xvnew, yvnew, pch = 19)
        d <- round(dbinom(x, n, p), digits)
        qs <- sort(c(x, n - x))
        if (show.d == TRUE & show.p == FALSE) 
            legend("topright", legend = bquote(paste(italic(f), 
                "(", .(x), ") = ", .(d), sep = "")), bty = "n", 
                cex = legend.cex)
        pv <- round(pbinom(x, n, p, lower.tail = TRUE), digits) * 
            2
        if (show.p == TRUE & show.d == FALSE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", italic(X) <= .(qs[1]), " and ", italic(X) >= 
                  .(qs[2]), ")=", .(pv), sep = "")), bty = "n", 
                cex = legend.cex)
        if (show.dist == TRUE) 
            legend("topleft", legend = bquote(paste(italic(X), 
                "  ~ ", italic(BIN), "(", .(n), ", ", .(p), ")", sep = "")), 
                bty = "n", cex = legend.cex)
    }
    if (tail == "middle") {
        xvnew <- c(xv[xv >= from & xv <= to])
        yvnew <- dbinom(xvnew, n, p)
        segments(xvnew, rep(0, length(xvnew)), xvnew, yvnew, 
            lwd = 2)
        points(xvnew, yvnew, pch = 19)
        pv <- round(pbinom(to, n, p) - pbinom(from, n, p), digits)
        if (show.p == TRUE) 
            legend("topright", legend = bquote(paste(italic(P), 
                "(", .(from) <= "", italic(X) <= .(to), ") = ", 
                .(pv), sep = "")), bty = "n", cex = legend.cex)
        if (show.dist == TRUE) 
            legend("topleft", legend = bquote(paste(italic(X), 
                "  ~ ", italic(BIN), "(", .(n), ", ", .(p), ")", sep = "")), 
                bty = "n", cex = legend.cex)
    }
}


#-------------------------Poisson distribution------------------------------#

shade.poi<-function(x=NULL,from=NULL,to=NULL,lambda=5,tail="X=x",show.p=TRUE,show.dist=TRUE,show.d=FALSE,legend.cex=.9,digits=5,...){
xv<-seq(0,round(lambda+2*lambda,0),length=round(lambda+2*lambda+1,0))
yv<-dpois(xv,lambda)
plot(xv,yv,type="h",xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")), col=gray(.4),lty=2,...)
points(xv,yv,col=gray(.4),pch=19)

if(tail=="X=x"){
d<-round(dpois(x,lambda),digits)
segments(x,0,x,d,lwd=2)
points(x,d,pch=19)
if(show.d==TRUE&show.p ==FALSE)legend("topright",legend=bquote(paste(italic(f),"(",.(x),") = ",.(d), sep ="")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)," = ",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)," = ",.(x),") = ",.(d), ",  ", italic(f), "(",.(x),") = ", .(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X), " ~ ", italic(POI),"(",.(lambda),")")),bty="n",cex=legend.cex)
}

if(tail=="lower"){
xvnew<-xv[xv<=x];yvnew<-dpois(xvnew,lambda)
segments(xvnew,rep(0,length(xvnew)),xvnew,yvnew,lwd=2)
points(xvnew,yvnew,pch=19)
d<-round(dpois(x,lambda),digits)
pv<-round(ppois(x,lambda,lower.tail=TRUE),digits)
if(show.d==TRUE&show.p ==FALSE)legend("topright",legend=bquote(paste(italic(f),"(",.(x),") = ",.(d), sep ="")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(pv), sep = "")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)<=.(x),") = ",.(pv), ",  ", italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X), " ~ ", italic(POI),"(",.(lambda),")")),bty="n",cex=legend.cex)
}

if(tail=="upper"){
xvnew<-xv[xv>=x];yvnew<-dpois(xvnew,lambda)
segments(xvnew,rep(0,length(xvnew)),xvnew,yvnew,lwd=2)
points(xvnew,yvnew,pch=19)
d<-round(dpois(x,lambda),digits)
pv<-round(ppois(x,lambda,lower.tail=FALSE),digits)
if(show.d==TRUE&show.p ==FALSE)legend("topright",legend=bquote(paste(italic(f),"(",.(x),") = ",.(d), sep ="")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==FALSE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(pv), sep = "")),bty="n",cex=legend.cex)
if(show.p==TRUE&show.d==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",italic(X)>=.(x),") = ",.(pv), ",  ", italic(f),"(",.(x),") = ",.(d), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X), " ~ ", italic(POI),"(",.(lambda),")")),bty="n",cex=legend.cex)
}

if(tail=="two")stop("two-tailed probability not implemented in 'shade.poi'")

if(tail=="middle"){
xvnew<-c(xv[xv>=from&xv<=to]);  yvnew<-dpois(xvnew,lambda)
segments(xvnew,rep(0,length(xvnew)),xvnew,yvnew,lwd=2)
points(xvnew,yvnew,pch=19)
pv<-round(ppois(to,lambda)-ppois(from,lambda),digits)
if(show.p==TRUE)legend("topright",legend=bquote(paste(italic(P),"(",.(from)<="",italic(X)<=.(to),") = ",.(pv), sep = "")),bty="n",cex=legend.cex)
if(show.dist==TRUE)legend("topleft",legend=bquote(paste(italic(X), " ~ ", italic(POI),"(",.(lambda),")")),bty="n",cex=legend.cex)
}
}
