# sim by pastor
sim.pastor=function(seed) {
    set.seed(seed)
    n=2000
    x=rbeta(n,4,4)
    x.expit=1/(exp(60*(x-.5))+1)
    x.bin.med=ifelse(x>median(x),1,0)
    eta = .5*(1-2*x.expit)
    y=rbern(n, expit(eta))
    ind = c(which(y==1)[1:500], which(y==0)[1:500])
    dat=data.frame(y,x.star=x,x.star.expit=x.expit)[ind,]
    dat$x.bin.med=ifelse(dat$x.star>median(dat$x.star),1,0)
    dat
}   


sim.my = function (n, seed, label, alpha, beta, e.=NULL, b.=NULL, tr.=NULL) {
    set.seed(seed)
    
    if (endsWith(label,"5")) mu=5.7 else if (endsWith(label,"4")) mu=4.7 else if (endsWith(label,"3")) mu=3.7 else stop("label not recognized")
    sd.x=1.6
    
    if(startsWith(label,"sigmoidimb")) { # imbalance
        x.star=c(rnorm(n-round(n/3), mu, sd=sd.x), mu-abs(rnorm(round(n/3), 0, sd=sd.x)))
    } else if(startsWith(label,"sigmoidlin")) { # unif
        x.star=runif(n)*4*sd.x + mu-2*sd.x
    } else if(startsWith(label,"sigmoidmix")) { # mixture
        x.star=c(rnorm(n*.6, mu, sd=sd.x), rep(mu-2*sd.x, n*.4))
    } else if(startsWith(label,"sigmoidgam")) { # gamma
        x.star=sd.x*scale(rgamma(n=n, 2.5, 1))+mu
    } else {
        x.star=rnorm(n, mu, sd=sd.x)
    } 
    
    if (startsWith(label, "elbow")) {
        x.star.tr = ifelse(x.star>tr., x.star-tr., 0)  
        X=cbind(1, x.star.tr)
    } else if (startsWith(label, "sigmoid")) {
        x.star.expit = expit.2pl(x.star, e=e., b=b.)        
        X=cbind(1, x.star.expit)
    } else stop ("label not supported: "%+%label)
    
    y=rbern(n, expit (X %*% c(alpha, beta)))
    
    sd.err=ifelse(x.star>3.5, 0.05, 1.5-0.414*x.star)            
    x=rnorm(n, x.star, sd.err)
    
 
    
    
    
    dat=data.frame (
        y=y,
        x.star=x.star,
        x.bin.med=ifelse(x.star>median(x.star), 1, 0),
        x.tri = ifelse(x.star>quantile(x.star,2/3),"High",ifelse(x.star>quantile(x.star,1/3),"Medium","Low")),
    
        x=x,
        
        x.tr.1=ifelse(x>log(100), x, 0) ,
        x.tr.2=ifelse(x>log(100), x, log(100)) ,
        x.tr.3=ifelse(x>3.5, x, 0) ,
        x.tr.4=ifelse(x>3.5, x, 3.5), 
        x.tr.5=ifelse(x>3.5, x, 3.5/2), 
        x.tr.6=ifelse(x>5, x, 5) ,
        x.ind =ifelse(x>3.5, 0, 1), 
        
        x.bin.35=ifelse(x>3.5, 1, 0), 
        x.bin.6=ifelse(x>6, 1, 0) ,
        x.bin.log100=ifelse(x>log(100), 1, 0) 
        
    )        
    dat$x.tri=as.factor(dat$x.tri)
    
    if (startsWith(label, "elbow")) {
        dat[["x.star.tr"]]=x.star.tr
    } else if (startsWith(label, "sigmoid")) {
        dat[["x.star.expit"]]=x.star.expit
    } else stop ("label not supported: "%+%label)
    
    dat
}
