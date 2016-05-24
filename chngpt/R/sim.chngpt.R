expit.2pl=function(t,e,b) sapply(t, function(t) 1/(1+exp(b*(t-e))))

sim.chngpt = function (label, 
    type=c("step","hinge","segmented","stegmented"),
    n, seed, 
    beta, coef.z=log(1.4), 
    x.distr="norm", 
    e., b.=-Inf,
    alpha.candidate=NULL, verbose=FALSE) 
{
    
    set.seed(seed)
    if (!requireNamespace("mvtnorm")) {print("mvtnorm does not load successfully"); return (NULL) }
    
    if (missing(type)) stop("type mssing")
    type<-match.arg(type)    
    
    mu=4.7 
    sd.x=1.6
    
    # generate covariates
    if(x.distr=="imb") { # imbalance
        x=c(rnorm(n-round(n/3), mu, sd=sd.x), mu-abs(rnorm(round(n/3), 0, sd=sd.x)))
        z=rep(1,n)
    } else if(x.distr=="lin") { # unif
        x=runif(n)*4*sd.x + mu-2*sd.x
        z=rep(1,n)
    } else if(x.distr=="mix") { # mixture
        x=c(rnorm(n*.6, mu, sd=sd.x), rep(mu-2*sd.x, n*.4))
        z=rep(1,n)
    } else if(x.distr=="gam") { # gamma
        x=sd.x*scale(rgamma(n=n, 2.5, 1))+mu
        z=scale(rgamma(n=n, 2.5, 1))        
    } else if(startsWith(x.distr,"norm")) { # normal
        if (x.distr=="norm") {
             rho=0
        } else if (x.distr=="norm3") {
            rho=0.3 
        } else if (x.distr=="norm6") {
            rho=0.6
        } else {
            stop("x.distr not supported: "%+%x.distr)
        }        
        tmp=mvtnorm::rmvnorm(n, mean = c(mu,0), sigma = matrix(c(sd.x^2,sd.x*rho,sd.x*rho,1),2)) # use mvtnorm
        x=tmp[,1]
        z=tmp[,2]    
    } else stop("x.distr not supported: "%+%x.distr)    
    
    x.star = expit.2pl(x, e=e., b=b.)  
    # get alpha
    alpha=chngpt::sim.alphas[[ifelse(label=="sigmoid5","sigmoid2",label)%+%"_"%+%x.distr]][e.%+%"", beta%+%""]
    if(is.null(alpha)) stop("alpha not found")
    
    # make design matrix and coefficients
    X=cbind(1,     z,        x,   x.star,   x.star*(x-e.), z*x,   z*x.star,   z*x.star*(x-e.))
    coef.=c(alpha, z=coef.z, x=0, x.star=0, x.hinge=0,     z.x=0, z.x.star=0, z.x.hinge=0)

    if (label=="sigmoid1") { 
    # intercept only
        coef.["x.star"]=beta
        coef.["z"]=0
    
    } else if (label %in% c("sigmoid2") ) { 
    # intercept + main effect
        if (type=="step") {
            coef.[1:5]=c(alpha, coef.z,          0,    beta,     0) 
        } else if (type=="hinge") {
            coef.[1:5]=c(alpha, coef.z,          0,       0,  beta) 
        } else if (type=="segmented") {
            coef.[1:5]=c(alpha, coef.z,  -log(.67),       0,  beta) 
        } else if (type=="stegmented") {
            coef.[1:5]=c(2,     coef.z,   log(.67), log(.67), beta) # all effects of x in the same direction, subject to perfect separation, though that does not seem to be the main problem
            #coef.[1:5]=c(0, coef.z, -log(.67), log(.67), beta) # effects of x and x.star in different direction
        }
        
    } else if (label %in% c("sigmoid3","sigmoid4","sigmoid5")) { 
    # intercept + main effect + interaction
        beta.var.name=switch(type,step="x.star",hinge="x.hinge",segmented="x.hinge",stegmented="x.hinge")
        coef.[beta.var.name]=switch(label,sigmoid3=log(.67),sigmoid4=-log(.67),sigmoid5=0)
        coef.["z."%+%beta.var.name]=beta
        if (type=="segmented") {coef.["x"]=tmp; coef.["z.x"]=log(.67) } 
        
    } else if (label=="sigmoid6") { # special treatment, model misspecification
        coef.=c(alpha, coef.z, log(.67),  beta)
        X=cbind(1, z, x.star, x.star*z^3)
    
    } else stop("label not supported: "%+%label) 
    
       
    if (verbose) {
        myprint(coef.)
        if (verbose==2) {
            str(X)
            print(colnames(X))
        }
    }
            
    y=rbern(n, expit(X %*% coef.))

    data.frame (
        y=y,
        z=z,         
        
        x=x,
        x.star=x.star,
        x.hinge=x.star*(x-e.),
        x.bin.med=ifelse(x>median(x), 1, 0),
        x.tri = factor(ifelse(x>quantile(x,2/3),"High",ifelse(x>quantile(x,1/3),"Medium","Low")), levels=c("Low","Medium","High")),
        
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
}
