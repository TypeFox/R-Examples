`ss.aipe.sc.ancova.sensitivity` <-
function(true.psi=NULL, estimated.psi=NULL, c.weights, 
desired.width=NULL, selected.n=NULL, mu.x=0, sigma.x=1, rho, divisor="s.ancova",
assurance=NULL, conf.level=.95, G=10000, print.iter=TRUE, detail=TRUE, ...)
{options(warn=-1)
if(divisor!="s.ancova" && divisor!="s.anova") stop("The argument 'divisor' must be either 's.ancova' or 's.anova'")

if(divisor=="s.ancova")
    {
    ##  standardized ANCOVA contrast using s.ancova as divisor
    
    if (is.null(estimated.psi) & is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the sample size per group).", call.=FALSE)
    if (!is.null(estimated.psi) & !is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the sample size per group), but not both.", call.=FALSE)
    
    if (sum(c.weights)!=0) stop("The sum of the contrast weights must be zero")
    if (sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")
    
    J<- length(c.weights)
    width<-desired.width
    sigma.ancova<-2
    sigma.anova<-sqrt(sigma.ancova^2/(1-rho^2))
    mu.y<-rep(NA, J)
    
    if (!is.null(estimated.psi))
        { n <- ss.aipe.sc.ancova(psi=estimated.psi, c.weights=c.weights, width=width, conf.level=conf.level, 
        assurance=assurance, ...) }
    else 
        { n <- selected.n}
    
    group<-gl(J, n)
    y.bar<-rep(NA, J)
    y.bar.adj<-rep(NA, J)
    cov.means<-rep(NA, J)
    
    G <- G
    Full.Width <- rep(NA, G)
    Width.from.psi.obs.Lower <- rep(NA, G)
    Width.from.psi.obs.Upper <- rep(NA, G)
    Type.I.Error.Upper <- rep(NA, G)
    Type.I.Error.Lower <- rep(NA, G)
    Type.I.Error <- rep(NA, G)
    Lower.Limit  <- rep(NA, G)
    Upper.Limit <- rep(NA, G)
    psi.obs <- rep(NA, G)
    psi.pop<- true.psi
    
    # create mu.y, the vector that contains the population mean of each group 
    # J-1 out of the J means are free to be any value, and are set to 1 here in this function
    # The mean for the 1st group whose contrast weight is not zero is determined by the rest J-1 means
    Psi<- psi.pop*sigma.ancova
    signal<-1
    for(j in 1: J)
        {if (c.weights[j]!=0 & signal==1) 
            {mu.y[j]<-0
            signal<-0
            k<-j
            j<-j+1
            }
        else
            {mu.y[j]<-1
            j<-j+1
            } 
        }
    mu.y.k<- (Psi - sum(c.weights*mu.y))/c.weights[k]
    mu.y[k]<-mu.y.k
        
    for (g in 1:G)
        {if(print.iter) cat(c(g), "\n")
        
        # generate random samples
        random.data<- ancova.random.data(mu.y=mu.y, mu.x=mu.x, sigma.y=sigma.anova, sigma.x=sigma.x, rho=rho, J=J, n=n)
        
        y<-random.data[, 1:J]
        x<-random.data[, (J+1):(J*2)]
        
        y.vector<-array(random.data[, 1:J])
        x.vector<-array(random.data[, (J+1):(J*2)])
        
        # fit the ANCOVA model with the random data
        fitted.model<- lm(y.vector~ x.vector+group)
        
        # extract the error variance from ANCOVA table
        error.var.ancova<-anova(fitted.model)[3,3]
        s.ancova<-sqrt(error.var.ancova)
        
        # obtain the regression coefficient and the adjusted means
        beta<- fitted.model$coefficients[2]
                
        for(j in 1:J)
            {y.bar[j]<-mean(y[,j])
            cov.means[j]<- mean(x[,j])
            y.bar.adj[j]<- mean(y[,j]) - beta*(mean(x[,j])-mean(x))
            }
        
        # obtain SSwithin.x from ANOVA on the covariate
        x.anova<-lm(x.vector ~ group)
        SSwithin.x<- anova(x.anova)[2,2]
        
        # calculate the observed psi and CI for psi
        psi.obs[g]<-sum(c.weights*y.bar.adj) / s.ancova    
        
        ci.psi<- ci.sc.ancova(adj.means=y.bar.adj, s.ancova=s.ancova, c.weights=c.weights, n=n, cov.means=cov.means, 
        SSwithin.x=SSwithin.x, conf.level=conf.level)
        
        psi.limit.lower<-ci.psi$psi.limit.lower
        psi.limit.upper<-ci.psi$psi.limit.upper
        
        # evaluate the obtained confidence interval
        Full.Width[g] <- abs(psi.limit.upper - psi.limit.lower)
        
        Width.from.psi.obs.Lower[g] <- psi.obs[g] - psi.limit.lower
        Width.from.psi.obs.Upper[g] <- psi.limit.upper - psi.obs[g]
        
        Type.I.Error.Upper[g] <- psi.pop > psi.limit.upper
        Type.I.Error.Lower[g] <- psi.pop < psi.limit.lower
        Type.I.Error[g] <- Type.I.Error.Upper[g] | Type.I.Error.Lower[g]
        
        Lower.Limit[g] <- psi.limit.lower
        Upper.Limit[g] <- psi.limit.upper
        }
    
    if(detail)
        { result<- list(
            Results=list(psi.obs=psi.obs, Full.Width=Full.Width, Width.from.psi.obs.Lower=Width.from.psi.obs.Lower,
        Width.from.psi.obs.Upper=Width.from.psi.obs.Upper, Type.I.Error.Upper=Type.I.Error.Upper,
        Type.I.Error.Lower=Type.I.Error.Lower, Type.I.Error=Type.I.Error, Lower.Limit=Lower.Limit,
        Upper.Limit=Upper.Limit), 
            
            Specifications=list(replications=G, True.psi=psi.pop, Estimated.psi=estimated.psi, Desired.Width=desired.width, 
        Confidence.Level=conf.level, assurance=assurance, Sample.Size.per.Group=n, Number.of.Groups=J, mu.x=mu.x, 
        sigma.x=sigma.x, rho=rho, divisor=divisor), 
            
            Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
        sd.full.width=sd(Full.Width), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width), mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower),
        mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper), Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
        Type.I.Error.Lower=mean(Type.I.Error.Lower)*100, Type.I.Error=mean(Type.I.Error)*100))
        }
    
    else {result<-list(
           Specifications=list(replications=G, True.psi=psi.pop, Estimated.psi=estimated.psi, Desired.Width=desired.width, 
           Confidence.Level=conf.level, assurance=assurance, Sample.Size.per.Group=n, 
           Number.of.Groups=J, mu.x=mu.x, sigma.x=sigma.x, rho=rho, divisor=divisor), 
    
            Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
        sd.full.width=sd(Full.Width), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width), mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower), mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper), Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
        Type.I.Error.Lower=mean(Type.I.Error.Lower)*100, Type.I.Error=mean(Type.I.Error)*100))
        }
    return(result)
    }
#################################################################################################
#################################################################################################

if(divisor=="s.anova")
## use s.anova to standardize the ANCOVA contrast

    {if(is.null(estimated.psi) & is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the sample size per group ).", call.=FALSE)
    if(!is.null(estimated.psi) & !is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the sample size per group), but not both.", call.=FALSE)
    
    if(sum(c.weights)!=0) stop("The sum of the contrast weights must be zero")
    if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")
    
    width<-desired.width
    sigma.anova<-2
    sigma.ancova<-sqrt(sigma.anova^2*(1-rho^2))
    
    if(!is.null(estimated.psi))
        { n <- ss.aipe.sc.ancova(Psi=estimated.psi*sigma.anova, sigma.ancova=sigma.ancova, sigma.anova=sigma.anova, c.weights=c.weights, 
        divisor="s.anova", width=width, conf.level=conf.level, assurance=assurance, ...) }
    else 
        { n <- selected.n}
    
    J<- length(c.weights)
    group<-gl(J, n)
    y.bar<-rep(NA, J)
    y.bar.adj<-rep(NA, J)
    cov.means<-rep(NA, J)
    mu.y<-rep(NA, J)
    
    G <- G
    Full.Width <- rep(NA, G)
    Width.from.psi.obs.Lower <- rep(NA, G)
    Width.from.psi.obs.Upper <- rep(NA, G)
    Type.I.Error.Upper <- rep(NA, G)
    Type.I.Error.Lower <- rep(NA, G)
    Type.I.Error <- rep(NA, G)
    Lower.Limit  <- rep(NA, G)
    Upper.Limit <- rep(NA, G)
    psi.obs <- rep(NA, G)
    psi.pop<- true.psi
    
    # create mu.y, the vector that contains the population mean of each group 
    # J-1 out of the J means are free to be any value, and are set to 1 here in this function
    # The mean for the 1st group whose contrast weight is not zero is determined by the rest J-1 means
    Psi<-sigma.anova*psi.pop
    signal<-1
    for(j in 1: J)
        {if (c.weights[j]!=0 & signal==1) 
            {mu.y[j]<-0
            signal<-0
            k<-j
            j<-j+1
            }
        else
            {mu.y[j]<-1
            j<-j+1
            } 
        }
    mu.y.k<- (Psi - sum(c.weights*mu.y))/c.weights[k]
    mu.y[k]<-mu.y.k
        
    for (g in 1:G)
        {if(print.iter) cat(c(g), "\n")
        
        # generate random samples
        random.data<- ancova.random.data(mu.y=mu.y, mu.x=mu.x, sigma.y=sigma.anova, sigma.x=sigma.x, rho=rho, J=J, n=n)
        
        y<-random.data[, 1:J]
        x<-random.data[, (J+1):(J*2)]
        
        y.vector<-array(random.data[, 1:J])
        x.vector<-array(random.data[, (J+1):(J*2)])
        
        # fit the ANCOVA model with the random data
        fitted.model<- lm(y.vector~ x.vector+group)
        
        # extract the error variance from ANCOVA table
        error.var.ancova<-anova(fitted.model)[3,3]
        s.ancova<-sqrt(error.var.ancova)
        
        # obtain the regression coefficient and the adjusted means
        beta<- fitted.model$coefficients[2]
                
        for(j in 1:J)
            {y.bar[j]<-mean(y[,j])
            cov.means[j]<- mean(x[,j])
            y.bar.adj[j]<- mean(y[,j]) - beta*(mean(x[,j])-mean(x))
            }
        
        # extract the error variance from the ANOVA on y
        anova.y<- lm(y.vector ~ group)
        error.var.anova<- anova(anova.y)[2,3]
        s.anova<-sqrt(error.var.anova)
        
        # obtain SSwithin.x from ANOVA on the covariate
        x.anova<-lm(x.vector ~ group)
        SSwithin.x<- anova(x.anova)[2,2]
        
        # calculate the observed psi and CI for psi based on s.anova
        psi.obs[g]<-sum(c.weights*y.bar.adj) / s.anova    
        
        ci.psi<-ci.sc.ancova(Psi=sum(c.weights*y.bar.adj), s.anova =s.anova, s.ancova=s.ancova, standardizer = "s.anova", 
        c.weights=c.weights, n=n, cov.means=cov.means, SSwithin.x=SSwithin.x, conf.level = conf.level)
        
        psi.limit.lower<-ci.psi$psi.limit.lower
        psi.limit.upper<-ci.psi$psi.limit.upper
        
        # evaluate the obtained confidence interval
        Full.Width[g] <- abs(psi.limit.upper - psi.limit.lower)
        
        Width.from.psi.obs.Lower[g] <- psi.obs[g] - psi.limit.lower
        Width.from.psi.obs.Upper[g] <- psi.limit.upper - psi.obs[g]
        
        Type.I.Error.Upper[g] <- psi.pop > psi.limit.upper
        Type.I.Error.Lower[g] <- psi.pop < psi.limit.lower
        Type.I.Error[g] <- Type.I.Error.Upper[g] | Type.I.Error.Lower[g]
        
        Lower.Limit[g] <- psi.limit.lower
        Upper.Limit[g] <- psi.limit.upper
        }
    
    if(detail)
        { result<-list(
            Results=list(psi.obs=psi.obs, Full.Width=Full.Width, Width.from.psi.obs.Lower=Width.from.psi.obs.Lower,
        Width.from.psi.obs.Upper=Width.from.psi.obs.Upper, Type.I.Error.Upper=Type.I.Error.Upper,
        Type.I.Error.Lower=Type.I.Error.Lower, Type.I.Error=Type.I.Error, Lower.Limit=Lower.Limit,
        Upper.Limit=Upper.Limit), 
            Specifications=list(replications=G, True.psi=true.psi, Estimated.psi=estimated.psi, rho=rho, 
        Desired.Width=desired.width, assurance=assurance, Sample.Size.per.Group=n, Number.of.Groups=J,
        mu.x=mu.x, sigma.x=sigma.x, divisor=divisor), 
            Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
        sd.full.width=sd(Full.Width), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
        mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower),
        mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
        Type.I.Error.Lower=mean(Type.I.Error.Lower)*100, Type.I.Error=mean(Type.I.Error)*100 ))
        }
    
    else {result<-list(
        Specifications=list(replications=G, True.psi=true.psi, Estimated.psi=estimated.psi, Desired.Width=desired.width, 
    assurance=assurance, Sample.Size.per.Group=n, Number.of.Groups=J, divisor=divisor), 
        Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
    sd.full.width=sd(Full.Width), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
    mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower),
    mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
    Type.I.Error.Lower=mean(Type.I.Error.Lower)*100, Type.I.Error=mean(Type.I.Error)*100 ) )
        }
    return(result)
    }
  
}

