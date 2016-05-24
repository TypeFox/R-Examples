ss.aipe.c.ancova.sensitivity<- function(true.error.var.ancova=NULL, est.error.var.ancova=NULL, true.error.var.anova=NULL, 
est.error.var.anova=NULL, rho, est.rho=NULL, G=10000, mu.y, sigma.y, mu.x, sigma.x, c.weights, width, conf.level=.95, 
assurance=NULL, certainty=NULL)
{
if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  

#######################################################################
ancova.random.data<- function(mu.y, mu.x, sigma.y, sigma.x, rho, J, n, randomized=TRUE)
    {if(length(mu.y)!=J) stop("'mu.y' should be a J by 1 vector that contains the means of the response variable in each group")
    
    if(randomized==FALSE)
        {if (length(mu.x)!=J) stop ("'mu.x' should be a J by 1 vector that contains the means of the covariate in each group")}
    if(randomized)
        {if (length(mu.x)!=1) stop ("In randomized design, the population means of the covariate in each group should be the same")}
    
    x<-matrix(NA, n,J)
    y<-matrix(NA, n,J)
    cov.matrix<- matrix(c(sigma.y^2, rho*sigma.y*sigma.x, rho*sigma.y*sigma.x, sigma.x^2), 2) 
    
    for(j in 1:J)
        {data<-MASS::mvrnorm(n, mu=c(mu.y[j], mu.x), Sigma=cov.matrix)
        y[,j]<-data[,1]
        x[,j]<-data[,2]
        }
    
    cbind(y,x)
    }
#######################################################################

if(!is.null(sigma.y) & !is.null(true.error.var.anova))
    {if(sigma.y!=sqrt(true.error.var.anova)) stop ("'sigma.y' and 'true.error.var.anova' should be the same")}

if(!is.null(assurance)& !is.null(certainty))
    {if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value")}
if(!is.null(certainty)) assurance<- certainty 

if(is.null(est.error.var.ancova)) 
    {if(is.null(est.error.var.anova)| is.null(est.rho) ) stop ("Please specify the either 'est.error.var.ancova', or both 'est.error.var.anova' and 'est.rho'")
    else est.error.var.ancova<- est.error.var.anova*(1-est.rho^2)
    }
if(length(mu.y)!=length(c.weights)) stop("'mu.y' must be a vector that contains the mean of each group")

n<- ss.aipe.c.ancova(error.var.ancova=est.error.var.ancova, c.weights=c.weights, width=width, conf.level=conf.level, 
assurance=assurance)
J<- length(c.weights)
Psi.pop<- sum(c.weights*mu.y)

Type.I.Error<- rep(NA, G)
Type.I.Error.Upper<- rep(NA, G)
Type.I.Error.Lower<- rep(NA, G)
se.Psi<- rep(NA, G)
se.Psi.res<-rep(NA, G)
se.res.vs.se.full<-rep(NA, G)
Psi.obs<-rep(NA, G)
width.obs<- rep(NA, G)
w.vs.omega<-rep(NA, G)

for(g in 1:G)
    {random.data<-ancova.random.data(mu.y=mu.y, mu.x=mu.x, sigma.y=sigma.y, sigma.x=sigma.x, rho=rho, J=J, n=n)
    
    y<-random.data[, 1:J]
    x<-random.data[, (J+1):(J*2)]
    
    y.vector<-array(random.data[, 1:J])
    x.vector<-array(random.data[, (J+1):(J*2)])
    group<- gl(J, n)
    
    # fit the ANCOVA model with the random data
    fitted.model<- lm(y.vector~ x.vector+group)
    
    # calculate adjusted Y means
    beta<- fitted.model$coefficients[2]
    y.bar.adj<- rep(NA, J)
    for (j in 1:J)
        {y.bar.adj[j]<- mean(y[,j]) - beta*(mean(x[,j])-mean(x))}
    
    # calculate the ANCOVA contrast
    Psi.obs[g]<- sum(c.weights*y.bar.adj)
    
    # extract the MSwithin from the ANCOVA table
    error.var.ancova<- anova(fitted.model)[3,3]
    
    # extract SSwithin.x from ANOVA on the covariate
    x.anova<-lm(x.vector ~ group)
    SSwithin.x<- anova(x.anova)[2,2]
    
    # calculate the standard error of the ANCOVA contrast
    x.bar<- rep(NA, J)
    for (j in 1:J)
        {x.bar[j]<- mean(x[,j]) }
        
    f.x.numerater<- (sum(c.weights*x.bar))^2
    f.x.denominator<- SSwithin.x
    sample.size.weighted<- sum(c.weights^2) / n 
    
    se.Psi2<- error.var.ancova*(sample.size.weighted + f.x.numerater/f.x.denominator)
    se.Psi[g]<- (sqrt(se.Psi2))
    se.Psi.res[g]<- sqrt(error.var.ancova*sample.size.weighted)
    se.res.vs.se.full[g]<- se.Psi.res[g]/se.Psi[g] 
    
    # calculate the confidence interval
    alpha<- 1-conf.level
    nu<- n*J-J-1    
    t.value<- qt(1-alpha/2, df=nu)
    
    ci.Psi<- list(upper=Psi.obs[g]+ t.value*se.Psi[g], lower=Psi.obs[g]- t.value*se.Psi[g])
    width.obs[g]<- ci.Psi$upper - ci.Psi$lower
    
    Type.I.Error.Upper[g] <- Psi.pop > ci.Psi$upper
    Type.I.Error.Lower[g] <- Psi.pop < ci.Psi$lower
    Type.I.Error[g] <- Type.I.Error.Upper[g] | Type.I.Error.Lower[g]
      
    if(width.obs[g]<width) w.vs.omega[g]<-TRUE
    else w.vs.omega[g]<- FALSE
    }

list(results=list(Psi.obs=Psi.obs, se.Psi=se.Psi, se.Psi.restricted=se.Psi.res, se.res.over.se.full=se.res.vs.se.full,
        width.obs=width.obs, Type.I.Error=Type.I.Error, Type.I.Error.Upper=Type.I.Error.Upper, 
        Type.I.Error.Lower=Type.I.Error.Lower), 
    summary=list(Type.I.Error=mean(Type.I.Error)*100, Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
        Type.I.Error.Lower=mean(Type.I.Error.Lower)*100, width.NARROWER.than.desired=mean(w.vs.omega)*100,
        Mean.width.obs=mean(width.obs), Median.width.obs=median(width.obs), 
        Mean.se.res.vs.se.full=mean(se.res.vs.se.full) ), 
    specifications=list(Psi.pop=Psi.pop, Contrast.Weights=c.weights, mu.y=mu.y, mu.x=mu.x, sigma.x=sigma.x,
        Sample.Size.per.Group=n, conf.level=conf.level,
        assurance=assurance, rho=rho, est.rho=est.rho, true.error.var.ANOVA=true.error.var.anova, est.error.var.ANOVA=est.error.var.anova))

}
