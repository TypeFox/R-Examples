`ss.aipe.sc.sensitivity` <-
function(true.psi=NULL, estimated.psi=NULL, c.weights, desired.width=NULL, selected.n=NULL, 
assurance=NULL, certainty=NULL, conf.level=.95, G=10000, print.iter=TRUE, detail=TRUE, ...)
{ if(is.null(estimated.psi) & is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the sample size per group ).", call.=FALSE)
if(!is.null(estimated.psi) & !is.null(selected.n)) stop("You must specify either \'estimated.psi\' or \'selected.n\' (i.e., the per group sample size), but not both.", call.=FALSE)

if(!is.null(assurance)& !is.null(certainty))
    {if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value")}
if(!is.null(certainty)) assurance<- certainty 

if(sum(c.weights)!=0) stop("The sum of the coefficients must be zero")
if(sum(c.weights[c.weights>0])>1) stop("Please use fractions to specify the contrast weights")

if(!is.null(estimated.psi))
    { n <- ss.aipe.sc(psi=estimated.psi, conf.level=conf.level, c.weights=c.weights, width=desired.width, assurance=assurance, Tolerance=1e-7)}
else 
    { n <- selected.n}
    
J<- length(c.weights)

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

# Create population group means for group 1 to group J 
signal=TRUE
mu.c <-rep(NA, J)
for(cw in 1:J)
    {if(c.weights[cw] > 0 & signal)
        { mu.c[cw] <- psi.pop/c.weights[cw]
        signal=FALSE
        }
        
    else{mu.c[cw] <- 0}
    }

for (i in 1:G)
    {if(print.iter) cat(c(i), "\n")

# generate random samples
    group.data<- array(NA, dim=c(n,J))
    x.bar<-rep(NA, J)
    sd.group<- rep(NA, J)
    
    for(p in 1:J)
        {group.data[,p]<- rnorm(n, mean=mu.c[p], sd=1)
        x.bar[p]<- mean(group.data[,p])
        sd.group[p]<- sd(group.data[,p])
        }

# calculate the observed psi and CI for psi
    s.pooled<- mean(sd.group)
    psi.obs[i]<- sum(c.weights*x.bar)/s.pooled
    lambda <- psi.obs[i]/sqrt( sum(c.weights^2) / n )
    lambda.limits <- conf.limits.nct(ncp=lambda, df=n*J-J, conf.level=conf.level) 
    psi.limit.upper <- lambda.limits$Upper.Limit*sqrt(sum(c.weights^2) / n)
    psi.limit.lower <- lambda.limits$Lower.Limit*sqrt(sum(c.weights^2) / n) 
     
    Full.Width[i] <- psi.limit.upper - psi.limit.lower
    Width.from.psi.obs.Lower[i] <- psi.obs[i] - psi.limit.lower
    Width.from.psi.obs.Upper[i] <- psi.limit.upper - psi.obs[i]
    
    Type.I.Error.Upper[i] <- psi.pop > psi.limit.upper
    Type.I.Error.Lower[i] <- psi.pop < psi.limit.lower
    Type.I.Error[i] <- Type.I.Error.Upper[i] | Type.I.Error.Lower[i]
    
    Lower.Limit[i] <- psi.limit.lower
    Upper.Limit[i] <- psi.limit.upper
    }
    
if(detail)
    { list(
        Results=list(psi.obs=psi.obs, Full.Width=Full.Width, Width.from.psi.obs.Lower=Width.from.psi.obs.Lower,
    Width.from.psi.obs.Upper=Width.from.psi.obs.Upper, Type.I.Error.Upper=Type.I.Error.Upper,
    Type.I.Error.Lower=Type.I.Error.Lower, Type.I.Error=Type.I.Error, Lower.Limit=Lower.Limit,
    Upper.Limit=Upper.Limit), 
        Specifications=list(replications=G, True.psi=psi.pop, Estimated.psi=estimated.psi, Desired.Width=desired.width, 
    assurance=assurance, Sample.Size.per.Group=n, Number.of.Groups=J), 
        Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
    sd.full.width=sqrt(var(Full.Width)), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
    mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower),
    mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
    Type.I.Error.Lower=mean(Type.I.Error.Lower)*100))
    }

else {list(
    Specifications=list(replications=G, True.psi=psi.pop, Estimated.psi=estimated.psi, Desired.Width=desired.width, 
assurance=assurance, Sample.Size.per.Group=n, Number.of.Groups=J), 
    Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
sd.full.width=sqrt(var(Full.Width)), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
mean.Width.from.psi.obs.Lower=mean(Width.from.psi.obs.Lower),
mean.Width.from.psi.obs.Upper=mean(Width.from.psi.obs.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
Type.I.Error.Lower=mean(Type.I.Error.Lower)*100))
    }
}
