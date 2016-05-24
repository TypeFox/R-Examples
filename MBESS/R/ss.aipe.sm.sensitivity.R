ss.aipe.sm.sensitivity <- function(true.sm=NULL, estimated.sm=NULL, desired.width=NULL, selected.n=NULL, 
assurance=NULL, certainty=NULL, conf.level=.95, G=10000, print.iter=TRUE, detail=TRUE, ...)
{ options(warn=-1)

if(is.null(estimated.sm) & is.null(selected.n)) stop("You must specify either 'estimated.sm' or 'selected.n' (i.e., the sample size per group ).", call.=FALSE)
if(!is.null(estimated.sm) & !is.null(selected.n)) stop("You must specify either 'estimated.sm' or 'selected.n' (i.e., the per group sample size), but not both.", call.=FALSE)

if(!is.null(assurance)& !is.null(certainty))
	{if(assurance!=certainty) stop("'assurance' and 'certainty' must have the same value") }
if(!is.null(certainty)) assurance<- certainty

if(!is.null(estimated.sm))
    { n <- ss.aipe.sm(sm=estimated.sm, conf.level=conf.level, width=desired.width, assurance=assurance, Tolerance=1e-7)}
else 
    { n <- selected.n}

G <- G
Full.Width <- rep(NA, G)
Width.from.sm.obs.Lower <- rep(NA, G)
Width.from.sm.obs.Upper <- rep(NA, G)
Type.I.Error.Upper <- rep(NA, G)
Type.I.Error.Lower <- rep(NA, G)
Type.I.Error <- rep(NA, G)
Lower.Limit  <- rep(NA, G)
Upper.Limit <- rep(NA, G)
sm.obs <- rep(NA, G)
sm.pop<- true.sm

for (i in 1:G)
    {if(print.iter) cat(c(i), "\n")
    sample.data<- rnorm(n, mean=sm.pop, sd=1)
    x.bar<- mean(sample.data)
    sm.obs[i]<- x.bar / sd(sample.data)
    
    lambda<- sm.obs[i]*sqrt(n)
    lambda.limits<- conf.limits.nct(ncp=lambda, df=n-1, conf.level=conf.level)
    sm.limit.upper <- lambda.limits$Upper.Limit/sqrt(n)
    sm.limit.lower <- lambda.limits$Lower.Limit/sqrt(n)  
    
    Full.Width[i] <- sm.limit.upper - sm.limit.lower
    Width.from.sm.obs.Lower[i] <- sm.obs[i] - sm.limit.lower
    Width.from.sm.obs.Upper[i] <- sm.limit.upper - sm.obs[i]
    
    Type.I.Error.Upper[i] <- sm.pop > sm.limit.upper
    Type.I.Error.Lower[i] <- sm.pop < sm.limit.lower
    Type.I.Error[i] <- Type.I.Error.Upper[i] | Type.I.Error.Lower[i]
    
    Lower.Limit[i] <- sm.limit.lower
    Upper.Limit[i] <- sm.limit.upper
    }
    
if(detail)
    { list(
        Results=list(sm.obs=sm.obs, Full.Width=Full.Width, Width.from.sm.obs.Lower=Width.from.sm.obs.Lower,
    Width.from.sm.obs.Upper=Width.from.sm.obs.Upper, Type.I.Error.Upper=Type.I.Error.Upper,
    Type.I.Error.Lower=Type.I.Error.Lower, Type.I.Error=Type.I.Error, Lower.Limit=Lower.Limit,
    Upper.Limit=Upper.Limit), 
        Specifications=list(replications=G, True.sm=sm.pop, Estimated.sm=estimated.sm, Desired.Width=desired.width, 
    assurance=assurance, Sample.Size=n), 
        Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
    sd.full.width=sqrt(var(Full.Width)), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
    mean.Width.from.sm.obs.Lower=mean(Width.from.sm.obs.Lower),
    mean.Width.from.sm.obs.Upper=mean(Width.from.sm.obs.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
    Type.I.Error.Lower=mean(Type.I.Error.Lower)*100))
    }
    
else {list(
    Specifications=list(replications=G, True.sm=sm.pop, Estimated.sm=estimated.sm, Desired.Width=desired.width, 
assurance=assurance, Sample.Size=n), 
    Summary=list(mean.full.width=mean(Full.Width), median.full.width=median(Full.Width),
sd.full.width=sqrt(var(Full.Width)), Pct.Width.obs.NARROWER.than.desired=mean(Full.Width<=desired.width),
mean.Width.from.sm.obs.Lower=mean(Width.from.sm.obs.Lower),
mean.Width.from.sm.obs.Upper=mean(Width.from.sm.obs.Upper), Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
Type.I.Error.Lower=mean(Type.I.Error.Lower)*100))
    }
}
