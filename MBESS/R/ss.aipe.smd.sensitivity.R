ss.aipe.smd.sensitivity <- function(true.delta=NULL, estimated.delta=NULL, desired.width=NULL, selected.n=NULL, 
assurance=NULL, certainty=NULL, conf.level=.95, G=10000, print.iter=TRUE, ...)
{
if (is.null(assurance) && !is.null (certainty)) assurance <-certainty
if (!is.null(assurance) && is.null (certainty)) assurance -> certainty

if(!is.null(assurance) && !is.null (certainty) )
	{if (assurance!=certainty)  stop("The arguments 'assurance' and 'certainty' must have the same value.")}

if(is.null(estimated.delta) & is.null(selected.n)) stop("You must specify either \'estimated.delta\' or \'selected.n\' (i.e., the per group sample size).", call.=FALSE)
if(!is.null(estimated.delta) & !is.null(selected.n)) stop("You must specify either \'estimated.delta\' or \'selected.n\' (i.e., the per group sample size), but not both.", call.=FALSE)

if(!is.null(estimated.delta))
{
n <- ss.aipe.smd(delta=estimated.delta, conf.level=conf.level, width=desired.width, which.width="Full", 
degree.of.certainty=certainty, Tolerance=1e-7)
}
else
{
n <- selected.n
}

G <- G

Full.Width <- rep(NA, G)
Width.from.d.Lower <- rep(NA, G)
Width.from.d.Upper <- rep(NA, G)
Type.I.Error.Upper <- rep(NA, G)
Type.I.Error.Lower <- rep(NA, G)
Type.I.Error <- rep(NA, G)
Low.Limit  <- rep(NA, G)
Upper.Limit <- rep(NA, G)
d <- rep(NA, G)

delta <- true.delta
for (i in 1:G)
{
if(print.iter==TRUE) cat(c(i),"\n")
d[i] <- smd(Group.1=rnorm(n,delta,1), Group.2 =rnorm(n,0,1))

ci.delta <-  ci.smd(smd = d[i], n.1 = n, n.2 = n, conf.level = conf.level)

Full.Width[i] <- ci.delta$Upper-ci.delta$Lower

Width.from.d.Lower[i] <- d[i]-ci.delta$Lower
Width.from.d.Upper[i] <- ci.delta$Upper-d[i]

Type.I.Error.Upper[i] <- delta > ci.delta$Upper
Type.I.Error.Lower[i] <- delta < ci.delta$Lower
Type.I.Error[i] <- Type.I.Error.Upper[i] | Type.I.Error.Lower[i]

Low.Limit[i] <- ci.delta$Lower
Upper.Limit[i] <- ci.delta$Upper

}

list(Results=list(d=d, Full.Width=Full.Width, Width.from.d.Lower=Width.from.d.Lower,
Width.from.d.Upper=Width.from.d.Upper, Type.I.Error.Upper=Type.I.Error.Upper,
Type.I.Error.Lower=Type.I.Error.Lower, Type.I.Error=Type.I.Error, Low.Limit=Low.Limit,
Upper.Limit=Upper.Limit), Specifications=list(replications=G, true.delta=delta, estimated.delta=estimated.delta, desired.width=desired.width, 
certainty=certainty, n.j=n), Summary=list(mean.full.width=mean(Full.Width),median.full.width=median(Full.Width),
sd.full.width=sqrt(var(Full.Width)), Pct.Less.Desired=mean(Full.Width<=desired.width),
mean.Width.from.d.Lower=mean(Width.from.d.Lower),
mean.Width.from.d.Upper=mean(Width.from.d.Upper),Type.I.Error.Upper=mean(Type.I.Error.Upper)*100,
Type.I.Error.Lower=mean(Type.I.Error.Lower)*100))
}
