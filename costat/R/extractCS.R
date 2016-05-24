extractCS <-
function(object, slot=c("startpar", "endpar", "convergence", "minvar",
	"pvals", "lcts"),
	coeftype=c("all", "alpha", "beta", "alphafunc", "betafunc"), solno, ...)
{

if (class(object) != "csFSS")
	stop("Input object not of correct type. Should be of class csFSS")

if (missing(solno))
	solno <- 1:nrow(object$startpar)

slot <- match.arg(slot)
coeftype <- match.arg(coeftype)

dcoeftofn <- function(v, n=256, filter.number=1, family="DaubExPhase",
	type="alpha")	{
	nv <- length(v)
	alpha <- v[1:(nv/2)]
	beta <- v[(nv/2+1):nv]
	ans <- coeftofn(alpha=alpha, beta=beta, n=n,
		filter.number=filter.number, family=family)
	if (type=="alpha")
		return(ans$alpha)
	else if (type=="beta")
		return(ans$beta)
	else
		stop("Unknown type")
	}

if (slot=="startpar" || slot=="endpar")	{

	if (slot=="startpar")
		sp <- object$startpar[solno,]
	else
		sp <- object$endpar[solno,]

	if (coeftype=="all")
		return(sp)

	else if (coeftype=="alpha")	{
		nc <- dim(sp)[2]
		sp <- sp[,1:(nc/2)]
		return(sp)
		}
	else if (coeftype=="beta")	{
		nc <- dim(sp)[2]
		sp <- sp[,(1+nc/2):nc]
		return(sp)
		}
	else if (coeftype=="alphafunc")	{
		ans <- apply(sp, 1, dcoeftofn, type="alpha", ...)
		return(t(ans))
		}
	else if (coeftype=="betafunc")	{
		ans <- apply(sp, 1, dcoeftofn, type="beta", ...)
		return(t(ans))
		}
	}

else if (slot=="convergence")	
	return(object$convergence[solno])

else if (slot=="minvar")	
	return(object$minvar[solno])

else if (slot=="pvals")	
	return(object$pvals[solno])

else if (slot=="lcts")	{
	sp <- object$endpar[solno,]
	ans <- apply(sp, 1, prodcomb, tsx=object$tsx, tsy=object$tsy,
		filter.number=object$filter.number,
		family=object$family)
	return(ans)
	}
}
