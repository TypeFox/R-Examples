"ss.aipe.smd" <- function(delta, conf.level=.95, width, which.width="Full", degree.of.certainty=NULL, 
assurance=NULL, certainty=NULL, ...)
{options(warn=-10)
if(!is.null(certainty)& is.null(degree.of.certainty)&is.null(assurance)) degree.of.certainty<-certainty
if (is.null(assurance) && !is.null (degree.of.certainty)& is.null(certainty)) assurance <-degree.of.certainty
if (!is.null(assurance) && is.null (degree.of.certainty)& is.null(certainty)) assurance -> degree.of.certainty

if(!is.null(assurance) && !is.null (degree.of.certainty) && assurance!=degree.of.certainty) 
stop("The arguments 'assurance' and 'degree.of.certainty' must have the same value.")

if(!is.null(assurance) && !is.null (certainty) && assurance!=certainty) 
stop("The arguments 'assurance' and 'certainty' must have the same value.")

if(!is.null(degree.of.certainty) && !is.null (certainty) && degree.of.certainty!=certainty) 
stop("The arguments 'degree.of.certainty' and 'certainty' must have the same value.")

char.expand(which.width, c("Full", "Lower", "Upper"), nomatch = stop("Problems with 'which.width' specification. You must choose either 'Full', 'Lower', or 'Upper'.", call.=FALSE))
alpha <- 1-conf.level

if(which.width=="Lower") stop("At this time \'Lower\' is not implemented; use \'Full\' width")
if(which.width=="Upper.") stop("At this time \'Upper\' is not implemented; use \'Full\' width")

delta.sign <- delta
delta <- abs(delta)

#if(warn==FALSE) options(warn=-10)

if(is.null(degree.of.certainty))
{
if(which.width == "Full") n <- ss.aipe.smd.full(delta, conf.level, width)
if(which.width == "Lower") n <- ss.aipe.smd.lower(delta, conf.level, width)
if(which.width == "Upper") n <- ss.aipe.smd.upper(delta, conf.level, width)
return(n)
}

if(!is.null(degree.of.certainty))
{
if((degree.of.certainty <= 0) | (degree.of.certainty >= 1)) stop("The 'degree.of.certainty' must either be NULL or some value greater than zero and less than unity.", call.=FALSE)
if(degree.of.certainty <= .50) stop("The 'degree.of.certainty' should be > .5 (but less than 1).", call.=FALSE)

if(which.width == "Full") 
{
n0 <- ss.aipe.smd(delta=delta, conf.level=conf.level, width=width, which.width="Full", ...)

Limit.2.Sided <- ci.smd(smd=delta, n.1=n0, n.2=n0, conf.level=NULL, alpha.lower=(1-degree.of.certainty)/2, 
     alpha.upper=(1-degree.of.certainty)/2, ...)$Upper

Limit.1.Sided <- ci.smd(smd=delta, n.1=n0, n.2=n0, conf.level=NULL, alpha.lower=0, alpha.upper=1-degree.of.certainty, ...)$Upper

determine.limit <- function(current.delta.limit=current.delta.limit, samp.size=n0, delta=delta, degree.of.certainty=degree.of.certainty)
{
Less <- pt(q=delta2lambda(delta=-current.delta.limit, n.1=samp.size, n.2=samp.size), df=2*samp.size-2, ncp=delta2lambda(delta=delta, n.1=samp.size, n.2=samp.size))
Greater <- 1 - pt(q=delta2lambda(delta=current.delta.limit, n.1=samp.size, n.2=samp.size), df=2*samp.size-2, ncp=delta2lambda(delta=delta, n.1=samp.size, n.2=samp.size))
Expected.Widths.Too.Large <- Less + Greater
return((Expected.Widths.Too.Large - (1-degree.of.certainty))^2)
}

Optimize.Result <- optimize(f=determine.limit, interval=c(Limit.1.Sided, Limit.2.Sided), delta=delta, degree.of.certainty=degree.of.certainty)
n <- ss.aipe.smd(delta=Optimize.Result$minimum, conf.level=1-alpha, width=width, which.width="Full", degree.of.certainty=NULL, ...)

}


if(which.width == "Lower") 
{
n0 <- ss.aipe.smd(delta, conf.level, width, which.width="Lower")
lims <- ci.smd(smd=delta, n.1=n0, n.2=n0, conf.level=degree.of.certainty)
Expected.Max <- max(abs(lims$Lower), abs(lims$Upper))
n <- ss.aipe.smd(delta=Expected.Max, conf.level=conf.level, width=width, which.width="Lower")
}

if(which.width == "Upper")
{
n0 <- ss.aipe.smd(delta, conf.level, width, which.width="Upper")
lims <- ci.smd(smd=delta, n.1=n0, n.2=n0, conf.level=degree.of.certainty)
Expected.Max <- max(abs(lims$Lower), abs(lims$Upper))
n <- ss.aipe.smd(delta=Expected.Max, conf.level=conf.level, width=width, which.width="Upper")
}

#if(warn==FALSE) options(warn=1)

return(n)
}
}
