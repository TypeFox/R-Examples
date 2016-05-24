sim.stpp <- function(class="poisson", s.region, t.region, npoints=NULL, nsim=1, ...)
{

if (sum(class==c("poisson","cluster","cox","infectious","inhibition"))==0) 
stop("class must be chosen among \"poisson\", \"cluster\", \"cox\", \"inhibition\", \"infectious\" ")

if (class=="poisson") res=rpp(s.region=s.region, t.region=t.region, npoints=npoints, nsim=nsim, ...)
if (class=="cluster") res=rpcp(s.region=s.region, t.region=t.region, npoints=npoints, nsim=nsim, ...)
if (class=="cox") res=rlgcp(s.region=s.region, t.region=t.region, npoints=npoints, nsim=nsim, ...)
if (class=="infectious") res=rinfec(s.region=s.region, t.region=t.region, npoints=npoints, nsim=nsim, ...)
if (class=="inhibition") res=rinter(s.region=s.region, t.region=t.region, npoints=npoints, nsim=nsim, ...)

invisible(return(res))
}

