#========================================================
#
#  Wrapper function for SABC algorithms
#
#========================================================

## Wrapper function for SABC algorithms. Runs function SABC.noninf() or 
## SABC.inf() depending on value of "method" argument. See documentation for 
## details of those functions.
## Furthermore, the distance functino f.dist is determined.
##
## ARGUMENTS:
## method:  Which algorithm to use. Available algorithms are "noninformative" 
##          and "informative", for the case of a non-informative or an 
##          informative prior.
## other arguments are passed to SABC.noninf() or SABC.inf()
##
## VALUES:
## Output of function SABC.noninf() or SABC.inf()

SABC <- function(r.model, r.prior, d.prior, 
                 n.sample, eps.init, iter.max, 
                 v=ifelse(method=="informative",0.4,1.2), 
                 beta=0.8,  
                 delta=0.1, resample=5*n.sample, 
                 verbose=n.sample, 
                 method="noninformative", adaptjump=TRUE, 
                 summarystats=FALSE, y=NULL, f.summarystats=NULL, ...)
{
  f.dist <- r.model
  
  if( !summarystats & !is.null(y) )
  {
    if(!all(y !=0))
      f.dist <- function(par,...)
      {
        x <- r.model(par,...)
        return(sum((x-y)^2))
      }
    else
      f.dist <- function(par,...)
      {
        x <- r.model(par,...)
        return(sum((x-y)/y)^2)
      }
  }
  
  if( !summarystats & is.null(y) & length( r.model( r.prior(),... ) )>1 )
  {
      stop("r.model needs to return distance to data or data vector y needs to be provided!")
  }
  
  if(summarystats & is.null(y)) stop("data vector y needs to be provided!")
  if(summarystats & is.null(f.summarystats) )
    f.summarystats <- function(x) return(c(x,x^2,x^3))
  
  if (method == "noninformative" | method == "noninf" | method == "non") {
    SABC.noninf(f.dist=f.dist, d.prior=d.prior, r.prior=r.prior, 
                n.sample=n.sample, eps.init=eps.init, iter.max=iter.max, 
                v=v, beta=beta, delta=delta, resample=resample, 
                verbose=verbose, adaptjump=adaptjump, 
                summarystats=summarystats, y=y, f.summarystats=f.summarystats, ...)
  }
  else if (method == "informative" | method == "inf") {
    SABC.inf(f.dist=f.dist, d.prior=d.prior, r.prior=r.prior, 
             n.sample=n.sample, eps.init=eps.init, iter.max=iter.max, 
             v=v, beta=beta, delta=delta, resample=resample, 
             verbose=verbose, adaptjump=adaptjump, 
             summarystats=summarystats, y=y, f.summarystats=f.summarystats, ...)
  }
  else {
    print("error in SABC(): method name not recognized")
  }
}
