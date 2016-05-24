## =============================================================================
##
## Several useful properties to be used with setup.prop
##
## =============================================================================

## exponential decline

 p.exp <- function(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)
       return(y.inf + (y.0-y.inf)*exp(-pmax(x-x.L,0)/x.att))

## linear decline

 p.lin <- function(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)
        return(pmin(y.0,pmax(y.inf,y.0-(y.0-y.inf)*(x-x.L)/x.att)))

## sigmoid decline

 p.sig <- function(x, y.0=1, y.inf=0.5, x.L=0, x.att=1)
        return(y.inf + (y.0-y.inf)*exp(-(x-x.L)/(0.25*x.att))/(1+exp(-(x-x.L)/(0.25*x.att))))

