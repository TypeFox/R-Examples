
########### exponential of identity functions   ##################

make.exp = function()
{
  explist = list()

  ##############################################################################

  explist$fn <- function(t,x,pars,more)
  {
     return(exp(x))
  }

  ##############################################################################

  explist$dfdx <- function(t,x,pars,more)
  {
     g = array(0,c(dim(x),dim(x)[2]))
     for(i in 1:dim(x)[2]){
    	 g[,i,i] = exp(x[,i])
     }
     return(g)
  }

  ##############################################################################

  explist$dfdp <- function(t,x,pars,more)
  {
     return( array(0,c(dim(x),length(pars))) )
  }

  ##############################################################################

  explist$d2fdx2 <- function(t,x,pars,more)
  {
    g =  array(0,c(dim(x),dim(x)[2],dim(x)[2]))
    for(i in 1:dim(x)[2]){
    	g[,i,i,i] = exp(x[,i])
    }
    return(g)
  }

  ##############################################################################

  explist$d2fdxdp <- function(t,x,pars,more)
  {
   return( array(0,c(dim(x),dim(x)[2],length(pars))) )
  }

  return(explist)
}

