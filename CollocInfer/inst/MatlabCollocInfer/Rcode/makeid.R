
########### identity functions   ##################

make.id = function()
{

  id <- list()

  ##############################################################################

  id$fn <- function(t,x,pars,more)
  {
     return(x)
  }

  ##############################################################################

  id$dfdx <- function(t,x,pars,more)
  {
     g = array(0,c(dim(x),dim(x)[2]))	
     for(i in 1:dim(x)[2]){
    	 g[,i,i] = 1
     }
     return(g)
  }

  ##############################################################################

  id$dfdp <- function(t,x,pars,more)
  {
     return( array(0,c(dim(x),length(pars))) )
  }

  ##############################################################################

  id$d2fdx2 <- function(t,x,pars,more)
  {
    return( array(0,c(dim(x),dim(x)[2],dim(x)[2])) )
  }

  ##############################################################################

  id$d2fdxdp <- function(t,x,pars,more)
  {
   return( array(0,c(dim(x),dim(x)[2],length(pars))) )
  }

  return(id)
}

