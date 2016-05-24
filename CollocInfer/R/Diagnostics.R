make.diagnostics = function()
{
    diagnostics.fn = function(t,x,p,more)
    {
        force = array(0,dim(x))

        force[,more$which] = more$psi %*% matrix(p,ncol(more$psi),length(more$which))

        return(more$fn(t,x,more$p,more$more) + force)
    }


    diagnostics.dfdx = function(t,x,p,more)
    {
      return(return(more$dfdx(t,x,more$p,more$more)))

    }

    diagnostics.dfdp = function(t,x,p,more)
    {
        force = array(0,c(dim(x),length(p)))

        k = ncol(more$psi)
        whichp = 0
        for(i in 1:length(more$which)){
          force[,more$which[i],whichp+(1:k)] = more$psi
          whichp = whichp + k
        }
       return(force)
    }


    diagnostics.d2fdx2 = function(t,x,p,more)
    {
      return(return(more$d2fdx2(t,x,more$p,more$more)))

    }

    diagnostics.d2fdxdp = function(t,x,p,more)
    {
      return(array(0,c(nrow(x),ncol(x),ncol(x),length(p))))
    }

    return( list( fn = diagnostics.fn,
                  dfdx = diagnostics.dfdx,
                  dfdp = diagnostics.dfdp,
                  d2fdx2 = diagnostics.d2fdx2,
                  d2fdxdp = diagnostics.d2fdxdp,
                  more = NULL
                )
          )
}
