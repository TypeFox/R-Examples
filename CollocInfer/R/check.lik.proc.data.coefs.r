check.lik.proc.data.coefs = function(lik=NULL,proc=NULL,data=NULL,times=NULL,coefs=NULL)
{
  if(!is.null(data)){
    if(!is.null(times)){
      if( length(times)!= dim(data)[1]){
        stop('Vector of observation times does not matich dimension of data')
      }
    }
    else{ times = 1:dim(data[1]) }
  }

  if(!is.null(lik) & !is.null(times)){
    if(dim(lik$bvals)[1] != length(times)){
      stop('Size of lik$bvals does not match number of observations')
    }
    if(ncol(lik$bvals) != nrow(coefs) ){
      stop('Number of basis functions does not match number of coefficients.')
    }
  }

  if(!is.null(proc)){
    if(is.list(proc$bvals)){
      if( !all( dim(proc$bvals$bvals) == dim(proc$vals$dbvals)) ){
        stop('proc bvals and dbvals dimensions do not match.')
      }
      if( !is.null(proc$more$qpts) ){
        if( nrow(proc$bvals$bvals) != length(proc$more$qpts)){
          stop('proc bvals object dimensions do not correspond to number of quadrature points')
        }
      }
    } else{
       if( nrow(proc$bvals) != length(proc$more$qpts)+1){
          stop('proc bvals object dimensions do not correspond to number of quadrature points')
        }
    }
    if(!is.null(coefs)){
      if(is.list(proc$bvals)){
        if( ncol(proc$bvals$bvals) != nrow(coefs) ){
          stop('Number of basis functions does not match number of coefficients.')
        }
      } else{
       if( ncol(proc$bvals) != nrow(coefs) ){
          stop('Number of basis functions does not match number of coefficients.')
        }
      }
      if(!is.null(proc$more$names)){
        if( length(proc$more$names) != ncol(coefs) ){
          stop('dimension of coefficients does not match length of state variable names')
        }
      }
    }
  }
}