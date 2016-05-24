
#'@export
par.to.sw.base = function(par, coeff=0.473){
  
  sw = par* coeff
}


#'@name par.to.sw
#'@aliases par.to.sw.base
#'
#'@title Convert PAR to shortwave
#'@description 
#'Returns incoming shortwave radiation by converting PAR measuremt.
#'
#'@usage
#'par.to.sw.base(par, coeff=0.473)
#'
#'par.to.sw(data, par.col='par', coeff=0.473)
#'
#'@param data Object of class data.frame with column name 'par' (units umol/m^2/sec)
#'@param par.col String of alternative name for PAR column 
#'@param par Numeric vector of PAR values (umol/m^2/sec)
#'@param coeff 
#'Numerical coefficient to convert PAR (umol/m^2/sec) to SW (W/m^2). 
#'Defaults to value from Britton and Dodd (1976).
#'
#'
#'@return 
#'#For par.to.sw
#'
#'Object of class data.frame with column name 'sw' and other values from \code{ts.data}
#'
#'#For par.to.sw.base
#'
#'Numeric vector of shortwave values with units W/m^2
#'
#'@keywords methods math
#'@references
#'Britton, C. M., and J. D. Dodd. \emph{Relationships of photosynthetically active radiation and shortwave irradiance.} 
#'Agricultural Meteorology 17, no. 1 (1976): 1-7.
#'@author
#'LakeMetabolizer
#'@seealso \link{sw.to.par}
#'@examples 
#'par <- 800
#'par.to.sw.base(par)
#'@export
par.to.sw = function(data, par.col='par', coeff=0.473){
  
  output = data
  
  indx = var.indx(data, par.col)
  
  data[,indx] = par.to.sw.base(data[,indx], coeff=coeff)
  names(data)[indx] = 'sw' #rename to par
  
  return(data)
}


