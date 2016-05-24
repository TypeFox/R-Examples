 sixty=function(scalar, trailsign = F) {

  if(length(scalar)!=1 ){
    stop('first parameter must contain 1 element')
  }	
  ss=abs(3600.0e0*scalar)
  mm=abs(60.0e0*scalar) 
  dd=abs(scalar) 
  result = numeric(3)
  result[1]= as.integer(dd) 
  result[2]= as.integer(mm-60.0e0*result[1])
  result[3]= ss - 3600.e0*result[1] - 60.0e0*result[2]

  if(scalar<0.0e0 ){ 
    if(trailsign) {
      result[1] = -result[1]
    }
    else {
      if(result[1]!=0 )
        result[1] = -result[1]
      else if(result[2]!=0 )
        result[2] = -result[2]
      else 
        result[3] = -result[3]
    } 
  }
  return(result)
}
