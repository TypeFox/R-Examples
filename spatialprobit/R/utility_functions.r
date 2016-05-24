###################
# DISTANCE RELATED
###################

################## ----------------------------- ########################
#@desc  convertes degrees in radians. similar to the python radians() function
radians <- function( degrees ){
  return( degrees * pi / 180 )
}

################## ----------------------------- ########################
#@desc  
#@source  http://www.movable-type.co.uk/scripts/latlong.html 
p2pdistance <- function(start, end){
  R = 6378 # earth radius in km
  lat1=as.numeric(start[1]); lon1=as.numeric(start[2])
  lat2=as.numeric(  end[1]); lon2=as.numeric(  end[2])
  dlat = radians(lat2-lat1)
  dlon = radians(lon2-lon1)
  a = sin(dlat/2)*sin(dlat/2)+cos(radians(lat1))*cos(radians(lat2))*sin(radians(dlon/2))*sin(radians(dlon/2))
  c = 2* atan2(sqrt(a), sqrt(1-a))
 return( R * c )
}

###################
# MATRIX PROCESSING
###################

################## ----------------------------- ########################
#@desc  
#@source  
matrix_is_na <- function( m ){
  for( i in 1:nrow( m ) )
    for( j in 1:nrow( m ) )
      if( ! is.na( m[i,j] ) ){      
       return( F )
      }
  return( T )
}


##########################
# SNA FUNCTIONS UNASIGNED 
##########################

################## ----------------------------- ########################
#@desc    function that outputs the compound of two networks without including
#         the diagonal.
#@source  90105 Fundamentals of social network homeworks

CR <- function( m1 , m2 , binary=T){
  diag(m1) <- diag(m2)<-0
  m3 <- m1 %*% m2
  diag( m3 )
  if( binary ){
    m3[ m3>0 ] <- 1
  }
  return( m3 )
}

