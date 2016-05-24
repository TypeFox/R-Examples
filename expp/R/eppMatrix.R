

 
eppMatrix <- function(data,  pairs = ~ male + female) {
	
	m = as.character(pairs[[2]][2])
	f  = as.character(pairs[[2]][3])
  
  # TODO: remove repeated lines, warn
  
  new('eppMatrix', male = as.character(data[, m]), female = as.character(data[, f]) )
	
  }













