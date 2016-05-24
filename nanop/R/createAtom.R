createAtom <- function(name, base=NA, sigma=0.01, scatterLength=NA, scatterFactor=NA){
  
  if(is.na(base[1]))
    base <- getBase(name)
  if(is.na(scatterLength[1]))
    scatterLength <- getScatterLength(name)
  if(is.na(scatterFactor[1]))
    scatterFactor <- getScatterFactor(name)
	
  list(name=name, sigma=sigma, scatterLength=scatterLength, scatterFactor=scatterFactor, base=base)	
}

getBase <- function(name){
  
  s <- pi/3
  if(name == "Na" || name == "Cu" || name == "Fe" || name == "Ca")
    base = matrix(c(0,0,0), ncol=3, byrow=TRUE)  #Na, Cu, Fe, Ca
  else if(name == "Cl")
	base = matrix(c(0.5,0,0), ncol=3, byrow=TRUE) #Cl
  else if(name == "Ti")	
	base = matrix(c(0.5,0.5,0.5), ncol=3, byrow=TRUE)  #Ti
  else if(name == "Zn" || name == "Mg"){     
    base <- matrix(c(c(cos(s),  (1/3)*sin(s),         0),   #Zn, Mg
	       		     c(cos(s), -(1/3)*sin(s),   cos(s))),
		           ncol=3, byrow=TRUE )
  }
  else if(name == "S"){
    base <- matrix(c(c(cos(s),  (1/3)*sin(s),        0.375),  #S
	                 c(cos(s), -(1/3)*sin(s),  (0.375-0.5))),
		           ncol=3, byrow=TRUE)
  }
  else if(name == "O3"){
    base <- matrix(c(c(0.5, 0.5, 0),     #O3
	                 c(0.5, 0,   0.5),  
			         c(0,   0.5, 0.5)),  
				   ncol=3, byrow=TRUE)
  }
  else
    stop("unknown atom name \n", immediate. = TRUE)
  
  base				 
}