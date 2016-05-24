specifyregion <-
function(dim, coord, radius=NULL, form=c("cube", "sphere", "manual"), fading=0){
	
	if((length(dim) < 2) || (length(dim) > 3)){
		stop("Dimensions should represent a 2D or 3D image.")
	}
	if(missing(form)){
		form <- "cube"
	}
	if((form=="cube") || (form=="sphere")){
		if(is.null(radius)){
			stop("Size of the activated region should be provided.")
		} 
		if(length(coord)!=length(dim)){
			stop("Mismatch between dimensions of image space and coordinates.")
		}
	}
	if(form=="manual"){
		if(is.matrix(coord)){
			if(ncol(coord)!=length(dim)){
				stop("Mismatch between dimensions of image space and coordinates.")
			}
		} else {
			stop("The coordinates should be a matrix.")
		}
	}
	
	act <- array(0, dim=dim)
	if(length(dim)==2){
          if(form=="cube"){
	    nv <- as.numeric(radius) + 1
	    for(i in (coord[1]-nv):(coord[1]+nv)){
	      if((i >= 1) && (i <= dim[1])){
		for(j in (coord[2]-nv):(coord[2]+nv)){
	          if((j >= 1) && (j <= dim[2])){
		    if(fading!=0){
		      act[i,j] <- (2*exp(-(((i-coord[1])^2+(j-coord[2])^2))*fading)+2)/4  
		    } else {
		      act[i,j] <- 1
		    }						
		  }
		}
	      }
	    }
	  }
	  if(form=="sphere"){
	    nv <- as.numeric(radius) + 1
	    for(i in (coord[1]-nv):(coord[1]+nv)){
	      if((i>=1)&&(i<=dim[1])){
		for(j in (coord[2]-nv):(coord[2]+nv)){
		  if((j>=1)&&(j<=dim[2])){
		    if(((i-coord[1])^2+(j-coord[2])^2)<= nv^2){
		      if(fading!=0){
		        act[i,j] <- (2*exp(-(((i-coord[1])^2+(j-coord[2])^2))*fading)+2)/4
		      } else {
			act[i,j] <- 1
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  if(form=="manual"){
	    if(fading!=0){
	      stop("Fading activation for manual coordinates is not implemented yet.")  
	    } else {
	      for(i in 1:nrow(coord)){
		  i.1 <- coord[i, 1]
                  i.2 <- coord[i, 2]
                  act[i.1, i.2] <- 1
	      }
	    }
	  }
	} else {
	  if(form=="cube"){
	    nv <- as.numeric(radius) + 1
	    for(i in (coord[1]-nv):(coord[1]+nv)){
	      if((i>=1)&&(i<=dim[1])){
	        for(j in (coord[2]-nv):(coord[2]+nv)){
	     	  if((j>=1)&&(j<=dim[2])){
		    for(k in (coord[3]-nv):(coord[3]+nv)){
		      if((k>=1)&&(k<=dim[3])){
		        if(fading!=0){
			  act[i,j,k] <- (3*exp(-((i-coord[1])^2+(j-coord[2])^2+(k-coord[3])^2)*fading)+3)/6
		        } else {
			  act[i,j,k] <- 1
		        }
		      }
		    }
		  }
	        }
	      }
	    }
	  }
	  if(form=="sphere"){
	    nv <- as.numeric(radius) + 1
            for(i in (coord[1]-nv):(coord[1]+nv)){
              if((i>=1)&&(i<=dim[1])){
                for(j in (coord[2]-nv):(coord[2]+nv)){
                  if((j>=1)&&(j<=dim[2])){
                    for(k in (coord[3]-nv):(coord[3]+nv)){
                      if((k>=1)&&(k<=dim[3])){
                        if(((i-coord[1])^2+(j-coord[2])^2 + (k-coord[3])^2)<= nv^2){
                          if(fading!=0){
                            act[i,j,k] <- (3*exp(-((i-coord[1])^2+(j-coord[2])^2+(k-coord[3])^2)*fading)+3)/6
                          } else {
                            act[i,j,k] <- 1
			  }
                        }
                      }
                    }
                  }
                }
              }
            }
	  }
	  if(form=="manual"){
	    if(fading!=0){
	      stop("Fading activation for manual coordinates is not implemented yet.")
	    } else {
	      for(i in 1:nrow(coord)){
                  i.1 <- coord[i, 1]
                  i.2 <- coord[i, 2]
                  i.3 <- coord[i, 3]
                  act[i.1, i.2, i.3] <- 1 
	      }
	    }
	  }
	}
	return(act)
}

