expo.scale <- function(DATA,center=TRUE,scale=TRUE){
	
	column.names <- colnames(DATA)
	DATA_dims <- dim(DATA)
	
	######THIS BLOCK INTENDED TO CREATE CENTERS AND SCALES BASED ON REQUESTS.
	if(class(scale)=="character"){
		if(tolower(scale)=="ss1"){ ##if you want to get SS1
			if(is.logical(center) && center){
				center <- apply(DATA, 2, mean, na.rm = TRUE)
			}else if(is.logical(center) && !center){
				center <- rep(0,DATA_dims[2])
			}##else, you are on your own. I perform rudimentary checks later for lengths and whatnot.
			scale <- apply(DATA, 2, sd, na.rm = TRUE)*(sqrt(DATA_dims[1]-1))
		}
		else if(tolower(scale)=="sd"){ ##if you want to get sd-norm; no center.
			center <- rep(0,DATA_dims[2])
			scale <- apply(DATA, 2, sd, na.rm = TRUE)
		}
		else if(tolower(scale)=="rms"){ ##I probably don't need this, but will be consistent.
			if(is.logical(center) && !center){
				center<-FALSE #=rep(0,DATA_dims[2])
				scale<-TRUE #=sqrt(colSums(DATA^2)/(nrow(DATA)-1))
			}else if(is.logical(center) && center){ ##you just wanted a z score then...
				center<-TRUE #=apply(DATA, 2, mean, na.rm = TRUE)	
				scale<-TRUE #=apply(DATA, 2, sd, na.rm = TRUE)
			}##else, you are on your own. I perform rudimentary checks later for lengths and whatnot.
		}		
		else if(tolower(scale)=="z"){ ##both z and SS1 need center/scale
			center<-TRUE #=apply(DATA, 2, mean, na.rm = TRUE)	
			scale<-TRUE #=apply(DATA, 2, sd, na.rm = TRUE)
		}else{ ## you made a booboo
			center<-TRUE
			scale<-TRUE 		
			print("Something is wrong with 'center' and 'scale'. 'center' and 'scale' both set to TRUE.")					
		}
		##will include other normalization schemes in the future. The Mu-methods will need row norms and other norms.
	}


	######THIS BLOCK INTENDED TO PERFORM A SET OF CHECKS
	if((!is.logical(center)) && (!(class(center)=="numeric" && length(center)==DATA_dims[2]))){
		center <- TRUE
		print("Something is wrong with 'center'. 'center' set to TRUE.")
	}
	if((!is.logical(scale)) && (!(class(scale)=="numeric" && length(scale)==DATA_dims[2]))){
		scale <- TRUE
		print("Something is wrong with 'scale'. 'scale' set to TRUE.")
	}	

	###NOW PERFORM THE ACTUAL NORMS.
	scale.info <- scale(DATA,center=center,scale=scale)	
		
	#center checks
	center.out <- attributes(scale.info)$`scaled:center`
	if(is.null(center.out)){
		center.out <- rep(0,DATA_dims[2]) ##create a 0 center
	}
	if(is.null(names(center.out))){
		names(center.out) <- column.names
	}
	#scale checks
	scale.out <- attributes(scale.info)$`scaled:scale`
	if(is.null(scale.out)){
		scale.out <- rep(1,DATA_dims[2]) ##create a 1 scale
	}
	if(is.null(names(scale.out))){
		names(scale.out) <- column.names
	}
	
	#this forces every data matrix to pass through here to have a center and a scale attribute.		
	attributes(scale.info)$`scaled:center` <- center.out
	attributes(scale.info)$`scaled:scale` <- scale.out	
		#recall DATA - 0 center * 1 scale = DATA.
	return(scale.info)
}