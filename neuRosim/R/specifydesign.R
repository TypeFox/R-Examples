specifydesign <-
function(onsets, durations, totaltime, TR, effectsize, accuracy=0.1, conv=c("none", "gamma", "double-gamma", "Balloon"), cond.names=NULL, param=NULL){

	if(is.list(onsets)){
	  ncond <- length(onsets)  
	} else {
	  if(is.numeric(onsets)){
	    onsets <- list(onsets)
	    ncond <- 1
	  }else{
	    stop("wrong format of onsets: should be list or vector")
	  }
	}
        if(is.list(durations)){
                if(length(durations)!=ncond){
                        stop("Mismatch between number of conditions and durations.")
                }
        } else {
	  if(is.numeric(durations)){
	    durations <- list(durations)
	  }else{
	    stop("wrong format of durations: should be a list or vector")
	  }
        }
       if(is.list(effectsize)){
                if(length(effectsize)!=ncond){
                        stop("Mismatch between number of conditions and effectsize.")
                }
        } else {
	  if(is.numeric(effectsize)){
	    effectsize <- list(effectsize)
	  }else{
	    stop("wrong format of effectsize: should be a list or vector")
	  }
        }
	if(!is.null(param)){
	     if(is.list(param)){
		  if(length(param)!=ncond){
                	        stop("Mismatch between number of conditions and param.")
		  }
	    } else {
		  stop("param should be a list")
	    }
	}

	if(missing(conv)){
		conv <- "none"
	}
	if(is.null(cond.names)){
		cond.names <- c(paste("C", 1:ncond, sep=""))
	}

	design.matrix <- matrix(0,(totaltime/TR),ncond)
	ix <- seq(1,totaltime/accuracy,TR/accuracy)
	for(i in 1:ncond){
		if(conv=="none"){
	                design.matrix[,i] <- stimfunction(totaltime, onsets[[i]], durations[[i]], accuracy)[ix]
 		}
		if(conv=="gamma"){
			s <- stimfunction(totaltime, onsets[[i]], durations[[i]], accuracy)
		    if(!is.null(param)){
			s.conv <- convolve(gammaHRF(seq(accuracy,totaltime,accuracy), param[[i]], verbose=FALSE), rev(s))
		    }else{
			s.conv <- convolve(gammaHRF(seq(accuracy,totaltime,accuracy), verbose=FALSE), rev(s))
		    }
			s.conv <- s.conv/max(s.conv)
			design.matrix[,i] <- effectsize[[i]]*s.conv[ix]
		}
		if(conv=="double-gamma"){
                        s <- stimfunction(totaltime, onsets[[i]], durations[[i]], accuracy)
                        s.conv <- convolve(canonicalHRF(seq(accuracy,totaltime,accuracy), param[[i]], verbose=FALSE), rev(s))
			s.conv <- s.conv/max(s.conv)
                        design.matrix[,i] <- effectsize[[i]]*s.conv[ix]
 		}
		if(conv=="Balloon"){
			s <- stimfunction(totaltime, onsets[[i]], durations[[i]], accuracy)
			s.conv <- balloon(s, totaltime, accuracy)
			s.conv <- s.conv/max(s.conv)
			design.matrix[,i] <- effectsize[[i]]*s.conv[ix]
		}
	}
	colnames(design.matrix) <- cond.names
	return(design.matrix)
}

