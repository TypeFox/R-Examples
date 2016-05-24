setClass(
		Class = "NG_Visualization2d",
		representation = representation(
				mat = "matrix",      ## transition matrix, should be transMat
				axisMat = "matrix",  ## History of which variable was at which axis
				transitionKind = "numeric"
		),
		contains = c("NG_Visualization","VIRTUAL")
)


ng_2dRotationMatrix <- function(viz,ngEnv) {
	
	mat <- viz@mat*0
	text <- ''
	state <- ngEnv$bulletState
	graph <- ngEnv$graph
	
	ifrom <- match(unlist(strsplit(state$from,split = graph@sep, fixed = TRUE)),viz@varList)
	
	if(state$to != ''){
		ito <- match(unlist(strsplit(state$to,split = graph@sep, fixed = TRUE)),viz@varList)
		
		cosT <- cos(state$percentage*pi/2)
		sinT <- sin(state$percentage*pi/2)
		
	}else{
		ito <- c(-1,-1)             
	}
	
	
	## We have 4 situations:
	## - leave node: add new col to axisMat
	## - same transition: leave axisMat
	## - arrive at node
	## - switch node
	
	## check if state is on the same transition
	if(all(c(state$from == viz@from,state$to == viz@to))){
		#text <- paste(text,'same transition.\n')
		
	}else if(all(viz@axisMat[,2] == 0) & (state$to != '')){ ## leave node
		#text <- paste(text,'leave node.\n')
		
		if((viz@axisMat[1,1] == ito[1]) | (viz@axisMat[2,1] == ito[2])){
			#text <- paste(text,'LN: 1-1 or 2-2 same.\n')
			viz@axisMat <- cbind(viz@axisMat[,1],ito)
			
		}else if((viz@axisMat[2,1] == ito[1]) | (viz@axisMat[1,1] == ito[2])){
			#text <- paste(text,'LN: 2-1 or 1-2 same.\n')
			viz@axisMat <- cbind(viz@axisMat[,1],ito[c(2,1)])
			
		}else{
			#text <- paste(text,'LN: 4D transition.\n')
			viz@axisMat <- cbind(viz@axisMat[,1],ito)
			
		}
		
		## check which axis stays the same
		if(viz@axisMat[1,1] == viz@axisMat[1,2] ){ ## x-axis the same
			viz@transitionKind <- 1
		}else if(viz@axisMat[2,1] == viz@axisMat[2,2]){ ## y-axis the same
			viz@transitionKind <- 2
		}else{ ## 4D transition
			viz@transitionKind <- 3
		}
		
	}else if((state$to == '') & all(range(viz@axisMat[,2]) == range(ifrom))){ ## arrive at new node
		#text <- paste(text,'arrive at new node.\n')
		viz@axisMat <- cbind(viz@axisMat[,2],c(0,0))
		viz@transitionKind <- 0
		
	}else if((state$to == '') & all(range(viz@axisMat[,1]) == range(ifrom))){ ## if going back to the same node as before
		#text <- paste(text,'arrive at old node.\n')
		viz@axisMat <- cbind(viz@axisMat[,1],c(0,0))
		viz@transitionKind <- 0
		
	}else{ ## switch node
		#text <- paste(text,'switch node.\n')
		viz@axisMat <- cbind(ifrom,c(0,0))
		viz@transitionKind <- 0
	}              
	
	
	## make transition matrix
	## x-axis
	
	if(viz@transitionKind == 0){ ## if at node
		
		mat[viz@axisMat[1,1],1] <- 1
		mat[viz@axisMat[2,1],2] <- 1
		
	}else if(viz@transitionKind == 1){ ## 3D transition x-axis same
		
		mat[viz@axisMat[1,1],1] <- 1
		
		mat[viz@axisMat[2,1],2] <- cosT
		mat[viz@axisMat[2,2],2] <- sinT
		
	}else if(viz@transitionKind == 2){ ## 3D transition y-axis same
		mat[viz@axisMat[2,1],2] <- 1
		
		mat[viz@axisMat[1,1],1] <- cosT
		mat[viz@axisMat[1,2],1] <- sinT
		
		
	}else if(viz@transitionKind == 3){## 4D transition
		
		mat[viz@axisMat[1,1],1] <- cosT
		mat[viz@axisMat[1,2],1] <- sinT
		
		mat[viz@axisMat[2,1],2] <- cosT
		mat[viz@axisMat[2,2],2] <- sinT
		
	}else{
		cat('Transition Fault.\n')
		mat[1,1] <- 1
		mat[2,2] <- 1
		
	}
	
	
	viz@from <- state$from
	viz@to <- state$to
	
	viz@mat <- mat
	
	
	## Debug output:	
#	cat(paste('\n~~~~~~~~~~~~~~~~~~~~~~~~~~\n',
#					'node: from:',state$from,', to:',state$to,'\n',
#					'from1:',ifrom[1], ', to1:', ito[1], ',  ',
#					viz@axisMat[1,1],viz@axisMat[1,2],
#					',   theta:',  round(state$percentage,3),'\n',
#					'from2:',ifrom[2], ', to2:', ito[2],',  ',
#					viz@axisMat[2,1],viz@axisMat[2,2],'\n'))
#	
#	
#	cat(paste(text,'\n\n'))
#	
#	cat(paste('transitionKind:',viz@transitionKind,'\n'))
#	
#	
#	cat('\n')
#	
#	cat(c(round(mat[,1],2),'\n'))
#	cat(c(round(mat[,2],2),'\n'))
#	
#	cat('\n')
#	
	
	return(viz)
	
}


ng_2d_xcoord <- function(viz,ngEnv) {
	## there are only two non zero elements per column
	ii <- which(viz@mat[,1]!=0)
	n <- length(ii)
	
	if (viz@scaled) {
		if(n == 1) {
			x <- ngEnv$scaledData[[viz@data]][,ii]
		}else if(n == 2) {
			x <- (ngEnv$scaledData[[viz@data]][,ii[1]]*viz@mat[ii[1],1]+ ngEnv$scaledData[[viz@data]][,ii[2]]*viz@mat[ii[2],1])					
		}else {
			stop("[2dAxis] neither a 3d rotation nor a 4d transition")
		}
	} else {
		if(n == 1) {
			x <- ngEnv$dataList[[viz@data]]@data[,ii]
		}else if(n == 2) {
			x <- (ngEnv$dataList[[viz@data]]@data[,ii[1]]*viz@mat[ii[1],1]+ ngEnv$dataList[[viz@data]]@data[,ii[2]]*viz@mat[ii[2],1])					
		}else {
			stop("[2dAxis] neither a 3d rotation nor a 4d transition")
		}
	}
	return(x)
}


ng_2d_ycoord <- function(viz,ngEnv) {
	ii <- which(viz@mat[,2]!=0)
	n <- length(ii)
	if (viz@scaled) {
		if(n == 1) {
			y <- ngEnv$scaledData[[viz@data]][,ii]
		}else if(n == 2) {
			y <- (ngEnv$scaledData[[viz@data]][,ii[1]]*viz@mat[ii[1],2]+ ngEnv$scaledData[[viz@data]][,ii[2]]*viz@mat[ii[2],2])					
		}else {
			stop("[2dAxis] not a 3 or 4 dimensional rotation")
		}	
	} else {
		if(n == 1) {
			y <- ngEnv$dataList[[viz@data]]@data[,ii]
		}else if(n == 2) {
			y <- (ngEnv$dataList[[viz@data]]@data[,ii[1]]*viz@mat[ii[1],2]+ ngEnv$dataList[[viz@data]]@data[,ii[2]]*viz@mat[ii[2],2])					
		}else {
			stop("[2dAxis] not a 3 or 4 dimensional rotation")
		}	 
	}
	return(y)
}


## ng_2d_distance: distance of point to plane for 3d rotation
## TODO: cases phi = 0 and phi = pi/2
ng_2d_dist <- function(viz,ngEnv) {
	ii_x <- which(viz@mat[,1]!=0)
	n_x <- length(ii_x)
	
	ii_y <- which(viz@mat[,2]!=0)
	n_y <- length(ii_y)
	
	## if necessary, norm with 	sqrt((sin/cos)^2 +1)
	if((n_x == 2) && (n_y == 1)) {
		if (abs(viz@mat[ii_x[1],1] - cos(ngEnv$bulletState$percentage*pi/2)) < 0.000001) {
			d <- (ngEnv$scaledData[[viz@data]][,ii_x[1]]*viz@mat[ii_x[2],1]/viz@mat[ii_x[1],1] - ngEnv$scaledData[[viz@data]][,ii_x[2]])
		} else {
			d <- (ngEnv$scaledData[[viz@data]][,ii_x[2]]*viz@mat[ii_x[1],1]/viz@mat[ii_x[2],1] - ngEnv$scaledData[[viz@data]][,ii_x[1]])
		}
	} else if ((n_x == 1) && (n_y == 2)) {
		if (abs(viz@mat[ii_y[1],2] - cos(ngEnv$bulletState$percentage*pi/2)) < 0.000001) {
			d <- (ngEnv$scaledData[[viz@data]][,ii_y[1]]*viz@mat[ii_y[2],2]/viz@mat[ii_y[1],2] - ngEnv$scaledData[[viz@data]][,ii_y[2]])
		} else {
			d <- (ngEnv$scaledData[[viz@data]][,ii_y[2]]*viz@mat[ii_y[1],2]/viz@mat[ii_y[2],2] - ngEnv$scaledData[[viz@data]][,ii_y[1]])
		}
	} else if ((n_x == 1) && (n_y == 1)) {
		d = "scatterplot"
	}else {
		d = NULL
	}
	return(d)
}



