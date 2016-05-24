pruneindicators<-function(x, At=0, Bt=0, sqrtIVt=0, max.indicators=4, verbose=FALSE) {

nonnested <- function (x, selection=NULL, verbose=FALSE) {
	if(is.null(selection)) selection = rep(TRUE, nrow(x$C))
	c = x$C[selection,]
	xc = x$XC[, selection]
	combs = row.names(c)
  	keep = rep(TRUE, ncol(xc))
  	for(c1 in 1:ncol(xc)) {
  		ccx1 = xc[,c1]>0
  	  	for(c2 in 1:ncol(xc)) { 	  	
	  		ccx2 = xc[,c2]>0
  			if(c1!=c2 && keep[c2]) {
  				if(sum(ccx1 & ccx2)==sum(ccx2) && sum(ccx1)>sum(ccx2)) {
  					keep[c2] = FALSE
  					if(verbose) cat(paste(combs[c2],"nested in",combs[c1],"\n"))
  				}
  				else if(sum(ccx1 & ccx2)==sum(ccx2) && sum(ccx1 & ccx2)==sum(ccx1)) {
  					if(verbose) cat(paste(combs[c2],"equal to",combs[c1],"\n"))
  					if(sum(c[c2,])> sum(c[c1,])) keep[c2] = FALSE
  					else if(sum(c[c1,])> sum(c[c2,])) keep[c1] = FALSE
  				}
  			}
      	}
 	}	
  	return(combs[keep])
}

	initCoverage<-coverage(x)
	if(verbose) cat(paste("Coverage of initial set of ",nrow(x$C)," indicators: ", round(initCoverage*100, digits=1),"%\n", sep=""))
	
	if(length(dim(x$A))==2) {
		selection<- x$A$lowerCI>=At & x$B$lowerCI>=Bt & x$sqrtIV$lowerCI>=sqrtIVt
	} else {
		selection<- x$A>=At & x$B>=Bt & x$sqrtIV>=sqrtIVt
	}
    if(sum(selection)==0) {
    	if(verbose) cat(paste("No indicator is valid using the given thresholds."))
    	return()
    }
	validCoverage<-coverage(x, selection=selection)
	if(verbose) cat(paste("Coverage of valid set of ",sum(selection)," indicators: ", round(validCoverage*100, digits=1),"%\n", sep=""))
	
	if(sum(selection)>1) {
	  NN <-nonnested(x, selection=selection, verbose=FALSE)
		selection <- row.names(x$C) %in% NN
		nnCoverage<-coverage(x, selection=selection)
		if(verbose) cat(paste("Coverage of valid set of ",sum(selection)," nonnested indicators: ", round(nnCoverage*100, digits=1),"%\n", sep=""))
	
	
		c = x$C[selection,]
		group.vec = x$group.vec
		xc = x$XC[, selection]

		#Preliminaries	
  		spnames = names(c)
  		spplist = names(c)

		indnames<-row.names(c)
		k <- length(indnames)
		j=1
		continue = (!is.null(max.indicators))
    	selmodFinal = selection
   		while(continue) {
      		co <- combn(k,j) #Generate subsets of indicators
   			if(verbose) cat(paste("Checking ",ncol(co)," subsets of ", j," indicator(s)", sep=""))
	      	keep2 = rep(FALSE,ncol(co))
   	   		maxcov = 0
	      	for(coi in 1:ncol(co)) { #Check coverage of subsets of indicators
	      		if(ncol(co)>10) if(coi%%round(ncol(co)/10)==0 && verbose) cat(".")
	      		selmod = selection
	      		selmod[selection]=FALSE
	      		selmod[selection][co[,coi]]=TRUE
	      		coicov = coverage(x, selection=selmod)
	      		if(coicov>maxcov) {
	      		   bestAtPoint= coi
	      		   maxcov = coicov
	      		}
		  		keep2[coi]=coicov== nnCoverage
	      	}
	        if(verbose) cat(paste(" maximum coverage: ", round(maxcov*100, digits=1),"%\n",sep=""))
	      	if(sum(keep2)>0) { #If at least one subset has the appropriate coverage keep it
	      		best = which(keep2)[1]
	      		selmodFinal[selection]=FALSE
	      		selmodFinal[selection][co[,best]]=TRUE
	      		finalCoverage<-coverage(x, selection=selmodFinal)
				if(verbose) cat(paste("Coverage of final set of ", j, " indicators: ",round(finalCoverage*100,digits=1),"%\n", sep=""))
	      		continue = FALSE
	      	} else {
		      	if(j==max.indicators) {
	    	  		selmodFinal[selection]=FALSE
		      		selmodFinal[selection][co[,bestAtPoint]]=TRUE
					if(verbose) cat(paste("\nCoverage maximum allowed set of ", j, " indicators: ",round(maxcov*100,digits=1),"%\n", sep=""))
		      	}
	      	}
	      	if(j<min(max.indicators,k)) j= j+1
  	    	else continue = FALSE
	    }
	} else {
		if(verbose) cat(paste("One valid indicator only. Stopping.\n", sep=""))
		selmodFinal = selection
	}
    indicators2 = x
    indicators2$C = as.data.frame(subset(indicators2$C, subset=selmodFinal))
    indicators2$XC = indicators2$XC[,selmodFinal]
    if(length(dim(indicators2$A))==2) {
    	indicators2$A = indicators2$A[selmodFinal,]
	    indicators2$B = indicators2$B[selmodFinal,]
    	indicators2$sqrtIV = indicators2$sqrtIV[selmodFinal,]
    } else {
    	indicators2$A = indicators2$A[selmodFinal]
	    indicators2$B = indicators2$B[selmodFinal]
    	indicators2$sqrtIV = indicators2$sqrtIV[selmodFinal]
    }

  	selSpp = colSums(indicators2$C)>0
  	indicators2$C = as.data.frame(indicators2$C[,selSpp])
  	row.names(indicators2$C)<-row.names(x$C)[selmodFinal]
  	names(indicators2$C)<-names(x$C)[selSpp]
  	
    #Select rows that contain the species or the group
    if(sum(selmodFinal)>1) {
	  	selRows = rowSums(indicators2$XC)>0 | indicators2$group.vec
	  	indicators2$group.vec = indicators2$group.vec[selRows]
	  	indicators2$XC = indicators2$XC[selRows,]
	} else {
	  	selRows = sum(indicators2$XC) | indicators2$group.vec
	  	indicators2$group.vec = indicators2$group.vec[selRows]
	  	indicators2$XC = indicators2$XC[selRows]
	}
    
	return(indicators2)
}



