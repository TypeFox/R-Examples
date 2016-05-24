broadPDF <- function(pdfob, sigma=0,  delta=NA, n=NA, nAtomTypes=1){
 
 
################################################
# function body
  
  if(is.na(delta))
    delta <- 0 
  
  if(delta==0)
    n <- 0 
	

  if (delta != 0 & is.na(n))
    n <- 2
  
  r <- pdfob$r
  dr <- r[2]-r[1]
  
  if(!pdfob$facs) { # uniform 
    gr <- pdfob$gr 
    nz <- which(gr != 0)
		 		 
    np <- .C("broadPDF",
            res = as.double(rep(0, length(r))),
            pdf = as.double(gr[nz]),
            r = as.double(r),
            rp = as.double(r[nz]),
            len = as.integer(length(r)),
            lenp = as.integer(length(nz)),
            sigi = as.double(sigma[1]),
			sigj = as.double(sigma[1]),
            delta = as.double(delta),
            n = as.integer(n)
          )$res
    
    np[which(is.na(np))] <- 0  ## fixes 0 issue
  
    pdfob$gr <- np  
  }
  else { # core-shell
    gr <- pdfob$gr
    grCCSS <- pdfob$gr_CCSS 
    grCS <- pdfob$gr_CS 
#	a1<-system.time({
	npCCSS <- list()
	for(i in 1:nAtomTypes){
	  nzCCSS <- which(grCCSS[[i]] != 0)
	  if(length(nzCCSS)>0) 
		npCCSS[[i]] <- .C("broadPDF",
			res = as.double(rep(0, length(r))),
			pdf = as.double(grCCSS[[i]][nzCCSS]),
			r = as.double(r),
			rp = as.double(r[nzCCSS]),
			len = as.integer(length(r)),
			lenp = as.integer(length(nzCCSS)),
			sigi = as.double(sigma[i]),
			sigj = as.double(sigma[i]),
			delta = as.double(delta),
			n = as.integer(n)
			)$res   					   
	  else
		npCCSS[[i]] <- rep(0, length(gr))
	}
#	})
#	cat("CCSS", a1, "\n")
	
	
	npCS <- list()
	
#	a2<-system.time({
	for(i in 2:nAtomTypes){
	  for(j in 1:(i-1)){
		k <- ceiling(i*(i-3)/2+j+1)
		nzCS <- which(grCS[[k]] != 0)
		if(length(nzCS)>0) 
		  npCS[[k]] <- .C("broadPDF",
			res = as.double(rep(0, length(r))),
			pdf = as.double(grCS[[k]][nzCS]),
			r = as.double(r),
			rp = as.double(r[nzCS]),
			len = as.integer(length(r)),
			lenp = as.integer(length(nzCS)),
			sigi = as.double(sigma[j]),
			sigj = as.double(sigma[i]),
			delta = as.double(delta),
			n = as.integer(n)
			)$res  
		 else 
		   npCS[[k]] <- rep(0, length(gr))	
	  }
	}
#	})
#	cat("CS", a2, "\n")

 #   a3<-system.time({
	for(i in 1:nAtomTypes){
	  npCCSS[[i]][which(is.na(npCCSS[[i]]))] <- 0
	}
    for(i in 2:nAtomTypes){
	  for(j in 1:(i-1)){
	    k <- ceiling(i*(i-3)/2+j+1)
	    npCS[[k]][which(is.na(npCS[[k]]))] <- 0
	  }
	}
	
	pdfob$gr  <- rep(0,length(r))
	for(i in 1:nAtomTypes){
	  pdfob$gr <- pdfob$gr + npCCSS[[i]]
	}
	
	for(i in 2:nAtomTypes){
      for(j in 1:(i-1)){
	    k <- ceiling(i*(i-3)/2+j+1)	    
		pdfob$gr <- pdfob$gr + npCS[[k]]
	  }
	}
	
    pdfob$gr_CCSS <- npCCSS 
    pdfob$gr_CS <- npCS
	
#	})
#	cat("the rest", a3, "\n")

  }
  
  pdfob
  
  
}
