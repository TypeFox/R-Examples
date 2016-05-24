# pid = min pid
# maxpass = max iteration nb
## ret=1 returns distance matrix
## ret=2 returns cost matrix
seqmathenikoff <- function(seqdata, with.miss=FALSE, full.matrix=TRUE, logoddmode=0) {
	n <- nrow(seqdata)
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize,
		" distinct events/states (", paste(alphabet,collapse="/"),")")

	## Gaps in sequences
	if (with.miss) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	else 
		if (any(seqdata==attr(seqdata,"nr")))
			stop("found missing values in sequences, please set 'with.miss' option to nevertheless compute distances")

	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.missing=with.miss)


#	dseq <- unique(seqdata)
	dseq <- seqdata
	mcorr <- match(seqconc(seqdata), seqconc(dseq))

	nd <- seqdim(dseq)[1]
	message(" [>] ", nd," distinct sequences")

	l <- ncol(dseq)
	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.missing=with.miss)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	debut <- Sys.time()

	message(" [>] computing new costs")

	## Function and arguments

 	
		retmat <- .Call("henikoff",
			as.integer(dseq),
			as.integer(dim(dseq)),
			as.integer(slength), 
			as.integer(alphsize),
			as.integer(logoddmode),
			PACKAGE="TraMineR")

	fin <- Sys.time()
	message(" (",round(difftime(fin,debut,units="mins"),2)," minutes)")

#	if (full.matrix && inherits(distances, "dist")) 
#		return(dist2matrix(distances))
#	else 
	
		rclab <- paste(alphabet,"->",sep="")
		
		
		if(with.miss)  { 
			retmat <- matrix(retmat, ncol=length(alphabet(seqdata))+1)
			dimnames(retmat) <- list(rclab,rclab)
			return(retmat)
		}
		
		retmat <- matrix(retmat, ncol=length(alphabet(seqdata)))
		dimnames(retmat) <- list(rclab,rclab)
		return(retmat)
	
}


seqmatsaltt <- function(seqdata, norm=FALSE, indel=2, sm=NULL,
	with.miss=FALSE, full.matrix=TRUE, optim=TRUE, pid=0.6, maxpass=100, ret=1, logoddmode=0) {



	## ======
	## CHECKS
	## ======
  if(is.null(sm)) {
	  sm <- seqsubm(seqdata, method="CONSTANT", cval=2)
  }	
  ## Taking care of correct normalization settings
  if(is.logical(norm)) {
    if (norm) {
        norm <- 1
      }
    
    else {
      norm <- 0
    }
  }
  else if (is.character(norm)) {

    normIndice <- match(norm, c("none", "maxlength", "gmean", "maxdist")) -1
   
    if (is.na(normIndice)) {  ##Not found
      stop("Unknow distance normalization method ", norm)
    }
    norm <- normIndice
    
  }else{
    stop("Unknow distance normalization method ", norm)
  }
    
	n <- nrow(seqdata)
	alphabet <- attr(seqdata,"alphabet")
	alphsize <- length(alphabet)
	message(" [>] ",n," sequences with ", alphsize,
		" distinct events/states (", paste(alphabet,collapse="/"),")")

	## Gaps in sequences
	if (with.miss) {
		alphabet <- c(alphabet,attr(seqdata,"nr"))
		alphsize <- length(alphabet)
		message(" [>] including missing value as additional state" )
	}
	else 
		if (any(seqdata==attr(seqdata,"nr")))
			stop("found missing values in sequences, please set 'with.miss' option to nevertheless compute distances")

	## ==============
	## Preparing data
	## ==============
	seqdata <- seqnum(seqdata, with.missing=with.miss)


#	dseq <- unique(seqdata)
	dseq <- seqdata
	mcorr <- match(seqconc(seqdata), seqconc(dseq))

	nd <- seqdim(dseq)[1]
	message(" [>] ", nd," distinct sequences")

	l <- ncol(dseq)
	slength <- seqlength(dseq)

	dseq <- seqasnum(dseq, with.missing=with.miss)

	message(" [>] min/max sequence length: ",min(slength),"/",max(slength))

	debut <- Sys.time()

	message(" [>] computing new costs")

	## Function and arguments

 	
		retmat <- .Call("saltt",
			as.integer(dseq),
			as.integer(dim(dseq)),
			as.integer(slength), 
			as.double(indel),
			as.integer(alphsize),
			as.double(sm),
			as.integer(norm),
			as.integer(optim),
			as.double(pid),
			as.integer(maxpass),
			as.integer(ret),
			as.integer(logoddmode),
			PACKAGE="TraMineR")

	fin <- Sys.time()
	message(" (",round(difftime(fin,debut,units="mins"),2)," minutes)")

#	if (full.matrix && inherits(distances, "dist")) 
#		return(dist2matrix(distances))
#	else 
	
	if(ret==1) {	
		return(matrix(retmat, ncol=dim(seqdata)[1]))
	}
	else {
		
		rclab <- paste(alphabet,"->",sep="")
		
		
		if(with.miss)  { 
			retmat <- matrix(retmat, ncol=length(alphabet(seqdata))+1)
			dimnames(retmat) <- list(rclab,rclab)
			return(retmat)
		}
		
		retmat <- matrix(retmat, ncol=length(alphabet(seqdata)))
		dimnames(retmat) <- list(rclab,rclab)
		return(retmat)
	}
}

