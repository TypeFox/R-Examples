## modified 3 Feb 2013
## this version allows for a more general form of inputs

ibdhap.make.calls <-function( qibd.filename = NULL, ids.filename = NULL,
                              qibd.file = NULL, ids.file = NULL, cutoff = .8){

  ## information form qibd file
  if( !is.null(qibd.file)){ ## specified the qibd file already loaded
    
    rawdat <- qibd.file
    print("qibd matrix accepted")
    n.col.rawdat <- apply( rawdat, 1, length )
    
  }else if(is.null(qibd.file)){ ## specified a file location to read
    
    n.col.rawdat <- count.fields(qibd.filename) ## number of coulmns of each row in rawdat file
    rawdat <- read.table(qibd.filename, fill =TRUE, colClasses = "numeric", col.names=1:max(n.col.rawdat))
    print("qibd matrix read from file")
    
  }else{ stop("Invalid qibd input") }
  
    
  ## information from ids file
  if( !is.null(ids.file)){ ## specified the ids file already loaded

    par.file <- ids.file
    print("ids matrix accepted")
    
  }else if(is.null(ids.file)){ ## specified a file location to read

    n.col.par <- max(count.fields(ids.filename)) ## max number of columns in the input file
    par.file <- read.table(ids.filename, fill =TRUE, colClasses = "numeric", col.names=1:n.col.par)
    print("ids matrix read from file")

  }else{ stop("Invalid ids input") }

  
  name.sets <- par.file[,1] ## names of the sets
  n.snps <- dim( rawdat[name.sets==name.sets[1],] )[1] ## number of snps
    
    
  ##geno.indicator, n.names etc were defined in original but in the
  ##end not used.  intended to get the names of the
  ##genotypes/haplotypes.
  
  ## empty matrix that will be the output (assuming all sets have the same number of markers
  ibd.states <- matrix(NA, ncol=length(name.sets), nrow=n.snps)
  colnames(ibd.states) <- name.sets
  
  
  count <- 1
  for( iset in name.sets){
    
    ## subset just the data for the markers in the set.
    xind <- which( rawdat[,1]==iset)
    n.col.i <- n.col.rawdat[xind][1] ## number of columns in the set i
    ibd.dat <- t(rawdat[xind, 2:n.col.i]) ## remove the set name (aready subsetted by set name)
    
    ## truncate based on "cutoff" : put 1 if cutoff is met, 0 otherwise
    ibd.trunc <- apply(ibd.dat[ 3:(n.col.i-1), ], c(1,2), meet.cutoff, cutoff)

    ## for each column (state probs for a marker) list ibd.state, or 0 if no state > cutoff
    ibd.states.temp<-apply(ibd.trunc, 2, return.ibd.val)

    ## add the marker names -- was not implemented in original

    ## save in the final matrix
    ibd.states[,count]=ibd.states.temp
    ##print(paste("pair", iset, "of", n.sets, "complete...", sep=" "))
    
    count <- count + 1
  }
  return(as.data.frame(ibd.states))   
}
