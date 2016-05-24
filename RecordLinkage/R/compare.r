# compare.r: Functions to group records to pairs.

# helper function, for internal use
makeBlockingPairs <- function(id_vec)
{
 # if list is empty, return empty matrix
 if (length(id_vec)==0)
    return (matrix(nrow=0, ncol=2))
 nPairs=sum(sapply(id_vec, function(x) length(x)* (length(x)-1) * 0.5))
 ret <- .C("makeBlockingPairs",id_vec, length(id_vec), 
   pairs=matrix(0L,nrow=nPairs, ncol=2), as.integer(nPairs),
   PACKAGE="RecordLinkage")
 return(ret$pairs)
}

compare.dedup <- function(dataset, blockfld=FALSE, phonetic=FALSE,
                    phonfun=pho_h, strcmp=FALSE,strcmpfun=jarowinkler, exclude=FALSE, 
                    identity=NA, n_match=NA, n_non_match=NA)
{
    # various catching of erronous input
    if (!is.data.frame(dataset) && !is.matrix(dataset))
        stop ("Illegal format of dataset")
    ndata=nrow(dataset) # number of records
    nfields=ncol(dataset)
    if (ndata<2) 
        stop ("dataset must contain at least two records")
    if (is.character(strcmp))
      strcmp <- match(strcmp, colnames(dataset))
    if (!is.numeric(strcmp) && !is.logical(strcmp))
        stop ("strcmp must be numeric or a single logical value")
    if (is.character(phonetic))
      phonetic <- match(phonetic, colnames(dataset))
    if (!is.numeric(phonetic) && !is.logical(phonetic))
        stop ("phonetic must be numeric or a single logical value")
    if (is.character(exclude))
      exclude <- match(exclude, colnames(dataset))
    if (!is.numeric(exclude) && !isFALSE(exclude))
        stop ("exclude must be numeric or FALSE")
    if (!isFALSE(strcmp) && any(is.na(strcmp) | strcmp <= 0 | strcmp > nfields))
        stop ("phonetic contains out of bounds index")
    if (!isFALSE(phonetic) && any(is.na(phonetic) | phonetic <= 0 | phonetic > nfields))
        stop ("phonetic contains out of bounds index")
    if (!isFALSE(exclude) && any(is.na(exclude) | exclude <= 0 | exclude > nfields))
        stop ("phonetic contains out of bounds index")
    if (!is.na(n_match) && !is.numeric(n_match))
      stop ("Illegal type for n_match!")
    if (!is.na(n_non_match) && !is.numeric(n_non_match))
      stop ("Illegal type for n_match!")

    if(!identical(blockfld, FALSE))
    {
      if (!is.list(blockfld) && !is.null(blockfld)) blockfld <- list(blockfld)
      if (!all(sapply(blockfld, function(x) class(x) %in% c("character", "integer", "numeric"))))
        stop("blockfld has wrong format!")
      blockfld <- lapply(blockfld, 
       function(x) {if (is.character(x)) match(x, colnames(dataset)) else (x)})
      if(any(unlist(blockfld) <= 0 | unlist(blockfld) > nfields))
        stop("blockfld countains out-of-bounds value!")
    }

    if(!identical(identity,NA))
    {
      if(length(identity)!=nrow(dataset))
      {
        stop("Length of identity vector does not match number of records!")
      }
    }

    # ensure dataset is a data frame and has column names
    ret=list()  # return object
    ret$data=as.data.frame(dataset)
    dataset=as.matrix(dataset, rownames.force=FALSE)
    dataset[dataset==""]=NA # label missing values



    # keep phonetics for blocking fields
    if (is.numeric(phonetic))
    {
        phonetic_block=intersect(phonetic,unlist(blockfld))
    }
    if (is.numeric(exclude))
    {        
        # adjust indices to list of included fields
        if (is.numeric(phonetic)) 
        {
            phonetic=setdiff(phonetic,exclude)
            phonetic=sapply(phonetic,function(x) return (x-length(which(exclude<x))))
        }
        if (is.numeric(strcmp))
        {
            strcmp=setdiff(strcmp,exclude)
            strcmp=sapply(strcmp,function(x) return (x-length(which(exclude<x))))       
        }
    }
    # issue a warning if both phonetics and string metric are used on one field
    if ((length(intersect(phonetic,strcmp))>0 && !isFALSE(strcmp) && !isFALSE(phonetic)) ||
         (isTRUE(strcmp) && !isFALSE(phonetic)) ||
         (isTRUE(phonetic) && !isFALSE(strcmp)))
    {
        warning(sprintf("Both phonetics and string metric are used on some fields",length(intersect(phonetic,strcmp))))
    }

    if (!is.function(phonfun))
    {
      stop("phonfun is not a function!")
    }

    
    if (!is.function(strcmpfun))
    {
      stop("strcmpfun is not a function!")
    }

        
# print("blocking beginnt")
    # Pair_ids collects ids of record pairs. It is a matrix because the following
    # rbind() calls are much faster than with a data.frame
   pair_ids=matrix(as.integer(0),nrow=0,ncol=2) # each row holds indices of one record pair
   if (isFALSE(blockfld))
   {
    if (is.na(n_match) || is.na(n_non_match))
   	{
	   pair_ids=t(unorderedPairs(ndata))
	 } else
	 {
		tempdat=data.frame(id=1:ndata,identity=identity)
		
		# Determine matches by join with identity vector
		pairs=merge(x=tempdat,y=tempdat,by=2)
		# limit to unordered pairs 
		match_ids=as.matrix(pairs[as.integer(pairs[,2])<as.integer(pairs[,3]),2:3], rownames.force=FALSE)
    n_available_matches=nrow(match_ids)
    n_available_non_matches=ndata*(ndata-1)/2 - n_available_matches
    if (n_available_matches < n_match && n_available_non_matches < n_non_match)
    {
      warning(sprintf("Only %d matches and %d non-matches!",
        n_available_matches, n_available_non_matches))
      # return all pairs
      pair_ids=t(unorderedPairs(ndata))
    } else 
    {
  		# draw required number
  		if (n_match>n_available_matches)
  		{
  			warning(sprintf("Only %d matches!",n_available_matches))
  			# in this case no sampling from match_ids
  		} else
  		{
  			s=sample(nrow(match_ids),n_match)
  			match_ids=match_ids[s,]
  		}

  		if (n_non_match > n_available_non_matches)
  		{
  			warning(sprintf("Only %d non-matches!",n_available_non_matches))
        all_pairs=t(unorderedPairs(ndata))
        is_match=(identity[all_pairs[,1]]==identity[all_pairs[,2]])
        non_match_ids=all_pairs[!is_match,]
        pair_ids=rbind(match_ids,non_match_ids)
      } else
     {
   	   # pairs already drawn are marked in A
   	   # drawing unordered pairs follows: Steven S. Skiena: The algorithm
   	   # design manual, 250-251.
  		A=list()
  		for (i in 1:n_non_match)
  		{
  			repeat
  		 	{
  		  		d1=sample(ndata,1)
  		  		d2=sample(ndata,1)
  				# If a match has been drawn, try again
  				if (identical(identity[d1],identity[d2]))
  				{
  				 	next
  				}
 		  		if (d1==d2)
 		  		  next
	  		  if (d1 > d2)
 		      {
 		        sorted_id=c(d2,d1)
	        } else
	        {
	          sorted_id=c(d1,d2)
           }
  				# If the pair has been drawn already, try again 
  				if (!is.null(A[[paste(sorted_id,collapse=" ")]]))
  			  	{
  			    	next
  			  	}
  			  	# mark pair as already drawn
  			  	A[[paste(sorted_id,collapse=" ")]]=sorted_id
  			  	break
  			}
  		}
   		non_match_ids=matrix(unlist(A),ncol=2,nrow=n_non_match,byrow=TRUE)
  		pair_ids=rbind(match_ids,non_match_ids)
  	 	rm(match_ids,non_match_ids,A)
     }
   }
	 }
   } else # with blocking
   { 
     for (blockelem in blockfld) # loop over blocking definitions
     {
      if (isTRUE(phonetic))
      {
        	block_data=phonfun(dataset)
      } else if (is.numeric(phonetic))
      {
        block_data=dataset
        block_data[,phonetic_block]=phonfun(dataset[,phonetic_block])
      } else
      {
        block_data=dataset
      }
      # for each record, concatenate values in blocking fields
      # do.call is faster than a former apply solution
      blockstr <- do.call(paste, as.data.frame(block_data[,blockelem]))
	  # exclude pairs with NA in blocking variable
	  # (paste just converts to "NA")
      for (i in blockelem)
      {
        is.na(blockstr)=is.na(block_data[,i])
      }
      rm(block_data)
     id_vec=tapply(1:ndata,blockstr,function(x) if(length(x)>1) return(x))
     id_vec=deleteNULLs(id_vec)
     id_vec=makeBlockingPairs(id_vec)
      rm(blockstr)
      # reshape vector and attach to matrix of record pairs
      if (nrow(id_vec)>0)
       pair_ids=rbind(pair_ids,id_vec)
      rm(id_vec)
    }
    # return empty data frame if no pairs are obtained
    if (length(pair_ids)==0)
    {
      stop("No pairs generated. Check blocking criteria.")
    }
    pair_ids <- data.table(pair_ids)
    setkeyv(pair_ids, names(pair_ids))
#    browser()
    pair_ids <- as.matrix(pair_ids[,1, by=names(pair_ids)], rownames.force=FALSE)[,1:2, drop=FALSE]
#    pair_ids=as.matrix(unique(as.data.frame(pair_ids)), rownames.force=FALSE)  # runs faster with data frame
   } # end else
    

    # remove excluded fields
    if (is.numeric(exclude))
    {        
        dataset=dataset[,-exclude, drop = FALSE]  # remove excluded columns
    }  

    # apply phonetic code                                                         
    if (!isFALSE(phonetic)) # true, if phonetic is TRUE or not a logical value
    {
        if (isTRUE(phonetic)) # true, if phonetic is a logical value and TRUE
        {    
            dataset=pho_h(dataset)
        } else # phonetic is not a logical value
        dataset[,phonetic]=pho_h(dataset[,phonetic])
    }

    left <- dataset[pair_ids[,1],,drop=FALSE]
    right <- dataset[pair_ids[,2],,drop=FALSE]
    # matrix to hold comparison patterns
    patterns=matrix(0,ncol=ncol(left),nrow=nrow(left)) 
    if (isTRUE(strcmp))
    {
        patterns=strcmpfun(as.matrix(left, rownames.force=FALSE),as.matrix(right, rownames.force=FALSE))
    } else if (is.numeric(strcmp)) 
    {
        patterns[,-strcmp]=(left[,-strcmp]==right[,-strcmp])*1
        patterns[,strcmp]=strcmpfun(left[,strcmp],right[,strcmp]) #*1
    } else
    {
       patterns=(left==right)*1
    }
    rm(left)
    rm(right)

    is_match=identity[pair_ids[,1]]==identity[pair_ids[,2]] # Matching status of pairs
    ret$pairs=as.data.frame(
					cbind(pair_ids,
                    patterns,
                    is_match)) # Matches

    if (is.numeric(exclude))
    {
     colnames(ret$pairs)=c("id1","id2",colnames(ret$data)[-exclude],"is_match")
    } else
    {
     colnames(ret$pairs)=c("id1","id2",colnames(ret$data),"is_match")
    }
    rownames(ret$pairs)=NULL

    ret$frequencies=apply(dataset,2,function(x) 1/length(unique(x)))
    ret$type="deduplication"
    class(ret)="RecLinkData"
    return(ret)
}


# Version for two data sets

# Requires that both have the same format

compare.linkage <- function(dataset1, dataset2, blockfld=FALSE, phonetic=FALSE,
                    phonfun=pho_h, strcmp=FALSE,strcmpfun=jarowinkler, exclude=FALSE, 
                    identity1=NA, identity2=NA, n_match=NA, n_non_match=NA)
{
    # various catching of erronous input
    if (!is.data.frame(dataset1) && !is.matrix(dataset1))
        stop ("Illegal format of dataset1")
    if (!is.data.frame(dataset2) && !is.matrix(dataset2))
        stop ("Illegal format of dataset2")
    if (ncol(dataset1) != ncol(dataset2))
        stop ("Data sets have different format")
    ndata1=nrow(dataset1) # number of records
    ndata2=nrow(dataset2)
    nfields=ncol(dataset1)
    if (ndata1<1 || ndata2<1) 
        stop ("empty data set")

    if (is.character(strcmp))
      strcmp <- match(strcmp, colnames(dataset1))
    if (!is.numeric(strcmp) && !is.logical(strcmp))
        stop ("strcmp must be numeric, character or a single logical value")
    if (!isFALSE(strcmp) && any(is.na(strcmp) | strcmp <= 0 | strcmp > nfields))
        stop ("strcmp contains out of bounds index")

    if (is.character(phonetic))
      phonetic <- match(phonetic, colnames(dataset1))
    if (!is.numeric(phonetic) && !is.logical(phonetic))
        stop ("phonetic must be numeric, character or a single logical value")
    if (!isFALSE(phonetic) && any(is.na(phonetic) | phonetic <= 0 | phonetic > nfields))
        stop ("phonetic contains out of bounds index")

    if (is.character(exclude))
      exclude <- match(exclude, colnames(dataset1))
    if (!is.numeric(exclude) && !is.logical(exclude))
        stop ("exclude must be numeric, character or a single logical value")
    if (!isFALSE(exclude) && any(is.na(exclude) | exclude <= 0 | exclude > nfields))
        stop ("exclude contains out of bounds index")

    if (!is.na(n_match) && !is.numeric(n_match))
      stop ("Illegal type for n_match!")
    if (!is.na(n_non_match) && !is.numeric(n_non_match))
      stop ("Illegal type for n_match!")

    if(!identical(blockfld, FALSE))
    {
      if (!is.list(blockfld) && !is.null(blockfld)) blockfld <- list(blockfld)
      if (!all(sapply(blockfld, function(x) class(x) %in% c("character", "integer", "numeric"))))
        stop("blockfld has wrong format!")
      blockfld <- lapply(blockfld, 
       function(x) {if (is.character(x)) match(x, colnames(dataset1)) else (x)})
      if(any(unlist(blockfld) <= 0 | unlist(blockfld) > nfields))
        stop("blockfld countains out-of-bounds value!")
    }

    if(!identical(identity1,NA))
    {
      if(length(identity1)!=nrow(dataset1))
      {
        stop("Length of identity1 does not match number of records!")
      }
    }

    if(!identical(identity2,NA))
    {
      if(length(identity2)!=nrow(dataset2))
      {
        stop("Length of identity2 does not match number of records!")
      }
    }

    dataset1=as.data.frame(dataset1)
    dataset2=as.data.frame(dataset2)
    ret=list()  # return object
    ret$data1=dataset1
    ret$data2=dataset2
    full_data1=as.matrix(dataset1, rownames.force=FALSE)
    full_data2=as.matrix(dataset2, rownames.force=FALSE)



    # keep phonetics for blocking fields
    if (is.numeric(phonetic))
    {
        phonetic_block=intersect(phonetic,unlist(blockfld))
    }
    if (is.numeric(exclude))
    {        
        dataset1=dataset1[,-exclude, drop = FALSE]  # remove excluded columns
        dataset2=dataset2[,-exclude, drop = FALSE]  # remove excluded columns
        # adjust indices to list of included fields
        if (is.numeric(phonetic)) 
        {
            phonetic=setdiff(phonetic,exclude)
            phonetic=sapply(phonetic,function(x) return (x-length(which(exclude<x))))
        }
        if (is.numeric(strcmp))
        {
            strcmp=setdiff(strcmp,exclude)
            strcmp=sapply(strcmp,function(x) return (x-length(which(exclude<x))))       
        }
    }
    # issue a warning if both phonetics and string metric are used on one field
    if ((length(intersect(phonetic,strcmp))>0 && !isFALSE(strcmp) && !isFALSE(phonetic)) ||
         (isTRUE(strcmp) && !isFALSE(phonetic)) ||
         (isTRUE(phonetic) && !isFALSE(strcmp)))
    {
        warning(sprintf("Both phonetics and string metric are used on some fields",length(intersect(phonetic,strcmp))))
    }
    dataset1[dataset1==""]=NA # label missing values
    dataset2[dataset2==""]=NA # label missing values
    full_data1[full_data1==""]=NA # label missing values
    full_data2[full_data2==""]=NA # label missing values
    dataset1=as.matrix(dataset1, rownames.force=FALSE)
    dataset2=as.matrix(dataset2, rownames.force=FALSE)

    if (!is.function(phonfun))
    {
      stop("phonfun is not a function!")
    }

    if (!isFALSE(phonetic)) # true, if phonetic is TRUE or not a logical value
    {
        if (isTRUE(phonetic)) # true, if phonetic is a logical value and TRUE
        {    
            dataset1=pho_h(dataset1)
            dataset2=pho_h(dataset2)
        } else # phonetic is not a logical value
        dataset1[,phonetic]=pho_h(dataset1[,phonetic])
        dataset2[,phonetic]=pho_h(dataset2[,phonetic])
    }
    
    if (!is.function(strcmpfun))
    {
      stop("strcmpfun is not a function!")
    }

    # Pair_ids collects ids of record pairs. It is a matrix because the following
    # rbind() calls are much faster than with a data.frame
   pair_ids=matrix(as.integer(0),nrow=0,ncol=2) # each row holds indices of one record pair
   if (isFALSE(blockfld))
   { 
     if (is.na(n_match) || is.na(n_non_match))
   	 {
        # full outer join
    	 pair_ids=merge(1:nrow(dataset1),1:nrow(dataset2),all=TRUE)
    	 # sort to enforce particular order
       pair_ids=pair_ids[order(pair_ids[,1],pair_ids[,2]),]
	}   else
	 {
		tempdat1=data.frame(id=1:ndata1,identity=identity1)
		tempdat2=data.frame(id=1:ndata2,identity=identity2)
		
		# Determine matches by join on identity vector 
		pairs=merge(x=tempdat1,y=tempdat2,by=2)
		match_ids=as.matrix(pairs[,2:3], rownames.force=FALSE)
    n_available_matches=nrow(match_ids)
    n_available_non_matches=ndata1*ndata2 - n_available_matches
    if (n_available_matches < n_match && n_available_non_matches < n_non_match)
    {
      warning(sprintf("Only %d matches and %d non-matches!",
        n_available_matches, n_available_non_matches))
      # return all pairs
      pair_ids=merge(1:nrow(dataset1),1:nrow(dataset2),all=TRUE)
    } else 
    {
  		# draw required number
  		if (n_match>n_available_matches)
  		{
  			warning(sprintf("Only %d matches!",n_available_matches))
  			# in this case no sampling from match_ids
  		} else
  		{
  			s=sample(nrow(match_ids),n_match)
  			match_ids=match_ids[s,]
  		}

  		if (n_non_match > n_available_non_matches)
  		{
  			warning(sprintf("Only %d non-matches!",n_available_non_matches))
        all_pairs=merge(1:nrow(dataset1),1:nrow(dataset2),all=TRUE)
        is_match=(identity1[all_pairs[,1]]==identity2[all_pairs[,2]])
        non_match_ids=all_pairs[!is_match,]
        pair_ids=rbind(match_ids,non_match_ids)
      } else
     {
   	   # mark pairs already drawn in A
  		A=list()
      for (i in 1:n_non_match)
      {
        repeat
        {
          d1=sample(ndata1,1)
          d2=sample(ndata2,1)
          # If a match has been drawn, try again
          if (identical(identity1[d1],identity2[d2]))
          {
            next
          }
          # If the pairs has been drawn already, try again
          if (!is.null(A[[paste(d1,d2)]]))
          {
            next
          }
          # Mark the pair as drawn
          A[[paste(d1,d2)]]=c(d1,d2)
          break
        }
      }
      non_match_ids=matrix(unlist(A),ncol=2,nrow=n_non_match,byrow=TRUE)
      pair_ids=rbind(match_ids,non_match_ids)
      rm(match_ids,non_match_ids,A)
     }
   }
	 }
   } else  # branch for blocking
   {
    if (!is.list(blockfld)) blockfld=list(blockfld)
    for (blockelem in blockfld) # loop over blocking definitions
    {
      if (isTRUE(phonetic))
      {
        block_data1=phonfun(full_data1)
        block_data2=phonfun(full_data2)
      } else if (is.numeric(phonetic))
      {
        block_data1=full_data1
        block_data1[,phonetic_block]=phonfun(full_data1[,phonetic_block])
        block_data2=full_data2
        block_data2[,phonetic_block]=phonfun(full_data2[,phonetic_block])
      } else
      {
        block_data1=full_data1
        block_data2=full_data2
      }
      # for each record, concatenate values in blocking fields
      # do.call is faster than a former apply solution
      blockstr1 <- do.call(paste, as.data.frame(block_data1[,blockelem]))
      blockstr2 <- do.call(paste, as.data.frame(block_data2[,blockelem]))
  	  # exclude pairs with NA in blocking variable
      for (i in blockelem)
      {
        is.na(blockstr1)=is.na(block_data1[,i])
        is.na(blockstr2)=is.na(block_data2[,i])
      }
      rm(block_data1)
      rm(block_data2)

    id_vec=merge(data.frame(id1=1:ndata1,blockstr=blockstr1),
                   data.frame(id2=1:ndata2,blockstr=blockstr2),
                   incomparables=NA)[,-1]

      rm(blockstr1)
      rm(blockstr2)
      # reshape vector and attach to matrix of record pairs
      if (nrow(id_vec)>0)
        pair_ids=rbind(pair_ids,id_vec)
      rm(id_vec)
    }
  if (length(pair_ids)==0)
  {
      stop("No pairs generated. Check blocking criteria.")
  }
  
    pair_ids=unique(as.data.frame(pair_ids))  # runs faster with data frame
  } # end else
    
  rm(full_data1,full_data2)
    left=dataset1[pair_ids[,1],,drop=FALSE]
    right=dataset2[pair_ids[,2],,drop=FALSE]
    # matrix to hold comparison patterns
    patterns=matrix(0,ncol=ncol(left),nrow=nrow(left)) 
    if (isTRUE(strcmp))
    {
        patterns=strcmpfun(as.matrix(left, rownames.force=FALSE),as.matrix(right, rownames.force=FALSE))
    } else if (is.numeric(strcmp)) 
    {
        patterns[,-strcmp]=(as.matrix(left[,-strcmp], rownames.force=FALSE)==as.matrix(right[,-strcmp], rownames.force=FALSE))*1
        patterns[,strcmp]=strcmpfun(as.matrix(left[,strcmp], rownames.force=FALSE),
          as.matrix(right[,strcmp], rownames.force=FALSE)) #*1
    } else
    {
       patterns=(left==right)*1
    }
    rm(left)
    rm(right)

    is_match=as.numeric(identity1[pair_ids[,1]]==identity2[pair_ids[,2]]) # match status of pairs
    ret$pairs=as.data.frame(cbind(pair_ids, patterns, is_match)) # Matches


    colnames(ret$pairs)=c("id1","id2",colnames(dataset1),"is_match")
    rownames(ret$pairs)=NULL

    ret$frequencies=apply(rbind(dataset1,dataset2),2,
      function(x) 1/length(unique(x)))
    ret$type="linkage"
    class(ret)="RecLinkData"
    return(ret)
}

 
