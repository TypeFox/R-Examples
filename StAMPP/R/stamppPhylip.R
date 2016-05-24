#######################################################
#
# Convert genetic distance matrix to PHYLIP text file
# 
# Luke Pembleton
# luke.pembleton@ecodev.vic.gov.au
#
#######################################################

stamppPhylip <-
function(distance.mat, file=""){

  rnames <- row.names(distance.mat) #individual or population names from Neis genetic distance matrix
  
  for(x in 1:(length(rnames))){ #if name is longer than 11 characters, truncate to 11 characters
    
    if(nchar(rnames[x])>11) rnames[x]=substr(rnames[x], 1, 11)
    
  }
  
  max.id <- max(nchar(rnames)) #length of the longest name
  
  for(i in 1:(length(rnames))){ #add spaces to the end of each name so that all names are the length of the longest name plus 4 spaces
    
    short <- (max.id+4)-nchar(rnames[i])
    
    for(c in 1:short){
      
      rnames[i] <- paste(rnames[i], "")
      
    }
    
  }
    
  rnames <- c(paste("  ",length(distance.mat[1,])), rnames) #insert the number of individuals or populations at the top of the PHYLIP file
  
  phylip.mat <- rbind("", distance.mat)
  
  row.names(phylip.mat) <- rnames #add the fomrated names to the PHYLIP file
  
  for (i in 1:(length(distance.mat[1,]))){
    
    c=i+1
    
    vec <- distance.mat[i,]
    
    vec <- sprintf("%.5f", vec) #format all genetic distances to 5 decimal places
    
    phylip.mat[c,]=vec
    
  }
  
  phylip.mat[phylip.mat=="-0.00000"]="0.00000"
  
  write.table(phylip.mat, file=file, sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE) #save the PHYLIP format distane matrix as a text file
  
}
