############ S-function: dup.matrix  ############

# For determining whether any rows in a matrix
# are duplicated. This is from the FUNFITS module.

# Last changed: 07 JUL 2001

dup.matrix <- function(mat)
{
    # Convert mat to a matrix

    mat <- as.matrix(mat)
    nc <- ncol(mat)
    temp <- matrix(match(c(mat),unique(c(mat))),ncol = nc)
    temp2 <- format(temp[, 1])
    if(nc > 1) 
       for(k in 2:nc) 
          temp2 <- paste(temp2,temp[,k],sep="X")
            
   return(dup(temp2))
}

########### End of dup.matrix #############
