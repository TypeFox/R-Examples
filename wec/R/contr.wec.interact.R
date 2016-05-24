contr.wec.interact <-
function(x1,x2)
{
    
    # Set up the baseline model matrix
    mm <- model.matrix(~x1*x2)
    #Remove the intercept from the model matrix
    mm <- mm[,-1]
    
    # Separate main from interaction effects
    n.int <- (length(levels(x1)) -1) * (length(levels(x2)) -1)
    n.main <- (length(levels(x1)) -1) + (length(levels(x2)) -1)
    mm.int <- mm[,(dim(mm)[2]-n.int+1):dim(mm)[2]]
    mm.main <- mm[,1:n.main]
    
    # Create a reference table, with in the cells names of variables+labels+interactions
    # This helpful in finding the correct number of observations in cells, as identified by the names of the interaction-categories 
    ref.x1 <- paste(names(attr(table(x1,x2), "dimnames"))[1], names(table(x1)), sep="")
    ref.x2 <- paste(names(attr(table(x1,x2), "dimnames"))[2], names(table(x2)), sep="")
    
    
    reftable <- matrix(data=NA, nrow=length(table(x1)), ncol=length(table(x2)))
    for (i in 1:dim(reftable)[1]){
        for(j in 1:dim(reftable)[2]) {
            reftable[i,j] <- paste(ref.x1[i], ref.x2[j], sep=":")
        }
    }
    
    # Figure out what the reference category is
    refcat.x1 <- setdiff(rownames(contrasts(x1)), colnames(contrasts(x1)))
    refcat.x2 <- setdiff(rownames(contrasts(x2)), colnames(contrasts(x2)))
    # Number of observations in reference category
    refcat.n <- table(x1,x2)[refcat.x1, refcat.x2]
    
    
    # CONDITION 1:
    # If an observation has negative values on ALL dummies, the coding is not the product but a POSITIVE fraction.
    # It is the number of obs. in the cell for which the two product terms are indicators divided by the number of obs. in the omitted categories
    cond1 <- apply(mm.main < 0, 1, sum) == n.main  
    
    ## Re-calculate contrasts for observations meeting condition 1 
    
    values <- NA
    for (i in 1:ncol(mm.int)){
        values[i] <- table(x1,x2)[which(reftable==colnames(mm.int)[i], arr.ind=TRUE)] / refcat.n
    }
    
    ## Update model matrix for interactions with new values
    mm.int[cond1,] <- rep(values, each = sum(cond1))
    
    
    # CONDITION 2:
    # Whenever the product of two dummies is <0, a NEGATIVE fraction has to be inserted
    # It is the number of observations in the cell for which the two product terms are indicators divided by 
    # the number of observations in the cell to which the particular observation belongs
    
    cond2 <- which(mm.int<0, arr.ind=TRUE)
    
    for (i in 1:nrow(cond2)) {
        
        num <-    table(x1,x2)[which(reftable==colnames(mm.int)[cond2[i,2]])]
        denom <-  table(x1,x2)[x1[cond2[i,1]], x2[cond2[i,1]]]
        
        mm.int[cond2[i,1], cond2[i,2]] <- -1 * num / denom
        
    }
    
    return(mm.int)
    
}
