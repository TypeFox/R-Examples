`check.ordiscores` <- 
function(x, ord, check.species=TRUE) {
    sitescores <- scores(ord,display="sites")
    if(nrow(x)!=nrow(sitescores)){
        cat("Warning: community data set and ordination result have different number of sites\n")
    }else{
        if(any(rownames(x)!=rownames(sitescores))){
            cat("Warning: names for sites are different in community data set and ordination result\n")
        }
    }
    if (check.species==T){
        specscores <- scores(ord,display="species")
        if(ncol(x)!=nrow(specscores)){
            cat("Warning: community data set and ordination result have different number of species\n")
        }else{
            if(any(colnames(x)!=rownames(specscores))){
                cat("Warning: names for species are different in community data set and ordination result\n")
            }
        }
    }
}

