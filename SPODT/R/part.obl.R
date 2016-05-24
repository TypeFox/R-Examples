part.obl <-
function(data, data.sp, perm, min.fils, ponderer)
{
    n.loc <- nrow(data.sp)
    x <- data$x[match(data.sp$loc, data$loc)]
    y <- data$y[match(data.sp$loc, data$loc)]
    ind <- ind.crsp(data.sp$loc, perm$loc1, perm$loc2)
    
    coeff <- rep(0,2)
    vic <- 0
    partition <- rep(0, n.loc)
    
    res <- .C("partObl",
              as.integer(min.fils), as.integer(ponderer),
              as.integer(n.loc), as.double(x), as.double(y), as.double(data.sp$z), as.integer(data.sp$n),
              as.double(cumsum(data.sp$z)), as.integer(cumsum(data.sp$n)), as.integer(length(which(duplicated(perm$pente)==TRUE))+1),
              as.double(perm$pente), as.integer(ind$i1), as.integer(ind$i2),
              as.double(coeff), as.double(vic), as.integer(partition)
              
             ) #ajout ,PACKAGE="SPODT"
             
    return(list(coeff=res[[14]], vic=res[[15]], part=res[[16]]))
}
