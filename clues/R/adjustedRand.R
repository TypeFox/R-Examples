# v0.2.4 on Feb. 3, 2009 by Weiliang Qiu
#  (1) moved some code out of the for loop
#
# cl1 --- partition 1 of the data set
# cl2 --- partition 2 of the data set
#
# flag = 1 --- Rand index
# flag = 2 --- Hubert and Arabie's adjusted Rand index
# flag = 3 --- Morey and Agresti's adjusted Rand index
# flag = 4 --- Fowlkes and Mallows's index
# flag = 5 --- Jaccard index
adjustedRand <- function(cl1, cl2,
    randMethod = c("Rand","HA", "MA", "FM", "Jaccard"))
{
    if(!is.vector(cl1)){
        stop("cl1 is not a vector!\n");
    }
    if(!is.vector(cl2)){
        stop("cl2 is not a vector!\n");
    }
    if(length(cl1) != length(cl2)){
        stop("two vectors have different lengths!\n");
    }
 
    len <- length(randMethod)
    if(len == 0)
    { stop("The argument 'randMethod' is empty!\n") }
 
    # unique values of elements in 'cl1'
    cl1u <- unique(cl1)
    # number of clusters in partition 1
    m1 <- length(cl1u)
    
    # unique values of elements in 'cl2'
    cl2u <- unique(cl2)
    # number of clusters in partition 2
    m2 <- length(cl2u)
  
    n <- length(cl1)
    randVec <- rep(0, len) 
    names(randVec) <- randMethod
    for(i in 1:len)
    {
        randMethod[i] <- match.arg(arg = randMethod[i], 
            choices = c("Rand","HA", "MA", "FM", "Jaccard"))
       
        flag <- match(randMethod[i], 
            c("Rand","HA", "MA", "FM", "Jaccard"))
     
        c.res <- .C("adjustedRand", 
            as.integer(cl1), 
            as.integer(cl1u), 
            as.integer(cl2), 
            as.integer(cl2u), 
            as.integer(m1), 
            as.integer(m2), 
            as.integer(n), 
            as.integer(flag),
            r = as.double(0)) 
        randVec[i] <- c.res$r
    }
    return(randVec)
}

