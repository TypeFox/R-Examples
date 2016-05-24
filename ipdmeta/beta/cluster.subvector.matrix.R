cluster.subvector.matrix <- function(
                         b,
                         n.groups)
{
		     #FUNCTION TO STACK EFFECTS IF MULTIPLE
    n = length(b)
    i = n/n.groups

    if(i==1){
    return(matrix(b,ncol=1))
        }
    else{

    B = matrix(0,n.groups,i)

    for(j in 1:i){
    B[,j] <- b[(n.groups*(j-1)+1):(n.groups*j)]
        }

    return(B)
    }
}

combined.subvector.matrix <- function(
                          b1,
                          b2,
                          n1,
                          n2,
                          as.matrix=TRUE)
{

	b1.matrix <- cluster.subvector.matrix(b1,n1)
	b2.matrix <- cluster.subvector.matrix(b2,n2)
	
	B <- rbind(b1.matrix,b2.matrix)

	if(!as.matrix) B = as.vector(B)

return(B)
}
