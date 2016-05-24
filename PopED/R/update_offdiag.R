update_offdiag <- function(mat,iOffIndex,add_value){
#Update a symmetric off diagnol matrix, where iOffIndex is the lower
#off diagonal triangular matrix index. This is updated (added) with the new_value

k=0
for(i in 1:size(mat,1)){
    for(j in 1:size(mat,1)){
        if((i<j)){
            k=k+1
            if((k==iOffIndex)){
                mat[i,j] = mat[i,j]+add_value
                mat[j,i] = mat[j,i]+add_value
                return
            }
        }
    }
}
return( mat ) 
}
