get_cov_matrix_index <- function(M,i){
#Get the row and col index of the Matrix M lower triangular number i (i=1
#to num off diag elements in lower triangular of M, starting with col1 then col2
#until coln
k=0
col = -1
row = -1
for(n in 1:size(M,2)){
    for(m in 1:size(M,1)){
        if((n<m)){
            k=k+1
        }
        if((i==k)){
            row = m
            col=n
            return
        }

    }
}
return(list( row= row,col =col )) 
}
