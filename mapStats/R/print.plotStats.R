print.plotStats <- function(x, 
                            horizontal.fill = TRUE,
                            num.row =1,
                            num.col =1,
                            ...) {

    #create vectors for row and column indexes

    
    indices <- data.frame(k=1:length(x), n_row = num.row, n_col=num.col)

    if ( horizontal.fill == TRUE ) {
    
        indices$rowind <- (ceiling(indices$k / indices$n_col)) %% indices$n_row
	indices$rowind[ indices$rowind == 0 ] <- indices$n_row[  indices$rowind == 0 ]    
	
	indices$colind <- indices$k %% indices$n_col
	indices$colind[ indices$colind == 0 ] <- indices$n_col[  indices$colind == 0 ]


     }

    else {
  
        indices$colind <- (ceiling(indices$k / indices$n_row)) %% indices$n_col
	indices$colind[ indices$colind == 0 ] <- indices$n_col[  indices$colind == 0 ]    
	
	indices$rowind <- indices$k %% indices$n_row
	indices$rowind[ indices$rowind == 0 ] <- indices$n_row[  indices$rowind == 0 ] 

    }

   indices$notlast <- (indices$k %% (indices$n_row * indices$n_col)) !=0  &  (indices$k != nrow(indices))
   indices <- indices[ , c("colind","rowind","n_col","n_row", "notlast") ]


   #print the objects

   
   for (k in 1:length(x)) { print(x[[ k ]], split=as.numeric(as.vector( indices[k, 1:4])), more = indices$notlast[k])  }

   x

}
   



