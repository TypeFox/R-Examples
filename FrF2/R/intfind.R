`intfind` <-
function(i,j,mat){
     ## function for finding the correct interaction coefficient 
     ## to be used in function iap
     sp <- NULL
     for (k in 1:ncol(mat)){
          if (!all(mat[,k]==c(i,j))) next 
          else sp <- k
          break
          }
     sp
    }

