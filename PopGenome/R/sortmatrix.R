#Cai JJ (2008) PGEToolbox: A Matlab toolbox for population genetics and evolution
#Journal of Heredity Jul-Aug;99(4):438-40. doi:10.1093/jhered/esm127
#modified
########################################################
################### sortmatrix #########################
########################################################
# A Matlab Function PGE Toolbox: unused

sortmatrix <- function (matrix){

 m <- dim(matrix)[1]
 n <- dim(matrix)[2]

if(m==1 || n==1){
  return(sort(matrix))

}else{
     
     rownames(matrix) <- 1:m   
     
     
     for(xx in 1:n){
         mat <- unique(matrix[,xx])
         if(length(mat)>1){
          break;
         }
     }
     
     sortmatrix <- sort(matrix[,xx])
     sortmatrix <- matrix[names(sortmatrix),]
}

return(sortmatrix)
}
