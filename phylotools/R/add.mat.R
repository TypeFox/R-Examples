#### Function add.mat as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


add.mat <-
function(mat1, mat2){
    nams <- c(mat1[,1], setdiff(mat2[,1], mat1[,1]))
    col2 <- rep("", length(nams))
    result <- data.frame(mames = nams, col2) 
    result[1:nrow(mat1), 2] <- mat1[,2] 
    result <- appendchar(result, "?")
    for(i in 1:nrow(result)){
        for(j in 1:nrow(mat2)){
   	        if( mat2[j, 1] == result[i, 1] ){
   		        result[i, 2] <- paste(result[i, 2], mat2[j, 2], sep = "")
   		    }
   		}
    }
    result <- appendchar(result, "?")
    return(result)
}

