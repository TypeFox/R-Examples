deletecodongaps <- function(matrix){

if(length(which(matrix==5))==0 & length(which(matrix==6))==0){return(matrix)}

erg <- apply(matrix,2,del)

id <- which(erg==1)
if(length(id)==0){return(matrix)}

mod     <- rep(-1,length(id))
codonid <- matrix(-1,length(id),3)

for(xx in 1:length(id)){

# x x gap
 if(id[xx]%%3==0){  
  codonid[xx,1] <- id[xx] -2
  codonid[xx,2] <- id[xx] -1
  codonid[xx,3] <- id[xx] 
 }
 if(id[xx]%%3==2){  
  codonid[xx,1] <- id[xx] -1
  codonid[xx,2] <- id[xx] 
  codonid[xx,3] <- id[xx] +1
 }
 if(id[xx]%%3==1){  
  codonid[xx,1] <- id[xx] 
  codonid[xx,2] <- id[xx] +1
  codonid[xx,3] <- id[xx] +2
 }

}
d <- as.vector(codonid)
d <- unique(d)
matrix <- matrix[,-d]

return(matrix)
}

################
##Function######
################

del <- function(vek){
if(is.element(5,vek)|is.element(6,vek)){ # gaps or unknown_positions
return(1)
}else{return(0)}
}



