my_read.nexus <- function(file){

#LISTE  <- read.nexus.data(file)
LISTE  <- APE_read.nexus.data(file) #Copyright APE package on CRAN  !
size   <- length(LISTE)

MAT <- NULL

for(xx in 1:size){

MAT <- rbind(MAT,LISTE[[xx]])

} 

if(is.element(".",MAT)){
 value  <- apply(MAT,2,function(x){
   x[x=="."] <- x[1]
   return(x)
 })
 MAT <- value
}

rownames(MAT) <- names(LISTE)
return(MAT) 

}
