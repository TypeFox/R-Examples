read.big.ms.output <- function(filename){

handle <- file(filename, "r")

listofMatrices <- list()
count <- 1
matrix <- NULL

BEG <- FALSE
while(1){

line       <- readLines(handle,n=1)
if(length(line)==0){
# the last matrix
if(length(matrix)==0){matrix <- 1}
listofMatrices[[count]] <- ff(matrix,dim=dim(matrix))
close(listofMatrices[[count]])
break
}
firstchar  <- substr(line,1,1)
secondchar <- substr(line,2,2)
        
	if((firstchar=="0" || firstchar=="1") && (secondchar=="0" || secondchar=="1" || secondchar=="") ){
	BEG    <- TRUE
        vec    <- strsplit(line,split="")[[1]]
	vec    <- as.numeric(vec)
        matrix <- rbind(matrix,vec)
        }
	

if(firstchar=="s" && BEG){
 if(length(matrix)==0){matrix <- 1}
 listofMatrices[[count]] <- ff(matrix,dim=dim(matrix))
 close(listofMatrices[[count]])
 matrix <- NULL
 count  <- count + 1 
}	



}
close(handle)
segsites  <- sapply(listofMatrices,function(x){if(length(x)==1){return(0)}else{return(dim(x)[2])}})
positions <- as.list(rep(NaN,length(segsites)))
return(list(gametes=listofMatrices,segsites=segsites,positions=positions))

}# End of function
