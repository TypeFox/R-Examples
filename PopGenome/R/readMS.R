readMS <- function(file, big.data=FALSE){

if(!big.data){
out     <- read.ms.output(file.ms.output=file)
gametes <- out$gametes
dir.create("SwapMS")

for(xx in 1:length(gametes)){

 d <- gametes[[xx]]
 d <- list(matrix=d,positions=NaN)
 samplename <- paste("ms_sample_",xx,".RD",sep="")
 save(d,file= file.path ("SwapMS",samplename) )

}
test <- readData("SwapMS", SNP.DATA=TRUE, FAST=TRUE, format="RData", big.data=big.data)
unlink("SwapMS",recursive=TRUE)
return(test)
}# end of !big.data

if(big.data){
out     <- read.big.ms.output(file)
gametes <- out$gametes
dir.create("SwapMS")

for(xx in 1:length(gametes)){
 open(gametes[[xx]])
 d <- gametes[[xx]][,]
 close(gametes[[xx]])
 d <- list(matrix=d,positions=NaN)
 samplename <- paste("ms_sample_",xx,".RD",sep="")
 save(d,file= file.path ("SwapMS",samplename) )

}
test <- readData("SwapMS", SNP.DATA=TRUE, FAST=TRUE, format="RData", big.data=big.data)
unlink("SwapMS",recursive=TRUE)
return(test)
}#end of big.data




}
