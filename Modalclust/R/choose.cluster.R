
##### Choose clusters according to user specified mode locations

#### Using knn
choose.cluster <- function(hmacobj,x=NULL,level=NULL,n.cluster=NULL){
### Works only for two dimensions
if(!dim(as.matrix(hmacobj$dat))[2]==2) stop("Choosing cluster not possible with more than two dimensions \n")


if(is.null(level) & is.null(n.cluster)){
 stop("provide either the level or the number of clusters")}
if(is.null(level) & !is.null(n.cluster)){  
 levels=hmacobj$level[which(hmacobj$n.cluster>=n.cluster)]
 
 
#### Possibility of adjustment of number of clusters get from MG's code
 if(length(levels)==0){
   stop("Cannot find a level with ", n.cluster," clusters \n Choose the number of clusters between ",
        min(hmacobj$n.cluster), " and ", max(hmacobj$n.cluster))}


 else {
   level=max(levels)
#   cat(level,hmacobj$n.cluster[level])
   if(hmacobj$n.cluster[min(which(hmacobj$level==level))]==n.cluster)
     cat("The level at which there are",n.cluster,"clusters is", level,"\n")
   else{
     cat("There are no levels with exactly",n.cluster,"clusters\n")
     cat("The level at which there are",hmacobj$n.cluster[min(which(hmacobj$level==level))],"clusters is", level,"\n")
   }   
}}
if(level> max(hmacobj$level) ) stop ("Provide a level not greater than ", max(hmacobj$level))


    plot(hmacobj$data,col="gray",pch=16)

if(!is.null(x)) {
	myclass=knn(hmacobj$data,x,hmacobj$membership[[level]],k=10)
	points(hmacobj$data[hmacobj$membership[[level]]==myclass,],col="red",pch=16)
	points(matrix(x,ncol=ncol(hmacobj$data)),col="green",pch=15)
	}

else{
more=TRUE;k=1
while(more){
	title(sub="Choose a point by clicking on the graph",line=3)
	title(sub="Stop the program by click anywhere outside the graph",line=4)
 	chosen.x=locator(1)
	x=c(chosen.x$x,chosen.x$y)


	if(x[1]<min(hmacobj$data[,1]) | x[1]> max(hmacobj$data[,1]) | x[2]<min(hmacobj$data[,2]) | x[2]> max(hmacobj$data[,2]) )  break
 
	myclass=knn(hmacobj$data,x,hmacobj$membership[[level]],k=10)
	points(hmacobj$data[hmacobj$membership[[level]]==myclass,],col=k+1,pch=16)

k=k+1
	}
	}
}
