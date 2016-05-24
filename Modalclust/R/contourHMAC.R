contour.hmac=function(x, n.cluster=NULL,level=NULL,prob=NULL,smoothplot=FALSE,...){
  #### Hard clustering of HMAC object
  ### 2d plot of hmac objects

 hmacobj=x
 
if(!dim(as.matrix(hmacobj$dat))[2]==2) stop("Contour plot not possible with more than two dimensions \n")


  require(MASS)
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

  member=hmacobj$membership[[level]]
  n.cluster=hmacobj$n.cluster[min(which(hmacobj$level==level))]

#### Choosing different color amps
palletes=c("Blues","Greens","Oranges","Reds")


if(smoothplot==TRUE) smoothScatter(hmacobj$dat)
else plot(hmacobj$dat)
mycols=rep("",length(member))
  	for(j in 1:n.cluster){
	points(hmacobj$dat[member==j,],col=j,pch=j)
#### Change bandwidth to take care of sd
  	}


d=kde2d(hmacobj$data[,1],hmacobj$data[,2])

if(is.null(prob)) contour(d, add=TRUE,...)
else contour(d,level=prob,add=TRUE,...)
title("Contour plot")

}


