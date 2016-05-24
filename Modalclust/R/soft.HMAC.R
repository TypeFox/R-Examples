#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  %%%%%%%  Soft Clustering from HMAC output
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

soft.hmac <- function(hmacobj,n.cluster=NULL,level=NULL,boundlevel=0.4,plot=TRUE){
if(class(hmacobj)!="hmac") stop (" Your object is not of class hmac")
if(is.null(level)+is.null(n.cluster)+(plot!="TRUE")==3) level=1
if(is.null(level) & is.null(n.cluster)){ 
par(ask=TRUE)
cat("As no specific level is specified plots are provided for all levels",unique(hmacobj$level),"\n\n")
for(i in 1:max(hmacobj$level)){
cat("Level",i,"..")
soft.hmac(hmacobj,level=i)
}
par(ask=FALSE)
}

else{
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


dat=as.matrix(hmacobj$data)


n=dim(dat)[1];
m=dim(dat)[2];


  #scale.dat=sd(dat)
  scale.dat=apply(dat,2,sd)
  Data=scale(dat,scale=scale.dat,center=FALSE)

  n.cluster=hmacobj$n.cluster[min(which(hmacobj$level==level))]
  member=hmacobj$membership[[level]]
  postmat=matrix(0,n,n.cluster)   
  sigma=hmacobj$sigmas[min(which(hmacobj$level==level))]; #select j-th bandwidth
 
  dens=matrix(0,n,n)
           for (i in 1:n)
                  dens[i,]=mydmvnorm(Data,Data[i,],sigma^2)
 
 post.prob=matrix(0,n,n.cluster)
  
      for(j in 1:n.cluster) ### Un normalized
           post.prob[,j]= apply(dens[,member==j,drop=FALSE],1,sum)       
   post.prob=post.prob/apply(post.prob,1,sum)


if(plot){
if(!dim(as.matrix(hmacobj$dat))[2]==2) stop("Soft clustering plot not available for this dimension \n Try using soft.hmac with plot=F to get the posterior probability values")
  if(n.cluster==2)
    mycol=rgb(post.prob[,1],post.prob[,2],blue=0)
 else if(n.cluster==1)
    mycol=rgb(post.prob[,1],green=0,blue=0)
  else
    mycol=rgb(post.prob[,1],post.prob[,2],post.prob[,3])
  plot(dat,type="n")

  points(dat[apply(post.prob>boundlevel,1,sum)>1,],col="black",pch=16)
boundary=dat[apply(post.prob>boundlevel,1,sum)>1,,drop=FALSE]
points(dat,col=mycol)


if(nrow(boundary))
legend("topleft",legend=paste("boundary points"),col="black",pch=16) 
title("Soft clustering produced by HMAC \n Colors represent different clusters",cex.main=1)
}

	boundary=dat[apply(post.prob>boundlevel,1,sum)>1,]
	output=list()
	output[["post.prob"]]=post.prob
	output[["boundary.point"]]=boundary 
	
if(plot==FALSE) return(output)
}
}
