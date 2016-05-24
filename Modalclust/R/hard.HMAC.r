
hard.hmac<- function(hmacobj,level=NULL,n.cluster=NULL,plot=TRUE,colors=1:6,...){
  #### Hard clustering of HMAC object
  if(class(hmacobj)!="hmac") stop (" Your object is not of class hmac")

if(is.null(level)+is.null(n.cluster)+(plot!="TRUE")==3) level=1
if(is.null(level) & is.null(n.cluster)){
par(ask=TRUE)
cat("As no specific level is specified plots are provided for all levels",unique(hmacobj$level),"\n\n")
for(i in 1:max(hmacobj$level)){
cat("Level",i,"..")
hard.hmac(hmacobj,level=i)
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

  member=hmacobj$membership[[level]]

hmacobj$dat=as.matrix(hmacobj$dat)
if(plot==TRUE){
if(dim(hmacobj$dat)[2]==1){
#plot(density(hmacobj$dat,bw=sd(hmacobj$dat)*min(hmacobj$sigmas[which(level==hmacobj$level)])),main="",xlab="Data")
plot(density(hmacobj$dat,bw=apply(hmacobj$dat,2,sd)*min(hmacobj$sigmas[which(level==hmacobj$level)])),main="",xlab="Data")


for(k in 1:max(unique(member))){
	points(hmacobj$dat[member==k],rep(0,sum(member==k)),pch="|",col=k,main=NULL)
	}
}

else plot(data.frame(hmacobj$dat),col=colors[member],...)


title(main="Hard clustering produced by HMAC \n Colors represent different clusters",cex.main=1,sub=paste("Level",level,"with number of clusters",length(unique(member))) )
}
else return(member)
}  
}



