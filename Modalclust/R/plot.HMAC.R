plot.hmac <- function(x,mycol=1:6,level=1,n.cluster=NULL,userclus=NULL,sep=.1,...){
### Start traces the tree from top to the level start
### Level draws the whole tree but colors it according to clsutering of the level specfied by level.
### One can also use user color to see how the two classications match
hmacobj=x

if(level==1 & !is.null(n.cluster)){  
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

if(level>max(hmacobj$level)) stop("Choose the level not greater than ", max(hmacobj$level))

orig.levels=max(hmacobj$level)

### put an if statement
nlevels=orig.levels+1

n=length(hmacobj$membership[[1]]) ### Number of datapoints
original=1:n
membmat=matrix(1,nlevels,n)

for(j in 1:orig.levels){
  membmat[j,]=hmacobj$membership[[j]]
}

calllist <- list()
for(j in 1:(nlevels-level+1))
  
  calllist[[j]]=membmat[(nlevels-j+1),]


ordered=do.call("order", calllist)
membmat=membmat[,ordered] #### ordered membership matrix

### If user does not provides a clustering scheme 
startx=1:n
x=startx
for(i in 2:n)
  if(membmat[level,i]!=membmat[level,i-1]) x[i:n]=x[i:n]+round(sep*n)


#cat(startx,"\n",x)


#x=1:n
plot(1,1,type="n",xlim=range(x),ylim=c((level-0.5),nlevels),xaxt="n",yaxt="n",bty="n",xlab="",ylab="clustering level",main="HMAC Dendrogram")
yxis=c(hmacobj$level,max(hmacobj$level)+1)
axis(2, at=yxis,labels=yxis, col.axis="black", las=2)

if(is.null(userclus)){
points(x,y=rep(level,n),col=mycol[membmat[level,]],type="h")
mtext(row.names(hmacobj$data)[ordered],at=x,side=1,las=2,col=mycol[membmat[level,]],cex=100/n)
}

### Else
else{
if(length(userclus)==n){
  mycol=rainbow(max(userclus))

points(x,y=rep(1,n),col=mycol[membmat[level,]],type="h")
#for(i in 1:n) 
#lines(rep(x[i],2),c((level-.7),level),col=mycol[membmat[level,i]])
mtext(row.names(hmacobj$data)[ordered],at=x,side=1,las=2,col=mycol[userclus])
           
title("\n \n Node names colored using user provided  groups")
}
else cat("Number of Datapoints does not match")
}

### Horizontal bar for lowest level

split.x=split(x,membmat[level,])

for(k in 1:length(split.x)){
  lines(x=range(split.x[[k]]),y=rep(level,2),col=mycol[k])
}


#for( i  in 2:n){    
#  if(membmat[level,i]==membmat[level,i-1]){
#lines(c(x[i],x[i-1]),rep(level,2),col=mycol[membmat[level,i]])
#}
#}

#cat(x,"\n")


for(j in (level):(nlevels-1)){

mid=sort(findmid(x,membmat[j,]))

for(k in 1:length(mid)) lines(rep(mid[k],2),c(j,j+1))

if(length(mid)>1){
  for(k in 1:(length(mid)-1)){
 #   cat(mid[k],mid[k+1],"...")
#   cat( findmemb(x,mid[k],membmat[j+1,]), findmemb(x,mid[k],membmat[j+1,]),"\n")
  #  if(membmat[(j+1),which(x==floor(mid[k]))] == membmat[(j+1),which(x==floor(mid[k+1]))])  {
      if(knn(x,mid[k],membmat[j+1,]) == knn(x,mid[k+1],membmat[j+1,])){
      lines(c(mid[k],mid[k+1]),c(j+1,j+1))}
  }
}


}

cat("\n")

}

findmid <- function(x,memb){
aa=split(x,memb)
bb=sapply(aa,mean)
return(bb)
}


