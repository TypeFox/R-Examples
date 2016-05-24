imageQT <-function (image,test=TOS2D,minsize=64,alpha=0.05,...){

# this function implements a quadtree decomposition of an image,
# splitting an image down to a minimum testing size. 

# Define internal convenience functions:

### 1. generic checking function: is.wholenumber ###

is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#####

### 2. index labelling function: next.index ###

next.index <-function(l=NULL) {

if(is.null(l)){
        labs<-as.character(0:3)
}
else{
        labs<-as.vector(outer(as.character(0:3),l, function(x,y){paste(y,x,sep="")}))
}

return(labs)
}

#####

### 3. list concatenating function newlist ###

newlist <-function(oldlist){

# assumes a list of length at least 1

ll<-length(oldlist)

newL<-NULL
for(i in 1:ll){
        newL<-c(newL,oldlist[[i]])
}

return(newL)

}

### 4. image splitting function: splitim ###

splitim <-function(im){

# a function which splits an image into its quadrants, whose output is a list
# next.index gives the indices of the quadrants according to:
# 
# 0: top-left; 1: bottom-left; 2: top-right; 3: bottom-right.

d<-dim(im)[1]
dh<-d/2

quads<-list()

quads[[1]]<-im[1:dh,1:dh]                 
quads[[2]]<-im[(dh+1):d,1:dh]             
quads[[3]]<-im[1:dh,(dh+1):d]             
quads[[4]]<-im[(dh+1):d,(dh+1):d]         

names(quads)<-next.index(names(quads))

return(quads)
}

### 5. decision function wrapper: statdecide ###

statdec <-function(X,test=TOS2D,alpha=0.05,...){

p<-test(X,...)$p.value

decision <- (p > alpha)

return(decision)
}


### start of main function ###

oa<-list(...)

si <-"smooth"%in%names(oa)

dname<-deparse(substitute(image))

depthlim<-ceiling(log2(minsize))

# check image dimension and define maxdepth:

r<-nrow(image)
c<-ncol(image)

if((r!=c)|!is.wholenumber(log2(r))|!is.wholenumber(log2(c))){
	stop("Please supply valid image!!\n\n")
}
else{
	depth<-log2(r)
}

indl<-resl<-list()

####################

# general procedure:
# 1. split image
# 2. test each subimage for stationarity -> four element vector of 0/1 (stationary is 1, nonstationary is 0)
# remove stationary ones from list of (sub)images left to split, recording subimage "tree index"
# repeat until minimum testing size is reached.

iml<-list()

iml[[1]]<-image
res<-FALSE			# always split original image
testdim<-dim(image)[1]
resl<-c(resl,res)
indl[[1]]<-NULL

imS<-indS<-list()

while((prod(res)==0)&(testdim>minsize)){
	cat("testing dimension is:",testdim,"\n") 		

	tmp<-next.index(names(iml))	#get total possible next names
	iml<-newlist(lapply(iml,splitim))
	testdim<-testdim/2
	if(si){
		res<-unlist(lapply(iml,statdec,alpha=alpha,...),use.names=F)
	}
	else{
		res<-unlist(lapply(iml,statdec,alpha=alpha,smooth=FALSE,...),use.names=F)
	}
	resl<-c(resl,list(res))

	imS<-c(imS,iml[which(res==1)])
	indS<-c(indS,tmp[which(res==1)])

	iml<-iml[which(res==0)]
	names(iml)<-tmp[which(res==0)]

	indl<-c(indl,list(tmp[which(res==0)]))		#record splitting names
}

l<-list(data.name=dname,indl=indl,resl=resl,imsize=r,imS=imS,indS=indS,minsize=minsize)

class(l)<-"imageQT"

return(l)

}

