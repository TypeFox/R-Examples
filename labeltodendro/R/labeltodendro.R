relabel<-function(x)
{
if (!is.vector(x)) {stop("x should be a vector")}
outrelabel <- .C("RM_Relabel", PACKAGE="labeltodendro", S = as.integer(x), Pat_L = as.integer(length(x)), N_Iter = as.integer(1), result = as.integer(rep(0,length(x))))
return(outrelabel$result)
}



tabletodendro<-function(labmat,freq,labels=NULL)
{
################FUNCTION
treemaker<-function(labmat,freqvec,labels=NULL,treeorder=1:ncol(labmat))
{
## Function
mergemaker<-function(mflab)
{
origlab<- -1:-length(mflab)
j<-0
firstlab<-origlab[1]
merge<-c()
labchange<-c()
for (i in 1:(length(mflab)-1))
  {
  secondlab<-origlab[i+1]
  if (mflab[i]==mflab[i+1])
    {
    merge<-rbind(merge,c(firstlab,secondlab))
    j<-j+1
    firstlab<-j
    if (i==(length(mflab)-1) ) {labchange<-c(labchange,j)}
    } else {
           labchange<-c(labchange,j)
           firstlab<-origlab[i+1]
          }
  }
label<-as.numeric(row.names(table(mflab)))[table(mflab)>1]
labchange<-labchange[!(labchange==0)]
labchange<-as.numeric(row.names(table(labchange)))
return(list(merge=merge,labchangerow=labchange,lab=label))
}
mergemakeroftwo<-function(mergematdv, heightdv, labchangerowdv, labdv, origlabdv, origlabag, heightag)
{
### Function
firstunmachedlab<-function(labdv,labag)
{
exist<-FALSE
from<-NULL
to<-NULL
for (i in 1:length(labdv))
   {
   if (!(labdv[i]==labag[i])) {from<-labdv[i];to<-labag[i];exist<-TRUE;break} 
   }
return(list(exist=exist,from=from,to=to))
}
### Function

workinglabdv<-origlabdv
workinglabag<-origlabag
object<-firstunmachedlab(origlabdv,origlabag)
changedmergemat<-mergematdv
changedheight<-heightdv
while (object$exist)
       {
       if (sum(labdv==object$from)>0) 
             {
             secondadd<-labchangerowdv[labdv==object$from]
             } else {secondadd<- - which(workinglabdv==object$from) }
       if (sum(labdv==object$to)>0)
             {
             firstadd<-labchangerowdv[labdv==object$to]
             } else {firstadd<- - which(workinglabdv==object$to)} #it was  - which(origlabdv==object$to)
       if ( (firstadd>0) & (secondadd>0)) 
            {
               labchangerowdv<-labchangerowdv[-which(labdv==object$from)]
               labdv<-labdv[-which(labdv==object$from)]
            }
       changedmergemat<-rbind(changedmergemat,c(firstadd,secondadd))
       changedheight<-c(changedheight,heightag)
       workinglabdv[workinglabdv==object$from]<-object$to
       workinglabdv[workinglabdv>object$from]<-workinglabdv[workinglabdv>object$from]-1
       labdv[labdv==object$from]<-object$to # it is added recently I am not sure
       labdv[labdv>object$from]<-labdv[labdv>object$from]-1
       if (sum(labdv==object$to)>0) 
               { 
                labchangerowdv[labdv==object$to]<-nrow(changedmergemat)  
                } else {
                        labdv<-c(labdv,object$to)
                        labchangerowdv<-c(labchangerowdv,nrow(changedmergemat))
                        }
       object<-firstunmachedlab(workinglabdv,origlabag)
       }
return(list(mergemat=changedmergemat, height=changedheight,lab=labdv,labchangerow=labchangerowdv))
}

## Function


maxfreq<-max(freqvec)
if (length(maxfreq)>1) {warning("the tree has equal maximas, the first will be used")}
maxindex<-which(freqvec==maxfreq)[1]
if (is.unsorted(freqvec))# if freq is an increasing vector, do not touch it 
{
  freqsum<-0
  for (i in (maxindex):1)
    {
    freqsum<-freqsum+freqvec[i]
    freqvec[i]<-freqsum
    }
}
origlabdv<-labmat[nrow(labmat),]
obj<-mergemaker(origlabdv)
mergematdv<-obj$merge
labdv<-obj$lab
labchangerowdv<-obj$labchangerow
if (is.null(obj$merge)) {heightdv<-c()} else {heightdv<-rep(freqvec[nrow(labmat)] ,nrow(mergematdv))}

for (i in (nrow(labmat)-1):1)
   {
   origlabag<-labmat[i,]
   heightag<-freqvec[i]
   obj<-mergemakeroftwo(mergematdv, heightdv, labchangerowdv, labdv, origlabdv, origlabag, heightag)
   mergematdv<-obj$mergemat
   heightdv<-obj$height
   labchangerowdv<-obj$labchangerow
   labdv<-obj$lab
   origlabdv<-origlabag
   }
if (is.unsorted(heightdv)) {warning("there are local maximas, the dendrogram is ill defined")}
myobj<-list(merge=mergematdv,height=heightdv,order=treeorder,
hcut=maxfreq,labels=labels,labmat=labmat,freq=freqvec)
oldClass(myobj)<-c("labclust","hclust")
return(myobj)
}

	
is.agg<-function(labmat)
{
istart<-1
checkvec<-TRUE
if (!(sum(labmat[1,]==1)<ncol(labmat))) {istart<-2}
for (i in istart:(nrow(labmat)-1))
    {
    checkvec<-c(checkvec,is.aggvec(labmat[i,],labmat[i+1,]))
    }
if(sum(checkvec)<length(checkvec)) {return(FALSE)} else {return(TRUE)}
}

################FUNCTION

x<-labmat
if (!is.matrix(x)) {stop("labmat should be a matrix")}
if (!is.vector(freq)) {stop("freq should be a vector")}

if (nrow(x)<2) {stop("the table consists of just the most frequent labelling, dendrogram is useless")}

if (!(length(freq)==nrow(x))) {stop("labmat and freq do not match")}
if (missing(freq)) {stop("frequency of labels is not specified")}
y<-x
clustorder<-1:ncol(x)
for (i in 1:(nrow(x)-1))
   {
   y[i,]<-relabel(y[i,]) # insists increasing order of labels
   clustorder<-clustorder[order(y[i,])] # finds approrpiate order of labels
   y<-y[,order(y[i,])]
   y[i+1,]<-relabel(y[i+1,]) # insists increasing order of labels
   }
i<-i+1
y[i,]<-relabel(y[i,]) # insists increasing order of labels
clustorder<-clustorder[order(y[i,])] # finds approrpiate order of labels

if (sum(y[1,]==rep(1,ncol(y)))<ncol(y)) {y<-rbind(rep(1,ncol(y)),y); freq<-c(0,freq)} # adds all data in one group
if (sum(y[nrow(y),]==1:ncol(y))<ncol(y)) {y<-rbind(y,1:ncol(y)); freq<-c(freq,0)} # adds all data are separate

labels<-labels[clustorder]
for (i in nrow(y):1)
   {
   y[i,]<-relabel(y[i,]) # insists increasing order of labels
   }
if(!is.agg(y)) {stop("x is not in agglomerative order, use selectlabels(.)")}
if (sum(freq<0)>0) {warning("freq contanins negative values, they will be shifted ");freq<-freq-min(freq)}
return(treemaker(y,freq,labels=labels,treeorder=clustorder))
}




colorplot<-function(x,h,horiz=FALSE,...)
{
# modified version, from HeatPlus package
	

dendroploth <- function(x, h, cluscol, leaflab= "none", lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, horiz=TRUE,...))
    }
    # Not nice, but necessary
    pn  = stats:::plotNode
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
    x = cut(x, h)
    plot(x[[1]],horiz=TRUE,leaflab="none", yaxs="i",xaxs="i",...)    
    x = x[[2]]
    K = length(x)
    if (missing(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=TRUE)
        x1 = x2 + 1
   }
abline(v=h,col="gray",lwd=lwd,lty=2) 
}


dendroplotv <- function(x, h, cluscol, leaflab= "none", horiz=TRUE, lwd=1, ...)
{
    if (missing(h)) {
        return(plot(x, leaflab=leaflab, ...))
    }
    # Not nice, but necessary
    pn  = stats:::plotNode
    opar = par()[c("col","lwd")]
    on.exit(par(opar))
    par(lwd=lwd)
    x = cut(x, h)
    plot(x[[1]],horiz=FALSE,leaflab="none", yaxs="i",xaxs="i",...)    
    x = x[[2]]
    K = length(x)
    if (missing(cluscol)) {
       cluscol = rainbow(K)
    }
    x1 = 1
    for (k in 1:K) {
        x2 = x1 + attr(x[[k]],"members")-1
        par(col=cluscol[k])
        pn(x1,x2, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(), horiz=FALSE)
        x1 = x2 + 1
   }
abline(h=h,col="gray",lwd=lwd,lty=2)   
}

if (sum(x$labmat[nrow(x$labmat),]==(1:ncol(x$labmat)))==ncol(x$labmat)) {ylim<-c(x$freq[nrow(x$labmat)],max(x$freq))} else {ylim<-c(0,max(x$height))}
if (!(class(x)[1]=="labclust")) {stop("x should be a labclust object")}
if (is.null(x$hcut) & missing(h) ) {stop("no guide of where to cut the tree use plot(x)")}
if(!is.null(x$labels)) {warning("to see labels use plot(x)")}
y<-x
y$labels<-NULL
if (missing(h))  {
      {if (horiz){dendroploth(as.dendrogram(y),h=y$hcut,xlim=ylim,...)} else 
      {dendroplotv(as.dendrogram(y),h=y$hcut,ylim=ylim,...)}}
   } else
   {
    if (horiz){dendroploth(as.dendrogram(y),h=h,xlim=ylim)} else 
    {dendroplotv(as.dendrogram(y),h=h,ylim=ylim,...)}
    }
}


plot.labclust<-function(x,...)
{
plot(as.dendrogram(x),...)
}



is.aggvec<-function(labagg,labdiv)
  {
is.increasing<-function(vec)
{
!is.unsorted(vec)
}
changepoints<-function(lab)
{
if (is.null(lab)) {return(NULL)} else {
change<-c()
for (i in 1:(length(lab)-1))
    {
    firstlab<-lab[i]
    secondlab<-lab[i+1]
    if (!(firstlab==secondlab)) {change<-c(change,i)}
    }
return(change)
}
}


  mylabdiv<-labdiv[order(labagg)]
  mylabagg<-labagg[order(labagg)]
  mylabdiv<-relabel(mylabdiv)
  mylabagg<-mylabagg[order(mylabdiv)]
  mylabdiv<-mylabdiv[order(mylabdiv)]
  if(is.increasing(mylabagg)&is.increasing(mylabdiv))
     {
     aindex<-changepoints(mylabagg)
     dindex<-changepoints(mylabdiv)
     if (!is.null(aindex))
     {
     for (i in 1:length(aindex))
       {
       if (! (sum(aindex[i]==dindex)>0)) {return(FALSE)}
       }
     return(TRUE)
     }else
     return(TRUE)
     }
return(FALSE)
  }

 
labeltotable<-function(x)
{

  ##################### FUNCTION
 


aggfinder<-function(labmat,freq,mostlab)
{
### Function	
firstagg<-function(labmat,freq,mostlab)
{
agglab<-c()
aggfreq<-c()
for (i in 1:nrow(labmat))
   {
   if (is.aggvec(labmat[i,],mostlab))
      { 
      return(list(lab=matrix(labmat[i,],ncol=ncol(labmat)),freq=freq[i]))
      }
   }
return(NULL)
} 

### Function	

if(nrow(labmat)<2) {stop("label matrix is too short")}
if (sum(apply(labmat,1,max)<max(mostlab))<1) 
{return(NULL)} else 
{
ordlab<-matrix(labmat[(apply(labmat,1,max)<max(mostlab)),],ncol=ncol(labmat))
ordfreq<-freq[(apply(labmat,1,max)<max(mostlab))]
myorder<-order(apply(ordlab,1,max),ordfreq,decreasing=TRUE)
ordlab<-matrix(ordlab[myorder,],ncol=ncol(labmat))
ordfreq<-ordfreq[myorder]
agg<-firstagg(ordlab,ordfreq,mostlab)
agglab<-agg$lab
aggfreq<-agg$freq
 while( (!is.null(agg)))
   {
    if (sum(apply(ordlab,1,max)<max(agg$lab))<1) {agg<-NULL} else{
    ordfreq<-ordfreq[(apply(ordlab,1,max)<max(agg$lab))]
     ordlab<-matrix(ordlab[(apply(ordlab,1,max)<max(agg$lab)),],ncol=ncol(labmat))
     agg<-firstagg(ordlab,ordfreq,agg$lab)
     agglab<-rbind(agglab,agg$lab)
     aggfreq<-c(aggfreq,agg$freq)}
     }
  rownames(agglab)<-NULL
  return(list(lab=agglab,freq=aggfreq))
  }
}



divfinder<-function(labmat,freq,mostlab)
{

### Function	
firstdiv<-function(labmat,freq,mostlab)
{
divlab<-c()
divfreq<-c()
for (i in 1:nrow(labmat))
   {
   if (is.aggvec(mostlab,labmat[i,]))
      { 
      return(list(lab=matrix(labmat[i,],ncol=ncol(labmat)),freq=freq[i]))
      }
   }
return(NULL)
} 

### Function	
	
if(nrow(labmat)<2) {stop("label matrix is too short")}
if (sum(apply(labmat,1,max)>max(mostlab))<1) 
{return(NULL)} else 
 {
ordlab<-matrix(labmat[(apply(labmat,1,max)>max(mostlab)),],ncol=ncol(labmat))
ordfreq<-freq[(apply(labmat,1,max)>max(mostlab))]
myorder<-rev(order(apply(ordlab,1,max),ordfreq,decreasing=TRUE))
ordlab<-matrix(ordlab[myorder,],ncol=ncol(labmat))
ordfreq<-ordfreq[myorder]
div<-firstdiv(ordlab,ordfreq,mostlab)
divlab<-div$lab
divfreq<-div$freq
 while( (!is.null(div)))
   {
     if (sum(apply(ordlab,1,max)>max(div$lab))<1) {div<-NULL} else{
     ordfreq<-ordfreq[(apply(ordlab,1,max)>max(div$lab))]
     ordlab<-matrix(ordlab[(apply(ordlab,1,max)>max(div$lab)),],ncol=ncol(labmat))
     div<-firstdiv(ordlab,ordfreq,div$lab)
     divlab<-rbind(divlab,div$lab)
     divfreq<-c(divfreq,div$freq)}
     }
  rownames(divlab)<-NULL
  return(list(lab=divlab,freq=divfreq))
  }
}


  ##################### FUNCTION

  
  
  
  
  if (!is.matrix(x)) {stop("x should be a matrix")}
S <- as.vector(t(x))
Pat_L <- ncol(x)
N_Iter <- nrow(x)

outrelabel <- .C("RM_Relabel", PACKAGE="labeltodendro",S = as.integer(S), Pat_L = as.integer(Pat_L), N_Iter = as.integer(N_Iter), result = as.integer(rep(0,length(S))))
matrix(outrelabel$result,ncol=ncol(x),byrow=TRUE)
S<-outrelabel$result

outcount <- .C("Rcount", PACKAGE="labeltodendro", S = as.integer(S), Pat_L = as.integer(Pat_L), N_Iter = as.integer(N_Iter), Pattern = as.integer(rep(0,N_Iter*Pat_L)), nrowD=as.integer(0),
Frequency = as.integer(rep(0,N_Iter)))

labmat<-matrix(outcount$Pattern[1:(outcount$nrowD*Pat_L)],ncol=ncol(x),byrow=TRUE)
freq<-outcount$Frequency[1:(outcount$nrowD)]
if (max(freq)<2) {stop("the most frequent label has frequency 1!")}
return(selectlabels(labmat,freq))
}

relabel.matrix<-function(x)
{
if (!is.matrix(x)) {stop("x should be a matrix")}
S <- as.vector(t(x))
Pat_L <- ncol(x)
N_Iter <- nrow(x)
outrelabel <- .C("RM_Relabel", PACKAGE="labeltodendro", S = as.integer(S), Pat_L = as.integer(Pat_L), N_Iter = as.integer(N_Iter), result = as.integer(rep(0,length(S))))
return(matrix(outrelabel$result,ncol=ncol(x),byrow=TRUE))
}

selectlabels<-function(labmat,freq)
{
############################# FUNCTIONS

aggfinder<-function(labmat,freq,mostlab)
{
### Function	
firstagg<-function(labmat,freq,mostlab)
{
agglab<-c()
aggfreq<-c()
for (i in 1:nrow(labmat))
   {
   if (is.aggvec(labmat[i,],mostlab))
      { 
      return(list(lab=matrix(labmat[i,],ncol=ncol(labmat)),freq=freq[i]))
      }
   }
return(NULL)
} 

### Function	

if(nrow(labmat)<2) {stop("label matrix is too short")}
if (sum(apply(labmat,1,max)<max(mostlab))<1) 
{return(NULL)} else 
{
ordlab<-matrix(labmat[(apply(labmat,1,max)<max(mostlab)),],ncol=ncol(labmat))
ordfreq<-freq[(apply(labmat,1,max)<max(mostlab))]
myorder<-order(apply(ordlab,1,max),ordfreq,decreasing=TRUE)
ordlab<-matrix(ordlab[myorder,],ncol=ncol(labmat))
ordfreq<-ordfreq[myorder]
agg<-firstagg(ordlab,ordfreq,mostlab)
agglab<-agg$lab
aggfreq<-agg$freq
 while( (!is.null(agg)))
   {
    if (sum(apply(ordlab,1,max)<max(agg$lab))<1) {agg<-NULL} else{
    ordfreq<-ordfreq[(apply(ordlab,1,max)<max(agg$lab))]
     ordlab<-matrix(ordlab[(apply(ordlab,1,max)<max(agg$lab)),],ncol=ncol(labmat))
     agg<-firstagg(ordlab,ordfreq,agg$lab)
     agglab<-rbind(agglab,agg$lab)
     aggfreq<-c(aggfreq,agg$freq)}
     }
  rownames(agglab)<-NULL
  return(list(lab=agglab,freq=aggfreq))
  }
}



divfinder<-function(labmat,freq,mostlab)
{

### Function	
firstdiv<-function(labmat,freq,mostlab)
{
divlab<-c()
divfreq<-c()
for (i in 1:nrow(labmat))
   {
   if (is.aggvec(mostlab,labmat[i,]))
      { 
      return(list(lab=matrix(labmat[i,],ncol=ncol(labmat)),freq=freq[i]))
      }
   }
return(NULL)
} 

### Function	
	
if(nrow(labmat)<2) {stop("label matrix is too short")}
if (sum(apply(labmat,1,max)>max(mostlab))<1) 
{return(NULL)} else 
 {
ordlab<-matrix(labmat[(apply(labmat,1,max)>max(mostlab)),],ncol=ncol(labmat))
ordfreq<-freq[(apply(labmat,1,max)>max(mostlab))]
myorder<-rev(order(apply(ordlab,1,max),ordfreq,decreasing=TRUE))
ordlab<-matrix(ordlab[myorder,],ncol=ncol(labmat))
ordfreq<-ordfreq[myorder]
div<-firstdiv(ordlab,ordfreq,mostlab)
divlab<-div$lab
divfreq<-div$freq
 while( (!is.null(div)))
   {
     if (sum(apply(ordlab,1,max)>max(div$lab))<1) {div<-NULL} else{
     ordfreq<-ordfreq[(apply(ordlab,1,max)>max(div$lab))]
     ordlab<-matrix(ordlab[(apply(ordlab,1,max)>max(div$lab)),],ncol=ncol(labmat))
     div<-firstdiv(ordlab,ordfreq,div$lab)
     divlab<-rbind(divlab,div$lab)
     divfreq<-c(divfreq,div$freq)}
     }
  rownames(divlab)<-NULL
  return(list(lab=divlab,freq=divfreq))
  }
}
	
	
##################### FUNCTIONS
	

mostlab<-labmat[which(freq==max(freq))[1],]
mostfreq<-max(freq)
if (mostfreq<2) {stop("the most frequent label has frequency 1!")}
if (length(which(freq==mostfreq))>1) {warning("labels have more than one maximum")}

agg<-aggfinder(labmat,freq,mostlab)
div<-divfinder(labmat,freq,mostlab)

aggdivlab<-c()
aggdivfreq<-c()

if (is.null(agg) ) 
     {warning("no agglomeration found")} else
     {
      aggdivlab<-rbind(aggdivlab,agg$lab[rev(1:nrow(agg$lab)),])
      aggdivfreq<-c(aggdivfreq,agg$freq[rev(1:length(agg$freq))])
     }

aggdivlab<-rbind(aggdivlab,mostlab)
aggdivfreq<-c(aggdivfreq,mostfreq)

if (is.null(div) ) 
     {warning("no divisions found")} else
     {
      aggdivlab<-rbind(aggdivlab,div$lab)
      aggdivfreq<-c(aggdivfreq,div$freq)
     }

rownames(aggdivlab)<-NULL
return(list(labmat=aggdivlab,freq=aggdivfreq))
}

labeltodendro<-function(x,labels=NULL)
{
if (!is.null(labels)) {if (!(length(labels)==ncol(x))) {stop("label length do not match with columns of x")}}
tx<-labeltotable(x)
return(tabletodendro(tx$lab,tx$freq,labels))
}
