nt.connect<-function(nt)
{
  returnECt=1# The *th connected part of the network.
  nel<-readshpnw(nt)
  nodelist<-nel[[2]]
  Nn<-length(nodelist)
  edgelist<-nel[[3]]
  Ne<-dim(edgelist)[1]
  #intialize the status of each node is 0, note: 0--unvisited, 1---visited
  # all the edges are intialized with a 0 value
  nst<-rep(0, length=Nn)
  econCt<-rep(0, length=Ne)
  
  ect<-1               
  vting<-matrix(edgelist[1,],ncol=3)
  econCt[1]<-ect
  unvted<-Ne-1
  unvtedEL<-edgelist[-1, ]
  fromn<-unvtedEL[,2]
  ton<-unvtedEL[,3]
  while (unvted>0)
  {
    #Suppose the 1st visiting node is 1st node
    n<-dim(vting)[1]
    vtingEidxs<-c()
    for (i in 1:n)
      {
        idxs1<-c(which(fromn==vting[i,2]),which(ton==vting[i,2]),which(fromn==vting[i,3]),which(ton==vting[i,3]))
        vtingEidxs<-c(vtingEidxs, idxs1)
      }
    if (length(vtingEidxs)>0)
      {
        eidxs<-norep(vtingEidxs)
        vtingE<-unvtedEL[eidxs, 1]
        econCt[vtingE]<-ect
        vting<-matrix(unvtedEL[eidxs,],ncol=3)
        unvtedEL<-matrix(unvtedEL[-eidxs,],ncol=3)
        unvted<-unvted-length(eidxs)
        fromn<-unvtedEL[,2]
        ton<-unvtedEL[,3]
      }
    else
      {
        vting<-matrix(unvtedEL[1,],ncol=3)
        ect<-ect+1
      }
  }
  bnt<-bbox(nt)
  plot.new()
  plot.window(xlim=c(bnt[1,1],bnt[1,2]),ylim=c(bnt[2,1],bnt[2,2]))
  ects<-norep(econCt)
  numConected<-length(ects)
  main<-paste(paste("There are ", as.character(numConected)), " self-connected parts in this data set")
  title(xlab="",ylab="",main=main)
  cols<-rainbow(numConected)
  Elns<-slot(nt, "lines")
  SLns<-as.SpatialLines.SLDF(nt)
  edf<-slot(nt, "data")
  for (i in 1:Ne)
  {
    lines(Elns[[i]], col=cols[which(econCt[i]==ects)])
  }
  #################return 
  idx<-MajEinV(econCt, ects)
  if (numConected==1) res<-nt
  else
  {
    ect<-ects[idx] 
    idxs<-which(econCt==ect)
    sldf<-SpatialLinesDataFrame(SLns[idxs], edf[idxs, ], match.ID=F)
    res<-sldf 
  }
  #res<-econCt  
  res
}

##get rid of all the repeated number in a vector
norep<-function(v)
{
  if (is.null(v))
  stop("V can't be an empty vector")
  vsorted<-sort(v)
  n<-length(v)
  res<-c(vsorted[1])
  if (n>1)
  {
    for (i in 2:n)
    {
      if (vsorted[i-1]!=vsorted[i]) res<-c(res, vsorted[i])
    }
  }
  res
}

#Get the majority of elements in a vector
MajEinV<-function(v, elv)
{
  if (is.null(v))
  stop("V can't be an empty vector")
  nelv<-elv
  ne<-length(elv)
  for (i in 1:ne)
  {
    nelv[i]<-length(which(elv[i]==v))
  }
  res<-which(max(nelv)==nelv)[1]
  res
}