coordmatch=function(coordref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg", radunit='asec', sep = ":", kstart=10, ignoreexact=FALSE){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  if (radunit %in% c("deg", "amin", "asec", "rad") == FALSE) {
    stop("radunit must be one of deg, amin, asec or rad")
  }
  origrad=rad
  coordref=rbind(coordref)
  coordcompare=rbind(coordcompare)
  N=length(coordref[,1])
  kmax=length(coordcompare[,1])
  if(length(rad)==1){rad=rep(rad,N)}
  if(length(rad) != N){stop("Length of rad must either be 1 or the same length as rows in coordref")}
  if (inunitref == "sex") {
    coordref[,1] = hms2deg(coordref[,1], sep = sep)
    coordref[,2] = dms2deg(coordref[,2], sep = sep)
    coordref=cbind(as.numeric(coordref[,1]),as.numeric(coordref[,2]))
  }
  if (inunitref == "rad") {
    coordref[,1] = coordref[,1] * 180/pi
    coordref[,2] = coordref[,2] * 180/pi
  }
  if (inunitcompare == "sex") {
    coordcompare[,1] = hms2deg(coordcompare[,1], sep = sep)
    coordcompare[,2] = dms2deg(coordcompare[,2], sep = sep)
    coordcompare=cbind(as.numeric(coordcompare[,1]),as.numeric(coordcompare[,2]))
  }
  if (inunitcompare == "rad") {
    coordcompare[,1] = coordcompare[,1] * 180/pi
    coordcompare[,2] = coordcompare[,2] * 180/pi
  }
  if (radunit == "asec"){
    rad=rad*(pi/180)/3600
  }
  if (radunit == "amin"){
    rad=rad*(pi/180)/60
  }
  if (radunit == "deg"){
    rad=rad*(pi/180)
  }
  userad=max(rad,na.rm = TRUE)
  
  coordrefxyz=sph2car(coordref,deg=TRUE)
  coordcomparexyz=sph2car(coordcompare,deg=TRUE)
  
  ksuggest=kstart
  while(is.na(ksuggest)==FALSE){
    tempmatch=nn2(coordcomparexyz,coordrefxyz,searchtype='radius',radius=userad,k=ksuggest)
    ignore=tempmatch[[1]]==0
    tempmatch[[2]][ignore]=NA
    tempmatch[[2]]=2*asin(tempmatch[[2]]/2)
    if(ignoreexact){remove=which(!(tempmatch[[2]]<=rad & tempmatch[[2]]>0))}else{remove=which(!tempmatch[[2]]<=rad)}
    tempmatch[[1]][remove]=0
    tempmatch[[2]][remove]=NA
    if (all(is.na(tempmatch[[2]][,ksuggest]))){
      kendmin=NA
    }else{
      kendmin=min(tempmatch[[2]][,ksuggest],na.rm = TRUE)
    }
    if(is.na(kendmin)==FALSE & ksuggest<kmax){
      ksuggest=ceiling(ksuggest*max(rad/tempmatch[[2]][,ksuggest],na.rm=TRUE))
      ksuggest=min(ksuggest,kmax)
    }else{
      ksuggest=NA
    }
  }

  keepcols=which(apply(is.na(tempmatch[[2]]),2,sum)<N)
  tempmatch[[1]]=matrix(tempmatch[[1]][,keepcols],nrow=N,byrow=FALSE)
  tempmatch[[2]]=matrix(tempmatch[[2]][,keepcols],nrow=N,byrow=FALSE)
  
  if (radunit == "asec"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180)/3600)
  }
  if (radunit == "amin"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180)/60)
  }
  if (radunit == "deg"){
    tempmatch[[2]]=tempmatch[[2]]/((pi/180))
  }
  
  if(length(tempmatch[[1]])>0){
    Nmatch=rowSums(tempmatch[[1]]!=0)
  }else{
    Nmatch=NA
  }
  
#   if(findbest & is.na(Nmatch[1])==FALSE){
#     bestID=tempmatch[[1]]
#     bestsep=tempmatch[[2]]
#     
#     if(besttype=='absbest'){
#       checkunique=unique(bestID)
#       checkunique=checkunique[checkunique>0]
#       for(i in 1:length(checkunique)){
#         allIDs=which(bestID==checkunique[i])
#         minsep=which.min(bestsep[allIDs])
#         allIDs=allIDs[-minsep]
#         bestID[allIDs]=0
#         bestsep[allIDs]=NA
#       }
#     }
#     
#     arrdim=dim(bestID)
#     bestmatches={}
#     while(any(bestID>0)){
#       nextbestwhich=which.min(bestsep)
#       nextbestarrID=arrayInd(nextbestwhich,arrdim)
#       nextbestID=bestID[nextbestwhich]
#       bestmatches=rbind(bestmatches,c(nextbestarrID[1,1],nextbestID,bestsep[nextbestwhich]))
#       bestsep[bestID==nextbestID]=NA
#       bestID[bestID==nextbestID]=0
#       bestsep[nextbestarrID[1,1],]=NA
#       bestID[nextbestarrID[1,1],]=0
#     }
#     output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch, bestmatches=bestmatches)
#   }else{
#     output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch)
#   }
  if(is.na(Nmatch[1])==FALSE){
    bestID=tempmatch[[1]][,1]
    bestsep=tempmatch[[2]][,1]
    select=which(bestID>0)
    bestsep=bestsep[select]
    bestID=bestID[select]
    orderID=order(bestsep)
    keep=!duplicated(bestID[orderID])
    bestrefID=select[orderID][keep]
    bestcompareID=bestID[orderID][keep]
    bestsep=bestsep[orderID][keep]
    reorderID=order(bestrefID)
    bestrefID=bestrefID[reorderID]
    bestcompareID=bestcompareID[reorderID]
    bestsep=bestsep[reorderID]
    bestmatch=cbind(refID=bestrefID, compareID=bestcompareID, sep=bestsep)
    output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch, bestmatch=bestmatch)
  }else{
    output=list(ID=tempmatch[[1]], sep=tempmatch[[2]], Nmatch=Nmatch, bestmatch=NA)
  }
  return(output)
}

coordmatchsing=function(RAref,Decref, coordcompare, rad=2, inunitref = "deg", inunitcompare="deg", radunit='asec', sep = ":", ignoreexact=FALSE){
  if (inunitref %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitref must be one of deg, rad or sex")
  }
  if (inunitcompare %in% c("deg", "rad", "sex") == FALSE) {
    stop("inunitcompare must be one of deg, rad or sex")
  }
  origrad=rad
  coordref=rbind(c(RAref,Decref))
  coordcompare=rbind(coordcompare)
  N=length(coordref[,1])
  if(N>1){stop("Only 1 RAref and Decref value allowed!")}
  if (inunitref == "sex") {
    coordref[,1] = hms2deg(coordref[,1], sep = sep)
    coordref[,2] = dms2deg(coordref[,2], sep = sep)
    coordref=cbind(as.numeric(coordref[,1]),as.numeric(coordref[,2]))
  }
  if (inunitref == "rad") {
    coordref[,1] = coordref[,1] * 180/pi
    coordref[,2] = coordref[,2] * 180/pi
  }
  if (inunitcompare == "sex") {
    coordcompare[,1] = hms2deg(coordcompare[,1], sep = sep)
    coordcompare[,2] = dms2deg(coordcompare[,2], sep = sep)
    coordcompare=cbind(as.numeric(coordcompare[,1]),as.numeric(coordcompare[,2]))
  }
  if (inunitcompare == "rad") {
    coordcompare[,1] = coordcompare[,1] * 180/pi
    coordcompare[,2] = coordcompare[,2] * 180/pi
  }
  
  coordrefxyz=sph2car(coordref,deg=TRUE)
  coordcomparexyz=sph2car(coordcompare,deg=TRUE)
  
  dotprod=coordcomparexyz[,1]*coordrefxyz[1]+coordcomparexyz[,2]*coordrefxyz[2]+coordcomparexyz[,3]*coordrefxyz[3]
  dotprod[dotprod< -1]=-1;dotprod[dotprod>1]=1
  ang=acos(dotprod)
  if (radunit == "asec"){
    ang=ang/((pi/180)/3600)
  }
  if (radunit == "amin"){
    ang=ang/((pi/180)/60)
  }
  if (radunit == "deg"){
    ang=ang/(pi/180)
  }
  if(ignoreexact){select=which(ang<=rad & ang>0)}else{select=which(ang<=rad)}
  if(length(select)==0){
    output=list(ID=NA, sep=NA, Nmatch=NA, bestmatch=NA)
  }else{
    ID=select
    sep=ang[select]
    orderID=order(sep)
    ID=ID[orderID]
    sep=sep[orderID]
    bestmatch=c(compareID=select[which.min(ang[select])],sep=min(ang[select]))
    output=list(ID=ID, sep=sep, Nmatch=length(select), bestmatch=bestmatch)
  }
  return(output)
}