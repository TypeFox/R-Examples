`drawPedigree` <-
function(Ped,cohorts=NULL,sex=NULL,dat=NULL, dots='n',plotfull='y',writeCohortLabels='n',links='all',sexInd=c(0,1),dotSize=0.001,dataDots='n',dataDots.cex=2,cohortLabs.cex=1, retain='informative',focal=NULL, sexColours=c('red','blue'), ...) {

for(x in 1:3) Ped[,x]<-as.character(Ped[,x])

  "prune"<-function(pedigree, keep, make.base=FALSE){
     ind.keep<-keep
     nind<-length(ind.keep)+1
     while(length(ind.keep)!=nind){
       nind<-length(ind.keep)
       ind.keep<-union(na.omit(unlist(pedigree[,2:3][match(ind.keep,pedigree[,1]),])), ind.keep)
     }
     pedigree<-pedigree[sort(match(ind.keep, pedigree[,1])),]
     if(make.base){
       if(any(match(pedigree[,2], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){
                 stop("Dams appearing after their offspring: use fixPedigree")}
       if(any(match(pedigree[,3], pedigree[,1])>match(pedigree[,1], pedigree[,1]), na.rm=T)){
                 stop("Sires appearing after their offspring: use fixPedigree")}
       phenotyped<-pedigree[,1]%in%keep
       delete<-rep(FALSE, dim(pedigree)[1])
       for(i in 1:dim(pedigree)[1]){
         nlinks<-phenotyped[i]+sum(pedigree[,2]%in%pedigree[,1][i])+sum(pedigree[,3]%in%pedigree[,1][i])+sum(is.na(pedigree[i,][2:3])==FALSE)
         if(nlinks<2 & phenotyped[i]==FALSE){
           pedigree[,2][which(pedigree[,2]==pedigree[,1][i])]<-NA
           pedigree[,3][which(pedigree[,3]==pedigree[,1][i])]<-NA
           delete[i]<-TRUE
         }
       }
       if(any(delete)){
         pedigree<-pedigree[-which(delete),]
       }
     }
     pedigree
  }

  names(Ped)[1]<-"id"
  if(names(Ped)[2]!="dam"|names(Ped)[3]!="sire"){
    if(names(Ped)[3]%in%c("mum","mom","mother","Mum","Mmom","Dam","Mother","MUM","MOM","DAM","MOTHER")){
      cat(paste("'mum' appears to be in third column, reordering to 'id','dam','sire'")); flush.console();
      Ped<-Ped[,c(1,3,2)]
    }
    if(names(Ped)[2]%in%c("mum","mom","mother","Mum","Mom","Dam","Mother","MUM","MOM","DAM","MOTHER")){
      names(Ped)[2]<-"dam"
      names(Ped)[3]<-"sire"
    }
    if(names(Ped)[2]!="dam"|names(Ped)[3]!="sire"){
     stop("Unable to identify column names, expecting 'id','dam','sire'")
    }
  }
  


  if(is.null(cohorts)==FALSE&&length(Ped[,1])!=length(cohorts)) 
    stop("Pedigree and cohorts differ in length.")

  if(is.null(dat)==FALSE&&is.numeric(dat)&&length(Ped[,1])!=dim(as.data.frame(dat))[1]) 
    stop("Pedigree and available data differ in length.")

  if(is.null(sex)==FALSE&&length(Ped[,1])!=length(sex)) 
    stop("Pedigree and sex differ in length.")



  if(is.null(cohorts)) cohorts<-kindepth(Ped[,1],Ped[,2],Ped[,3])
  cohorts<-as.numeric(as.character(cohorts))
  Ped$cohorts<-cohorts
  plotWidth<-1
  plotHeight<-1
  cohortSizes<-table(Ped$cohorts)
  names(cohortSizes)<-names(table(Ped$cohorts))

  cohortIndex<-array(1,dim=length(cohortSizes))
  names(cohortIndex)<-names(table(Ped$cohorts))

  Ped$xlocs<-NA
  Ped$ylocs<-NA

  scaledCohorts=NULL

  for(x in 1:length(Ped[,1])) {
    Ped$xlocs[x]<-cohortIndex[as.character(Ped$cohorts[x])]/(plotWidth*cohortSizes[as.character(Ped$cohort[x])])-1/(2*plotWidth*cohortSizes[as.character(Ped$cohorts[x])])
    scaledCohorts<-cohorts-min(cohorts)
    scaledCohorts<-scaledCohorts/max(scaledCohorts)
    scaledCohorts<-(scaledCohorts/1.1)+0.05
    Ped$ylocs[x]<-plotHeight-scaledCohorts[x] 
    cohortIndex[as.character(Ped$cohorts[x])]<-cohortIndex[as.character(Ped$cohorts[x])]+1
  }
  Ped$xlocs<-Ped$xlocs*0.94+0.03
  if(writeCohortLabels=='y') Ped$xlocs<-Ped$xlocs*0.95+0.05

  Ped$matxlocs<-Ped$xlocs[match(Ped$dam,Ped$id)]
  Ped$matylocs<-Ped$ylocs[match(Ped$dam,Ped$id)]

  Ped$patxlocs<-Ped$xlocs[match(Ped$sire,Ped$id)]
  Ped$patylocs<-Ped$ylocs[match(Ped$sire,Ped$id)]

  if(dots=='y') {
    Ped$dots<-"black"
    if(is.null(sex)==FALSE) {
      for(x in 1:length(sex)) {
        if(sex[x]==sexInd[2]) Ped$dots[x]<-sexColours[1]
        if(sex[x]==sexInd[1]) Ped$dots[x]<-sexColours[2]
      }
    }
  }

  if(dataDots=='y') {
    if(is.numeric(dat)) dataPed<-subset(Ped, rowSums(as.data.frame(as.numeric(dat)))>0)
    if(is.character(dat)) dataPed<-subset(Ped, Ped[,1]%in%dat)
  }

  Ped.subset<-NULL

  if(is.null(dat)==FALSE) {

    if(is.numeric(dat)){
      avail<-rowSums(as.data.frame(dat))
      keep<-subset(Ped$id,avail>0)
    } else {
      keep<-dat
    }

    if(retain=='informative'){
      Ped.subset<-prune(Ped,keep,make.base=TRUE)
    }else{
      Ped.subset<-prune(Ped,keep,make.base=FALSE)
    }
 
    Ped.subset$xlocs<-Ped$xlocs[match(Ped.subset$id,Ped$id)]
    Ped.subset$ylocs<-Ped$ylocs[match(Ped.subset$id,Ped$id)]
    Ped.subset$matxlocs<-Ped.subset$xlocs[match(Ped.subset$dam,Ped.subset$id)]
    Ped.subset$matylocs<-Ped.subset$ylocs[match(Ped.subset$dam,Ped.subset$id)]
    Ped.subset$patxlocs<-Ped.subset$xlocs[match(Ped.subset$sire,Ped.subset$id)]
    Ped.subset$patylocs<-Ped.subset$ylocs[match(Ped.subset$sire,Ped.subset$id)]
    Ped.subset$dots<-Ped$dots[match(Ped.subset$id,Ped$id)]
  }

  Ped.focal<-NULL
  if(is.null(focal)==FALSE){
    Ped.focal<-subset(Ped,Ped$id==as.character(focal[1]))
    if(focal[2]=='offspring'){
      o<-subset(Ped,as.character(Ped$dam)==as.character(focal[1])|
                    as.character(Ped$sire)==as.character(focal[1]))
      Ped.focal<-rbind(Ped.focal,Ped[which(Ped$id %in% o$id),])
    }
    if(focal[2]=='parents'){
      f<-subset(Ped,Ped$id==as.character(focal[1]))
      if(is.na(f$sire)==FALSE) Ped.focal<-rbind(Ped.focal,subset(Ped,Ped$id==f$sire))
      if(is.na(f$dam)==FALSE) Ped.focal<-rbind(Ped.focal,subset(Ped,Ped$id==f$dam))
    }


    if(focal[2]=='ancestors'|focal[2]=='kin'){
      for(x in 1:10) {
        Ped.focal<-rbind(Ped.focal,Ped[which(Ped$id %in% Ped.focal$dam),])
        Ped.focal<-rbind(Ped.focal,Ped[which(Ped$id %in% Ped.focal$sire),])
      }
    }


    if(focal[2]=='descendents'|focal[2]=='kin'){
      for(x in 1:10) {
        Ped.focal<-rbind(Ped.focal,Ped[which(Ped$dam %in% Ped.focal$id),])
        Ped.focal<-rbind(Ped.focal,Ped[which(Ped$sire %in% Ped.focal$id),])
      }
    }




    Ped.focal$matxlocs.focal<-Ped.focal$xlocs[match(Ped.focal$dam,Ped.focal$id)]
    Ped.focal$matylocs.focal<-Ped.focal$ylocs[match(Ped.focal$dam,Ped.focal$id)]

    Ped.focal$patxlocs.focal<-Ped.focal$xlocs[match(Ped.focal$sire,Ped.focal$id)]
    Ped.focal$patylocs.focal<-Ped.focal$ylocs[match(Ped.focal$sire,Ped.focal$id)]
  }

  if(is.null(dat)==FALSE) {
    cat(paste("Individuals in full pedigree:",length(Ped[,1]),"\n")); flush.console();
    cat(paste("Individuals in informative pedigree subset:",length(Ped.subset[,1]),"\n")); flush.console();
  }

  if(is.null(dat)&is.null(focal)){
    if(links=='all'|links=='mums')
      grid.segments(Ped$xlocs,Ped$ylocs,Ped$matxlocs,Ped$matylocs,gp=gpar(col=sexColours[1]))
    if(links=='all'|links=='dads')
    grid.segments(Ped$xlocs,Ped$ylocs,Ped$patxlocs,Ped$patylocs,gp=gpar(col=sexColours[2]))
    if(dots=='y') for(x in 10:0) grid.circle(Ped$xlocs,Ped$ylocs,r=dotSize,gp=gpar(col=Ped$dots,fill=Ped$dots))
  }else{
    if(plotfull=='y'&sexColours[1]=='red') {
      if(links=='all'|links=='mums')
           grid.segments(Ped$xlocs,Ped$ylocs,Ped$matxlocs,Ped$matylocs,gp=gpar(col="gray"))
      if(links=='all'|links=='dads')
           grid.segments(Ped$xlocs,Ped$ylocs,Ped$patxlocs,Ped$patylocs,gp=gpar(col="gray"))
    }
    if(plotfull=='y'&sexColours[1]!='red') {
      if(links=='all'|links=='mums')
           grid.segments(Ped$xlocs,Ped$ylocs,Ped$matxlocs,Ped$matylocs,gp=gpar(col=colours()[354]))
      if(links=='all'|links=='dads')
           grid.segments(Ped$xlocs,Ped$ylocs,Ped$patxlocs,Ped$patylocs,gp=gpar(col=colours()[354]))
    }
    if(is.null(dat)==FALSE&(links=='all'|links=='mums'))
         grid.segments(Ped.subset$xlocs,Ped.subset$ylocs,Ped.subset$matxlocs,Ped.subset$matylocs,gp=gpar(col=sexColours[1]))
    if(is.null(dat)==FALSE&(links=='all'|links=='dads'))
         grid.segments(Ped.subset$xlocs,Ped.subset$ylocs,Ped.subset$patxlocs,Ped.subset$patylocs,gp=gpar(col=sexColours[2]))

    if(is.null(focal)==FALSE&(links=='all'|links=='mums'))
         grid.segments(Ped.focal$xlocs,Ped.focal$ylocs,Ped.focal$matxlocs.focal,Ped.focal$matylocs.focal,gp=gpar(col=sexColours[1]))
    if(is.null(focal)==FALSE&(links=='all'|links=='dads'))
         grid.segments(Ped.focal$xlocs,Ped.focal$ylocs,Ped.focal$patxlocs.focal,Ped.focal$patylocs.focal,gp=gpar(col=sexColours[2]))

    if(dots=='y') grid.circle(Ped$xlocs,Ped$ylocs,r=dotSize,gp=gpar(col=Ped$dots,fill=Ped$dots))
    if(dataDots=='y') grid.circle(dataPed$xlocs,dataPed$ylocs,r=dotSize*dataDots.cex,gp=gpar(col='black',fill='black'))
  }
  if(writeCohortLabels=='y'){
    grid.text(names(table(cohorts)),x=rep(0.05,length(table(cohorts))),y=plotHeight-as.numeric(names(table(scaledCohorts))),gp=gpar(cex=cohortLabs.cex))
  }
}

