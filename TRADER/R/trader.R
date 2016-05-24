##########################################################################################
##########################################################################################
##########################################################################################
# OVERAL FUNCTIONS

##########################################################################################
# Do separate analysis by Fraver & White called absolute increase.
# NEEDS data
# OPTIONAL ...
# RETURNS do figures and tables with prefix
absoluteIncreaseALL<- function(data,abs=NULL,abs.threshold=NULL,m1=10,m2=10,buffer=2,prefix="ai",
                               drawing=TRUE,gfun=mean, length=2, storedev=jpeg,...) {
  
  if ( is.null(abs) )
    abs<-absIncrease(data,m1,m2)
  if ( is.null(abs.threshold) )
    abs.threshold<- absTreshold(abs)
  releases <- absoluteIncrease(data,abs,abs.threshold,buffer=buffer,gfun=gfun,length=length)
  #figures
  if (drawing){
    for(i in 1:length(data)) {      
      storedev(paste(prefix,"_",gsub("/","_",names(data)[i]),".",deparse(substitute(storedev)),sep=""),pointsize = 10)
      plotRelease(data,abs,releases, i, method="FraverWhite",addHLines=c(abs.threshold),...)
      dev.off()
    }
  }
  #tables
  write.table ( abs, paste(prefix,"_change.csv", sep = ""), row.names=F,sep="\t")
  #write.table ( releases$releases, paste(prefix,"_releases_tops.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$all_releases, paste(prefix,"_releases_all.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), row.names=F,sep="\t")
  #return(releases)
}

##########################################################################################
# Do separate analysis by Nowacki & Abrams called growth averaging.
# NEEDS dplr data
# OPTIONAL ...
# RETURNS do figures and tables with prefix
growthAveragingALL<- function(data,releases=NULL,m1=10,m2=10, buffer=2, drawing=TRUE,
                               criteria=0.25, criteria2=0.5, prefix="ga",gfun=mean,
                              length=2, storedev=jpeg, ...) {
  if ( is.null(releases) )
    releases<-noblabrams(data,m1=m1,m2=m2,buffer=buffer,length=length,
                        criteria=criteria,criteria2=criteria2,black=FALSE,gfun=gfun)
  
  #figures
  if (drawing){
    for(i in 1:length(data)) {
      storedev(paste(prefix,"_",gsub("/","_",names(data)[i]),".",deparse(substitute(storedev)),sep=""))
      plotRelease(data,releases$change,releases, i, method="NowackiAbrams",
                   addHLines=c(criteria,criteria2),...)
      dev.off()
    }
  }
  #tables
  releases$onlyModerate$year=rownames(releases$onlyModerate)
  releases$onlyModerate = releases$onlyModerate[,c(ncol(releases$onlyModerate),1:(-1+ncol(releases$onlyModerate)))]
  releases$onlyMajor$year=rownames(releases$onlyMajor)
  releases$onlyMajor = releases$onlyMajor[,c(ncol(releases$onlyMajor),1:(-1+ncol(releases$onlyMajor)))]
  
  
  write.table ( releases$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releases$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
  #return(releases)
}

##########################################################################################
# Do separate analysis by Black & Abrams called boundary line.
# NEEDS dplr data
# OPTIONAL ...
# RETURNS do figures and tables with prefix
boundaryLineALL<- function(data,releases=NULL,m1=10,m2=10, boundary=NULL, buffer=2, 
                               criteria=0.2, criteria2=0.5, segment=0.5, segment2=0.5,
                           prefix="bl", drawing=TRUE,gfun=mean,length=2, 
                           notop=10,notop2=10,storedev=jpeg,...) {
  if ( is.null(releases) )
    releases<-noblabrams(data,m1=m1,m2=m2,boundary=boundary,buffer=buffer,gfun=gfun,length=length,
                        criteria=criteria,criteria2=criteria2,segment=segment,
                     segment2=segment2,storedev=storedev,notop=notop,notop2=notop2,black=TRUE)
  
  #figures
  if (drawing){
    for(i in 1:length(data)) {
      #print(paste(i))
      storedev(paste(prefix,"_",gsub("/","_",names(data)[i]),".",deparse(substitute(storedev)),sep=""))
      plotRelease(data,releases$change,releases, i, method="BlackAbrams",
                   addHLines=c(criteria, criteria2))
      dev.off()
    }
  }
  #tables
  releases$onlyModerate$year=rownames(releases$onlyModerate)
  releases$onlyModerate = releases$onlyModerate[,c(ncol(releases$onlyModerate),1:(-1+ncol(releases$onlyModerate)))]
  releases$onlyMajor$year=rownames(releases$onlyMajor)
  releases$onlyMajor = releases$onlyMajor[,c(ncol(releases$onlyMajor),1:(-1+ncol(releases$onlyMajor)))]
  
  
  write.table ( releases$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releases$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
  #return(releases)
}

##########################################################################################
# Do separate analysis by Splechtna.
# NEEDS dplr data
# OPTIONAL ...
# RETURNS do figures and tables with prefix
splechtnaALL<- function(data, releases=NULL,m1=10,m2=10, boundary=NULL, buffer=2, drawing=TRUE,
                            criteria=0.2, criteria2=0.5,segment=0.5, segment2=0.5,prefix="sp",
                        gfun=mean,length=2,notop=10,notop2=10,storedev=jpeg,...) {
  if ( is.null(releases) )
    releases<-splechtna(data,m1=m1,m2=m2,boundary=boundary,buffer=buffer,gfun=gfun,
                           criteria=criteria,criteria2=criteria2,segment=segment,
                        segment2=segment2,length=length,notop=notop,notop2=notop2,storedev=storedev)
  
  #figures
  if (drawing){
    for(i in 1:length(data)) {
      storedev(paste(prefix,"_",gsub("/","_",names(data)[i]),".",deparse(substitute(storedev)),sep=""))
      plotRelease(data,releases$change,releases, i, method="Splechtna",
                   addHLines=c(criteria, criteria2),...)
      dev.off()
    }
  }
  #tables
  releases$onlyModerate$year=rownames(releases$onlyModerate)
  releases$onlyModerate = releases$onlyModerate[,c(ncol(releases$onlyModerate),1:(-1+ncol(releases$onlyModerate)))]
  releases$onlyMajor$year=rownames(releases$onlyMajor)
  releases$onlyMajor = releases$onlyMajor[,c(ncol(releases$onlyMajor),1:(-1+ncol(releases$onlyMajor)))]
  
  
  write.table ( releases$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releases$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releases$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releases$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
  #return(releases)
}


##########################################################################################
# Do all analysis and produce 4 figures in one
# NEEDS dplr data
# RETURNS do figures and tables with prefix
doAll<- function(data,m1=10,m2=10, boundary=NULL, buffer=2, 
                 criteria=0.2, criteria2=0.5, segment=0.5,segment2=0.5,prefix="all", 
                 gfun=mean,length=2, notop=10,notop2=10,storedev=jpeg,
                 drawing=TRUE,...) {
  
  abs<-absIncrease(data,m1,m2)
  abs.threshold<- absTreshold(abs)
  releasesFW <- absoluteIncrease(data,abs,abs.threshold,m1=m1,m2=m2,buffer=buffer,gfun=gfun,
                                 length=length)
  releasesNA<-noblabrams(data,m1=m1,m2=m2,buffer=buffer,criteria=criteria+0.05,
                        criteria2=criteria2,black=FALSE,gfun=gfun,length=length)
  releasesBA<-noblabrams(data,m1=m1,m2=m2,boundary=boundary,buffer=buffer,criteria=criteria,
                        criteria2=criteria2,segment=segment,black=TRUE,gfun=gfun,length=length
                     ,segment2=segment2,notop=notop,notop2=notop2,storedev=storedev)
  releasesS<-splechtna(data,m1=m1,m2=m2,boundary=boundary,buffer=buffer,criteria=criteria,
                          criteria2=criteria2,segment=segment,gfun=gfun,length=length
                       ,segment2=segment2,notop=notop,notop2=notop2,storedev=storedev)
  
  if (drawing) {
    for(i in 1:length(data)) {      
      storedev(paste(prefix,"_",gsub("/","_",names(data)[i]),".",deparse(substitute(storedev)),sep=""))
      par(mfrow=c(2,2))      
      plotRelease(data,abs,releasesFW, i, method="FraverWhite",addHLines=c(abs.threshold),...)
      plotRelease(data,releasesNA$change,releasesNA, i, method="NowackiAbrams",
                   addHLines=c(criteria+0.05,criteria2),plotfirst=FALSE,...)
      plotRelease(data,releasesBA$change,releasesBA, i, method="BlackAbrams",
                   addHLines=c(criteria, criteria2),plotfirst=FALSE,...)
      plotRelease(data,releasesS$change,releasesS, i, method="Splechtna",
                   addHLines=c(criteria, criteria2),plotfirst=FALSE,...)
      dev.off()
    }
  }
  #tables
  
  prefix<-"ai"
  write.table ( abs, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releasesFW$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesFW$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesFW$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""), 
                row.names=F,sep="\t")
  write.table ( releasesFW$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), 
                row.names=F,sep="\t")
  prefix<-"ga"
  releasesNA$onlyModerate$year=rownames(releasesNA$onlyModerate)
  releasesNA$onlyModerate = releasesNA$onlyModerate[,c(ncol(releasesNA$onlyModerate),1:(-1+ncol(releasesNA$onlyModerate)))]
  releasesNA$onlyMajor$year=rownames(releasesNA$onlyMajor)
  releasesNA$onlyMajor = releasesNA$onlyMajor[,c(ncol(releasesNA$onlyMajor),1:(-1+ncol(releasesNA$onlyMajor)))]
  
  
  write.table ( releasesNA$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releasesNA$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesNA$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesNA$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""),
                row.names=F,sep="\t")
  write.table ( releasesNA$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), 
                row.names=F,sep="\t")
  write.table ( releasesNA$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releasesNA$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
  
  prefix<-"bl"
  releasesBA$onlyModerate$year=rownames(releasesBA$onlyModerate)
  releasesBA$onlyModerate = releasesBA$onlyModerate[,c(ncol(releasesBA$onlyModerate),1:(-1+ncol(releasesBA$onlyModerate)))]
  releasesBA$onlyMajor$year=rownames(releasesBA$onlyMajor)
  releasesBA$onlyMajor = releasesBA$onlyMajor[,c(ncol(releasesBA$onlyMajor),1:(-1+ncol(releasesBA$onlyMajor)))]
  
  
  write.table ( releasesBA$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releasesBA$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesBA$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)
  write.table ( releasesBA$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""),
                row.names=F,sep="\t")
  write.table ( releasesBA$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), 
                row.names=F,sep="\t")
  write.table ( releasesBA$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releasesBA$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
  
  prefix<-"sp"
  releasesS$onlyModerate$year=rownames(releasesS$onlyModerate)
  releasesS$onlyModerate = releasesS$onlyModerate[,c(ncol(releasesS$onlyModerate),1:(-1+ncol(releasesS$onlyModerate)))]
  releasesS$onlyMajor$year=rownames(releasesS$onlyMajor)
  releasesS$onlyMajor = releasesS$onlyMajor[,c(ncol(releasesS$onlyMajor),1:(-1+ncol(releasesS$onlyMajor)))]
  
  write.table ( releasesS$change, paste(prefix,"_change.csv", sep = ""), sep="\t",row.names=F)
  #write.table ( releasesS$releases, paste(prefix,"_releases_tops.csv", sep = ""), sep="\t",row.names=F)  
  write.table ( releasesS$all_releases, paste(prefix,"_releases_all.csv", sep = ""), sep="\t",row.names=F)  
  write.table ( releasesS$years_list_total, paste(prefix,"_releases_years_total.csv", sep = ""),
                row.names=F,sep="\t")
  write.table ( releasesS$pgc, paste(prefix,"_releases_values_total.csv", sep = ""), 
                row.names=F,sep="\t")
  write.table ( releasesS$onlyMajor, paste(prefix,"_releases_Only_Major.csv", sep = ""), row.names=F,sep="\t")
  write.table ( releasesS$onlyModerate, paste(prefix,"_releases_Only_Moderate.csv", sep = ""), row.names=F,sep="\t")
}

##########################################################################################
# Fraver and White 2005 analysis.
# NEEDS dplr data OR absIncrease, abs threshold
# OPTIONAL abs,abs threshold, m1, m2, buffer
# RETURNS release table and years
absoluteIncrease<- function(data, abs=NULL, abs.threshold=NULL, m1=10, m2=10, 
                            buffer=2,gfun=mean, length=2) {
#data<-mdata;m1=10;m2=10;buffer=4;gfun=mean;length=3; abs<-NULL;abs.threshold<-NULL
  print(paste("## Fraver & White analysis!"))
  if ( is.null(abs) )
    abs<-absIncrease(data,m1,m2,gfun=gfun)
  if ( is.null(abs.threshold) )
    abs.threshold<- absTreshold(abs)

  print(paste("Absolute threshold ",round(abs.threshold,3),"m1",m1,"m2",m2,
              "Buffer",buffer,"Length",length))
  
  releases1 <- abs[,1][-c(1)]
  for(t in 2:length(abs)){ #loop across series (trees)
    temp <- c()
    for(f in 2:length(abs[,t])){ #loop within a single series
      # test if all three points are not NA
      if( !is.na(abs[f,t]) & !is.na(abs[f-1,t]) & !is.na(abs[f+1,t])) { #TOP3
        # store only focal point if it is the highest
        if(abs[f,t] > abs.threshold & abs[f,t] > abs[f-1, t] & abs[f,t]>abs[f+1,t]){
          # store in rel only if
          rel <- abs[f,t]
        } else {
          rel <- NA
        }
      } else if ( ( (f==length(abs[,t])) | ( sum(!is.na(abs[f+1:length(abs[,t]),t]))==0 ) ) 
                  & !is.na(abs[f,t]) & !is.na(abs[f-1,t]) &
                    (abs[f,t] > abs.threshold) & (abs[f,t] > abs[f-1, t])  ) {
        rel <- abs[f,t] # tail
      } else if ( ( (f==2) | ( sum(!is.na(abs[1:(f-2),t]))==0 ) ) 
                  & !is.na(abs[f,t]) & !is.na(abs[f-1,t]) &
                    (abs[f-1,t] > abs.threshold) & (abs[f-1,t] > abs[f, t])  ) {
        rel <- abs[f-1,t] # head
      } else {
        rel <- NA
      }
      temp <- append(temp, rel)
    }
    releases1 <- as.data.frame(cbind(releases1, temp))
  }  
  names(releases1)[2:length(names(releases1))] <- names(abs)[-c(1)]
  
  #all releases above threshold
  releases1All <- abs[-c(1),]
  releases1AllInc <- releases1All > abs.threshold 
  for(t in 2:length(releases1All)){ #loop across series
    releases1All[,t]<- releases1AllInc[,t]*releases1All[,t]
    releases1All[,t]<- ifelse(releases1All[,t]==0, NA,releases1All[,t])
  }
  rm(releases1AllInc,t)
    
  release_list <- reduceByLB(releases=releases1,above=releases1All,
                             buffer=buffer,length=length,type=1)
  
  release_list_vals <- releases1
  for ( t in 2:length(releases1) ){
    release_list_vals[,t]<-ifelse(release_list_vals[,1] %in% release_list[[t-1]],
                                  release_list_vals[,t],NA)
  }
   
  rs<-writeReleaseStats(release_list,"Total number of releases is")
    
  norel<-plotNORelease(data,rs, criteria=round(abs.threshold,3),prefix="relai")
  
  return( list("releases" = releases1,"years"=release_list, 
               "years_list_total" =norel, "pgc"=release_list_vals,
               "all_releases"=releases1All) )    
}

##########################################################################################
# Nowacki and Abrams 1997, Black and Abrams 2003 or "pure boundary line"
# NEEDS dplr data
# OPTIONAL buffer, length, m1, m2, criteria, criteria2, black, prior, change
# RETURNS release table, years and pgc
# black=FALSE -> Nowacki and Abrams 1997
# black=TRUE -> Black and Abrams 2003
noblabrams<-function(data=NULL,prior=NULL,change=NULL,m1=10,m2=10,boundary=NULL,
                    buffer=2,criteria=0.25,criteria2=0.5,segment=0.5,segment2=0.5,
                 black=FALSE,gfun=mean, length=2,notop=10,notop2=10,storedev=jpeg){
  #data=mdata;prior=NULL;change=NULL;m1=10;m2=10;boundary=NULL;buffer=3;criteria=0.25;criteria2=0.5;segment=0.5;segment2=0.5;black=TRUE;gfun=mean;length=4;storedev=jpeg;gfun=mean;notop=10;notop2=10
  if ( black ) {
    print(paste("## Black & Abrams analysis!"))
    print(paste("Criteria",criteria,"Criteria2",criteria2, "m1",m1,"m2",m2, 
                "Buffer",buffer,"Length",length,"Segment",segment,"Segment2",segment2))
  } else{
    print(paste("## Nowacki & Abrams analysis!"))
    print(paste("Criteria",criteria,"Criteria2",criteria2, "m1",m1,"m2",m2, 
                "Buffer",buffer,"Length",length))
  }
    
  if ( is.null(prior) )
    prior <- priorGrowth(data, m1=m1, m2=m2, gfun=gfun) #prior growth   
  if ( is.null(change) )
    change<-PGC(data, m1=m1, m2=m2, gfun=gfun) #percent growth change
    
  names(change)[2:length(names(change))] <- names(data)  
  
  if ( black ) { # black=TRUE -> Black and Abrams 2003
    
    bo<-boundaryGet(data,prior=prior,change=change,m1=m1,m2=m2,segment=segment,
                    segment2=segment2,gfun=gfun,notop=notop,notop2=notop2)
    if ( is.null(boundary) ) {
      bofit<-boundaryFit(bo$bo,bo$x,bo$y)
      boundary<-bofit$fun
      rsq<-bofit$rsq
    } else {
      #boundaryFit(bo$bo,bo$x,bo$y,boundary=boundary)
      rsq<-NULL
    }
    plotBoundary(bo$bo,bo$x,bo$y,boundary=boundary,rsq=rsq,criteria=criteria,
                 criteria2=criteria2,storedev=storedev)
        
    scaled <- c(prior[1:length(change[,1]),1]) 
    for(i in 2:length(prior)){
      temp <- c()
      for(n in 1:length(prior[,1])){
        if(!is.na(prior[n,i])){
          mval <- boundary(prior[n,i])
          if( is.na(mval) | (mval == 0)){
            temp1<-NA
          }  else
            temp1 <- change[n,i]/mval
        } else {
          temp1 <- NA
        }
        temp <- append(temp, temp1)
      }
      scaled <- as.data.frame(cbind(scaled, temp))
    }
    names(scaled)[2:length(names(scaled))] <- names(data)
    #Count releases that fill the criteria given at the setup
    releases3 <- PGCreleases(scaled,criteria)
    
    releases3All <- scaled
    releases3AllInc <- releases3All >= criteria
    for(t in 2:length(releases3All)){ #loop across series
      releases3All[,t]<- releases3AllInc[,t]*releases3All[,t]
      releases3All[,t]<- ifelse(releases3All[,t]==0, NA,releases3All[,t])
    }
    
    #all releases above threshold    
    releases3All2 <- scaled
    releases3AllInc <- releases3All2 >= criteria2
    for(t in 2:length(releases3All2)){ #loop across series
      releases3All2[,t]<- releases3AllInc[,t]*releases3All2[,t]
      releases3All2[,t]<- ifelse(releases3All2[,t]==0, NA,releases3All2[,t])
    }
    
    rm(releases3AllInc,t)
  } else {
    #Count releases that fill the criteria given at the setup
    releases3 <- PGCreleases(change,criteria)
    
    releases3All <- change
    releases3AllInc <- releases3All >= criteria    
    for(t in 2:length(releases3All)){ #loop across series
      releases3All[,t]<- releases3AllInc[,t]*releases3All[,t]
      releases3All[,t]<- ifelse(releases3All[,t]==0, NA,releases3All[,t])
    }    
    
    releases3All2 <- change
    releases3AllInc <- releases3All2 >= criteria2
    for(t in 2:length(releases3All2)){ #loop across series
      releases3All2[,t]<- releases3AllInc[,t]*releases3All2[,t]
      releases3All2[,t]<- ifelse(releases3All2[,t]==0, NA,releases3All2[,t])
    }    
    
    rm(releases3AllInc,t)
  }
  
  names(releases3)[2:length(names(releases3))] <- names(data)
  
  #Drop consequtive releases for which the difference in years is < buffer
  release_list3 <- reduceByLB(releases=releases3,above=releases3All,
                              buffer=buffer,length=length,type=1)
  release_list32 <- reduceByLB(releases=releases3,above=releases3All2,
                              buffer=buffer,length=length,type=1)
  
  
  # remove major relese from the moderate
  all_release_list3 <- release_list3
  release_list3<-removeMajorFromModerate(release_list3,release_list32,0)  
  
  #List percent growth changes to separate moderate and major releases
  release_list3_pgc <- reduceByLB(releases=releases3,above=releases3All,
                                      buffer=buffer,length=length,type=2)
  release_list32_pgc <- reduceByLB(releases=releases3,above=releases3All2,
                                      buffer=buffer,length=length,type=2)
  
  # remove major relese from the moderate
  all_release_list3_pgc <- release_list3_pgc
  release_list3_pgc<-removeMajorFromModerate(release_list3,release_list32,0.0,
                                             release_list3_pgc)  
  
  rs<-writeReleaseStats(release_list3,paste("Total number of releases >=",criteria,
                                            "& <",criteria2,"is"))
  rs2<-writeReleaseStats(release_list32,paste("Total number of releases >=",criteria2,"is"))
                        
  if(black) {
    mpref<-"bl"
  } else
    mpref<-"ga"
  
  norel<-plotNORelease(data,rs,rs2, criteria=criteria, criteria2=criteria2,
                   prefix=paste("rel",mpref,sep=""))
  
  release_list_vals <- releases3
  for ( t in 2:length(releases3) ){
    release_list_vals[,t]<-ifelse(release_list_vals[,1] %in% all_release_list3[[t-1]],
                                  release_list_vals[,t],NA)
  }
  
  if ( black )
    return ( list("releases" = releases3,"years"=all_release_list3, "change" = scaled ,
                  "pgc"=release_list_vals,"years_list_total" =norel,
                  "all_releases"=releases3All,
                  "onlyModerate"=relListToDataFrame(release_list3,data),
                  "onlyMajor"=relListToDataFrame(release_list32,data) 
                  ) )
  else
    return ( list("releases" = releases3,"years"=all_release_list3, "change" = change ,
                  "pgc"=release_list_vals,"years_list_total" =norel,
                  "all_releases"=releases3All,
                  "onlyModerate"=relListToDataFrame(release_list3,data),
                  "onlyMajor"=relListToDataFrame(release_list32,data) 
                  ) )
}

##########################################################################################
# Splechtna et al. 2005 type of release analysis
# first filter out growth changes > splechtna (as defined in the setup)
# NEEDS dplr data, percent growth change, prior
# OPTIONAL buffer, m1,m2,criteria
# RETURNS release table, years and pgc
splechtna<-function(data,change=NULL,prior=NULL,m1=10,m2=10,boundary=NULL,buffer=2,
                       criteria=0.2,criteria2=0.5,segment=0.5,gfun=mean, length=2, 
                    segment2=0.5, notop=10,notop2=10,storedev=jpeg){
 #data=mdata;prior=NULL;change=NULL;m1=10;m2=10;boundary=NULL;buffer=2;criteria=0.25;criteria2=0.5;segment=0.2;gfun=mean;length=2
  print(paste("## Splechtna analysis!"))
  if ( is.null(prior) )
    prior <- priorGrowth(data, m1=m1, m2=m2,gfun=gfun) #prior growth   
  if ( is.null(change) )
    change<- PGC(data, m1=m1, m2=m2,gfun=gfun) #percent growth change
  
  print(paste("Criteria",criteria,"Criteria2",criteria2, "m1",m1,"m2",m2,
              "Buffer",buffer,"Length",length,"Segment",segment,"Segment2",segment2))
  
  releases3All <- change
  releases3AllInc <- releases3All >= criteria
  for(t in 2:length(releases3All)){ #loop across series
    releases3All[,t]<- releases3AllInc[,t]*releases3All[,t]
    releases3All[,t]<- ifelse(releases3All[,t]==0, NA,releases3All[,t])
  }
  releases3All2 <- change
  releases3AllInc <- releases3All2 >= criteria2
  for(t in 2:length(releases3All2)){ #loop across series
    releases3All2[,t]<- releases3AllInc[,t]*releases3All2[,t]
    releases3All2[,t]<- ifelse(releases3All2[,t]==0, NA,releases3All2[,t])
  }
  rm(releases3AllInc,t)
  
  releases_filter <- PGCreleasesSplechtna(change, criteria)
  filter_years <- reduceByLB(releases=releases_filter,above=releases3All,length=length,
                             buffer=buffer,type=1)#years(releases_filter,buffer)
  filter_pgc <- reduceByLB(releases=releases_filter,above=releases3All,length=length,
                           buffer=buffer,type=2)#PGCrel(releases_filter,buffer=buffer)
  filter_pg <- reduceByLB(releases=releases_filter,above=releases3All,length=length,
                          buffer=buffer,val=prior)#PG(releases_filter, prior, buffer)
  
  bo<-boundaryGet(data,prior=prior,change=change,m1=m1,m2=m2,segment=segment,segment2=segment2,
                  gfun=gfun,notop=notop,notop2=notop2)
  if ( is.null(boundary) ) {
    bofit<-boundaryFit(bo$bo,bo$x,bo$y)
    boundary<-bofit$fun
    rsq<-bofit$rsq
  } else {    
    rsq<-NULL
  }
  plotBoundary(bo$bo,bo$x,bo$y,boundary=boundary,rsq=rsq,criteria=criteria,
               criteria2=criteria2,storedev=storedev)
  
  scaled <- c(prior[1:length(change[,1]),1]) 
  for(i in 2:length(prior)){
    temp <- c()
    for(n in 1:length(prior[,1])){
      if(!is.na(prior[n,i])){
        mval<-boundary(prior[n,i])
        if( is.na(mval) | (mval == 0)){
          temp1<-NA
        }  else
        temp1 <- change[n,i]/mval
      } else {
        temp1 <- NA
      }
      temp <- append(temp, temp1)
    }
    scaled <- as.data.frame(cbind(scaled, temp))
  }
  names(scaled)[2:length(names(scaled))] <- names(data)
  names(change)[2:length(names(change))] <- names(data)
  
  #scale the candidate releases to the boundary line
  release_list4 <- as.list(rep(list(0),length(filter_years)))
  release_list4_pgc <- as.list(rep(list(0),length(filter_years)))
  
  for(i in 1:length(filter_years)){
    
    for(n in 3:length(filter_years[[i]])){
      mval<-boundary(filter_pg[[i]][n])
      if( is.na(mval) | (mval== 0) ){
        temp<-NA
      } else
        temp <- filter_pgc[[i]][n]/mval
      temp2 <- c()
      temp3 <- c()
      if(!is.na(temp)){
        if(temp > criteria){
          temp2 <- append(temp2, filter_years[[i]][n])
          mval<-boundary(filter_pg[[i]][n])
          if ( is.na(mval) | (mval ==0) ) {
            temp3 <- append(temp3,NA)
          }else
            temp3 <- append(temp3, filter_pgc[[i]][n]/mval)
        }
      }
      release_list4[[i]] <- append(release_list4[[i]], temp2)
      release_list4_pgc[[i]] <- append(release_list4_pgc[[i]], temp3)
    }
  }
  
  
  #rs<-writeReleaseStats(release_list4,release_list4_pgc,criteria=criteria, criteria2=criteria2,
  #                      buffer=buffer,length=length)
  
  release_list42<-release_list4
  for (i in 1:length(release_list4)){
    release_list42[[i]] <- release_list42[[i]] [release_list4_pgc[[i]] >= criteria2]
  }

  # remove major relese from the moderate
  all_release_list4 <- release_list4
  release_list4<-removeMajorFromModerate(release_list4,release_list42,0)    
  
  rs<-writeReleaseStats(release_list4,paste("Total number of releases >=",criteria,
                                            "& <",criteria2,"is"))
  rs2<-writeReleaseStats(release_list42,paste("Total number of releases >=",criteria2,"is"))
  
  norel<-plotNORelease(data,rs,rs2, criteria=criteria, criteria2=criteria2,prefix="relsp")
  
  inyears<-unlist(release_list4)
  inyears<-inyears[inyears!=0]
  
  release_list_vals <- change
  for ( t in 2:length(change) ){
    release_list_vals[,t]<-ifelse(release_list_vals[,1] %in% all_release_list4[[t-1]],
                                  release_list_vals[,t],NA)
  }
  
  return ( list("releases" = change,"years"=all_release_list4, "change" = scaled ,
                "pgc"=release_list_vals,"years_list_total" =norel,
                "all_releases"=releases3All,
                "onlyModerate"=relListToDataFrame(release_list4,data),
                "onlyMajor"=relListToDataFrame(release_list42,data) ) 
    )
}

##########################################################################################
##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################
##########################################################################################
# HELP FUNCTIONS

##########################################################################################
# compute absolute increases
# NEEDS dplr data
# RETURNS absIncrease
absIncrease <- function(data, m1=10, m2=10,gfun=mean){
  series <- names(data)
 
  years <- c( as.numeric(row.names(data)[ m1:(length(data[,1])-m2) ]) )  
  
  for(i in 1:length(data)){  
    inc.abs <- c()    
    for(n in m1:(length(data[,i])-m2)){ 
      if( is.na(data[(n+1-m1), i]) == FALSE & is.na(data[(n+m2),i])==FALSE){
        prior.abs <- gfun(data[(n+1-m1):n,i])
        post.abs <- gfun(data[(n+1):(n+m2),i])
        temp <- post.abs-prior.abs
      } else {
        temp <- NA
      }
      inc.abs <- append(inc.abs, temp)
    }
    years <- as.data.frame(cbind(years, inc.abs))
  }
  names(years)[2:length(names(years))]<-series
  return(years)
}

##########################################################################################
# "Blind" definition of the absolute-increase threshold of Fraver & White 2005 (1.25*standard deviation)
# NEEDS absIncrease
# RETURN absIncrease theshold
absTreshold<-function(abs,tvalue=1.25){
  return ( tvalue * 
             sd( unlist(abs[2:length(names(abs))] ), na.rm = TRUE))
}

##########################################################################################
# get priors
# NEEDS dplr data
# RETURNS prior
priorGrowth <- function(data, m1=10, m2=10,gfun=mean){
  
  series<-names(data)  
  years<-as.numeric(row.names(data)[m1+2:(length(data[,1])-m2-m1)])
  
  for(i in 1:length(data)){  # for all individuals
    prior <- c()
    
    for(n in m1+2:(length(data[,i])-m2-m1)){ 
      if( is.na(data[(n-(m1+1)), i]) == FALSE & is.na(data[(n+m2),i])==FALSE ){
        prior.growth <- gfun(data[(n-m1):n,i])
        temp <- prior.growth
      } else {
          temp <- NA
      }
      prior <- append(prior, temp)
    }
    years <- as.data.frame(cbind(years, prior))
  }
  return( years )
}


##########################################################################################
# get percentage growth change
# NEEDS dplr data
# RETURNS percent growth change
PGC <- function(data, m1=10, m2=10,gfun=mean){
  #data; m1=10; m2=10;gfun=mean
  series<-names(data)
  
  years<-as.numeric(row.names(data)[m1+2:(length(data[,1])-(m2+m1))])
  
  for(i in 1:length(data)){  
    pgc <- c()
    for(n in m1+2:(length(data[,i])-(m2+m1))){ 
      if(!is.na(data[(n+1-m1), i])& !is.na(data[(n+m2+1),i])){
        prior.growth <- gfun(data[(n+1-m1):n,i])
        post.growth <- gfun(data[(n+1):(n+1+m2),i])
        temp <-(post.growth-prior.growth)/prior.growth
      } else {
          temp <- NA
      }
      pgc <- append(pgc, temp)
    }
    years <- as.data.frame(cbind(years, pgc))    
  }
  return (years)
}

##########################################################################################
# NEEDS percent growth change
# RETURNS release years according criteria
PGCreleases <- function(change,criteria=0.2){
  
  releases2 <- c(change[2:length(change[,1]),1])
  
  for(t in 2:length(change)){
    temp <- c()
    
    for(f in 2:length(change[,t])){ 
      if(!is.na(change[f,t]) & !is.na(change[f-1,t]) & !is.na(change[f+1, t])){ #TOP3
        
        if(change[f,t] > criteria & change[f,t] > change[f-1, t] & change[f,t] > change[f+1,t]){
          rel <- change[f,t]
        } else {
          rel <- NA
        }
      } else if ( ( (f==length(change[,t])) | ( sum(!is.na(change[f+1:length(change[,t]),t]))==0 ) ) 
                  & !is.na(change[f,t]) & !is.na(change[f-1,t]) &
                    (change[f,t] > criteria) & (change[f,t] > change[f-1, t])  ) {
        rel <- change[f,t]  
      } else if ( ( (f==2) | ( sum(!is.na(change[1:(f-2),t]))==0 ) ) 
                  & !is.na(change[f,t]) & !is.na(change[f-1,t]) &
                    (change[f-1,t] > criteria) & (change[f-1,t] > change[f, t])  ) {
        rel <- change[f-1,t] # head
      } else {
        rel <- NA
      }
      temp <- append(temp, rel)
    }
    releases2 <- as.data.frame(cbind(releases2, temp))
  }  
  return (releases2)
}

##########################################################################################
# NEEDS percent growth change
# RETURNS release years according criteria
PGCreleasesSplechtna <- function(change,criteria=0.5){  
  return( PGCreleases(change,criteria) )
}

##########################################################################################
# reduce release list according buffer and length
# NEEDS releases
# OPTIONAL buffer, length
# RETURNS release_list
reduceByLB<-function(releases, above, buffer=2, type=1, length=2,val=NULL){
  #releases=releases1;above=releases1All; buffer=2
  #releases=releases3;above=releases3All
  
  if ( length <1 ){
    print(paste("Length must be 1 or more!"))
    return (NULL)
  }
  if ( buffer <1 ){
    print(paste("Buffer must be 1 or more!"))
    return (NULL)
  }
  
  release_list <- ( rep(list(0),length(releases)-1) )
  if ( is.null(val)){    
    if ( type==1 ){ # type =1 years
      for(i in 2:length(releases)){                       
        #temp<-append(c(0,0), releases[,1][!is.na(releases[,i])] )    
        #temp2<-temp[ temp-c(0,temp[-c(length(temp))]) >=buffer]        
        temp<-releases[,1][!is.na(releases[,i])]
        if (length(temp)>0) {
          temp2<-append(c(),temp[1])
          if ( length(temp) > 1){
            for( r in 2:length(temp) ){
              tdif<- temp[r] - temp2[length(temp2)]
              if ( tdif >= buffer ) { # difference big enough
                temp2<-append(temp2,temp[r])
              } else if ( sum ( !is.na(above[ 
                grep(temp2[length(temp2)],above[,1])[1]:grep(temp[r],above[,1])[1]
                  ,i ])) == (tdif+1) ) {# is continuous
                if ( above[grep(temp2[length(temp2)],above[,1])[1],i] 
                   < above[grep(temp[r],above[,1],above[,1])[1],i] ) {
                  temp2[length(temp2)] <- temp[r]
                }                 
              }
            }
          }
        } else
          temp2<-c()        
        
        temp3<-c()
        for (f in temp2 ){
          index<-grep(f,above[,1])[1]
          if (sum(!is.na(above[ max(1,index-length+1):min(dim(above)[1],index+length-1),i])) >=length)
            temp3<-append(temp3,f)
        }
        release_list[[i-1]] <- append(release_list[[i-1]], temp3 )      
      }
    } else { # type =2 values
        for(i in 2:length(releases)){                   
          #temp<-append(c(0,0), releases[,1][!is.na(releases[,i])] ) 
          #tempVal<-append(c(0,0), releases[,i][!is.na(releases[,i])] ) 
          #temp2<-temp[ temp-c(0,temp[-c(length(temp))]) >=buffer]
          #temp2Val<-tempVal[ temp-c(0,temp[-c(length(temp))]) >=buffer]
          
          temp<-releases[,1][!is.na(releases[,i])]
          tempVal<-releases[,i][!is.na(releases[,i])]
          if (length(temp)>0) {
            temp2<-append(c(),temp[1])
            temp2Val<-append(c(),tempVal[1])
            if ( length(temp) > 1){
              for( r in 2:length(temp) ){
                tdif<- temp[r] - temp2[length(temp2)]
                if ( tdif >= buffer ) { # difference big enough
                  temp2<-append(temp2,temp[r])
                  temp2Val<-append(temp2Val,tempVal[r])
                } else if ( sum ( !is.na(above[ 
                  grep(temp2[length(temp2)],above[,1])[1]:grep(temp[r],above[,1])[1]
                  ,i ])) == (tdif+1) ) {# is continuous
                  if ( above[grep(temp2[length(temp2)],above[,1])[1],i] 
                       < above[grep(temp[r],above[,1],above[,1])[1],i] ) { # further is biger
                    temp2[length(temp2)] <- temp[r]
                    temp2Val[length(temp2Val)] <- tempVal[r]
                  }                 
                }
              }
            }
          } else {
            temp2<-c()            
            temp2Val<-c()            
          }
                    
          temp3<-c()
          for (f in temp2 ){
            index<-grep(f,above[,1])[1]
            if (sum(!is.na(above[max(1,index-length+1):min(dim(above)[1],index+length-1),i])) >=length)
              temp3<-append(temp3,f)
          }
          
          release_list[[i-1]] <- append(release_list[[i-1]], temp2Val[temp2 %in% temp3] )
        }
    }    
  } else{ # get values from val
    for(i in 2:length(releases)){                   
      #temp<-append(c(0,0), releases[,1][!is.na(releases[,i])] ) 
      #tempVal<-append(c(0,0), val[,i][!is.na(releases[,i])] ) 
      #temp2<-temp[ temp-c(0,temp[-c(length(temp))]) >=buffer]
      #temp2Val<-tempVal[ temp-c(0,temp[-c(length(temp))]) >=buffer]
      
      temp<-releases[,1][!is.na(releases[,i])]
      tempVal<-val[,i][!is.na(releases[,i])]
      if (length(temp)>0) {
        temp2<-append(c(),temp[1])
        temp2Val<-append(c(),tempVal[1])
        if ( length(temp) > 1){
          for( r in 2:length(temp) ){
            tdif<- temp[r] - temp2[length(temp2)]
            if ( tdif >= buffer ) { # difference big enough
              temp2<-append(temp2,temp[r])
              temp2Val<-append(temp2Val,tempVal[r])
            } else if ( sum ( !is.na(above[ 
              grep(temp2[length(temp2)],above[,1])[1]:grep(temp[r],above[,1])[1]
              ,i ])) == (tdif+1) ) {# is continuous
              if ( above[grep(temp2[length(temp2)],above[,1])[1],i] 
                   < above[grep(temp[r],above[,1],above[,1])[1],i] ) { # further is biger
                temp2[length(temp2)] <- temp[r]
                temp2Val[length(temp2Val)] <- tempVal[r]
              }                 
            }
          }
        }
      } else {
        temp2<-c()            
        temp2Val<-c()            
      }
      
      temp3<-c()
      for (f in temp2 ){
        index<-grep(f,above[,1])[1]
        if (sum(!is.na(above[max(1,index-length+1):min(dim(above)[1],index+length-1),i])) >=length)
          temp3<-append(temp3,f)
      }
      release_list[[i-1]] <- append(release_list[[i-1]], temp2Val[temp2 %in% temp3] )
    }
  }
  return(release_list)
}

##########################################################################################
# determines points on boundary line
# NEEDS dplr data
# OPTIONAL prior,change,m1,m2,segment,notop
# RETURNS boundary points, and x and y
boundaryGet<-function(data,prior=NULL,change=NULL,m1=10,m2=10,segment=0.5,segment2=0.5,
                      notop=10,notop2=10,gfun=mean, bfun=mean){
  # bfun=mean; gfun=mean;notop=10
  if ( is.null(prior) )
    prior <- priorGrowth(data, m1, m2, bfun) #prior growth   
  if ( is.null(change) )
    change<-PGC(data, m1, m2, bfun) #percent growth change
  
  all.prior <- as.numeric(unlist(prior[,2:length(prior)]))
  all.change <- as.numeric(unlist(change[,2:length(prior)]))
  
  #prior growth segments
  temp <- as.data.frame(cbind(all.prior, all.change))
  
  #divided into segments 0-1cm by segment2, 1-max byt segment
  if (segment2 <= 1) {
    no1 <- round( 1/segment2,0)
  } else
    no1 <-0
  no2 <- round(( (max(all.prior, na.rm = TRUE))/segment), 0)
#  no<-no1+no2
  
  #the first value in each list is the upper boundary of each prior growth segment.
  segment_list <-as.list( 
    c( seq(segment2,segment2*no1,by=segment2) ,seq(1+segment,segment*no2,by=segment) )
    )
  
  no<-length(segment_list)
  
  #arrange percent growth change values into corresponding prior growth segments
  for(n in 1:no1) {    
    doon<- !is.na(all.prior)
    segment_list[[n]]<- append(segment_list[[n]], all.change[doon][(all.prior[doon] < segment_list[[n]][1]) & 
                        (all.prior[doon] >= segment_list[[n]][1] - segment2)])
  }
  for(n in (no1+1):no) {    
    doon<- !is.na(all.prior)
    segment_list[[n]]<- append(segment_list[[n]], all.change[doon][(all.prior[doon] < segment_list[[n]][1]) & 
                                                                     (all.prior[doon] >= segment_list[[n]][1] - segment)])
  }
  segment_mids <- c()
  segment_tops <- c()
  observations <- c()
  
  #the segments
  for(n in 1:no){
    if(length(segment_list[[n]]) > notop2){
      temp <- rev(sort(segment_list[[n]][-c(1)]))
      temp2 <- bfun(temp[1:notop],na.rm=T)      
    } else if ( length(segment_list[[n]])>1 ) {
      temp2 <- bfun(segment_list[[n]][-c(1)],na.rm=T)      
    } else
      temp2 <-0
    
    temp3 <- segment_list[[n]][1]
    if (n <= no1)
      temp3<-temp3-(segment2/2)
    else
      temp3<-temp3-(segment/2)
    temp4 <- length(segment_list[[n]])
    
    segment_mids <- append(segment_mids, temp3)
    segment_tops <- append(segment_tops, temp2)
    observations <- append(observations, temp4)
  }
  
  #observations <- as.vector(observations)
  segments <- c(as.vector(segment_mids))
  tops <- c(as.vector(segment_tops))
  
  boundaries <- cbind(segments, tops)
  
  x <- boundaries[1:dim(boundaries)[1],1] 
  y <- boundaries[1:dim(boundaries)[1],2]
  
  return( list( "bo"= as.data.frame(boundaries), "x" = all.prior, "y"=100*all.change) )
}


##########################################################################################
# returns boundary function or test given one
# NEEDS boundaries
# RETURNS boundary line as function, rsq and the best model
# morefun determines usage of more fitting function
boundaryFit<- function(boundaries,x,y,boundary=NULL,prefix="bo",store=TRUE,
                       storedev=pdf,initNLS=NULL){
  #x<-bo$x;y<-bo$y;boundaries<-bo$bo; store=FALSE
  
  doto<-which(boundaries$tops <0)[1]
  doto<-ifelse(is.na(doto),dim(boundaries)[1],doto)
  boundaries <- boundaries[1:(doto-1),]
  #boundaries$tops <- ifelse( boundaries$tops>0,boundaries$tops,0) 
  
  flm<-lm(tops~segments,data=boundaries)  
  print("--Summary of y=a+bx fit.")
  print(summary(flm))  
  fpoly<-lm(tops~segments+I(segments^2),data=boundaries)
  print("--Summary of y=a+bx+cx^2 fit.")
  print(summary(fpoly))  
  fr2<-c(summary(flm)$r.squared,summary(fpoly)$r.squared) 
  
  TSS<-sum((boundaries$tops-mean(boundaries$tops))^2)
  
  if (is.null(initNLS)){ 
    mya<-5; myb<--1
    } else { 
      mya<- initNLS$a; myb<- initNLS$b 
    }
  
  fr2inc<-length(fr2)
  tryCatch(
    {fexp0<-nls(tops~a*exp(b*segments),data=boundaries,start=list(a=mya,b=myb));
    fr2<-c(fr2,1-(deviance(fexp0)/TSS))}
  , error=function(e) print(paste("y=ae^bx","nls error:", e$message)))
  if (fr2inc == length(fr2) ) {fr2<-c(fr2,0);fexp0<-lm(boundaries$tops~1)}
  else{
    print("--Summary of y=ae^bx fit.")
    print(summary(fexp0))
  }
  
  if (is.null(initNLS)){ mya<-3; myb<--1; myc<-5} 
  else { mya<- initNLS$a; myb<- initNLS$b; myc<- initNLS$c;  }
  
  fr2inc<-length(fr2)
  tryCatch(
    {fexp1<-nls(tops~c+a*exp(b*segments),data=boundaries,start=list(a=mya,b=myb,c=myc));
    fr2<-c(fr2,1-(deviance(fexp1)/TSS)) }
  , error=function(e) print(paste("y=c+ae^bx","nls error:",e$message)))
  if (fr2inc == length(fr2) ) {fr2<-c(fr2,0);fexp1<-lm(boundaries$tops~1)}
  else{
    print("--Summary of y=c+ae^bx fit.")
    print(summary(fexp1))
  }
  
  if (is.null(initNLS)){ mya<-10; myb<--3; myc<-5; myd<--2} 
  else { mya<- initNLS$a; myb<- initNLS$b; myc<-initNLS$c;  myd<-initNLS$d; }
  
  fr2inc<-length(fr2)
  tryCatch(
    {fexp2<-nls(tops~c+d*segments+a*exp(b*segments),data=boundaries,start=list(a=mya,b=myb,c=myc,d=myd));
     fr2<-c(fr2,1-(deviance(fexp2)/TSS)) }
  , error=function(e) print(paste("y=c+dx+ae^bx","nls error:",e$message)))
  if (fr2inc == length(fr2) ) {fr2<-c(fr2,0);fexp2<-lm(boundaries$tops~1)}
  else{
    print("--Summary of y=c+dx+ae^bx fit.")
    print(summary(fexp2))
  }
  
  if (is.null(initNLS)){ mya<-2; myb<--1; myc<-20; myd<--3} 
  else { mya<- initNLS$a; myb<- initNLS$b; myc<-initNLS$c;  myd<-initNLS$d; }
  
  fr2inc<-length(fr2)
  tryCatch(
    {fexp3<-nls(tops~a*exp(b*segments)+c*exp(d*segments),data=boundaries,start=list(a=mya,b=myb,c=myc,d=myd))
    fr2<-c(fr2,1-(deviance(fexp3)/TSS)) }
  , error=function(e) print(paste("y=ae^bx+ce^dx","nls error:",e$message)))
  if (fr2inc == length(fr2) ) {fr2<-c(fr2,0);fexp3<-lm(boundaries$tops~1)}
  else{
    print("--Summary of y=ae^bx+ce^dx fit.")
    print(summary(fexp3))
  }
  
  if (is.null(initNLS)){ mya<-10; myb<--2} 
  else { mya<- initNLS$a; myb<- initNLS$b }
  
  fexp4<-lm(tops~segments+exp(segments)+segments:exp(segments),data=boundaries)
  fr2<-c(fr2,summary(fexp4)$r.squared)
  
  fr2inc<-length(fr2)
  tryCatch(
  { flog0<-nls(tops~a+b*log(segments),data=boundaries,start=list(a=mya,b=myb))
    fr2<-c(fr2,1-(deviance(flog0)/TSS)) }
  , error=function(e) print(paste("y=a+blog(x)","nls error:",e$message)))
  if (fr2inc == length(fr2) ) {fr2<-c(fr2,0);flog0<-lm(boundaries$tops~1)}
  else{
    print("--Summary of y=a+blog(x) fit.")
    print(summary(flog0))
  }
  
  flog1<-lm(tops~segments+log(segments)+segments:log(segments),data=boundaries)
  print("--Summary of y=a+bx+clog(x)+dxlog(x) fit.")
  print(summary(flog1))  
  fr2<-c(fr2,summary(flog1)$r.squared)
  
  # which R2 is max
  imax<-which.max(fr2)
  
  print(paste("The fitted boundary line summary!"))
  bbestmodel<-c()
  
  if (imax == 1) { #lm
    bline<-function(x) { pmax( as.numeric(flm$coefficients[1]) + 
                           as.numeric(flm$coefficients[2])*x,0) }
    print(paste("Linear model y=a+bx was the best!"))
    print(summary(flm))
    bbestmodel<-flm
  } else if (imax == 2) { #poly
    bline<-function(x) { pmax(as.numeric(fpoly$coefficients[1]) + 
                                as.numeric(fpoly$coefficients[2])*x +
                                as.numeric(fpoly$coefficients[3])*x*x,0) }
    print(paste("Polynomial model y=a+bx+cx^2 was the best!"))
    print(summary(fpoly))
    bbestmodel<-fpoly
  } else if (imax == 3) { #exp0  
    mys<-summary(fexp0)
    bline<-function(x) { pmax(mys$coefficients[1,1]*exp(mys$coefficients[2,1]*x),0) } 
    print(paste("Exponential model y=ae^bx was the best!"))
    print(mys)  
    bbestmodel<-fexp0
  } else if (imax == 4) { #exp1
    mys<-summary(fexp1)  
    bline<-function(x) { pmax(mys$coefficients[1,1]*exp(mys$coefficients[2,1]*x)
                              +mys$coefficients[3,1],0) }     
    print(paste("Exponential model y=c+ae^bx was the best!"))
    print(mys)
    bbestmodel<-fexp1
  } else if (imax == 5) {  #exp2
    mys<-summary(fexp2)  
    bline<-function(x) { pmax(mys$coefficients[1,1]*exp(mys$coefficients[2,1]*x)+
                                mys$coefficients[3,1]+ mys$coefficients[4,1]*x,0) }     
    print(paste("Exponential model y=c+dx+ae^bx was the best!"))
    print(mys)
    bbestmodel<-fexp2
  } else if (imax == 6) { # exp3
    mys<-summary(fexp3)  
    bline<-function(x) { pmax(mys$coefficients[1,1]*exp(mys$coefficients[2,1]*x)+
                                mys$coefficients[3,1]*exp(mys$coefficients[4,1]*x),0) }     
    print(paste("Exponential model y=ae^bx+ce^dx was the best!"))
    print(mys)
    bbestmodel<-fexp3
  } else if (imax == 7) { # exp4
    bline<-function(x) { pmax(as.numeric(fexp4$coefficients[1]) + 
                                as.numeric(fexp4$coefficients[2])*x +
                                as.numeric(fexp4$coefficients[3])*log(x) +
                                as.numeric(fexp4$coefficients[4])*log(x)*x,0) }
    print(paste("Logarithmic model y=a+bx+ce^x+dxe^x was the best!"))
    print(summary(fexp4))
    bbestmodel<-fexp4
  } else if (imax == 8) {# log0
    mys<-summary(flog0)  
    bline<-function(x) { pmax(mys$coefficients[1,1]+mys$coefficients[2,1]*log(x),0) }     
    print(paste("Logaritmic model y=c+ae^bx was the best!"))
    print(mys)
    bbestmodel<-flog0
  } else if (imax ==9) { #log1
    bline<-function(x) { pmax(as.numeric(flog1$coefficients[1]) + 
                                as.numeric(flog1$coefficients[2])*x +
                                as.numeric(flog1$coefficients[3])*log(x) +
                                as.numeric(flog1$coefficients[4])*log(x)*x,0) }
    print(paste("Logarithmic model y=a+bx+clog(x)+dxlog(x) was the best!"))
    print(summary(flog1))
    bbestmodel<-flog1
  }
  
  
  if(store){
    write.table(boundaries, "boundaries.csv", sep="\t",row.names=F)
    mycol<-rainbow(length(fr2))
    
    storedev(paste(prefix,"_boundaries.",deparse(substitute(storedev)),sep=""))
    plot(x, y, xlab = "prior growth [mm]", ylab = "growth change [%]", pch = 21, col = "black")
    points(boundaries$segments, boundaries$tops*100,col="darkblue",pch=15)
    abline(h=0,lty="dashed")
    lines(boundaries$segments,predict(flm)*100,col=mycol[1])
    lines(boundaries$segments,predict(fpoly)*100,col=mycol[2])    
    lines(boundaries$segments,predict(fexp0)*100,col=mycol[3])
    lines(boundaries$segments,predict(fexp1)*100,col=mycol[4])
    lines(boundaries$segments,predict(fexp2)*100,col=mycol[5])
    lines(boundaries$segments,predict(fexp3)*100,col=mycol[6])
    lines(boundaries$segments,predict(fexp4)*100,col=mycol[7])
    lines(boundaries$segments,predict(flog0)*100,col=mycol[8])
    lines(boundaries$segments,predict(flog1)*100,col=mycol[9])
    
    if ( ! is.null(boundary) ) {
      #lines(boundaries$segments,predict(fbou)*100,col="black")
      lines(boundaries$segments,boundary(boundaries$tops)*100,col="black")
      legend("topright",legend=paste(c("y=a+bx R2=","y=a+bx+cx^2 R2=","y=ae^bx R2=",
                                       "y=c+ae^bx R2=","y=c+dx+ae^bx R2=",
                                       "y=ae^bx+ce^dx R2=","y=a+bx+ce^x+dxe^x R2=",
                                       "y=a+blog(x) R2=","y=a+bx+clog(x)+dxlog(x) R2=","given"),
                                     c(round(fr2,3),"")),
             lty=rep(1,length(mycol)+1),col=c(mycol,"black"))
    } else
      legend("topright",legend=paste(c("y=a+bx R2=","y=a+bx+cx^2 R2=","y=ae^bx R2=",
                                       "y=c+ae^bx R2=","y=c+dx+ae^bx R2=",
                                       "y=ae^bx+ce^dx R2=","y=a+bx+ce^x+dxe^x R2=",
                                       "y=a+blog(x) R2=","y=a+bx+clog(x)+dxlog(x) R2="),
                                     round(fr2,3)),lty=rep(1,length(mycol)),col=mycol)
    dev.off()
  } 
  return(list("fun" = bline, "rsq" = fr2[imax], "bestModel" = bbestmodel))
}


##########################################################################################
writeReleaseStats<-function(release_list,mytext){
      
  inyears<-unlist(release_list) 
  inyears<-inyears[inyears!=0]
  release_list_count <- length(inyears)
  print(paste(mytext,release_list_count))  
  print(table(inyears))
  return(inyears)
}


##########################################################################################
# Tranform list of releases in data frame
relListToDataFrame<-function(release_list,data){
  relDF<-data
  relDF[!is.na(relDF)]<-0
  for(i in 1:length(release_list) ) { # i<-1
    tlist<-release_list[[i]]
    tlist<-tlist[-c(1)]
    relDF[,i]<-ifelse( (rownames(relDF) %in% tlist) & (relDF[,i] == 0 ), 1, NA )
  }  
  return(relDF)
}



##########################################################################################

removeMajorFromModerate<-function(mod,maj,zero=0,on=NULL){
  for (r in 1:length(mod)) {
    rll<-length(mod[[r]])    
    if ( rll >1 ){
      rll2<-length(maj[[r]])    
      if ( rll2 >1 ){
        if ( is.null(on) ) {
          mod[[r]] <- c(zero,
                      mod[[r]][2:rll][! mod[[r]][2:rll] %in% maj[[r]][2:rll2]  ])
        } else {
          on[[r]] <- c(zero,
                       on[[r]][2:rll][! mod[[r]][2:rll] %in% maj[[r]][2:rll2]  ])
        }
      }
    }    
  }
  if ( is.null(on) ) {
    return(mod)
  } else {
    return(on)
  }
}
##########################################################################################
##########################################################################################
##########################################################################################
# MAIN PLOT FUNCTIONS

# Plot data and releses
# NEEDS dplr data, abs, release data
# RETURNS plot 
# method = c("FraverWhite", "NowackiAbrams", "BlackAbrams","Splechtna")
plotRelease<-function(data, abs, rel, treeno=1, method = "FraverWhite", 
                         type="l", xlab=NULL, ylab = NULL, main=NULL, col = c("black","lightblue"),
                         addHLinesCol = c("olivedrab","red","darkblue"), addHLines = c(NULL, NULL, NULL),
                         addHLinesText = c("","",""), smallcex= 0.85, plotfirst=TRUE, plotpoints=FALSE,
                         ...) {
  if (method == "FraverWhite" ) {
    if ( is.null(ylab) )
      ylab <- "absolute increase [mm]"
    if ( is.null(main) )
      main <- paste("Fraver & White (2005)", sep = "")
    mytimes <- 1
  } else {    
    mytimes <- 100
  }
  
  if (method == "NowackiAbrams" ) {
    if ( is.null(main) )
      main <- paste("Nowacki & Abrams (1997)", sep = "") 
    if ( is.null(ylab) )
      ylab <- "growth change [%]"
  }
  if (method == "BlackAbrams" ) {
    if ( is.null(main) )
      main <- paste("Black & Abrams (2003)", sep = "")
    if ( is.null(ylab) )
      ylab <- "boundary line [%]"
  }
  if (method == "Splechtna" ) {
    if ( is.null(main) )
      main <- paste("Splechtna et al. (2005)", sep = "")
    if ( is.null(ylab) )
      ylab <- "boundary line [%]"
  }
  
  if ( is.null(xlab) )
    xlab <- "years"
  
  par(mar = c(5.1,4.1,4.1,5.1))
  # plot data
  cls<- ! is.na(data[,treeno ])
  plot ( rownames( data[cls,] ), data[cls,treeno ], 
         type = type, xlab=xlab, ylab = "ring widths [mm]", main = main, col = col[1], frame= F, ... )
  
  if (plotfirst)
    mtext(paste("first year ", min(rownames(data[!is.na(data[,treeno]) ,]),na.rm=T), " ", sep = ""), 
        side = 3, line = -1, cex = smallcex, adj=1)
  par(new=TRUE)
  
  cls2<- ! is.na(abs[,treeno+1])
  if ( sum(cls2,na.rm=T) >0 ) {
    #second y range
    myrange<-c(mytimes*range(abs[,treeno+1],na.rm=T))
    myrange[2] <- max(myrange[2],(addHLines[1]*mytimes)+5,na.rm=T)
    #empty plot
    plot ( rownames( data[cls,] ), data[cls,treeno ],type = "n", xlab="", ylab = "",
         main = "",xaxt="n",yaxt="n",frame=F,ylim=c(myrange) )
    # plot releases
    lines(abs[cls2,1], mytimes*abs[cls2,treeno+1], type=type, col=col[2], lwd=2, ...)
  
    axis(4)
    mtext(ylab,side=4,line=3, ...)
    # plot horizontal lines  
    for (i in 1:length(addHLines) ) {
      abline(h = mytimes*addHLines[i], col = addHLinesCol[i],lty = "dotted" )    
    }
    if ( length(addHLines) >0 )
      text(x=min( abs[cls2,1],na.rm=T ),y=mytimes*(addHLines)+((myrange[2]-myrange[1])/20),
       labels=round(addHLines*mytimes,2),col=addHLinesCol,cex = smallcex)
  
    # plot points
    if (plotpoints) {
      if ( (method == "BlackAbrams") | (method == "Splechtna") )
        points(rel$years[[treeno]][2:length(rel$pgc[[treeno]])], 
           mytimes*rel$pgc[[treeno]][2:length(rel$pgc[[treeno]])], pch = 19, cex = 0.7)
      else
        points(rel$releases[,1], mytimes*rel$releases[,treeno+1], pch = 19, cex = 0.7)
    }  
  
    # plot vertical lines
    for(t in 2:length(rel$years[[treeno]]) ){
      abline(v = rel$years[[treeno]][t], lty = "dotted")
      text(rel$years[[treeno]][t]-5, -0.3, rel$years[[treeno]][t], srt = 90, cex=smallcex)
    }
  } else
    print( paste("No release data for",names(data)[treeno]) )
}

# returns boundary function
# NEEDS boundaries
# RETURNS plot the best of given boundary line
plotBoundary<- function(boundaries,x,y,boundary,rsq=NULL,prefix="bo",criteria=0.2,criteria2=0.5,
                        store=TRUE,storedev=jpeg){
  #x<-bo$x;y<-bo$y; boundaries<-bo$bo;boundary=bo2$fun;criteria=0.2;criteria2=0.5;store=TRUE
  if (deparse(substitute(storedev)) =="storedev")
    storedev<-jpeg
  if (store)
    pdf(paste(prefix,"_boundary.pdf",sep=""))
  
  plot(x, y, xlab = "prior growth [mm]", ylab = "growth change [%]", type='n')
  #points(boundaries$segments, boundaries$tops*100,col="darkblue",pch=15)
  all<-boundary(x)*100
  
  # draw up to negative or the last value
  doto<-which(boundaries$tops <0)[1]
  doto<-ifelse(is.na(doto),dim(boundaries)[1],doto)
  
  tochoose<-(all*criteria > y) | (x>boundaries$segments[doto] )
  points(x[tochoose],y[tochoose],pch = 4, col = "gray75",)
  tochoose<-all <= y & (x<=boundaries$segments[doto] )
  points(x[tochoose],y[tochoose],pch = 20, col = "darkblue",)
  tochoose<-(all*criteria2 <= y) & (all > y) & (x<=boundaries$segments[doto])
  points(x[tochoose],y[tochoose],pch = 20, col = "darkblue")
  tochoose<-(all*criteria <= y) & (all*criteria2 > y)
  points(x[tochoose],y[tochoose],pch = 4, col = "black",)
  
  abline(h=0,lty="dashed")
  
  toline<-seq(from=0.01,to=boundaries$segments[doto],length=50)
  lines(toline,boundary(toline)*100,col="red")
  lines(toline,boundary(toline)*(criteria2*100),col="blue")
  lines(toline,boundary(toline)*(criteria*100),col="green")
  legend("topright",legend=c("boundary line",paste(criteria2,"boundary line"),
                             paste(criteria,"boundary line")),
         col=c("red","blue","green"),lty=c(1,1,1))
  if ( ! is.null(rsq) )
    mtext(paste("R2=",round(rsq,3)," ") ,side = 3, line = -5, adj=1)
  
  if (store)
    dev.off()
}
  
# returns plot of releases
# NEEDS data and releases
# RETURNS plot total releases and number of trees
plotNORelease<-function(data,inyears,in2years=NULL,criteria,criteria2=NULL,prefix="rel",
                        store=TRUE,storedev=jpeg){
  #data<-mdata;inyears<-rs;in2years<-rs2;criteria=0.2;criteria2<-0.5;store=F
  notrees<-rowSums(data>0,na.rm=T)
  
  a<-table(inyears) 
  if (! is.null(in2years) ) {
    b<-table(in2years)
    ab<-table( c(inyears, in2years) )
    btot<-notrees[names(notrees) %in% names(b)]
    bdif<-round((b*100)/btot,2)    
  } else {
    ab<-a
  }
  
  myxlim<-as.numeric( c( names(ab)[1], names(ab)[length(ab)] ))
  
  abtot<-notrees[names(notrees) %in% names(ab)]
  abtotADJ<-ifelse(abtot>ab,abtot,ab)
  abdif<-round((ab*100)/abtotADJ,2)
  # change solve problem with >100% abdif<-round((ab*100)/abtot,2)
  
  if (store)
    storedev(paste(prefix,"_inyears.",deparse(substitute(storedev)),sep=""))
  
  par(mar = c(5.1,4.1,4.1,5.1))
  plot(c(min(abdif),max(abdif))~c(myxlim[1],myxlim[2]),xlim=myxlim, xlab="years",ylab="% of trees with release",type='n')
  points(abdif,xlim=myxlim, xlab="",ylab="",col="blue")
  if (!  is.null(in2years) )
    points(bdif,xlim=myxlim, xlab="years",ylab="number of releases",col="red")
  par(new=TRUE)  
  plot( c(min(notrees),max(notrees))~c(myxlim[1],myxlim[2]),xlim=myxlim, xlab="",ylab="",
        xaxt='n',yaxt='n',type='n')
  lines(notrees~(rownames(data)),type="l",ylab="",xaxt='n',yaxt='n', 
     xlab="",xlim=myxlim)
  axis(4)
  mtext("number of trees",side=4,line=3)
  if (!  is.null(in2years) )
    legend("topleft",c(paste(">=",criteria," <",criteria2,sep=""),paste(">=",criteria2,sep=""),
                       "number of trees"), col=c("blue","red","black"),lty=c(1,1,1),
           lwd=c(2,2,2))
  else
    mtext(paste(" ",">",criteria," ",sep=""), side = 3, line = -1, adj=0)
  
  if (store)
    dev.off()
  
  if (!  is.null(in2years) ) {
    btotret<-ab
    btotret[! names(ab) %in% names(b)]<-0
    btotret[ names(ab) %in% names(b)]<-btotret[ names(ab) %in% names(b)]-b    
    bdifret<-round((btotret*100)/abtot,2)
    myret<-data.frame("AllReleasesFreq"=abdif, "NumberOfAllReleses"=as.numeric(ab),
                      "ModerateReleasesFreq"=as.numeric(bdifret),
                      "NumberOfModeraeReleases"=as.numeric(btotret), "NumberOfAllTrees"=abtot)
  } else
    myret<-data.frame("AllReleasesFreq"=abdif, "NumberOfAllReleses"=as.numeric(ab),"NumberOfAllTrees"=abtot)
  names(myret)[1:2] <-c("AllReleasesYear","AllReleasesFreq")
  return(myret)
}

# returns plot of growth of inidvidual tree
# NEEDS data 
# RETURNS plot growth and polynom
plotGrowth<-function(data=NULL,prefix="growth", polynom=4, store=TRUE,storedev=jpeg, ...){
  if(is.null(data)){
    return(paste("Data must be specified as argument!"))
  }
  
  for(i in 1:length(data)) {
    if (store)
      storedev(paste(prefix,"_tree",names(data)[i],".",deparse(substitute(storedev)),
                     sep=""),pointsize = 10)
    
    treena<-is.na(data[,i])
    mtree<-data[!treena,i]
    plot(mtree~rownames(data)[!treena],xlab="years",ylab="rings width",type="l",col="gray", ... )
    mtext(paste(names(data)[i]," ",sep=""), side = 3, line = -1, adj=1)
    fpoly<-lm(mtree~poly(as.numeric(rownames(data)[!treena]),polynom))  
    lines(as.numeric(rownames(data)[!treena]),predict(fpoly),lwd=2 )
    
    if (store)
      dev.off()
  }  
}

##########################################################################################
# Plot first years of trees
# NEEDS dplr data, misspith
# RETURNS list of first year for each tree
plotFirstYears <- function(data=NULL, misspith=NULL,store=TRUE,storedev=jpeg,
                           prefix="fy",...){
  # data<-mdata;store=FALSE
  if(is.null(data)){
    return(paste("Data must be specified as argument!"))
  }
  wdata<-data    
  if (! is.null(misspith) ) {    
    if ( dim(data)[2] == dim(misspith)[1]) {
      ymin<-as.numeric(rownames(wdata)[1])
      for(i in 1:dim(misspith)[1]) {
        if ( misspith[i,2]>0 ){
          imin<-rownames(wdata)[!is.na(wdata[, grep(misspith[i,1],colnames(wdata))[1]  ])][1]
          iwhich<-grep(imin,rownames(wdata))[1]
          if ( as.numeric(imin) - misspith[i,2] >= ymin ){
            wdata[ (iwhich - misspith[i,2]) :iwhich, grep(misspith[i,1],colnames(wdata))[1]]<-1
          } else
            wdata[ grep(ymin,rownames(wdata))[1] :iwhich, grep(misspith[i,1],colnames(wdata))[1]]<-1
        }
      }
    }else 
      print(paste("Data and misspith must have same length!"))
  }
  
  firsty<-data.frame(tree=names(wdata),firstyear=NA)
  for (i in 1:length(firsty$tree)) {
    firsty[i,2] <- (rownames(wdata) [ !is.na(wdata[,i])]) [1]
  }
  write.table(firsty,paste(prefix,"_firstyears.csv",sep=""),row.names=F)
  
  if(store) {
    storedev(paste(prefix,"firstyears.",deparse(substitute(storedev)),sep=""),pointsize = 10,...)
    write.table(data.frame(sum=rowSums(wdata>0,na.rm=T),year=rownames(wdata)),"firstyear.csv",
                row.names=F)
  }
  
  plot(rowSums(wdata>0,na.rm=T)~rownames(wdata),xlab="years",ylab="sample depth",type='l')
  if ( ! is.null(misspith)) {
    lines(rownames(data),rowSums(data>0,na.rm=T),lty=2)
    legend("topleft",c("with missing rings","based on first measured year"),lty=c(1,2))
  }
  
  if(store)
    dev.off()  
}