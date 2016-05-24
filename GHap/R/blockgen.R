#Function: ghap.blockgen
#License: GPLv3 or later
#Modification date: 2 Feb 2016
#Written by: Yuri Tani Utsunomiya
#Contact: ytutsunomiya@gmail.com
#Description: Generate blocks based on sliding windows

ghap.blockgen<-function(
  phase,
  windowsize=10,
  slide=5,
  unit="marker"
){
  
  #Check if phase is a GHap.phase object
  if(class(phase) != "GHap.phase"){
    stop("Argument phase must be a GHap.phase object.")
  }
  
  #Initialize vectors
  BLOCK <- rep(NA,times=phase$nmarkers.in)
  BP1 <- rep(NA,times=phase$nmarkers.in)
  BP2 <- rep(NA,times=phase$nmarkers.in)
  SIZE <- rep(NA,times=phase$nmarkers.in)
  CHR <- rep(NA,times=phase$nmarkers.in)
  NSNP <- rep(NA,times=phase$nmarkers.in)
  bp<-phase$bp[which(phase$marker.in)]
  
  if(unit == "kb"){
    
    #Transform windowsize to bp
    windowsize<-windowsize*1e+3
    slide=slide*1e+3
    
    id1<-seq(1,max(bp),by=slide);
    id2<-id1+windowsize;
    id1<-id1[id2<=max(bp)]
    id2<-id2[id2<=max(bp)]
    for(i in 1:length(id1)){
      slice <- which(bp >= id1[i] & bp <= id2[i])
      BLOCK[i]<-paste("CHR",phase$chr,"_B",i,sep="")
      BP1[i] <- id1[i]
      BP2[i] <- id2[i]
      SIZE[i] <- id2[i]-id1[i]
      CHR[i] <- phase$chr
      NSNP[i] <- length(slice)
    }
    
  }else if(unit == "marker"){
    
    id1<-seq(1,phase$nmarkers.in,by=slide);
    id2<-id1+(windowsize-1);
    id1<-id1[id2<=phase$nmarkers.in]
    id2<-id2[id2<=phase$nmarkers.in]
    for(i in 1:length(id1)){
      slice <- seq(id1[i],id2[i],by=1)
      BLOCK[i]<-paste("CHR",phase$chr,"_B",i,sep="")
      BP1[i] <- bp[id1[i]]
      BP2[i] <- bp[id2[i]]
      SIZE[i] <- bp[id2[i]]-bp[id1[i]]
      CHR[i] <- phase$chr
      NSNP[i] <- length(slice)
    }
    
  }
  
  BLOCK<-na.omit(BLOCK)
  CHR<-na.omit(CHR)
  BP1<-na.omit(BP1)
  BP2<-na.omit(BP2)
  SIZE<-na.omit(SIZE)
  NSNP<-na.omit(NSNP)
  
  
  return(data.frame(BLOCK,CHR,BP1,BP2,SIZE,NSNP,stringsAsFactors = FALSE));
  
}