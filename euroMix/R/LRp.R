#Computes p-value. Wrapper for  the function pvalue.machine. Takes as input
#matrices with sample data, victim profile and suspect profile
#SO FAR ONLY WORKS FOR ONE SUSPECT (AND ONE VICTIM)!!
#forensim must be added as required in euroMix description!!
LRp <- function(sampleData, victimData, suspectData, db, hp, hd, prD, prC ) {
  
  if(!inherits(suspectData,'list')) suspectData <- list(suspectData)
  
  Nmarkers <- nrow(sampleData)
  markers <- rownames(sampleData)
  
  Nsuspects <- length(suspectData) 
  
  shp <- grepl("S",hp,ignore.case=TRUE)
  vhp <- grepl("V",hp,ignore.case=TRUE)
  uhp <- grep("U",hp,ignore.case=TRUE)
  
  shd <- grepl("S",hd,ignore.case=TRUE)
  vhd <- grepl("V",hd,ignore.case=TRUE)
  uhd <- grep("U",hd,ignore.case=TRUE)
  
  isp <- ivp <- isd <- ivd <- numeric()
  isp2 <- isd2 <- 1:(Nsuspects*2); ivp2 <- ivd2 <- 1:2
  if(any(shp)) {isp <- 1:(Nsuspects*2); isp2 <- numeric()}
  if(any(vhp)) {ivp <- 1:2; ivp2 <- numeric()}
  if(any(shd)) {isd <- 1:(Nsuspects*2); isd2 <- numeric()}
  if(any(vhd)) {ivd <- 1:2; ivd2 <- numeric()}
  
  lr <- numeric()
  A <- numeric()
  for(i in 1:Nmarkers){
    r <- as.numeric(sampleData[i,which(!is.na(sampleData[i,]))])
    v <- as.numeric(victimData[i,])
    s <- sapply(suspectData, function(x) as.numeric(x[i,]))
    afreq <- db[db[,1]==markers[i],3]
    A[i] <- choose(length(afreq)+1,2)
    names(afreq) <- db[db[,1]==markers[i],2]
    
    pD <- rep(prD,(length(shp)+length(vhp)+length(uhp))*2)
    
    lrp <- likEvid(r, T=c(s[isp],v[ivp]), V=c(s[isp2],v[ivp2]), x=length(uhp), theta=0, 
                   prDHet=pD, prDHom=pD^2, prC=prC, freq=afreq)
    lrd <- likEvid(r, T=c(s[isd],v[ivd]), V=c(s[isd2],v[ivd2]), x=length(uhd), theta=0, 
                   prDHet=pD, prDHom=pD^2, prC=prC, freq=afreq)
    
    lr[i] <- lrp/lrd
  }
  LR.suspect <- prod(lr)
  
  
  #LR when replacing suspect with all possible genotypes (i.e. all random men)
  LR.table <- P.table <- matrix(0,Nmarkers,max(A))
  for(i in 1:Nmarkers){
    afreq <- db[db$Marker==markers[i],]$Frequency
    names(afreq) <- db[db$Marker==markers[i],]$Allele
    #All possible genotypes for the given marker
    gtAll <- as.matrix(expand.grid(names(afreq),names(afreq)))
    gtAll <- gtAll[gtAll[,1]<=gtAll[,2],]
    #likelihood ratios and corresponding genotype probabilities
    lrAll <- numeric()
    pAll <- numeric()
    for(j in 1:nrow(gtAll)) {
      rm <- gtAll[j,]
      r <- as.numeric(sampleData[i,which(!is.na(sampleData[i,]))])
      v <- as.numeric(victimData[i,])
      lrp <- likEvid(r, T=c(rm[isp],v[ivp]), V=c(rm[isp2],v[ivp2]), x=length(uhp), theta=0, 
                     prDHet=pD, prDHom=pD^2, prC=prC, freq=afreq)
      lrd <- likEvid(r, T=c(rm[isd],v[ivd]), V=c(rm[isd2],v[ivd2]), x=length(uhd), theta=0, 
                     prDHet=pD, prDHom=pD^2, prC=prC, freq=afreq)
      lrAll[j] <- lrp/lrd
      pAll[j] <- ifelse(rm[1]==rm[2], prod(afreq[as.character(rm)]), 2*prod(afreq[as.character(rm)]) )
    }
    LR.table[i,1:length(lrAll)] <- lrAll
    P.table[i,1:length(lrAll)] <- pAll
  }
  #The p-value
  pvalue <- pvalue.machine( LR.suspect, LR.table, P.table )
  
  return(list(LR=LR.suspect,pvalue=pvalue))
}
