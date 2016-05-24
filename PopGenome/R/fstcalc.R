fstcalc <- function(matrix_pol,populations,outgroup=FALSE,data){



### ACHTUNG !!! Hier kann man zwar eine Outgroup definieren, diese wird jedoch nur als weitere Population betrachtet !



# If there is only one population defined --------------------------
if(length(populations)<1){#&& outgroup[1]==FALSE){
FSTALL   <- NaN
PIA      <- as.matrix(NaN)
PIA2     <- as.matrix(NaN)
PIW      <- as.matrix(NaN)
FSTHKY   <- as.matrix(NaN)
FSTHKY2  <- as.matrix(NaN)
FSTPAIR  <- as.matrix(NaN)
FSTPAIR2 <- as.matrix(NaN)
FST1ALL  <- as.vector(NaN)
SV1      <- as.matrix(NaN)
SV2      <- as.matrix(NaN)

return(list(PIA=PIA,PIA2=PIA2,PIW=PIW,FSTALL=FSTALL,FSTHKY=FSTHKY,FSTPAIR=FSTPAIR,FST1ALL=FST1ALL,SV1=SV1,SV2=SV2,FSTHKY2=FSTHKY2,FSTPAIR2=FSTPAIR2))
}
# -------------------------------------------------------------------
### include outgroup  <-------------------------------------- OUTGROUP is INCLUDED !!!!!!!!!!!!!!!!
 
 npops <- length(populations)

 if(outgroup[1]!=FALSE){
 populationsoutgroup    <- c(populations,list(outgroup))
 npopsoutgroup          <- length(populations)
 }

###############################################################
# data is a list returned from get_data_test 
# stat_length:(NUMBER OF NUCLEOTIDES)
###############################################################

if(length(data$n.valid.sites)!=0){
valid <- TRUE
stat_length <- numeric(npops)
for (x in 1:npops){
    stat_length[x] <- sum(data$n.nucleotides[populations[[x]]])
}

# stat_total_length: NUMBER OF ALL POSITIONS #################
 stat_total_length <- data$n.valid.sites 
##############################################################
}

if(length(data$n.valid.sites)==0){
  valid <- FALSE
  samplesize <- numeric(npops)
  for(xx in 1:npops){
     samplesize[xx] <- length(populations[[xx]])
  }  

}
#####################################


###


flag <- 0                 # <-------------------------------------------- HKY CORRECTION

###########################
# initialize
############################

matrix_sv   <- data$transitions
SV1         <- matrix(0,npops,npops)
SV1[upper.tri(SV1,TRUE)]<- NA
SV2        <- SV1
PIA        <- SV1
apply_id   <- t(which(SV1==0,TRUE)) 

#print(apply_id)

# initialize
#########################

diagonal1 <- numeric(npops)
diagonal2 <- numeric(npops)
PIW       <- numeric(npops)
FREQ_LIST <- vector("list",npops)

#--------------------------------------------------------------

for(x in 1:npops){

if(length(populations[[x]])==0){next;}

erg           <- get_freq(populations[[x]],matrix_pol,matrix_sv) # Function
diagonal1[x]  <- erg$SV0
diagonal2[x]  <- erg$SV1
PIW[x]        <- erg$PIW

# If there is only one population defined -------------------------- berechne nur PIW
if(length(populations)==1){#&& outgroup[1]==FALSE){
FSTALL   <- NaN
PIA      <- as.matrix(NaN)
PIA2     <- as.matrix(NaN)
PIW           <- t(as.matrix(PIW))
colnames(PIW) <- paste("pop",1:npops)
FSTHKY   <- as.matrix(NaN)
FSTHKY2  <- as.matrix(NaN)
FSTPAIR  <- as.matrix(NaN)
FSTPAIR2 <- as.matrix(NaN)
FST1ALL  <- as.vector(NaN)
SV1      <- as.matrix(NaN)
SV2      <- as.matrix(NaN)

return(list(PIA=PIA,PIA2=PIA2,PIW=PIW,FSTALL=FSTALL,FSTHKY=FSTHKY,FSTPAIR=FSTPAIR,FST1ALL=FST1ALL,SV1=SV1,SV2=SV2,FSTHKY2=FSTHKY2,FSTPAIR2=FSTPAIR2))
}
# -------------------------------------------------------------------


FREQ_LIST[[x]]<- erg$FREQ

}

diag(SV1) <- diagonal1
diag(SV2) <- diagonal2

### PIW
### FREQ_LIST
### SV1
### SV2


TEST <- numeric(npops)

#APPLY ! ####  PIA, PIA1ALL ###################
#---------------------------------------------#

ap_erg <- apply(apply_id,2,function(pop){ # PERMUTATION OVER ALL POPULATION PAIRS

idd <- which((FREQ_LIST[[pop[1]]]["freq0",]*FREQ_LIST[[pop[2]]]["freq0",])>0)
pia <- (FREQ_LIST[[pop[1]]]["freq1",idd]*FREQ_LIST[[pop[2]]]["freq2",idd] + FREQ_LIST[[pop[1]]]["freq2",idd]*FREQ_LIST[[pop[2]]]["freq1",idd])/
       (FREQ_LIST[[pop[1]]]["freq0",idd]*FREQ_LIST[[pop[2]]]["freq0",idd])

PIA[pop[1],pop[2]] <<- sum(pia) # every pair of population !!!!

if(pop[2] < npops && length(populations[[pop[1]]])>=1 && length(populations[[pop[2]]])>=1){ 
TEST[pop[1]] <<- TEST[pop[1]] + sum(pia) # PIA1ALL
TEST[pop[2]] <<- TEST[pop[2]] + sum(pia) # PIA1ALL
}

id1 <-   which(matrix_sv[idd]==1) 
id2 <-   which(matrix_sv[idd]==2)


sv1 <-  (FREQ_LIST[[pop[1]]]["freq1",id1]*FREQ_LIST[[pop[2]]]["freq2",id1] + FREQ_LIST[[pop[1]]]["freq2",id1]*FREQ_LIST[[pop[2]]]["freq1",id1])/
                       (FREQ_LIST[[pop[1]]]["freq0",id1]*FREQ_LIST[[pop[2]]]["freq0",id1])
   
sv2 <- (FREQ_LIST[[pop[1]]]["freq1",id2]*FREQ_LIST[[pop[2]]]["freq2",id2] + FREQ_LIST[[pop[1]]]["freq2",id2]*FREQ_LIST[[pop[2]]]["freq1",id2])/
                      (FREQ_LIST[[pop[1]]]["freq0",id2]*FREQ_LIST[[pop[2]]]["freq0",id2])


SV1[pop[1],pop[2]]  <<- sum(sv1) # GLOBAL damit beim AUfruf das Veraenderte SV uebernommen wird
SV2[pop[1],pop[2]]  <<- sum(sv2) # GLOBAL


return()

})#End of Apply!
######################################

pop           <- paste("pop",1:npops)
colnames(SV1) <- pop
rownames(SV1) <- pop
colnames(SV2) <- pop
rownames(SV2) <- pop
colnames(PIA) <- pop
rownames(PIA) <- pop
K       <- PIA[npops,1:(npops-1)]  # the populations that are tested against the outgroup
PIA1ALL <- TEST

# Missing values statistics.K und pia1all

if(flag==1){  # FLAG #

   #TCGA for all populations 
   if(is.matrix(data$biallelic.compositions[populations[[1]],])){TCGA1 <- colSums(data$biallelic.compositions[populations[[1]],])}
   if(is.vector(data$biallelic.compositions[populations[[1]],])){TCGA1 <- sum(data$biallelic.compositions[populations[[1]],])}


    for(x in 2:npops){
      if(length(populations[[x]])==1){
         TCGA <- data$biallelic.compositions[populations[[x]],]
      }else{
         TCGA  <- colSums(data$biallelic.compositions[populations[[x]],])
      }   
         
         TCGA1 <- rbind(TCGA1,TCGA) 
    }

  TCGA <- TCGA1
  rownames(TCGA) <- pop


#### HKY #################################### Correction ------
HKY <- matrix(0,npops,npops)
# HKY for each population
for(x in 1:npops){ 
if(valid){n1 <- ceiling(stat_length[x]/stat_total_length)}else{n1 <- samplesize[x]}
gT <- TCGA[x,"tcga1"]
gC <- TCGA[x,"tcga2"]
gG <- TCGA[x,"tcga3"]
gA <- TCGA[x,"tcga4"]
P1 <- SV1[x,x] * (gA*gG/(gA*gG + gT*gC))
P2 <- SV1[x,x] * (gT*gC/(gA*gG + gT*gC)) 
Q <-  SV2[x,x] 
HKY[x,x]  <- tn93(gT,gC,gG,gA,P1,P2,Q,stat_total_length)
}# End of for 

#HKY for PERMUTATIONS TESTS
ap_erg2 <- apply(apply_id,2,function(pop){
gT <- TCGA[pop[1],"tcga1"] + TCGA[pop[2],"tcga1"]
gC <- TCGA[pop[1],"tcga2"] + TCGA[pop[2],"tcga2"]
gG <- TCGA[pop[1],"tcga3"] + TCGA[pop[2],"tcga3"]
gA <- TCGA[pop[1],"tcga4"] + TCGA[pop[2],"tcga4"]
P1 <- SV1[pop[1],pop[2]] * (gA*gG/(gA*gG + gT*gC))
P2 <- SV1[pop[1],pop[2]] * (gT*gC/(gA*gG + gT*gC)) 
Q <-  SV2[pop[1],pop[2]] 
HKY[pop[1],pop[2]]  <<- tn93(gT,gC,gG,gA,P1,P2,Q,stat_total_length)
})#END APPLY


}#END if flag

FSTHKY  <- matrix(NA,npops,npops) 
FSTPAIR <- matrix(NA,npops,npops)

### --------------------- FST PAIR-PAIR ---------------------------- ###

##APPLY ####
ap_erg3 <- apply(apply_id,2,function(pop){
if(PIA[pop[1],pop[2]] && length(populations[[pop[1]]])>= 1 && length(populations[[pop[2]]]) >=1){
  FSTPAIR[pop[1],pop[2]] <<- 1.0 - ((PIW[pop[1]] + PIW[pop[2]])/2.0)/(PIA[pop[1],pop[2]]) 

}else{
  FSTPAIR[pop[1],pop[2]] <<- NaN
}
  
if(flag==1){# HKY CORRECTION

     if(!is.na(HKY[pop[1],pop[2]])){ 
        if(HKY[pop[1],pop[2]] > 0 && length(populations[[pop[1]]])>= 1 && length(populations[[pop[2]]]) >=1){
           FSTHKY[pop[1],pop[2]] <<- 1.0 - ((diag(HKY)[pop[1]] + diag(HKY)[pop[2]])/2.0)/(HKY[pop[1],pop[2]])   
        }
     }else{
      FSTHKY[pop[1],pop[2]] <<- NaN 
     }  
  }#Ende if FLAG

}) # END of APPLY

colnames(FSTHKY)  <- pop
rownames(FSTHKY)  <- pop
colnames(FSTPAIR) <- pop
rownames(FSTPAIR) <- pop

#### FSTALL

if(npops>1){
 #perm is a permutation matrix which describes the FORS
 #perm <-sapply(s,function(x){return(s[-x])})

 #--------------FST1ALL -----------------------------#

FST1ALL <- rep(0,(npops-1))

# ONLY POPULATIONS !
#########################

if(outgroup[1]!=FALSE){TIL <- (npops-1)}
if(outgroup[1]==FALSE){TIL <-  npops   }
 
for(pop1 in 1:TIL) {
 piw <- 0
 ncw <- 0
   for(pop2 in 1:TIL) {
       if(pop2 != pop1 && length(populations[[pop2]]) >= 1) {
          piw <- piw + PIW[pop2]
          ncw <- ncw + 1
	 }
   }
	piw <- piw/ncw
	nca <- ncw
	if(PIA1ALL[pop1] && length(populations[[pop1]]) >= 1 && ncw > 0){
	   FST1ALL[pop1] <- 1 - ((PIW[pop1] + piw)/2)/(PIA1ALL[pop1]/nca)
	}else{FST1ALL[pop1] <- NaN}
}

apply_id2 <- rbind(apply_id[2,],apply_id[1,])

#APPLY #####------------FSTALL ----------------------------#####

if(outgroup[1]!=FALSE){ # Outgroup
  
  erg <- apply(apply_id2,2,function(pop){

  spiw <-0
  ncw <-0
  spia <-0
  nca <-0

if(length(populations[[pop[1]]]) >= 1) {
     spiw <- PIW[pop[1]]
     ncw  <- 1
}

# Only for Populations: no outgroup !!! look at if condition !!! #
#------------------------------------------------------------------#
if(length(populations[[pop[1]]]) >= 1 && length(populations[[pop[2]]]) >= 1 && pop[2]<npops){
	 spia <- PIA[pop[2],pop[1]]
	 nca <-1
}
vecerg <- c(spiw,ncw,spia,nca)
return(vecerg)

})
rownames(erg) <- c("spiw","ncw","spia","nca")

#### FSTALL ##################################################
# only for populations #
########################################

suni     <- unique(apply_id2[1,]) 
firstfor <- match(suni,apply_id2[1,])

spiw <- sum(erg[1,firstfor])
ncw  <- sum(erg[2,firstfor])
spia <- sum(erg[3,])
nca <-  sum(erg[4,])

 ####---------------------- FSTALL --------------------------####
  if(spia){ 
     FSTALL = 1 - (spiw/ncw)/(spia/nca)
  }else{FSTALL <- 0}

}#END if outgroup[1]!=FALSE

if(outgroup[1]==FALSE){ # No outgroup

erg <- apply(apply_id2,2,function(pop){

  spiw <-0
  ncw  <-0
  spia <-0
  nca  <-0

if(length(populations[[pop[1]]]) >= 1) {
     spiw <-PIW[pop[1]]
     ncw  <- 1
}

# ALL Populations
#------------------------------------------------------------------#
if(length(populations[[pop[1]]]) >= 1 && length(populations[[pop[2]]]) >= 1){
	 spia <- PIA[pop[2],pop[1]]
	 nca <-1
}
vecerg <- c(spiw,ncw,spia,nca)
return(vecerg)

})
rownames(erg) <- c("spiw","ncw","spia","nca")

#### FSTALL ##################################################
# only for populations #
########################################

suni     <-  unique(apply_id2[1,]) 
firstfor <-  match(suni,apply_id2[1,])
spiw     <-  sum(erg[1,firstfor]) + PIW[npops] #auch die letzte Population mit rein
ncw      <-  npops
spia     <-  sum(erg[3,])
nca      <-  sum(erg[4,])

 ####---------------------- FSTALL --------------------------####
  if(spia){ 
     FSTALL = 1 - (spiw/ncw)/(spia/nca)
  }else{FSTALL <- 0}
}


}#Ende if npops >1

#clean Workspace
#rm(SV1)
#rm(SV2)
FSTHKY[FSTHKY==-1000] <- NaN


##### ------------------- FSTHKY2/FSTPAIR2
### 
##### FSTHKY2: store FSTHKY in a comfortable style
##### useful for multi analysis
##################################################

pp <- combn(1:(npops),2)
nn <- paste("pop",pp[1,1],"/pop",pp[2,1],sep="")
if(dim(pp)[2]>1){ # more than 2 Populations
 for(xx in 2:dim(pp)[2]){
    m <- paste("pop",pp[1,xx],"/pop",pp[2,xx],sep="")
    nn <- c(nn,m)
 } 
}#END if


#### Konvert the results in a confortable multi locus approach
FSTHKY2  <- apply(pp,2,function(x){ 
                 return(FSTHKY[x[2],x[1]])})
FSTPAIR2 <- apply(pp,2,function(x){ 
                 return(FSTPAIR[x[2],x[1]])})
PIA2     <- apply(pp,2,function(x){ 
                 return(PIA[x[2],x[1]])})

FSTHKY2             <- as.matrix(FSTHKY2)
row.names(FSTHKY2)  <- nn  
FSTPAIR2            <- as.matrix(FSTPAIR2)
row.names(FSTPAIR2) <- nn 
PIA2                <- as.matrix(PIA2)
row.names(PIA2)     <- nn 
##############################################

diag(PIA)     <- PIW
PIW           <- t(as.matrix(PIW))
colnames(PIW) <- paste("pop",1:npops)

return(list(PIA=PIA,PIA2=PIA2,PIW=PIW,FSTALL=FSTALL,FSTHKY=FSTHKY,FSTPAIR=FSTPAIR,FST1ALL=FST1ALL,SV1=SV1,SV2=SV2,FSTHKY2=FSTHKY2,FSTPAIR2=FSTPAIR2))


}

#------------------------------------------------------------------#
#------------------------------------------------------------------#
#############################
# FUNCTIONS #################
#############################

get_freq <- function(pop,matrix,matrix_sv){

# no apply !!!
# PIW !!!

if(length(pop)==1){ # Only one gene
    mat   <- matrix[pop,]
    lm    <- length(mat)
    freq1 <- numeric(lm)
    freq2 <- numeric(lm)
    freq3 <- numeric(lm)
    
    freq1[mat==0]     <- 1   
    freq2[mat==1]     <- 1
    freq3[is.na(mat)] <- 1
     
    freq0 <- freq1 + freq2

}else{

  if(dim(matrix)[2]>1){ # biallelic matrix
  freq1    <- colSums(matrix[pop,] == 0,na.rm=TRUE)
  freq2    <- colSums(matrix[pop,] == 1,na.rm=TRUE)
  freq3    <- colSums(is.na(matrix[pop,]))
  freq0    <- freq1 + freq2
  }

  if(dim(matrix)[2]==1){ # biallelic vector
  freq1    <- sum(matrix[pop] == 0,na.rm=TRUE)
  freq2    <- sum(matrix[pop] == 1,na.rm=TRUE)
  freq3    <- sum(is.na(matrix[pop]))
  freq0    <- freq1 + freq2
  }

}

FREQ <- rbind(freq1,freq2,freq3,freq0)
row.names(FREQ) <- c("freq1","freq2","freq3","freq0")

id <- which(FREQ["freq0",]>1) # Diejenigen Spalten die mind. eine NULL oder EINS haben

#print(id)
PIW <- (FREQ["freq1",id]*FREQ["freq2",id])/(FREQ["freq0",id]*(FREQ["freq0",id ]-1)/2) 
PIW <- sum(PIW)

id1 <-   which(matrix_sv[id]==TRUE) 
id2 <-   which(matrix_sv[id]==FALSE)

SV0 <- (FREQ["freq1",id1]*FREQ["freq2",id1])/(FREQ["freq0",id1]*(FREQ["freq0",id1]-1)/2)
SV1 <- (FREQ["freq1",id2]*FREQ["freq2",id2])/(FREQ["freq0",id2]*(FREQ["freq0",id2]-1)/2)

SV0 <- sum(SV0)
SV1 <- sum(SV1)

return(list(SV0=SV0,SV1=SV1,PIW=PIW,FREQ=FREQ))

}

