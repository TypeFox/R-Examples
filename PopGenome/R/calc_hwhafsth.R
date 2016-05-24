
calc_hwhafsth <- function(matrix_pol,populations,outgroup=FALSE,simulation=FALSE,only.haplotype.counts=FALSE){


# Imprtant for Coalescent Simulation

if(simulation){
rownames(matrix_pol) <- 1:dim(matrix_pol)[1]
}

# Only one polymorphic site
if(is.vector(matrix_pol)){
   matrix_pol <- as.matrix(matrix_pol)
   warning("#---------------> Only one polymorphic site <-------------------#")
}

if(length(rownames(matrix_pol))==0){rownames(matrix_pol) <- 1:dim(matrix_pol)[1]}

### WICHTIG !!!! ##########################
# EINE DEFINIERTE OUTGROUP WIRD FUER HAPW und HAPA 
# ALS WEITERE POPULATION BETRACHTET !!!!!!!!!!!!!
# HAPW1ALL betrachtet nur die Populationen
######################################

popsize <- length(populations)
npops   <- popsize
init    <- numeric(npops)
init2   <- matrix(0,npops,npops)

#####################
# delete gap-columns#
# define matrix_hap #
#####################

rowcol <- which(is.na(matrix_pol),arr.ind=TRUE)
if(length(rowcol)!= 0){

 col         <- unique(rowcol[,"col"])
 matrix_hap  <- matrix_pol[,-col,drop=FALSE]
 numvec      <- as.numeric(matrix_hap)
 matrix_hap  <- matrix(numvec,ncol=dim(matrix_hap)[2],nrow=dim(matrix_hap)[1])

}else{

 matrix_hap <- matrix_pol
 matrix_hap <- as.numeric(matrix_hap)
 matrix_hap <- matrix(matrix_hap,dim(matrix_pol)[1] , dim(matrix_pol)[2])

}
 

# If gaps in every site # matrix_hap is empty
if(dim(matrix_hap)[2]==0){return(list(matrix_hap=as.matrix(NaN),hapw=0,nhgesamt=NaN,nh=NaN,hapa=as.matrix(NaN),hapamatrix=as.matrix(NaN),hapw1all=NaN,hapa1all=NaN,
fsth1all=NaN,fsth=as.matrix(NaN),fsthmatrix=as.matrix(NaN),fsthALL=NaN,sfreqh=as.matrix(NaN),Gst=as.matrix(NaN),Gstmatrix=as.matrix(NaN),GstAll=NaN,PIW_nei=0,
PIA_nei=as.matrix(NaN),HST=NaN,HSTpair=as.matrix(NaN),Gst_Hudson=NaN,Gst_Hudson_pair=as.matrix(NaN),KST=NaN,fstnALL= NaN))}
#-------------------------------------------


#### INIT hapwvek ######
hapwvek <- init # define hpwvek
nam     <- paste("pop", 1:npops)
names(hapwvek) <- nam
#######################
nh        <- vector(,npops)
names(nh) <- nam


## Calculate Frequencies ###################################
############################################################

rownames(matrix_hap) <- rownames(matrix_pol)
matrix_hap_sub       <- matrix_hap[unique(unlist(populations)),,drop=FALSE]#### Wegen sfreq !!!!

uniquematrix         <- unique(matrix_hap_sub) 
nhgesamt             <- dim(uniquematrix)[1]
sfreqh               <- matrix(NaN,npops,nhgesamt)
 
rownames(sfreqh)     <- nam
colnames(sfreqh)     <- rownames(uniquematrix)
rownames(matrix_hap) <- NULL

 for(xx in 1:npops){
  if(length(populations[[xx]])==0){next;}
  sfreqh[xx,] <- get_sfreqh(uniquematrix,matrix_hap[populations[[xx]],,drop=FALSE])  
  nh[xx]      <- length(which(sfreqh[xx,]!=0)) # number of haplotypes for each population
 
 }

# If there is only one Haplotype in the whole Population
if(nhgesamt==1||only.haplotype.counts){return(list(matrix_hap=as.matrix(NaN),hapw=0,nhgesamt=nhgesamt,nh=nh,hapa=as.matrix(0),hapamatrix=as.matrix(NaN),hapw1all=NaN,hapa1all=NaN,
fsth1all=NaN,fsth=as.matrix(NaN),fsthmatrix=as.matrix(NaN),fsthALL=NaN,sfreqh=sfreqh,Gst=as.matrix(NaN),Gstmatrix=as.matrix(NaN),GstAll=NaN,PIW_nei=0,
PIA_nei=as.matrix(NaN),HST=NaN,HSTpair=as.matrix(NaN),Gst_Hudson=NaN,Gst_Hudson_pair=as.matrix(NaN),KST=NaN,fstnALL= NaN))}
#-------------------------------------------



### Generate Hap Compare Matrix
### NA hier Substitutionsmodelle
#####################################################################

if(nhgesamt > 1){

SUBMATRIX <<- matrix(,nhgesamt,nhgesamt)
colnames(SUBMATRIX) <- rownames(uniquematrix)
rownames(SUBMATRIX) <- rownames(uniquematrix)

comb      <- combn(nhgesamt,2)
             apply(comb,2,function(pop){
             hap1 <- uniquematrix[pop[1],]
             hap2 <- uniquematrix[pop[2],]
             SUBMATRIX[pop[1],pop[2]] <<- sum(hap1!=hap2)
             SUBMATRIX[pop[2],pop[1]] <<- sum(hap1!=hap2)
             })
             
SUBMATRIX <- SUBMATRIX             

}
####################################################################


# count hapw for each population within diversities
# Within Diversities 
##################################################

# Create happairs ! perhaps there is al better internal R - Function !!!

happairsbetween <- NULL
vek      <- 1:dim(uniquematrix)[1]

for(xx in 1:dim(uniquematrix)[1]){

com <- vek[-xx]
val <- rep(xx,length(com))
mat <- rbind(val,com)

happairsbetween <- cbind(happairsbetween,mat) 

}


happairswithin      <- combn(dim(uniquematrix)[1],2)


Ki            <- init  # Hudson Paper :  A statistical Test for Detecting Geographic Subdivision
names(Ki)     <- nam
Kistar        <- init  # Hudson Paper :  A statistical Test for Detecting Geographic Subdivision K*
names(Kistar) <- nam

for(xx in 1:npops){
  
  if(length(populations[[xx]])==0){next;} # Wenn eine Population nicht vorkommt ! Multi Locus PopGenome Ansatz !  kann weg !
  
     div     <- apply(happairswithin,2,function(hap){
                freq       <- sfreqh[xx,]
                return(freq[hap[1]]*freq[hap[2]])
     }) 
     
     divsubX  <- apply(happairswithin,2,function(hap){ # Look at Distance
                freq       <- sfreqh[xx,]
                return(c(freq[hap[1]]*freq[hap[2]]*SUBMATRIX[hap[1],hap[2]],freq[hap[1]]*freq[hap[2]]*log(1+SUBMATRIX[hap[1],hap[2]])))
     }) 
               
  divsub      <- divsubX[1,]
  divsub2     <- divsubX[2,]
  hapwvek[xx] <- sum(div,na.rm=TRUE)
  p_size      <- length(populations[[xx]])
  Ki[xx]      <- sum(divsub,na.rm=TRUE)     # Hudson Paper 1992  within diversities Ki
  Kistar[xx]  <- sum(divsub2,na.rm=TRUE)    # Hudson Paper 1992  within diversities Ki*

  if(p_size>1){   
     hapwvek[xx]     <- hapwvek[xx]/((p_size*(p_size-1))/2)
     Ki[xx]          <- Ki[xx]/((p_size*(p_size-1))/2)     # Hudson Paper 1992 within diversities Ki
     Kistar[xx]      <- Kistar[xx]/((p_size*(p_size-1))/2) # Hudson Paper 1992 within diversities Ki*
 }
}


# Do the same for population pairs

if(npops > 1){

 pairs <- combn(npops,2)
 #--------------
 
 KT           <- matrix(NaN,npops,npops) # Hudson 92 between diversities KT
 hapamatrix   <- matrix(NaN,npops,npops)
 hapamatrixN  <- matrix(NaN,npops,npops)

 # Apply
 hapavek     <- apply(pairs,2,function(x){
     
     hapa <- NaN
  if(length(populations[[x[1]]])!=0 &  length(populations[[x[2]]])!=0){
     
      m1_size <- length(populations[[ x[1] ]]) # size of population 1
      m2_size <- length(populations[[ x[2] ]]) # size of population 2     
     
     freqpop1 <- sfreqh[x[1],]
     freqpop2 <- sfreqh[x[2],]
     
            div     <-  apply(happairsbetween,2,function(hap){
                        return(freqpop1[hap[1]]*freqpop2[hap[2]])
            })
            
            divsub  <-  apply(happairsbetween,2,function(hap){
                        return(freqpop1[hap[1]]*freqpop2[hap[2]]*SUBMATRIX[hap[1],hap[2]]) # Look at Distance

            })
     
     hapa  <- sum(div)
     kk    <- sum(divsub)
     hapaN <- kk/(m1_size*m2_size)
     hapa  <- hapa/(m1_size*m2_size)
     hapamatrix[x[1],x[2]]  <<- hapa
     hapamatrixN[x[1],x[2]] <<- hapaN
     KT[x[1],x[2]]  <<- (Ki[x[1]] + Ki[x[2]] + kk)/choose(m1_size+m2_size,2) # Hudson 1992 KT
  }  
     return(hapa)

}) 


  
#### --- Names of population pairs --- #### 
nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
if(dim(pairs)[2]>1){ # more than 2 Populations
 for(xx in 2:dim(pairs)[2]){
    m <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
 } 
}#END if
### --- ------------------------------ ####

hapavek           <- as.matrix(hapavek)
rownames(hapavek) <- nn
colnames(hapavek) <- "Diversity between populations(Haplotype)"

}else{hapavek <- as.matrix(NaN);hapamatrix <- NaN} # End of if (popsize > 1)

## -------------- #########################
##   What we have #####
# Hapw for population pairs ---> hapavek
#########################################

if(npops > 1){
 
  freqall  <- colSums(sfreqh,na.rm=TRUE) 
  div      <- apply(happairswithin,2,function(hap){
                   return(freqall[hap[1]]*freqall[hap[2]])
             }) 
  
  hapw1all <- sum(div,na.rm=TRUE)
  hapw1all <- hapw1all/(npops*(npops-1)/2)
  #### ENDE hapw1all #####################################
  
  
  ########### hapa1all/fsth/fsthmatrix ## only populations ########
  #################################################################
  
  fsthmatrix <- init2  
  hapa1all   <- init
  poppairs   <- combn(npops,2)
  
  fsthvek  <- apply(poppairs,2,function(x){         
     
            pop1 <- x[1]
            pop2 <- x[2]
            hapa1all[pop1] <<- hapa1all[pop1] + hapamatrix[pop1,pop2] 
            hapa1all[pop2] <<- hapa1all[pop2] + hapamatrix[pop1,pop2]
            fsth <- 1 - ((hapwvek[pop1]+hapwvek[pop2])/2)/hapamatrix[pop1,pop2]
            fsthmatrix[pop2,pop1] <<- fsth
            return(fsth)
  
  }) 
  
  fsthvek <- as.matrix(fsthvek)
  rownames(fsthvek) <-  nn
  colnames(fsthvek) <- "fsth"

}else{hapw1all <- NaN;hapa1all <- NaN;fsthvek <- as.matrix(NaN);fsthmatrix <- as.matrix(NaN)} # End of if(only_populations >1)



##################################################################### 
######calculate Gst (Nei 1973) for pair-pair comparisons ############
###### GST ##########################################################
#-------------------------------------------------------------------#
 
 
##### GST #############################################################

Dklmatrix <- init2 # Fuer Gstall wichtig !
rownames(Dklmatrix) <- nam
colnames(Dklmatrix) <- nam
Gstmatrix <- init2


if(npops>1){  
 
 Gst <-  apply(pairs,2,function(x){

    pop1 <- x[1]
    pop2 <- x[2]
    if(length(populations[[pop1]]) == 0 | length(populations[[pop2]])==0){
      Gstmatrix[pop2,pop1] <<- NaN
      Dklmatrix[pop1,pop2] <<- NaN
      return(NaN) 
    } 
    ss_1 <- sfreqh[pop1,]*sfreqh[pop1,]/(length(populations[[pop1]]))^2
    ss_2 <- sfreqh[pop2,]*sfreqh[pop2,]/(length(populations[[pop2]]))^2 
    Js_1 <- sum(ss_1,na.rm=TRUE)
    Js_2 <- sum(ss_2,na.rm=TRUE)
    Js   <- Js_1 + Js_2
     
    size_pop1 <- length(populations[[pop1]])
    size_pop2 <- length(populations[[pop2]])
       
    temp <- (sfreqh[pop1,]/size_pop1 - sfreqh[pop2,]/size_pop2)^2/2
    Dkl  <- sum(temp,na.rm=TRUE)
        
    Js  <- Js/2
    Hs  <- 1- Js
    Dm  <- Dkl/2
    Ht  <- Hs + Dm
    Gst <- Dm/Ht
    
    Gstmatrix[pop2,pop1] <<- Gst
    Dklmatrix[pop1,pop2] <<- Dkl # Fuer GstAll wichtig !
      
  return(Gst)  
     
 })
  Gst <- as.matrix(Gst)
  rownames(Gst) <- nn
  colnames(Gst) <- "Gst"
  
 }#End if pop_and_outgroup>1
 else{Gst<-as.matrix(NaN)} 
 
### FSTH1ALL ###############################################################
#--------------------------------------------------------------------------#
### only for Populationen ##################################################
#--------------------------------------------------------------------------#
### Across the rest !!!!! ################################################## <----- Across all the rest

if(npops > 1){

fsth1all <- init
hapwvek2 <- hapwvek[1:npops]
names(fsth1all) <- paste("pop",1:npops,"/rest",sep="")

 for (xx in 1:npops){
    nca  <- npops - 1  # Anzahl der uebrigen Populationen
    hapw <- sum(hapwvek2[-xx],na.rm=TRUE)/nca
    fsth1all[xx] <- 1 - ((hapwvek[xx]+hapw)/2)/(hapa1all[xx]/nca)
 }  

}else{fsth1all <- NaN} # End if(only_populations)

## ENDE FSTH1ALL ##########################################################  
####### FSTHALL #### ----- ONLY THE POPULATIONS ----------- !!!!!!!!!!!!!!!
  
if(npops > 1){
   
   shapw  <- sum(hapwvek[1:npops],na.rm=TRUE)
   shapwN <- sum(Ki[1:npops],na.rm=TRUE)
   ncw    <- npops
   
   
   sh   <- apply(poppairs,2,function(x){
         pop1 <- x[1]
         pop2 <- x[2]
         return(c(hapamatrix[pop1,pop2],hapamatrixN[pop1,pop2]))      
   })
   
   shapa  <- sum(sh[1,],na.rm=TRUE)
   shapaN <- sum(sh[2,],na.rm=TRUE)
   nca    <- length(sh[1,])
   
   if(shapa){
      fsthALL <- 1 - (shapw/ncw)/(shapa/nca)
      fstnALL <- 1 - (shapwN/ncw)/(shapaN/nca)
   }else{fsthALL <- NaN;fstnALL <- NaN}

}else{fsthALL <- NaN;fstnALL<- NaN} # End if(only_populations >1)

###################### END FSTHALL #########################################
  
###################### GstALL ##############################################
#################   calculate GstALL (G'st = Dm/H't) #######################
#################   !!!  Only Populations !!!        #######################
############################################################################


if(npops > 1){

nhp    <- 0
Jstemp <- init

for(xx in 1:npops){
   if(length(populations[[xx]])>0){nhp <- nhp + 1}
   size_pop   <- length(populations[[xx]])      
   Jstemp[xx] <- sum((sfreqh[xx,])^2/(size_pop)^2,na.rm=TRUE)
}

 Js <- sum(Jstemp,na.rm=TRUE)
 Js <- Js/nhp
 Hs <- 1-Js
 
# Dklvector <- init
# print(poppairs)
 
 Dklvector <- apply(poppairs,2,function(x){
     
            pop1 <- x[1]
            pop2 <- x[2]
            
            return(Dklmatrix[pop1,pop2])
            
           # if(!is.na(Dklmatrix[pop1,pop2])){
           #  Dklvector[pop1] <<- Dklvector[pop1] + Dklmatrix[pop1,pop2] 
           #  Dklvector[pop2] <<- Dklvector[pop2] + Dklmatrix[pop1,pop2]      
           # } 
 })
 
 Dm       <- sum(Dklvector,na.rm=TRUE)
 rm(Dklvector)
 Dm       <- Dm/((nhp*(nhp-1)))
 Ht       <- Hs + Dm
 GstAll   <- Dm/Ht


}else{GstAll <- NaN} # End if(only_populations > 1)

## calcPi (within) #########################################################

PIW_nei          <- rep(NaN,npops)
names(PIW_nei)   <- nam
H_within         <- rep(NaN,npops)
names(H_within)  <- nam
HAPIDS           <- colnames(sfreqh)

rownames(matrix_hap) <- rownames(matrix_pol)

for(xx in 1:npops){

   if(length(populations[[xx]])==0){next;}
   
   freqq    <- sfreqh[xx, ]/length(populations[[xx]]) 
   freqq    <- freqq[which(freqq!=0)]

   
   if(nh[xx]>1){vergl <-combn(nh[xx],2)}else{PIW_nei[xx]<-0;H_within[xx]<-0;next;}   
   
   res      <- apply(vergl,2,function(x){
   hap1     <- matrix_hap[names(freqq[x[1]]),]
   hap2     <- matrix_hap[names(freqq[x[2]]),]
   div      <- hap1!=hap2 # ohne substituitionsmodell
   div      <- sum(div)
   return(2*freqq[x[1]]*freqq[x[2]]*div)               
   })

n            <- length(populations[[xx]])
PIW_nei[xx]  <- (n/(n-1))*sum(res,na.rm=TRUE)
H_within[xx] <- (n/(n-1))*(1-sum((sfreqh[xx,]/n)^2,na.rm=TRUE)) # Hudson A statistical Test for Detecting Geographic Subdivision eq (4)

}

## calcPi (between)  #######################################################

if(npops>1){

 PIA_nei <- apply(pairs,2,function(x){         
      
      if(length(populations[[x[1]]])!=0 & length(populations[[x[2]]])!=0){ 
           
           n1        <- length(populations[[x[1]]]) 
           n2        <- length(populations[[x[2]]]) 
           freqq1    <- sfreqh[x[1], ]/n1 
           freqq1    <- freqq1[which(freqq1!=0)]
           freqq2    <- sfreqh[x[2], ]/n2 
           freqq2    <- freqq2[which(freqq2!=0)]
             
             pi <- 0
             for(xx in 1:length(freqq1)){
                 for(yy in 1:length(freqq2)){
                    hap1     <- matrix_hap[names(freqq1[xx]),]
                    hap2     <- matrix_hap[names(freqq2[yy]),]
                    div      <- hap1!=hap2
                    div      <- sum(div) 
                    pi       <- pi + freqq1[xx]*freqq2[yy]*div 
                 }
             }   
                                              
       return((n1+n2)/((n1+n2)-1)*pi)
 
      }else{return(NaN)} # if !=0
}) 

PIA_nei           <- as.matrix(PIA_nei)
rownames(PIA_nei) <- nn

}else{PIA_nei <- as.matrix(NaN)} # pop_and_outgroup

### PAIR - PAIR Comparison #######################
##################################################


##################################################################### 
######calculate Gst (Nei 1973) for pair-pair comparisons ############
###### GST2 #########################################################
#-------------------------------------------------------------------#
# Note this calculation is from the Hudson Paper !!!!!!!!!!!!!!!!!!!!
##### GST2 ##########################################################

Gst2matrix <- init2
HSTmatrix  <- init2 
H_between  <- init2
K_pairwise <- init2
Ks         <- init2

if(npops>1){  
 
 Gst2 <-  apply(pairs,2,function(x){

      pop1 <- x[1]
      pop2 <- x[2]
    
    if(length(populations[[pop1]]) == 0 | length(populations[[pop2]])==0){
      Gst2matrix[pop2,pop1] <<- NaN
      return(NaN) 
    } 
    
    size_pop1  <- length(populations[[pop1]])
    size_pop2  <- length(populations[[pop2]])
    
    ss_1 <- (sfreqh[pop1,]/size_pop1)^2 # frequencies p^2
    ss_2 <- (sfreqh[pop2,]/size_pop2)^2 # frequencies 
         
    totalsize  <- size_pop1 + size_pop2
    
    H1 <-  (size_pop1/(size_pop1-1)) * (1 - sum(ss_1,na.rm=TRUE))
    H2 <-  (size_pop2/(size_pop2-1)) * (1 - sum(ss_2,na.rm=TRUE))
    
    w1   <- size_pop1/totalsize
    w2   <- size_pop2/totalsize
    m1   <- (size_pop1-2)/(totalsize-4)
    m2   <- (size_pop2-2)/(totalsize-4)
 
    Hs                     <- w1*H1 + w2*H2
    Hs2                    <- m1*H1 + m2*H2
    H_between[pop1,pop2]   <<- Hs2
    Ks[pop1,pop2]          <<- (size_pop1/totalsize)*Ki[pop1] + (1-(size_pop1/totalsize))*Ki[pop2] 
    K_pairwise[pop1,pop2]  <<- 1 - Ks[pop1,pop2]/KT[pop1,pop2]
    
    totalfreq       <- colSums(sfreqh[c(pop1,pop2),])/totalsize
    Ht_dash         <- 1 - sum(totalfreq^2) + Hs/(2/(sum(1/c(size_pop1,size_pop2)))*2)   
    Ht              <- (totalsize/(totalsize-1))* (1-sum(totalfreq^2))
    HST             <- 1 - Hs2/Ht
    Gst2            <- 1 - Hs/Ht_dash
    
    Gst2matrix[pop2,pop1] <<- Gst2
    HSTmatrix[pop2,pop1]  <<- HST
    
    return(Gst2)
         
 })
 
  Ks          <- Ks
  H_between   <- H_between
  HSTmatrix   <- HSTmatrix
  K_pairwise  <- K_pairwise
  Gst2 <- as.matrix(Gst2)
  rownames(Gst2) <- nn
  colnames(Gst2) <- "Gst (Hudson)"
   
 }#End if pop_and_outgroup>1
 
 else{Gst2<-as.matrix(NaN)} 


############################################################################################
# ALL ###
# ###########################################################################################
# Hudson : GST/HST  #########################################################################
# A statistical Test for Detecting Geographic Subdivision
#############################################################################################

###########
##ALL######
###########

if(npops > 1){

# calculate weightings
nx  <- lapply(populations,function(pop){
       return(length(pop))
       })
nx  <- unlist(nx)       
w   <- nx/sum(nx)
# -------------------------   
   
Hs  <- w*H_within
Hs  <- sum(Hs) # within diversity

Hs2 <- sum(H_between,na.rm=TRUE)/choose(npops,2)

totalsize <- length(unlist(populations))
totalfreq <- colSums(sfreqh)/totalsize

Ht              <- (totalsize/(totalsize-1))* (1-sum(totalfreq^2))
 
#Ht              <- sum(H_between,na.rm=TRUE)/choose(npops,2)
Ht_dash         <- 1 - sum(totalfreq^2) + Hs/(length(nx)/(sum(1/nx))*npops)

KST             <- 1 - mean(Ks,na.rm=TRUE)/mean(KT,na.rm=TRUE) # have to be tested !!!
HST             <- 1 - Hs2/Ht          # Hudson eq (2)  works well !
Gst_Hudson      <- 1 - Hs/Ht_dash      # Hudson eq (5)  works well !    

}else{HST<- NaN;Gst_Hudson<-NaN; KST <- NaN}

rownames(Gstmatrix)  <- nam
colnames(Gstmatrix)  <- nam
rownames(HSTmatrix)  <- nam
colnames(HSTmatrix)  <- nam
rownames(fsthmatrix) <- nam
colnames(fsthmatrix) <- nam

if(npops>1){
 hapamatrix           <- t(hapamatrix)
 diag(hapamatrix)     <- hapwvek
 rownames(hapamatrix) <- nam
 colnames(hapamatrix) <- nam
}else{
 hapamatrix <- as.matrix(hapwvek)
}

# fsth ist haplotype.FST for population pairs

return(list(matrix_hap=matrix_hap,hapw=hapwvek,nhgesamt=nhgesamt,nh=nh,hapa=hapavek,hapamatrix=hapamatrix,hapw1all=hapw1all,hapa1all=hapa1all,
fsth1all=fsth1all,fsth=fsthvek,fsthmatrix=fsthmatrix,fsthALL=fsthALL,sfreqh=sfreqh,Gst=Gst,Gstmatrix=Gstmatrix,GstAll=GstAll,PIW_nei=PIW_nei,
PIA_nei=PIA_nei,HST=HST,HSTpair=HSTmatrix,Gst_Hudson=Gst_Hudson,Gst_Hudson_pair=Gst2,KST=KST,fstnALL=fstnALL))

} # End of Function hwhafsth 


###############################################
###### FUNCTIONS ##############################
###############################################

############## get_sfreqh ######################
get_sfreqh <- function(uniquematrix, matrix){ 
       
 matrix_a    <- uniquematrix
 matrix_b    <- matrix

 # if(is.vector(matrix_b)){matrix_b <- as.matrix(t(matrix_b))} # wenn eine Population nur eine Sequenz
 
apply(matrix_a,1,function(x){ 
     comp_a <- x
     freq   <- 0
     for (yy in 1:dim(matrix_b)[1]){
         comp_b <- matrix_b[yy,]
        # if(all(comp_a==comp_b)){
        #    freq <- freq + 1      
        # }
	  
 	  if(.Call("Ccompare",comp_a,comp_b)){
		freq <- freq + 1
	  }    
     }
 return(freq)
     }) 
     
}



