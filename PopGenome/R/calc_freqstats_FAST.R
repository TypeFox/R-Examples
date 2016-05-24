
calc_freqstats_FAST<- function(matrix,populations,outgroup=FALSE,data=NULL,methods=FALSE,Achaz=FALSE){


############################################
# data is a class returned from get_data_test
############################################

nsamtot <-  dim(matrix)[1]
npops   <-  length(populations)
sfreq   <-  paste("sfreq",0:nsamtot)

####################################################
# Following values we need from data
# LENGTH und TOTAL_LENGTH
####################################################

data_length <- numeric(npops)

if(length(data)!=0){
  #LENGTH (NUMBER OF NUCLEOTIDES)
  for (x in 1:npops){
    data_length[x] <- sum(data$n.nucleotides[populations[[x]]])
  }
  #NUMBER OF ALL POSITIONS
  data_total_length <- data$n.valid.sites
  
}else{ # no data defined
  samplesize <- numeric(npops)
  for(xx in 1:npops){
     #samplesize[xx] <- length(populations[[xx]])
     #correct in case there are NaNs in the biallelic.matrix (SNP data, sliding windows) 
      samplesize[xx] <- sum(!is.na(matrix[populations[[xx]],,drop=FALSE]))/dim(matrix)[2]
  }  
}

##################################################
##################################################

if(outgroup[1]==FALSE){ # <----------------------------------------------------------------------------- without Outgroup

# INIT
# popnames               <- paste("pop",1:npops)
FREQ_ERG               <- vector("list",npops)
THETA_ERG              <- matrix(0,6,npops) # 6 THETA VALUES
# rownames(THETA_ERG)    <- c("S","thetaS","thetaT","thetaFL","thetaSA","thetaTA")
# colnames(THETA_ERG)    <- popnames

init       <- rep(NaN,npops) 

Taj_D      <- init # Tajima
#names(Taj_D)  <- popnames
FuLi_F <- init # Fu Lis F*
#names(FuLi_F) <- popnames
FuLi_D <- init # Fu Lis D*
#names(FuLi_D) <- popnames
Yach        <-init
#names(Yach) <- popnames
HnFw        <- init
#names(HnFw) <- popnames
Ez          <- init
#names(Ez)   <- popnames


for(xx in 1:npops){
  
 # if(length(populations[[xx]])==0){next;}
 
  FREQ    <- vector(,(nsamtot+1))  
  matrix1 <- matrix[populations[[xx]],,drop=FALSE]
  if(is.matrix(matrix1)){nsamtot <- dim(matrix1)[1]}
  if(is.vector(matrix1)){nsamtot <- 1}
  if(is.matrix(matrix1)){ee  <-  .Call("ap_pop_C",matrix1)}
  #if(is.matrix(matrix1)){ee <-  ap_pop_FAST(matrix1)}
  #if(is.matrix(matrix1)){ee <-  apply(matrix1,2,ap_pop_FAST)} # APPLY
  if(is.vector(matrix1)){ee  <-  sapply(matrix1,ap_pop_FAST)}  # Only one sequence
  # row.names(ee) <- c("S","thetaS","thetaT","thetaFL","thetaSA","thetaTA",sfreq[1:(nsamtot+1)])
  # print(ee[,1:4])
  # stop("")
  erg                   <- rowSums(ee)
  rm(ee)
  # names(erg)            <- row.names(ee)
  FREQ_ERG[[xx]]        <- erg[7:length(erg)] # note if you want to add other theta values  
  THETA_ERG[,xx]        <- erg[1:6]
  # names(FREQ_ERG[[xx]]) <- sfreq[1:(nsamtot+1)]
  
  
  #### STATISTICS #####################
  ####TAJIMA ##########
  if(length(data)!=0){ n <- ceiling(data_length[xx]/data_total_length)}else{n <- samplesize[xx]}
  coef      <- init_coef(n)
  Taj_D[xx] <- tajima_d(THETA_ERG[3,xx],THETA_ERG[1,xx],coef)
  ########################################
  #####  Fu and Li's F* ##################
  FuLi_F[xx] <- fl_f2(n,THETA_ERG[4,xx],THETA_ERG[1,xx],THETA_ERG[3,xx],coef)
  ##### Fu and Li's D* ###################
  len <- length(FREQ_ERG[[xx]])
  FuLi_D[xx] <- fl_d2(n,FREQ_ERG[[xx]][2]+FREQ_ERG[[xx]][len-1],THETA_ERG[1,xx],coef) # use singletons
  ########################################
  ######## ACHAZ #########################
  
if(Achaz){

  init    <- vector(,floor(length(populations[[xx]])/2))
  w1      <- init
  w2      <- init
  sfreqn  <- init  
  popsize <- length(populations[[xx]])
  sfreq   <- FREQ_ERG[[xx]][2:popsize] # monomorph positions are not relevant # ACHTUNG bei sehr kleinen Populationen !!!
  
  for(x in 1:floor(popsize/2)) {
				if(x == 1) {w1[x] <- 0
				}else{
					if(x == popsize-x){w1[x] <- popsize/2.0
					}else{w1[x] <- popsize}
				}
				if(x == 1){w2[x] <- 0
				}else{
					if(x == popsize-x){ w2[x] <- popsize/(x*(popsize - x)*2.0)
					}else{ w2[x] <- popsize/(x*(popsize - x))}
				}
				if(x == popsize-x){ sfreqn[x] <- as.integer(((sfreq[x]+ sfreq[popsize-x])/2))
				}else{sfreqn[x] <- sfreq[x] + sfreq[popsize-x]
				}
  }
  Yach[[xx]] <- freqtestn_achaz(popsize,sfreqn,THETA_ERG["S",xx],w1,w2)
} # End if Achaz  
 
 }# End for ueber alle populationen

  #----------------------------------------------------------------
  THETA_ERG             <- rbind(THETA_ERG,matrix(NaN,2,npops))
  # rownames(THETA_ERG)   <- c("S","thetaS","thetaT","thetaFL","thetaSA","thetaTA","thetaFW","thetaL")
  # ---------------------------------------------------------------
 FREQ_ERG <- as.matrix(FREQ_ERG)
 return(list(FREQ=FREQ_ERG,THETA=THETA_ERG,taj_D=Taj_D,FuLi_F=FuLi_F,FuLi_D=FuLi_D,Yach=Yach,HnFw=HnFw,Ez=Ez))

}#END: IF NO OUTGROUP


if(outgroup[1]!=FALSE){ #<------------------------------------------------------------------------------------ Outgroup included

    outgroup_size <- length(outgroup)
    outgroup      <- matrix[outgroup,,drop=FALSE]
   
   # Get ANCESTRAL  
   if(outgroup_size>1){
      ancestral <- apply(outgroup,2,get_outgroup)
   }else{
      ancestral <- sapply(outgroup,get_outgroup)
   }
   #------------ 
  
  FREQ_OUTGROUP1 <- ancestral[-1,]        # FREQ ANCESTRAL
  ancestral      <- ancestral[1, ]        # ANCESTRAL
 
 if(is.matrix(FREQ_OUTGROUP1)){FREQ_OUTGROUP        <- rowSums(FREQ_OUTGROUP1)}
 if(is.vector(FREQ_OUTGROUP1)){FREQ_OUTGROUP        <- sum(FREQ_OUTGROUP1)}
 if(is.matrix(FREQ_OUTGROUP1)){names(FREQ_OUTGROUP) <- sfreq[1:(outgroup_size+1)]}
 
 #------------------------------------
 #INIT# ------------------------------
 
  FREQ_ERG               <- vector("list",npops+1)# outgroup included
  FREQ_ERG               <- t(as.matrix(FREQ_ERG))
  popnames               <- paste("pop",1:npops)
  colnames(FREQ_ERG)     <- c(paste("pop",1:npops),"outgroup")
  THETA_ERG              <- matrix(0,8,npops)      # 6 THETA VALUES
  rownames(THETA_ERG)    <- c("S","thetaS","thetaT","thetaFL","thetaSA","thetaTA","thetaFW","thetaL")
  colnames(THETA_ERG)    <- paste("pop",1:npops)
  FREQ_ERG[[npops+1]]    <- FREQ_OUTGROUP
  
  init          <- rep(NaN,npops)
  Taj_D         <- init
  #names(Taj_D)  <- popnames
  FuLi_D <- init
  #names(FuLi_D) <- popnames
  FuLi_F        <- init
  #names(FuLi_F) <- popnames
  HnFw        <- init   # Fay Wu normalized
  #names(HnFw) <- popnames
  Ez          <- init   #  E Zeng
  #names(Ez)   <- popnames
  Yach        <- init
  #names(Yach) <- popnames

for(xx in 1:npops){

       m <- matrix[populations[[xx]],]
       if(is.vector(m)){length <- 1 } # only one gene in the population
       if(is.matrix(m)){length <- dim(m)[1] }
       m <- rbind(m,ancestral)
       # if(is.matrix(m)){erg <- apply(m,2,ap_pop_ancestral)}
       if(is.matrix(m)){erg <- .Call("ap_pop_ancestral_C",m)}
       if(is.vector(m)){erg <- apply(m,1,ap_pop_ancestral)}
       row.names(erg) <- c("S","thetaS","thetaT","thetaFL","thetaSA","thetaTA","thetaFW","thetaL",sfreq[1:(length+1)])
       erg <- rowSums(erg)
       THETA_ERG[,xx]  <- erg[1:8]
       FREQ_ERG[[xx]]  <- erg[9:length(erg)]       

## ----- STATISTICS ----------------------------------------------------
  ####TAJIMA ####################################
  if(length(data)!=0){ n <- ceiling(data_length[xx]/data_total_length)}else{n <- samplesize[xx]}
  coef       <- init_coef(n)
  Taj_D[xx]  <- tajima_d(THETA_ERG["thetaT",xx],THETA_ERG["S",xx],coef)
  ############################################### 
  ### FU LI D* ##################################
  FuLi_D[xx] <- fl_d(n,FREQ_ERG[[xx]]["sfreq 1"],THETA_ERG["S",xx],coef)
  ### FU LI F* ##################################
  FuLi_F[xx] <- fl_f(n,FREQ_ERG[[xx]]["sfreq 1"],THETA_ERG["S",xx],THETA_ERG["thetaT",xx],coef) 
  ### Fay Wu normalized #########################
  HnFw[xx]   <- fay_wu_normalized2(n,THETA_ERG["thetaL",xx],THETA_ERG["thetaS",xx],THETA_ERG["S",xx],coef,THETA_ERG["thetaT",xx])
  ###  E Zeng ###################################
  Ez[xx]     <- E_zeng(n,THETA_ERG["thetaL",xx],THETA_ERG["thetaS",xx],THETA_ERG["S",xx],coef)
  ############ Achaz ############################
   if(Achaz){ 
    w1 <- c(0, seq(length(populations[[xx]])-2,1)) # ?? 0 3 2 1 # Population size 5
    w2 <- c(0,1/(3:length(populations[[xx]])-1))   # ?? 0 1/2 1/3 1/4 # Population size 5
    Yach[xx]   <- freqtesto_achaz(length(populations[[xx]]),FREQ_ERG[[xx]],THETA_ERG["S",xx],w1,w2)
   }
  ################################################
  
 }# End of for

 
 
}# End of outgroup =TRUE


return(list(FREQ=FREQ_ERG,THETA=THETA_ERG,taj_D=Taj_D,FuLi_D=FuLi_D,FuLi_F=FuLi_F,HnFw=HnFw,Ez=Ez,Yach=Yach)) #return a list


}#END of function calc_freqstats2

#----------- FUNCTIONS ------------------------------#

###############################################
# FUNCTIONS return theta values for each column
# (APPLY)
###############################################

ap_pop_ancestral <- function(mat){

#POP
       S       <-0
       thetaS  <-0
       thetaT  <-0
       thetaFL <-0
       thetaSA <-0
       thetaTA <-0
       thetaFW <- 0
       thetaL  <- 0
       ende    <- length(mat)
       freq    <- numeric(4)
       sfreq   <- numeric(length(mat)) # EINSEN, NULLEN
       freq[2] <- sum(mat[1:(ende-1)]==0,na.rm=TRUE)
       freq[3] <- sum(mat[1:(ende-1)]==1,na.rm=TRUE)
       freq[4] <- sum(is.na(mat[1:(ende-1)]))
       freq[1] <- freq[2] + freq[3]
     
     if(freq[1]){
         if(mat[ende]==0){    # ancestral is  0
           sfreq[freq[3]+1]<- 1 # Anzahl der Einsen (MUTATIONEN) # number of mutations
           if((freq[3]>0) && (freq[3]<freq[1])){S <- 1}
         }
         if(mat[ende]==1){     #ancestral is  1
           sfreq[freq[2]+1] <- 1 # ANZAHL der Nullen (MUTATIONEN) # number of mutations
           if((freq[2]>0) && (freq[2]<freq[1])){
              S<-1
           }
           fr<- freq[2]
           freq[2]<-freq[3]
           freq[3] <- fr
         }
      }

     # For missing Values

     if((freq[3]>0) && (freq[3]<freq[1])){

        an <- sum(1/1:(freq[1]-1)) # modified

        thetaS <- 1/an
        thetaT <- (freq[2]*freq[3])/(freq[1]*(freq[1]-1)/2)

        if(mat[ende]==1 || mat[ende]==0){
           if(freq[3]==1){thetaFL <- 1}
           thetaFW <-(freq[3]*freq[3])/(freq[1]*(freq[1]-1)/2)
           thetaL <- freq[3]/(freq[1]-1)
           if(freq[3]>1){
              thetaSA <- (1/(an-1))
              thetaTA <- (freq[3]*freq[2]/(freq[1]*(freq[1]-1)/2)*freq[1]/(freq[1]-2))
           }
        }
     }

   matrixerg    <- numeric(8)
   matrixerg[1] <- S
   matrixerg[2] <- thetaS
   matrixerg[3] <- thetaT
   matrixerg[4] <- thetaFL
   matrixerg[5] <- thetaSA
   matrixerg[6] <- thetaTA
   matrixerg[7] <- thetaFW
   matrixerg[8] <- thetaL
   
   matrixerg <- c(matrixerg,sfreq)

return(matrixerg)

}#End: function ap_pop_ancestral

#----- FUNCTION: ap_pop -----------------------------------------#
   
ap_pop_FAST <- function(matrix){

   
       freq       <- numeric(3) 
       matrixerg  <- numeric(6) # + dim(matrix)[1]
       S          <-0
       thetaS     <-0
       thetaT     <-0
       thetaFL    <-0
       thetaSA    <-0
       thetaTA    <-0

  res <- apply(matrix,2,function(mat){

       le      <- length(mat) + 1
       sfreq   <- numeric(le) # EINSEN, NULLEN
       # freq    <- numeric(4)
        
       # modified
       freq[2] <- sum(mat==0,na.rm=TRUE)
       freq[3] <- sum(mat==1,na.rm=TRUE)
       
       # freq[4] <- sum(is.na(mat))
       freq[1] <- freq[2] + freq[3]
       
       if(freq[1]){
          sfreq[freq[3]+1] <- 1 # Anzahl der Einsen in der Spalte
          if((freq[3]>0) && (freq[3]< freq[1])){
             S <- 1
          } # S !!!
       }

       # For missing Values // or all ONES

       if((freq[3] > 0) && (freq[3]<freq[1])){

          an <- sum(1/1:(freq[1]-1)) # modified

          thetaS <-  1/an
          thetaT <- (freq[2]*freq[3])/((freq[1]*(freq[1]-1))/2)
          if((freq[3]==1) || (freq[2] ==1)){thetaFL <-  1*(freq[1]-1)/freq[1]}
          if((freq[3] > 1) && (freq[2] >1)){
             thetaSA <- 1/(an - freq[1]/(freq[1]-1))
             thetaTA <- freq[3]*freq[2]/(freq[1]*(freq[1]-1)/2) * (freq[1]-1)/(freq[1]-3)
          }
       }

   # matrixerg    <- numeric(6)
   matrixerg[1] <- S
   matrixerg[2] <- thetaS
   matrixerg[3] <- thetaT
   matrixerg[4] <- thetaFL
   matrixerg[5] <- thetaSA
   matrixerg[6] <- thetaTA

matrixerg2 <- c(matrixerg,sfreq)

return(matrixerg2) })
return(res)

}

#------------- FUNCTION: get_outgroup ----------------------------#

get_outgroup <- function(outgroup){

       freqo <- numeric(4)
       sfreq <- numeric(length(outgroup)+1)
       #modified
       freqo[2] <-sum(outgroup==0,na.rm=TRUE)
       freqo[3] <-sum(outgroup==1,na.rm=TRUE)
       freqo[4] <-sum(is.na(outgroup))
       freqo[1] <- freqo[2] + freqo[3]

      # Get ancestral
      if(freqo[1]){
         if((freqo[2]!=freqo[1]) && (freqo[2]!=0)){
          sfreq[freqo[2]+1]<-1 # 0 is the MUTATION - Number of values==0
          ancestral<- 2 # not included in test
         }else{
             if(freqo[2]==freqo[1]){ancestral<-0}
             if(freqo[3]==freqo[1]){ancestral<-1}
         }
      }else{ancestral<- 2 } # not included in test

return(c(ancestral,sfreq))

}

### ADDITIONAL STATISTIC FUNCTIONS

# Fu Li D*
fl_d <- function(sample_size,fr1,S,coef){

n <- sample_size
if(S==0 || coef[1] < 1.5){return(NaN)}

re <- fr1
an <- coef[1]
vd <- coef[9]
ud <- coef[10]

D <-  (S - an*re)/sqrt(ud*S +vd*S*S)

return(D)
}

#TAJIMA ##############################################

tajima_d <- function (k,S,coef){

 S_D <- NaN

if(S==0 || coef[1] < 1.51 ){return(NaN)}

an <- coef[1]
ut <- coef[4]
vt <- coef[3]

S_D <- (k-(S/an))/(sqrt((ut*S)+ (vt*S*(S))))

if(is.na(S_D)){return(NaN)}
if(abs(S_D) < 1.0e-15){return(NaN)}

return (S_D)

}


#########################################
# Fu and Li's D* */

 fl_d2 <- function(sample_size,fr1w,S,coef) # NO outgroup */
{


	 D2 <- NaN

	if(S == 0 || coef[1] < 1.51){return(NaN)}

	rs <- fr1w
	n <- sample_size
	an  <- coef[1]
	vd2 <- coef[5]
	ud2 <- coef[6]
	D2  <- (S/an - rs*((n-1)/n)) / sqrt(ud2*S + vd2*S*S)

	return(D2)
}

################################################

fl_f <- function(sample_size,fr1,S,piw,coef){   # amb outgroup

  n <- sample_size

	if(S == 0 || coef[1] < 1.5) return(NaN);

	re <- fr1;
	vf <- coef[11]
	uf <- coef[12]

	F  = (piw - re) / sqrt(uf*S + vf*S*S);

	return(F)
}


################################################

#/* Fu and Li's F* */

fl_f2 <- function(sample_size,fr1w,S,pii,coef) #/* NO outgroup */

{
 if(S == 0 || coef[1] < 1.51){return(NaN)}

	rs  <- fr1w  # Divpop hingegen macht as.integer()
	n   <- sample_size
	vf2 <- coef[7]
	uf2 <- coef[8]
	
	F2  <- (pii -(((n-1)/n)*rs)) / sqrt(uf2*S + vf2*S*S)
 
	return(F2)
}

#################################################
fay_wu_normalized2 <- function(n,thetaL,thetaw,S,coef,pi) # eq 11 and 12 from Fay and Wu H NORMALIZED (Zeng et al. Genetics 2006 174: 1431-9)*/
{

    if(pi == 0 || n < 2) {return(NaN)}
          an <- coef[1]
	  bn <- coef[2]

	varpiTL = thetaw * ((n-2))/(6*((n-1))) + S*(S-1)/(an*an+bn) / (9*(n*(n-1)*(n-1))) * (18*n*n*(3*n+2)*(bn+1/(n*n)) - (88*n*n*n + 9*n*n - 13*n + 6))

	H <- (pi - thetaL)/sqrt(varpiTL)

    return(H)
    
}
###################################################
E_zeng <- function (n,thetaL,thetaw,S,coef) # (eq 13 and 14 from Zeng et al. Genetics 2006 174: 1431-9)*/

{
    if(thetaw == 0 || n < 2) {return(NaN)}
    an <- coef[1]
   	bn <- coef[2]
	  varLW <- thetaw * (n/(2*(n-1)) - 1/an) + S*(S-1)/(an*an+bn) * (bn/(an*an) + 2*bn*(n/(n-1))*(n/(n-1)) - 2*(n*bn-n+1)/((n-1)*an) - (3*n+1)/((n-1)))

	E <- (thetaL - thetaw)/sqrt(varLW)

    return(E)
}
######################################################









