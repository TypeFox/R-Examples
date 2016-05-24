permute <- function (obj,iter){

# obj contains an obeject of class data [[1]] and stat[[2]]
data <- obj[[1]]
stat <- obj[[2]]

populations <- data@populations
npops       <- length(populations)
popnames    <- paste("pop",1:npops)
outgroup    <- data@outgroup
matrix_pol  <- data@matrix_pol
#---------------------------------------------------------#

nsamtot    <- dim(matrix_pol)[1]
length_seq <- dim(matrix_pol)[2]

### Create permutation indices ####
### -------------------------------
###################################

perm <- matrix(,iter,nsamtot)
for(xx in 1:iter){
   perm[xx,] <- sample(1:nsamtot) # permute all gensequences
}

### Create permutationsmatrices of matrix_pol
matrix_perm <- vector("list",iter)
for(xx in 1:iter){
  matrix_perm[[xx]] <- matrix_pol[perm[xx,],]
}

# Calculate statistics for permutation matrix_pol
#-------------------------------------------------#

#INIT#
erg           <- vector("list",iter) # store the permutation results
perm_fstall   <- rep(NA,iter)
perm_fsthall  <- perm_fstall
perm_fst1all  <- matrix(NA,npops,iter)
perm_fsth1all <- perm_fst1all
rownames(perm_fst1all) <- popnames
rownames(perm_fsth1all)<- popnames
perm_gst     <- rep(NA,iter)


### Count, when perm values better than original values
#------------------------------------------------------
iall         <- 0  # FSTALL  Nucleotides
niteriall    <- 0  # FSTALL  Nucleotides
ihall        <- 0  # FSTHALL Haplotypes
niterihall   <- 0  # FSTHALL Haplotypes
ighall       <- 0  # GstALL
niterighall  <- 0  # GstALL
                         
i1           <- matrix(,npops,iter)  # fst1all
niteri1      <- matrix(,npops,iter)  # fst1all
ih1          <- matrix(,npops,iter)  # fsth1all
niterih1     <- matrix(,npops,iter)  # fsth1all

### The definition of populations is the same
for(xx in 1:iter){

    #-- fstcalc ---- FUNCTION
    erg[[xx]]          <- fstcalc(matrix_perm[[xx]],populations,outgroup,data) # FST for nucleotides
    #--- fstall
    perm_fstall[xx]    <- erg[[xx]]$FSTALL  # Fst
    # print(perm_fstall[xx]);print(stat@FSTALL)
    if(perm_fstall[xx]!=-1000 & stat@FSTALL!=-1000){if(stat@FSTALL<=perm_fstall[xx]){iall<-iall+1} ;niteriall <- niteriall + 1}
    
    #---fst1all 
    # print(perm_fst1all);print(erg[[xx]]$FST1ALL)
    perm_fst1all[,xx]  <- erg[[xx]]$FST1ALL  # pop against rest nucleotides
    
    id1                <- which((perm_fst1all[,xx]!=-1000 & stat@FST1ALL!=-1000) & stat@FST1ALL<=perm_fst1all[,xx])
    i1[id1,xx]         <- 1 # perm is better
    id2                <- which((perm_fst1all[,xx]!=-1000 & stat@FST1ALL!=-1000))
    niteri1[id2,xx]    <- 1 # Count the compares    

    # - calc_hwhafsth ---- FUNCTION   
    erg[[xx]]          <- calc_hwhafsth(matrix_perm[[xx]],populations,outgroup)# FST for haplotypes
    #--- Gst
    perm_gst[xx]       <- erg[[xx]]$GstAll   # GST 
    if(perm_gst[xx]!="NaN" & stat@GstAll!="NaN"){if(stat@GstAll<=perm_gst[xx]){ighall<-ighall+1}; niterighall <- niterighall + 1}    
    #---fsthall
    perm_fsthall[xx]   <- erg[[xx]]$fsthALL  # Fst for haplotypes
    if(perm_fsthall[xx]!="NaN" & stat@FSTHALL!="NaN"){if(stat@FSTHALL<=perm_fsthall[xx]){ihall <- ihall+1};niterihall <- niterihall +1}
        
    #--- fsth1all
    perm_fsth1all[,xx] <- erg[[xx]]$fsth1all # Fst pop against rest haplotypes
    id1                <- which((perm_fsth1all[,xx]!="NaN" & stat@FSTH1ALL!="NaN") & stat@FSTH1ALL<=perm_fsth1all[,xx])
    ih1[id1,xx]        <- 1 # perm is better
    id2                <- which((perm_fsth1all[,xx]!="NaN" & stat@FSTH1ALL!="NaN"))
    niterih1[id2,xx]   <- 1 # Count the compares
    print("perm_fsth1all")
    print(perm_fsth1all)
    print("FSTH1ALL (original)")
    print(stat@FSTH1ALL)
    print(ih1)


}

# Calculate P-Values
p_fstall  <- iall/niteriall         # P-Value FSTALL
p_gstall  <- ighall/niterighall     # P-Value GSTALL

p_fst1all   <- rep(0,npops)         # P-Value FST1ALL
p_fsth1all  <- rep(0,npops)         # P-Value FSTH1ALL

for(xx in 1:npops){
 p_fst1all  [xx] <- sum(i1[xx,],na.rm=TRUE)/sum(niteri1[xx,],na.rm=TRUE)  
 p_fsth1all [xx] <- sum(ih1[xx,],na.rm=TRUE)/sum(niterih1[xx,],na.rm=TRUE)
}


# --------------------------------------------------- #
# # permute pairs of pops #
# (die Populationen untereinander werden permutiert !?)
# --------------------------------------------------- #

 erg2 <- vector("list",iter)
 perm_fstpair <- matrix(,npops,npops)
 rownames(perm_fstpair) <- popnames
 colnames(perm_fstpair) <- popnames

 niteri    <- 0 # FSTPAIR
 i            <- matrix(,npops,npops) # FSTPAIR
 iall         <- matrix(,npops,npops) # FSTPAIR
 
 icount    <- 0
 niterih   <- 0 #FSTH
 ih           <- matrix(,npops,npops) # FSTH
 ihall        <- matrix(,npops,npops) # FSTH
 ihcount   <- 0
 niterigh  <- 0 # Gst
 igh          <- matrix(,npops,npops) # Gst
 ighall       <- matrix(,npops,npops)
 ighcount  <- 0
 
 poppairs <- combn(npops,2)
  
 res <- apply(poppairs,2,function(x){
     pop1           <- unlist(populations[[x[1]]])
     pop2           <- unlist(populations[[x[2]]])
     pop12          <- c(pop1,pop2)      
     pop1_id        <- x[1]
     pop2_id        <- x[2]
     
     
    for(xx in 1:iter){
        
       perm     <- sample(pop12)
       newpops  <- list(perm[1:length(pop1)],perm[(length(pop1)+1):length(pop12)]) # Permute the sequences in the popualtions
       
       ## --------------------------
       ## Make Tests for Perm Matrix
       ## FSTPAIR Nucleotides 
       #############################
       
       erg2         <- fstcalc(matrix_pol,newpops,outgroup=FALSE,data) # FST for nucleotides # Only the populations 
       perm_fstpair <- erg2$FSTALL
          
       if((stat@FSTPAIR[pop2_id,pop1_id]!=-1000 & perm_fstpair!=-1000) & (stat@FSTPAIR[pop2_id,pop1_id] <= perm_fstpair)){icount <- icount + 1}        
       if((stat@FSTPAIR[pop2_id,pop1_id]!=-1000 & perm_fstpair!=-1000)){niteri <- niteri + 1}
       
       # calc FSTPAIR Haplotypes ----------------------------------------------------------------------------------
       erg2          <- calc_hwhafsth(matrix_pol,newpops,outgroup=FALSE)
       perm_fsthpair <- erg2$fsthALL
       
       if((stat@FSTHMATRIX[pop2_id,pop1_id]!="NaN" & perm_fsthpair!="NaN") & (stat@FSTHMATRIX[pop2_id,pop1_id]<=perm_fsthpair)){ihcount <- ihcount + 1}
       if((stat@FSTHMATRIX[pop2_id,pop1_id]!="NaN" & perm_fsthpair!="NaN")){niterih <- niterih + 1}
        
        
        # Gst ------------------------------------------------------------------------------------------------------
        perm_Gst <- erg2$GstAll
       
       if((stat@GstMATRIX[pop2_id,pop1_id]!="NaN" & perm_Gst!="NaN") & (stat@GstMATRIX[pop2_id,pop1_id]<=perm_Gst)){ighcount <- ighcount + 1}
       if((stat@GstMATRIX[pop2_id,pop1_id]!="NaN" & perm_Gst!="NaN")){niterigh <- niterigh + 1}
    
    }
       i  [x[2],x[1]]   <<- icount      # perm was better # FSTPAIR
       ih [x[2],x[1]]   <<- ihcount     # perm was better # FSTHPAIR
       igh[x[2],x[1]]   <<- ighcount    # perm was better # GSTPAIR
       
       iall  [x[2],x[1]] <<- niteri      # all compares   # FSTPAIR
       ihall [x[2],x[1]] <<- niterih     # all compares   # FSTHPAIR
       ighall[x[2],x[1]] <<- niterigh    # all compares   # GSTPAIR
 })
 
 # Generate P-values
 #-----------------------#
 p_fstpair  <- i/iall      # P-Value FSTPAIR
 p_fsthpair <- ih/ihall    # P-Value FSTHPAIR
 p_gstpair  <- igh/ighall  # P-Vaue  GSTPAIR
 
 return(list(p_fstall=p_fstall,p_gstall=p_gstall,p_fst1all=p_fst1all,p_fsth1all=p_fsth1all,p_fstpair=p_fstpair,p_fsthpair=p_fsthpair,p_gstpair=p_gstpair))
}
