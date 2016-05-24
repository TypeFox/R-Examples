jointfreqdist <- function(matrix_pol,populations,outgroup=FALSE){


npops    <- length(populations)
popnames <- paste("pop",1:npops)

if(outgroup[1]!=FALSE){    ######################### ---------- OUTGROUP
     outgmatrix <- matrix_pol[outgroup,,drop=FALSE]
   if(length(outgroup)>1){                      # Outgroup is a matrix
      erg <- apply(outgmatrix,2,get_monomorph)
   }else{                                       # Outgroup is a vector
      erg <- sapply(outgmatrix,get_monomorph)   
   }
   
 outgmonoid       <- !is.na(erg)
 MATRIX           <- matrix_pol[,outgmonoid] # the matrix where the outgroup is monomorph  
 cols             <- which(outgmonoid)       # the ids of cols where the outgroup is monomorph
 anc              <- erg[outgmonoid]         # The outgroup monomorph values

# Calculate the Frequencies of each population 
 
 # Init ---------------------------
 jfd           <- matrix(NaN,npops,length(anc))
 rownames(jfd) <- popnames
 colnames(jfd) <- cols
 #-----------------------------------

 for(xx in 1:npops){

#  if(length(populations)==0){next;}

  m        <- MATRIX[populations[[xx]],,drop=FALSE]
  m_anc    <- rbind(m,anc)
  fr       <- apply(m_anc,2,getFREQanc)
  jfd[xx,] <- fr  

 }


}# End of if outgroup[1]!=FALSE

if(outgroup[1]==FALSE){ ##########---------------------- NO OUTGROUP

# Calculate the Frequencies of each population 
 
 # INIT ------------------------------------------
  jfd           <- matrix(NaN,npops,dim(matrix_pol)[2])
  rownames(jfd) <- popnames
 # -----------------------------------------------
 
for(xx in 1:npops){

  # if(length(populations[[xx]])==0){next;}
   
    m              <- matrix_pol[populations[[xx]],,drop=FALSE]
    jfd[xx,]       <- apply (m,2,function(x){
                       eins <- sum(x==1,na.rm=TRUE)
                       null <- sum(x==0,na.rm=TRUE)
                       if(eins<=null){
		       f     <- eins/(null+eins)
                       }else{
 		       f     <- null/(eins+null)        	
                       }
                      return(f)
                      })
                           
}

#eins      <- colSums(m==1,na.rm=TRUE)
#null      <- colSums(m==0,na.rm=TRUE)

}# End of if outgroup[1]==FALSE


return(list(jfd=jfd))

}# End of Function


###########################################################
# SUBFUNCTIONS
###########################################################

getFREQanc <- function(m_anc){
# last element is the outgroup value
anc <- m_anc[length(m_anc)]
vek <- m_anc[1:(length(m_anc)-1)]
gapids <- !is.na(vek)            # delete gaps
vek    <- vek[gapids]            # delete gaps
mutations    <- sum(vek!=anc)
nonmutations <- sum(vek==anc) + mutations
f <- mutations/nonmutations
return(f)

}

get_monomorph <- function(vek){
  gapids <- !is.na(vek)
  vek    <- vek[gapids]
  ss     <- unique(vek) 
  lenss  <- length(ss) 
  if(lenss==1) {return(as.numeric(ss))}  #  monomorph ss is the monomorph value
  if(lenss==2) {return(NaN)}             #  biallelic  
  if(lenss==0) {return(NaN)}             #  all are gaps
}
