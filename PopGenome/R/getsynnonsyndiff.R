###############################################################################
## GETSYNNONSYNDIFF - Return matrix  of Syn- Nonsyn- differences between codons

# Syntax: [SynDif,AsynDif] = getsynnonsyndiff(icode)
#
# Inputs:
#     icode      - index of genetic code

# Outputs:
#    SynDif      - 64x64 matrix of syn- differences between codons
#    AsynDif     - 64x64 matrix of nonsyn- differences between codons
##############################################################################

getsynnonsyndiff_my<-function(...){
  l <- list(...) # Put the arguments in a (named) list
  
  
  if(length(l)==0){
       icode<-1
      }
  else{
       icode<-l[[1]]
      }
      
  if(icode>13)
  return(print("icode should between 1 to 13"))

########################   codontable    ##################################################
codontable<-list()
codontable[1]<- "FFLLSSSSCCW*YY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #   'Standard' # geaendert !
codontable[2]<- "FFLLSSSSCCWWYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMMTTTTSS**NNKK-"  #   'Vertebrate Mitochondrial' # geaendert !
codontable[3]<- "FFLLSSSSCCWWYY**TTTTPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMMTTTTSSRRNNKK-"  #   'Yeast Mitochondrial' # geaendert !
codontable[4]<- "FFLLSSSSCCWWYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #   'Mold, Protozoan, and CoelenterateMitochondrial and Mycoplasma/Spiroplasma'# geaendert !
codontable[5]<- "FFLLSSSSCCWWYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMMTTTTSSSSNNKK-"  #   'Invertebrate Mitochondrial # geaendert !
codontable[6]<- "FFLLSSSSCCW*YYQQLLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #   'Ciliate, Dasycladacean and Hexamita Nuclear'# geaendert !
codontable[7]<- "FFLLSSSSCCWWYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSSSNNKN-"  #   'Echinoderm Mitochondrial' # geaendert !
codontable[8]<- "FFLLSSSSCCWCYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #   'Euplotid Nuclear'# geaendert !
codontable[9]<- "FFLLSSSSCCW*YY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #    '"Bacterial"' # geaendert !
codontable[10]<-"FFLLSSSSCCW*YY**LLSLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #    'Alternative Yeast Nuclear' # geaendert !
codontable[11]<-"FFLLSSSSCCWWYY**LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMMTTTTSSGGNNKK-"  #    'Ascidian Mitochondrial' # geaendert !
codontable[12]<-"FFLLSSSSCCWWYY*YLLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSSSNNKN-"  #    'Flatworm Mitochondrial' # geaendert !					
codontable[13]<-"FFLLSSSSCCW*YYQ*LLLLPPPPRRRRHHQQVVVVAAAAGGGGDDEEIIMITTTTSSRRNNKK-"  #    'Blepharisma Nuclear' # geaendert !

######################## CODON  and   TABLE   ,syndDif  and  AsynDif  ##################################################
 s1<-c(rep(1,16),rep(2,16),rep(3,16),rep(4,16))
 s2<-c(rep(1,4),rep(2,4),rep(3,4),rep(4,4))
 s3<-c(1,2,3,4)
 s<-c(s1,rep(s2,4),rep(s3,16))
CODON<- matrix(s,ncol=3)
TABLE<-codontable

SynDif<-matrix(rep(0,times=4096),ncol=64)
AsynDif<-matrix(rep(0,times=4096),ncol=64)

########################## find the position for stopcodon "*"  in a icode #######################

stops=vector()
for(i in 1:64){
    if(substr(TABLE[[icode]],i,i)=="*") {
         stops<-c(stops,i)
        }
    }
    

#source("i_countSynNonsynDiff.r")    
##########################################################################################
for(i in 1:64)
for(j in 1:64){
    Diff<-i_countSynNonsynDiff(i,j,stops,CODON,icode,TABLE)
    SynDif[j,i]<-Diff$SynDif
    AsynDif[j,i]<-Diff$AsynDif
    }
#############################################################################
    m<-list(SynDif=SynDif,AsynDif=AsynDif)
    return(m)



}



#########################
#####  sub function #####
#########################
################################################################################
#######
#######
#############################################################################

i_countSynNonsynDiff<-function(codon1,codon2,stops,CODON,icode,TABLE)  {

SynDif<-0
AsynDif<-0

###################   ## there is stop codon in codon1 oder/and  codon2 ##################################
if(!is.na(match(codon1,stops))||!is.na(match(codon2,stops))){
    SynDif<--1
    AsynDif<--1
    Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
    return(Diff)
    }
    
##############################################################################

ndiff<-3-length(which(CODON[codon1,]==CODON[codon2,]))

#####################################################################
if(ndiff==0){
    SynDif<-0
    AsynDif<-0
    Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
    return(Diff)
    }
##########################################################
if(ndiff==1) {
    if(substr(TABLE[[icode]],codon1,codon1)==substr(TABLE[[icode]],codon2,codon2)){
        SynDif<-1
        AsynDif<-0
        }
    else{
         SynDif<-0
        AsynDif<-1
        }
    Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
    return(Diff)
    }
########################################################
if(ndiff==2){
    rep_codon1<-matrix(rep(as.vector(CODON[codon1,]),64), nrow=64,byrow=TRUE)
    rep_codon2<-matrix(rep(as.vector(CODON[codon2,]),64), nrow=64,byrow=TRUE)

    A1<-matrix(rep(NA,times=64*3),ncol=3)
    for (i in 1:64)
    for (j in 1:3){
        if(CODON[i,j]==rep_codon1[i,j])
        A1[i,j]<-0
        else
        A1[i,j]<-1
        }

    A2<-matrix(rep(NA,times=64*3),ncol=3)
    for (i in 1:64)
    for (j in 1:3){
        if(CODON[i,j]==rep_codon2[i,j])
        A2[i,j]<-0
        else
        A2[i,j]<-1
        }
        
   A12<-vector()
   A1_position<-vector()  ###the positon with sum(A1[i,])=1
   A2_position<-vector()
   for(i in 1:64){
        if(sum(A1[i,])==1){
            A1_position <-c(A1_position,i)
        }
        if(sum(A2[i,])==1){
            A2_position <-c(A2_position,i)
        }
    }
    
    for(i in 1:length(A1_position))
    for(j in 1:length(A2_position)){
        if(A1_position[i]==A2_position[j])
            A12<-c(A12,A1_position[i])
        }
  
   # for(i in 1:length(A1_position)){
        #if(!is.na(match(A1_position[i],A2_position)))
            #A12<-c(A12,A1_position[i])
       # }
    
    # setdiff(A12,stops)
    flag<-vector()
    for(i in 1:length(A12))
    for(j in 1:length(stops)){
        if(A12[i]==stops[j]){
           flag<-c(flag,i)
           break
        }
    }
    if(length(flag)>0)
        A12<-A12[-flag]
    

#############    [n,m]  #############
    if(is.vector(A12))              #
       m<-length(A12)               #
    if(is.matrix(A12))              #
       m<-dim(A12)[2]               #
#####################################
      if(m==0){      
          SynDif<-0
          AsynDif<-2
          Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
          return(Diff)
          }

      for(i in 1:m){
          if(substr(TABLE[[icode]],codon1,codon1)==substr(TABLE[[icode]],A12[i],A12[i]))
              SynDif<-SynDif+1
          else
              AsynDif<-AsynDif+1
          if(substr(TABLE[[icode]],codon2,codon2)==substr(TABLE[[icode]],A12[i],A12[i]))
              SynDif<-SynDif+1
          else
              AsynDif<-AsynDif+1
          }
      SynDif<-SynDif/m
     	AsynDif<-AsynDif/m

      Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
      return(Diff)
      }
###############################################################
if(ndiff==3){
    rep_codon1<-matrix(rep(as.vector(CODON[codon1,]),64), nrow=64,byrow=TRUE)
    rep_codon2<-matrix(rep(as.vector(CODON[codon2,]),64), nrow=64,byrow=TRUE)

    A1<-matrix(rep(NA,times=64*3),ncol=3)
    for (i in 1:64)
    for (j in 1:3){
        if(CODON[i,j]==rep_codon1[i,j])
        A1[i,j]<-0
        else
        A1[i,j]<-1
        }

    A2<-matrix(rep(NA,times=64*3),ncol=3)
    for (i in 1:64)
    for (j in 1:3){
        if(CODON[i,j]==rep_codon2[i,j])
        A2[i,j]<-0
        else
        A2[i,j]<-1
        }

    X1<-vector()  ###the positon with sum(A1[i,])=1
    X2<-vector()  ###the positon with sum(A2[i,])=1
    for(i in 1:64){
        if(sum(A1[i,])==1)
            X1 <-c(X1,i)
        }
        
    for(i in 1:64){
        if(sum(A2[i,])==1)
            X2 <-c(X2,i)
        }
    


    Y1<-vector()  ###the positon with sum(A1[i,])=2
    Y2<-vector()  ###the positon with sum(A2[i,])=2
    for(i in 1:64){
        if(sum(A1[i,])==2)
            Y1 <-c(Y1,i)
        }
    for(i in 1:64){
        if(sum(A2[i,])==2)
            Y2 <-c(Y2,i)
        }
    
    
    Z1<-vector()  # Z1 = intersect(X1,Y2);
    Z2<-vector()  # Z2 = intersect(X2,Y1);
    
    for(i in 1:length(X1))
    for(j in 1:length(Y2)){
        if(X1[i]==Y2[j])
            Z1<-c(Z1,X1[i])
        }
    
    for(i in 1:length(X2))
    for(j in 1:length(Y1)){
        if(X2[i]==Y1[j])
            Z2<-c(Z2,X2[i])
        }
        
    
########################################################

    Z1<-c(Z1[1],Z1[1],Z1[2],Z1[2],Z1[3],Z1[3])
    Z2<-c(Z2[1],Z2[2],Z2[1],Z2[3],Z2[2],Z2[3])
    SixPath<-matrix(c(rep(codon1,times=6),Z1,Z2,rep(codon2,times=6)),ncol=6,byrow=TRUE)
    
    
    SixPath_index<-matrix(rep(1,times=dim(SixPath)[1]*dim(SixPath)[2]),nrow=dim(SixPath)[1])
   
    for(i in 1:dim(SixPath)[1])
    for(j in 1:dim(SixPath)[2]){
        for(k in 1:length(stops)){
            if(SixPath[i,j]==stops[k]){ 
                SixPath_index[i,j]<-0
                break
                }
            }
        }


#############    [n,m]  #############################################
                                                                    #
                                                                    #
    picker<-vector()                                                #
    for(j in 1:dim(SixPath_index)[2]){                              #
        if(sum(SixPath_index[,j])==4)                               #
            picker<-c(picker,j)                                     #
        }
                                                                    #
    if(length(picker)>0){                                                                     
        new_SixPath<-vector()                                           #
        for(i in 1:length(picker)){                                     #
            new_SixPath<-cbind(new_SixPath,SixPath[,picker[i]])        #
            }                                                           #
                                                                    #
        if(is.vector(new_SixPath))                                      #
            m<-length(new_SixPath)                                       #
        if(is.matrix(new_SixPath))                                      #
            m<-dim(new_SixPath)[2]                                       #
          # n<-dim(new_SixPath)[1]
        }                                       #
        
    if(length(picker)==0) 
        m<-0
#####################################################################

    if(m==0){
        SynDif<-1
        AsynDif<-2
        Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
        return(Diff)
        }
        
    for(j in 1:m){
        #a1<-substr(TABLE[[icode]],new_SixPath[1,j],new_SixPath[1,j])
       # a2<-substr(TABLE[[icode]],new_SixPath[2,j],new_SixPath[2,j])
        #a3<-substr(TABLE[[icode]],new_SixPath[3,j],new_SixPath[3,j])
        #a4<-substr(TABLE[[icode]],new_SixPath[4,j],new_SixPath[4,j])
        
        if(substr(TABLE[[icode]],new_SixPath[1,j],new_SixPath[1,j])==substr(TABLE[[icode]],new_SixPath[2,j],new_SixPath[2,j]))
            SynDif<-SynDif+1 
        else
            AsynDif<-AsynDif+1
             
        if(substr(TABLE[[icode]],new_SixPath[2,j],new_SixPath[2,j])==substr(TABLE[[icode]],new_SixPath[3,j],new_SixPath[3,j]))
            SynDif<-SynDif+1
        else
            AsynDif<-AsynDif+1
             
        if(substr(TABLE[[icode]],new_SixPath[3,j],new_SixPath[3,j])==substr(TABLE[[icode]],new_SixPath[4,j],new_SixPath[4,j]))
            SynDif<-SynDif+1
        else
             AsynDif<-AsynDif+1

        }
        
    SynDif<-SynDif/m
    AsynDif<-AsynDif/m
    
    Diff<-list(SynDif=SynDif,AsynDif=AsynDif)
        return(Diff)
    }
    
##########################################################


}       #end for sub function


   

