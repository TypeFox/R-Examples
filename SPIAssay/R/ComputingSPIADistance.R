####################################################################
#
# ComputingSPIADistance
#  Input: - two vectors of SNPs, one for each cell line
#           SNPs are coded as 0,1,2,NA (AA,BB,AB,NoCall)
#         - significative digits (default 5)
#         - verbose mode
#  Output: a row with informations of about the distance
#           1. SPIA distance  
#           2. number of valid calls
#           3. number of total calls
#           4. number of calls where one of the two SPNs are not available
#           5. number of calls where both SNPs are not available
#           6. number of calls where SNP change from {AA,BB} to AB or from AB to {AA,BB}
#           7. number of calls where SNP change from AA to BB or from BB to AA
#
####################################################################

ComputingSPIADistance<-function(vector1, vector2, defaultDigists = 5, verbose = FALSE){
  
  ## Control that the two vectors has the same length 
  if(length(vector1)!=length(vector2)){
    #cat("SPIA: error in function ComputingSPIADistance input vectors must have the same length\n");
    return(-1)
    }  
  
   vector1 <- factor(vector1)
   vector2 <- factor(vector2)
   levels(vector1) <- c(0,1,2)
   levels(vector2) <- c(0,1,2)
  
  ##Build contingency table of vector1 and vector2
  SummaryTable<-table(vector1,vector2,exclude=NULL);
  
  lenVet<-length(vector1);

  ##COMPUTATION START
  if (verbose) {message("SPIA: computing distance:\n")}
       
  ##Both SNPs are not available
  counterBothNA <- SummaryTable[4,4]
  if (verbose) {message("SPIA: number of calls where both SNPs are not available:",counterBothNA,"\n") }
  
  ##One of the two SPNs are not available
  counterOneNA <- sum(SummaryTable[4,c(1:3)], SummaryTable[c(1:3),4])
  if (verbose) {message("SPIA: number of calls where one of the two SNPs are not available:",counterOneNA,"\n") }
  
  ##From AA to BB or from BB to AA
  counterDiffAvsBorBvsA <- SummaryTable["1","0"]+SummaryTable["0","1"]
  if (verbose) {message("SPIA: number of calls where AA becomes BB or BB becomes AA:",counterDiffAvsBorBvsA,"\n") }
  
  ##From {AA,BB} to AB or from AB to {AA,BB}
  counterDiffAorBvsABorviceversa <- SummaryTable["2","0"]+
                                    SummaryTable["2","1"]+
                                    SummaryTable["0","2"]+
                                    SummaryTable["1","2"]
  if (verbose) {message("SPIA: number of calls where AA or BB become AB or vice versa:",counterDiffAorBvsABorviceversa,"\n") }                                    

  ##From {AA,BB} to AB
  counterDiffABvsAorB <- SummaryTable["0","2"] + SummaryTable["1","2"]
  if (verbose) {message("SPIA: number of calls where AA or BB become AB:",counterDiffABvsAorB,"\n") }
  
  ##Both SNPs are Homozygous
  counterBothHomoz <- SummaryTable["0","0"] + SummaryTable["1","1"]
  if (verbose) {message("SPIA: number of calls where both SNPs are homozygous:",counterBothHomoz,"\n") }

  ##Both SNPs are Hetherozygous
  counterBothHeter<-SummaryTable["2","2"]
  if (verbose) {message("SPIA: number of calls where both SNPs are heterozygous:",counterBothHeter,"\n") }

  ##Computing SPIA Distance
  somma<-counterDiffAorBvsABorviceversa+counterDiffAvsBorBvsA
  counter<-counterDiffAorBvsABorviceversa+counterDiffAvsBorBvsA+counterBothHomoz+counterBothHeter;
  distance<-somma/counter;
  if (verbose) {message("SPIA: distance is ",distance,"\n") }

  ##Return    
  return(c(signif(distance,digits=defaultDigists),counter,lenVet,counterOneNA,counterBothNA,counterDiffAvsBorBvsA,counterDiffAorBvsABorviceversa,counterDiffABvsAorB,counterBothHomoz,counterBothHeter))
  
}
