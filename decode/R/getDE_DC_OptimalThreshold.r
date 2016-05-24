#####################################################
#' Perform chi-square optimization
#' @param t_result The t-statistics
#' @param MaxGene Number of genes in expression data
#' @param d_r DC measures
#' @param minSupport The minimum expected frequency in contingency table
#' @return The optimal threshold information
#' @export
###################################################### 
# DECODE (c) Copyright 2014 by The Hong Kong Polytechnic University, Department of Health Technology and Informatics 
# Written by Thomas Lui
# Permission is granted to copy and use this program provided no fee is
# charged for it and provided that this copyright notice is not removed.
##########################################################################
getDE_DC_OptimalThreshold = function(t_result, MaxGene, d_r, minSupport) {

  optimal=data.frame(chi2Value=double(), pValue=double(), tCutOff=double(), rCutOff=double(), 
                     obsA=double(), obsB=double(), obsC=double(), obsD=double(), 
					 expA=double(), expB=double(), expC=double(), expD=double())
					 
  rank_abs_t = rank(t_result$absTScore, ties.method= "min")
  # for every gene i
  for (i in 1:MaxGene) {
	optimal[i,"pValue"]=1
    optimal[i,"chi2Value"]=0
	optimal[i,"tCutOff"]=0
	optimal[i,"rCutOff"]=0
	optimal[i,"obsA"]=0
	optimal[i,"obsB"]=0
	optimal[i,"obsC"]=0
	optimal[i,"obsD"]=0
	optimal[i,"expA"]=0
	optimal[i,"expB"]=0
	optimal[i,"expC"]=0
	optimal[i,"expD"]=0
    rank_one_d_r = rank(d_r[i,], ties.method= "min")

	
    # consider each threshold candidate
    for (j in 1:MaxGene) {

      currentTcutoff = t_result$absTScore[j]
      currentRcutoff = d_r[i,j]
	  
      tempA = sum(rank_abs_t >=rank_abs_t[j] & rank_one_d_r >=rank_one_d_r[j])
	  tempB = MaxGene-rank_one_d_r[j]+1-tempA
	  tempC = MaxGene-rank_abs_t[j]+1-tempA
	  tempD = MaxGene-tempA -tempB-tempC

      expectedA = (tempA+tempB)*(tempA+tempC)/ MaxGene;
      expectedB = (tempA+tempB)*(tempB+tempD)/ MaxGene;
      expectedC = (tempA+tempC)*(tempC+tempD)/ MaxGene;
      expectedD = (tempD+tempB)*(tempD+tempC)/ MaxGene;

	  if (expectedA<minSupport || expectedB<minSupport|| expectedC<minSupport|| expectedD<minSupport) {
	    next
	  }	
	  expectedValues = c(expectedA, expectedB, expectedC, expectedD)

      # # # # # # # # # # # # # # # # # # # # # # # # # # 
      # chi-square
      # # # # # # # # # # # # # # # # # # # # # # # # # # 
      observedValues = matrix(c(tempA, tempB, tempC, tempD),nrow=2,ncol=2)
      chi2stat = sum((observedValues-expectedValues)^2 / expectedValues)

	  # record this data point if it is min p 
      if (chi2stat> optimal[i,"chi2Value"]) {
        pValue = pchisq(chi2stat,1,lower.tail=FALSE)
	    optimal[i,"pValue"]=pValue
	    optimal[i,"chi2Value"]=chi2stat
	    optimal[i,"tCutOff"]=currentTcutoff
	    optimal[i,"rCutOff"]=currentRcutoff
	    optimal[i,"obsA"]=tempA
	    optimal[i,"obsB"]=tempB
	    optimal[i,"obsC"]=tempC
	    optimal[i,"obsD"]=tempD
	    optimal[i,"expA"]=expectedA
	    optimal[i,"expB"]=expectedB
	    optimal[i,"expC"]=expectedC
	    optimal[i,"expD"]=expectedD
      }	  
    }
	if (i %%1==0) {
      print(sprintf("Gene id: %d",i))
	  #print(sprintf("%d   %f",i,optimal[i,"chi2Value"]))
    }
  }

  adjustedPValues = getBonferroniPValue(optimal$pValue)  
  adjustedPValues = getFDR(adjustedPValues)  
  optimal$pValue =adjustedPValues

  return(optimal)
}