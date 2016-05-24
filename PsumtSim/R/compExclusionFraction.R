#
# Function to compute fraction of cells with a response to category which are
# excluded by pre-testing, using Poisson simulation of counts in background 
# and response intervals.
#
compExclusionFraction<-function(bkg,resps,numTrialsPerCat,
																 pretestP,anovaP,
																 showProgress=FALSE, numCells=1000) {
	numCats<-length(resps)
	
  totPreselected<-0
	totCatSelective<-0
	totCatSelectivePreselected<-0
	
	totTrials<-numCats*numTrialsPerCat
	catLabels<-factor(rep(1:numCats,numTrialsPerCat))
	
	for (i in 1:numCells) {
	  if (showProgress) print(i)
	  
		bkgCnts<-rpois(totTrials,bkg)
		respCnts<-rpois(totTrials,resps)
		
		# check if any category is significant, and so cell is pre-selected
		presel<-FALSE
		for (j in 1:numCats) {
			catSel<-catLabels==levels(catLabels)[j]
			pVal<-t.test(bkgCnts,respCnts[catSel],alternative='two.sided')$p.value
			if (pVal<pretestP) {
				presel<-TRUE
				break
			}
		}
		
		# check for effect of category
		mod<-glm(respCnts~catLabels)
		pVal<-anova(mod,test='F')$"Pr(>F)"[2]
		
		# count results
		if (pVal<anovaP) totCatSelective<-totCatSelective+1
		if (presel) {
			totPreselected<-totPreselected+1
			if (pVal<anovaP) {
			  totCatSelectivePreselected<-totCatSelectivePreselected+1
			}
		}
	}
  
  list( exclusionFrac=1-(totPreselected/numCells),
        catSelectiveFrac=totCatSelective/numCells,
        catSelExclFrac=1-(totCatSelectivePreselected/totCatSelective))
}
