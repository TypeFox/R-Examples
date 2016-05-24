# 
# Function to compute excess rejection fraction of a sequential pre-test by category
# followed by an ANOVA for an effect of category using Poisson distributed counts.
#
compRejectionFraction<-function(bkgLevel,respLevel,numCats,pretestP,anovaP,
																 showProgress=FALSE,
																 numTrialsPerCat=10,numCells=1000) {

	totPreselected<-0
	totRejected<-0
	
	totTrials<-numCats*numTrialsPerCat
	catLabels<-factor(rep(1:numCats,numTrialsPerCat))
	
	for (i in 1:numCells) {
	  if (showProgress) print(i)
	  
		bkgCnts<-rpois(totTrials,bkgLevel)
		respCnts<-rpois(totTrials,respLevel)
		
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
		
		# check for effect of category if cell is pre-selected
		if (presel) {
			totPreselected<-totPreselected+1
			mod<-glm(respCnts~catLabels)
			pVal<-anova(mod,test='F')$"Pr(>F)"[2]
			if (pVal<anovaP) totRejected<-totRejected+1
		}
	}
  
  list( exclusionFrac=1-(totPreselected/numCells),
        rejectionFrac=totRejected/totPreselected )
}