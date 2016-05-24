computeCombinedScores <-
function (ScoresSTRING,hscores) 
{
	prior=0.063 
	priorscores=ScoresSTRING-prior 
	priorscores[which(priorscores<0)]=0 
	priorscores=priorscores/(1-prior) 
	priorscores[which(priorscores>1)]=1 
	priorscores[,which(colnames(ScoresSTRING)=="cooccurence")]=priorscores[,which(colnames(ScoresSTRING)=="cooccurence")]*(1-hscores) 
	priorscores[,which(colnames(ScoresSTRING)=="textmining")]=priorscores[,which(colnames(ScoresSTRING)=="textmining")]*(1-hscores)
	CombinedScores=1-apply(1-priorscores,1,prod) 
	CombinedScores=CombinedScores*(1-prior) 
	if (("cooccurence" %in% colnames(ScoresSTRING)) | ("textmining" %in% colnames(ScoresSTRING)))
	{
		CombinedScores=CombinedScores+prior 
	}
	return(CombinedScores)
}
