ShowAllSimulations<-function(sym,y,x)
{
 	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	MaxUtility<-caMaxUtility(sym,y,x)
	BTLmodel<-caBTL(sym,y,x)
	LogitModel<-caLogit(sym,y,x)
	TotalUtility<-totalsimutility(sym,y,x)
	allsimul<-cbind(TotalUtility,MaxUtility,BTLmodel,LogitModel)
	print(round(allsimul,2))
}