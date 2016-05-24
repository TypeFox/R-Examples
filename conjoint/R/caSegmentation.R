caSegmentation<-function(y,x,c=3)
{
	options(contrasts=c("contr.sum","contr.poly"))
	outdec<-options(OutDec="."); on.exit(options(outdec))
	options(OutDec=",")
	y<-m2v(y)
	Usi<-caTotalUtilities(y,x)
	segm<-kmeans(Usi,Usi[initial.Centers(Usi,c),])
	return(segm)
}