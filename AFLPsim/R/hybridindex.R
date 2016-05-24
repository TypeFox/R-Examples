hybridindex <-
function(data){
introdata<-sim2introgress(data)
dataf<-prepare.data(admix.gen=introdata$admix.gen, loci.data=introdata$loci,parental1=introdata$parental1, parental2=introdata$parental2,pop.id=F, ind.id=F)
indexh<-est.h(introgress.data=dataf, loci.data=introdata$loci, fixed=FALSE)
indexh}
