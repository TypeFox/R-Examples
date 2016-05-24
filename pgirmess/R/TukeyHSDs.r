
TukeyHSDs<-function(TukeyHSD.object){
# Patrick Giraudoux 19.1.2004
# Simplify the list of TuckeyHSD objet keeping
# only the significant differences.
res<-TukeyHSD.object
x<-TukeyHSD.object[[1]]
suppressWarnings(y<-!is.nan(sqrt(x[,2]*x[,3])))
res[[1]]<-x[y,,drop=FALSE]
res
}
