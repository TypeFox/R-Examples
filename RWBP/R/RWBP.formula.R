RWBP.formula <-
function(formula,data,...,nn_k=10,min.clusters=8,clusters.iterations=6,clusters.stepSize=2,alfa=0.5,dumping.factor=0.9){
m<- match.call(expand.dots=FALSE)
m$...<-NULL
m[[1]]<- as.name("model.frame")
m<- eval(m,parent.frame())
Terms<- attr(m,"terms")
Y<- model.extract(m,"response")
X<- m[,-attr(Terms,"response"),drop=FALSE]
RWBP(X,nn_k,min.clusters,clusters.iterations,alfa,dumping.factor,...)
}
