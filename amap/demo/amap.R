###################################################
### chunk number 1: 
###################################################
cat("\n------  CLUSTERING TOOLS -------\n")
cat("\n------  Hierarchical clustering  -------\n")

data(USArrests)
h = hcluster(USArrests)
plot(h)
readline("Next")


###################################################
### chunk number 2: 
###################################################
cat("\n------  Hierarchical clustering using function heatmap -------\n")
heatmap(as.matrix(USArrests),
        hclustfun=hcluster,
        distfun=function(u){u})
readline("Next")


###################################################
### chunk number 3: 
###################################################
cat("\n------  Parralelized Hierarchical clustering  -------\n")
h = hclusterpar(USArrests,nbproc=4)   
readline("Next")


###################################################
### chunk number 4: 
###################################################
cat("\n------  K-means clustering  -------\n")
Kmeans(USArrests,centers=3,method="correlation")
readline("Next")


###################################################
### chunk number 5: 
###################################################
cat("\n------  ROBUST TOOLS  -------\n")
cat("\n------  A robust variance computation  -------\n")
data(lubisch)
lubisch <- lubisch[,-c(1,8)]
varrob(scale(lubisch),h=1)
readline("Next")


###################################################
### chunk number 6: 
###################################################
cat("\n------  A robust principal component analysis  -------\n")
p <- acpgen(lubisch,h1=1,h2=1/sqrt(2))
plot(p)
readline("Next")

###################################################
### chunk number 6: 
###################################################
cat("\n------  Another robust principal component analysis  -------\n")
p <- acprob(lubisch,h=4)
plot(p)
readline("Next")


