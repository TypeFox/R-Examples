## ------------------------------------------------------------------------
library(ClustVarLV)

## ------------------------------------------------------------------------
data(apples_sh)
# 43 sensory attributes of 12 varieties of apple from southern hemisphere
senso<-apples_sh$senso
# Scores of liking given fy 60 consumers for each of the 12 varieties of apple
pref<-apples_sh$pref

## ---- results="hide"-----------------------------------------------------
resclv_senso <- CLV(X = senso, method = "directional", sX = TRUE)
# option sX=TRUE means that each attribute will be auto-scaled (standard deviation =1)

# Print of the 'clv' object 
print(resclv_senso)

## ---- fig.width=3, fig.height=3------------------------------------------
# Dendrogram of the CLV hierarchical clustering algorithm :
plot(resclv_senso,"dendrogram")
# Graph of the variation of the clustering criterion
plot(resclv_senso,"delta")

## ------------------------------------------------------------------------
# Summary the CLV results for a partition into 4 groups
summary(resclv_senso,K=4)

## ---- fig.width=4, fig.height=4------------------------------------------
# Representation of the group membership for a partition into 4 groups
plot_var(resclv_senso,K=4,label=T,cex.lab=0.8)

## ---- fig.width=6, fig.height=6------------------------------------------
plot_var(resclv_senso,K=4,beside=T)

## ---- results="hide"-----------------------------------------------------
# Extract the group membership of each variable
get_partition(resclv_senso,K=4,type="vector")
# or 
get_partition(resclv_senso,K=4,type="matrix")

# Extract the group latent variables 
get_comp(resclv_senso,K=4)

## ---- results="hide"-----------------------------------------------------
res.segext<- CLV(X = pref, Xr = senso, method = "local", sX=TRUE, sXr = TRUE)

print(res.segext)

## ---- fig.width=3, fig.height=3------------------------------------------
plot(res.segext,"dendrogram")
plot(res.segext,"delta") 

## ------------------------------------------------------------------------
table(get_partition(res.segext,K=2),get_partition(res.segext,K=3))

## ------------------------------------------------------------------------
get_load(res.segext,K=3)

## ------------------------------------------------------------------------
res.clvkm.rd<-CLV_kmeans(X = pref, Xr = senso, method = "local", sX=TRUE,
                         sXr = TRUE, clust=3, nstart=100)

## ------------------------------------------------------------------------
res.clvkm.hc<-CLV_kmeans(X = pref, Xr = senso, method = "local", sX=TRUE,
                        sXr = TRUE, clust=res.segext[[3]]$clusters[1,])

## ------------------------------------------------------------------------
table(get_partition(res.segext,K=3),get_partition(res.clvkm.hc,K=3)) 

## ------------------------------------------------------------------------
table(get_partition(res.segext,K=3),get_partition(res.clvkm.rd,K=3)) 

