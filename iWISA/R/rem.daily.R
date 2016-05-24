rem.daily <-
function(data1){
#==================================================================
n.station<-dim(data1)[2] # each column of the data matrix are records for each individual station
N<-dim(data1)[1]         # number of observations for each station
# matrix to store results
recon<-matrix(0, ncol = n.station, nrow = N)
#
# vector: average of all stations' records
mean.v<-numeric(N)
if(n.station==1) mean.v<-numeric(N) else mean.v<-apply(data1, 1, mean)
#apply(X, MARGIN, FUN, ...)	
#MARGIN: which the function will be applied over. For a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. 
#
# prepare to remove mean
data1.centr<-matrix(NA, ncol=n.station, nrow = N)
#
# remove pointwise mean of all stations
for(i in 1:n.station){data1.centr[,i]<-data1[,i] - mean.v}
#   
# remove periodic component
deco<-pca.sq(data1.centr)
#
# Pseudo-Sq component
SQ <- deco$sq
for(i in 1:n.station) recon[,i] <- deco$diff[,i] + mean.v
#
list(recon = recon, SQ = SQ)
}
