#'@import fda splines
#'@importFrom stats quantile median 

pca.sq <-
function(data) #replace the |largest| 5% of scores by median
{
period = 1440
n.station<-dim(data)[2] #each row of the data matrix are records for each individual station
N<-dim(data)[1]         # number of observations for each station
# how many subseries of length=period, here how many days
how.many<-floor(N/period)
# z is a temporary matrix
z<-matrix(data = NA, ncol = n.station, nrow = (how.many*period))
for(i in 1:n.station){z[, i]<-data[1:(how.many*period), i]}
# create an array how.many x period, for all stations
z.new<-array(NA, dim=c( period, how.many,n.station))
for(k in 1:n.station){ z.new[, , k]<- matrix(z[ ,k], nrow=period)}
# prepare to filter 1 PC scores
scores<-matrix(NA, ncol=n.station, nrow = how.many)
harm<-matrix(NA, ncol=n.station, nrow = period)
# convert each station records into functional object
# store the scores in the matrix
# create basis
d<-create.bspline.basis(rangeval=c(0, period), nbasis=159)
t<-seq(1, period, 1)
#convert data into a functional object
for(k in 1:n.station){
    fd<-Data2fd(t,z.new[,,k],basisobj=d)
    #compute the first PC
    pca<-pca.fd(fd, nharm = 1, centerfns = F)
    scores[, k]<-pca$scores[,1]
    harm[, k]<-eval.fd(c(1:period), pca$harmonics[1])
}
#replace the highest 5% of scores with the median score
q<-numeric(n.station)
for(i in 1:n.station) q[i]<-quantile(abs(scores[,i]), probs=0.9)
med<-numeric(n.station)
for(i in 1:n.station) med[i]<-median(scores[,i])
delta<-matrix(0, ncol=n.station, nrow = how.many)
for(j in 1:n.station){
    for(i in 1:how.many){
         if(abs(scores[i,j]) > q[j]) delta[i,j]<-1
                        }
}
ind<-apply(delta,1,sum)
for(j in 1:n.station){
    for(i in 1:how.many){
         if(ind[i]==(n.station)) scores[i,j]<-med[j]
                        }
}
# to extract pseudo-Sq
SQ<-matrix(data = NA, ncol = n.station, nrow = N)
#SQa<-matrix(data = NA, ncol = n.station, nrow = N)
#mean.fd<-eval.fd(c(1:1440),mean(fd))
for(j in 1:n.station){
    for(i in 1:how.many){
        s<-period*(i-1)+1
        e<-period*i
#d.sq<-(harm[,j]*scores[i,j]+mean.fd[j])
#SQ[s:e,j]<-d.sq
d.sq<-(harm[,j]*scores[i,j])
SQ[s:e,j]<-d.sq
                        }
}
list(sq=SQ, diff=(data-SQ))
}
