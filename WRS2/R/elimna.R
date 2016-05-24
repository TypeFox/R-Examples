elimna <-
function(m){
#
# remove any rows of data having missing values
#
if(is.null(dim(m)))m<-as.matrix(m)
ikeep<-c(1:nrow(m))
for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
elimna<-m[ikeep[ikeep>=1],]
elimna
}
