################################################################################
#  Below Not Used Yet                                                          #
#                                                                              #
################################################################################
#
#rknn.FNN<- function(data, y, k=5, r=500, mtry=trunc(sqrt(ncol(data))))
#{
#
#   n<- nrow(data);
#   index1<- which(y==0); n1<- length(index1)
#   index2<- which(y==1); n2<- length(index2)
#
#   p<- ncol(data);
#
#   selected<- matrix(0, nrow=r, ncol=mtry);
#   acc<- numeric(r);
#
#   for(i in 1:r){
#        dset<- c(sample(index1, n1/2), sample(index2, n2/2));
#
#        fset<- sample(p, mtry);
#        selected[i,]<- fset;
#        aknn<- get.knnx(data=data[dset, fset], query=data[-dset, fset], k=k);
#        signs<- (y[dset][aknn$nn.index] ==  y[-dset])*2-1;
#        #signs<- y[dset][aknn$nn.index] ==  y[-dset];
#        acc[i]<- sum(1/aknn$nn.dist^2*signs)
#   }
#
#  selected<- as.vector(selected);
#  features<- unique(selected);
#
#  mreal<- length(features);
#  support<- numeric(mreal);
#
#  for(j in 1:mreal){
#      acc.ind<- (which(selected==features[j])-1)%%r+1;
#      support[j]<- mean(acc[acc.ind]);
#  }
#
#  names(support)<- colnames(data)[features];
#
#  res<- sort(support, decreasing=T);
#  attr(res, "meanacc")<- mean(acc);
#
#  return(res);
#}
#
#rknn.dist<- function(data, r=500, k=1,
#                    mtry=trunc(sqrt(ncol(data))),
#                    y=NULL
#)
#{
#   n<- nrow(data);
#   p<- ncol(data);
#
#   res<- matrix(0, n, n);
#
#   for(i in 1:r){
#        fset<- sample(p, mtry);
#        nn<- get.knn(data[, fset], k=k);
#        for(j in 1:k){
#          index<- cbind(1:n, nn$nn.index[,j])
#          res[index]<- res[index]+nn$nn.dist[,j];
#        }
#   }
#
#  res<- (res+t(res))/2
#  #res<- sqrt(1- res/(r*k));
#  if(!is.null(y)) rownames(res)<- y;
#
#  return(as.dist(res))
#}
#
#################################################################################
