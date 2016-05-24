Gini <-
function(y,p=rep(1,length(y)))
{
if (length(p)==0) {p<-rep(1,length(y))} else {p};
dataset<-cbind(y,p)
ord_y<-order(y)
dataset_ord<-dataset[ord_y,]
y<-dataset_ord[,1]
p<-dataset_ord[,2]  
N<-sum(p)
yp<-y*p
C_i<-cumsum(p)
num_1<-sum(yp*C_i)
num_2<-sum(yp)
num_3<-sum(yp*p)
m_y<-(sum(yp)/N)
G_num<-(2/N^2)*num_1-(1/N)*num_2-(1/N^2)*num_3
G<-(1/m_y)*G_num
list(Gini=G)
}
