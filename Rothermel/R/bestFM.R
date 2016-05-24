bestFM <-
function (obs,m,u,slope) {

SFM_metric<-get (data (SFM_metric,
            envir = environment ()))

rmse <- function (pred,obs) {
  sqrt ((sum((obs-pred)^2))/length(obs))
}

pred<-matrix(rep(NA,53*length(obs)),nrow=length(obs))
pred<-as.data.frame(pred)
error<-rep(NA,53)
bias<-rep(NA,53)

for (i in 1:53) {
  pred[,i]<-ros(
    modeltype=SFM_metric[i,"Fuel_Model_Type"],
    w=SFM_metric[i,2:6],
    s=SFM_metric[i,7:11],
    delta=SFM_metric[i,"Fuel_Bed_Depth"],
    mx.dead=SFM_metric[i,"Mx_dead"],
    h=SFM_metric[i,14:18],
    m=m,
    u=u,
    slope=slope)[15]
  
  error[i]<-rmse(pred[,i],obs)
  bias[i]<-mean(pred[,i]-obs)
  
}

df<-data.frame(error,bias)
rownames(df)<-names(bias)<-rownames(SFM_metric)
df<-df[with(df, order(error)), ]

if (length(obs)>1) return(df)
if (length(obs)==1) return(bias)
}
