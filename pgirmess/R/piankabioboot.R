"piankabioboot" <-
function(dataframe1,dataframe2,B=1000,probs=c(0.025,0.975)){


    obs<-piankabio(dataframe1,dataframe2)
    bo<-rep(0,B)
    for (i in 1:B){
    x1<-sample(1:length(dataframe1[,1]),replace=TRUE)
    x2<-sample(1:length(dataframe2[,1]),replace=TRUE)
    bo[i]<-piankabio(dataframe1[x1,],dataframe2[x2,])
    }
    x<-t(c(obs_mean=obs,boot_mean=mean(bo),emp_CI=quantile(bo,probs)))
    row.names(x)<-c("Pianka's index")
    return(round(x,3))
}
