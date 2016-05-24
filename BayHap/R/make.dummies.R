`make.dummies` <-
function(covariate){
    
    ifelse((!is.factor(covariate)),f<-as.factor(covariate),f<-covariate)
    ref<-levels(f)[1]
    names.cat<-levels(f)
    dims<-length(levels(f))*length(f)
    m<-matrix(rep(0,dims),ncol=length(levels(f)),nrow=length(f))
    i<-1
    for (i in 1:length(f)){
        if(f[i]!=levels(f)[1]){
            m[i,names.cat==f[i]]<-1
        }
    }
    m<-as.data.frame(m[,-1])
    names(m)<-levels(f)[-1]
    m

}

