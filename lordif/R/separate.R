separate <-
function(resp,flag,gr) {
    resp.nodif<-resp[,!flag,drop=FALSE] 
    resp.dif<-resp[,flag,drop=FALSE]
    nobs<-length(gr) 
    gr.freq<-table(gr) 
    ng<-length(gr.freq) 
    gr.label<-names(gr.freq) 
    ndif<-sum(flag)
    dif.items<-names(resp.dif) 
    sparse.resp<-data.frame(matrix(NA,nobs,ndif*ng)) 
    colnames(sparse.resp)<-paste0(rep(dif.items,rep(ng,ndif)),".",rep(1:ng,ndif)) 
    for (i in 1:ndif) {
      for (j in 1:ng) {
        sparse.resp[gr==gr.label[j],ng*i-ng+j]<-resp.dif[gr==gr.label[j],i] 
      }
    }
    out<-data.frame(resp.nodif,sparse.resp)
    return(out)
  }
