calctheta <-
function(ipar,resp.data,theta,prior.mean=0.0,prior.sd=1.0,model="GRM") {
    if (!(model %in% c("GRM","GPCM"))) {
      warning("model must be either \"GRM\", or \"GPCM\"; will be reset to default")
      model<-"GRM"
    }
    prior<-dnorm((theta-prior.mean)/prior.sd) 
    pp<-calcprob(ipar,theta,model=model)
    nExaminees<-nrow(resp.data) 
    nq<-length(theta) 
    posterior<-matrix(rep(prior,nExaminees),nExaminees,nq,byrow=T)
    for (i in 1:nrow(ipar)) {
      resp<-matrix(resp.data[,i],nExaminees,1)
      prob<-t(pp[,i,resp])
      prob[is.na(prob)]<-1.0
      posterior<-posterior*prob
    }
    EAP<-posterior%*%theta/rowSums(posterior)
    SEM<-sqrt(rowSums(posterior*(matrix(theta,nExaminees,nq,byrow=T)-matrix(EAP,nExaminees,nq))^2)/rowSums(posterior))
    return(list(theta=as.vector(EAP),SE=SEM))
  }
