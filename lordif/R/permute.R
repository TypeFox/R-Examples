permute <-
function(obj,alpha=0.01,nr=100) {
    call<-match.call()
    if (class(obj)!="lordif") stop(paste(deparse(substitute(obj))," must be of class lordif"))
    if (alpha<0 || alpha>1) {
      warning("alpha must be a fraction between 0 and 1; will be reset to .01")
      alpha<-.01
    }
    if (nr<100) warning("number of replications is less than 100")
    if (nr<=0) {
      warning("number of replications is not a positive integer; will be reset to 100")
      nr<-100
    }
    options<-obj$options
    ipar<-obj$ipar
    group<-obj$group
    weights<-obj$weights
    ncat<-obj$ncat
    resp<-obj$recoded
    theta<-obj$calib$theta
    selection<-obj$selection
    nobs<-length(theta)
    ni<-nrow(ipar)
    chi12<-matrix(NA,nr,ni);rownames(chi12)<-paste0("Rep",1:nr);colnames(chi12)<-paste0("I",selection)
    chi13<-matrix(NA,nr,ni);rownames(chi13)<-paste0("Rep",1:nr);colnames(chi13)<-paste0("I",selection)
    chi23<-matrix(NA,nr,ni);rownames(chi23)<-paste0("Rep",1:nr);colnames(chi23)<-paste0("I",selection)
    pseudo12.CoxSnell<-matrix(NA,nr,ni);rownames(pseudo12.CoxSnell)<-paste0("Rep",1:nr);colnames(pseudo12.CoxSnell)<-paste0("I",selection)
    pseudo13.CoxSnell<-matrix(NA,nr,ni);rownames(pseudo13.CoxSnell)<-paste0("Rep",1:nr);colnames(pseudo13.CoxSnell)<-paste0("I",selection)
    pseudo23.CoxSnell<-matrix(NA,nr,ni);rownames(pseudo23.CoxSnell)<-paste0("Rep",1:nr);colnames(pseudo23.CoxSnell)<-paste0("I",selection)
    pseudo12.Nagelkerke<-matrix(NA,nr,ni);rownames(pseudo12.Nagelkerke)<-paste0("Rep",1:nr);colnames(pseudo12.Nagelkerke)<-paste0("I",selection)
    pseudo13.Nagelkerke<-matrix(NA,nr,ni);rownames(pseudo13.Nagelkerke)<-paste0("Rep",1:nr);colnames(pseudo13.Nagelkerke)<-paste0("I",selection)
    pseudo23.Nagelkerke<-matrix(NA,nr,ni);rownames(pseudo23.Nagelkerke)<-paste0("Rep",1:nr);colnames(pseudo23.Nagelkerke)<-paste0("I",selection)
    pseudo12.McFadden<-matrix(NA,nr,ni);rownames(pseudo12.McFadden)<-paste0("Rep",1:nr);colnames(pseudo12.McFadden)<-paste0("I",selection)
    pseudo13.McFadden<-matrix(NA,nr,ni);rownames(pseudo13.McFadden)<-paste0("Rep",1:nr);colnames(pseudo13.McFadden)<-paste0("I",selection)
    pseudo23.McFadden<-matrix(NA,nr,ni);rownames(pseudo23.McFadden)<-paste0("Rep",1:nr);colnames(pseudo23.McFadden)<-paste0("I",selection)
    beta12<-matrix(NA,nr,ni);rownames(beta12)<-paste0("Rep",1:nr);colnames(beta12)<-paste0("I",selection)
    cat(paste0("Start time: ",date(),"\n\n"))
    for (r in 1:nr) {
      group.random<-group[sample(nobs)]
      out<-rundif(1:ni,resp,theta,group.random,options$criterion,options$alpha,options$beta.change,options$pseudo.R2,options$R2.change,weights)
      chi12[r,]<-out$stats$chi12
      chi13[r,]<-out$stats$chi13
      chi23[r,]<-out$stats$chi23
      pseudo12.CoxSnell[r,]<-out$stats$pseudo12.CoxSnell
      pseudo13.CoxSnell[r,]<-out$stats$pseudo13.CoxSnell
      pseudo23.CoxSnell[r,]<-out$stats$pseudo23.CoxSnell
      pseudo12.Nagelkerke[r,]<-out$stats$pseudo12.Nagelkerke
      pseudo13.Nagelkerke[r,]<-out$stats$pseudo13.Nagelkerke
      pseudo23.Nagelkerke[r,]<-out$stats$pseudo23.Nagelkerke
      pseudo12.McFadden[r,]<-out$stats$pseudo12.McFadden
      pseudo13.McFadden[r,]<-out$stats$pseudo13.McFadden
      pseudo23.McFadden[r,]<-out$stats$pseudo23.McFadden
      beta12[r,]<-out$stats$beta
      cat(paste0(" Replication: ",r,"\n"))
    }
    cat(paste0("\nEnd time: ",date(),"\n"))
    stat<-c("chi12","chi13","chi23","pseudo12.CoxSnell","pseudo13.CoxSnell","pseudo23.CoxSnell","pseudo12.Nagelkerke","pseudo13.Nagelkerke","pseudo23.Nagelkerke","pseudo12.McFadden","pseudo13.McFadden","pseudo23.McFadden","beta12")
    cutoff<-matrix(NA,ni,length(stat))
    for (i in 1:ni) {
      cutoff[i,1]<-getcutoff(chi12[,i],alpha,F)
      cutoff[i,2]<-getcutoff(chi13[,i],alpha,F)
      cutoff[i,3]<-getcutoff(chi23[,i],alpha,F)
      cutoff[i,4]<-getcutoff(pseudo12.CoxSnell[,i],alpha,T)
      cutoff[i,5]<-getcutoff(pseudo13.CoxSnell[,i],alpha,T)
      cutoff[i,6]<-getcutoff(pseudo23.CoxSnell[,i],alpha,T)
      cutoff[i,7]<-getcutoff(pseudo12.Nagelkerke[,i],alpha,T)
      cutoff[i,8]<-getcutoff(pseudo13.Nagelkerke[,i],alpha,T)
      cutoff[i,9]<-getcutoff(pseudo23.Nagelkerke[,i],alpha,T)
      cutoff[i,10]<-getcutoff(pseudo12.McFadden[,i],alpha,T)
      cutoff[i,11]<-getcutoff(pseudo13.McFadden[,i],alpha,T)
      cutoff[i,12]<-getcutoff(pseudo23.McFadden[,i],alpha,T)
      cutoff[i,13]<-getcutoff(beta12[,i],alpha,T)
    }
    colnames(cutoff)<-stat
    rownames(cutoff)<-paste0("I",selection)
    out<-list(call=call,chi12=chi12,chi13=chi13,chi23=chi23,
              pseudo12.CoxSnell=pseudo12.CoxSnell,pseudo13.CoxSnell=pseudo13.CoxSnell,pseudo23.CoxSnell=pseudo23.CoxSnell,
              pseudo12.Nagelkerke=pseudo12.Nagelkerke,pseudo13.Nagelkerke=pseudo13.Nagelkerke,pseudo23.Nagelkerke=pseudo23.Nagelkerke,
              pseudo12.McFadden=pseudo12.McFadden,pseudo13.McFadden=pseudo13.McFadden,pseudo23.McFadden=pseudo23.McFadden,
              beta12=beta12,alpha=alpha,nr=nr,cutoff=data.frame(cutoff))
    class(out)<-"lordif.MC"
    return(out)
  }
