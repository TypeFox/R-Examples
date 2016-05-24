#' @export
#' @references Stephen J. Ganocy and Jiayang Sun (2014), "Heteroscedastic Change Point Analysis and Application to Footprint Data", J of Data Science, In Press.
hcp<-
  function(dataset,jlo,jhi,klo,khi,method=c("C-LLL","C-LQL","U-LLL","MMP"),variance = c("121", "Common", "Differ"),plot=c("FALSE","TRUE"),sigma21,r1,s1,sigma22,r2,s2,sigma23,r3,s3){
    jlow<-jlo
    jhigh<-jhi
    klow<-klo
    khigh<-khi
    sortdata<-dataset[order(dataset[,1]),]
    x<-sortdata[,1]
    y<-sortdata[,2]
    n<-length(x)
    jlo<-length(x[x<jlow])
    jhi<-length(x[x<jhigh])
    klo<-length(x[x<klow])
    khi<-length(x[x<khigh])
    method<-match.arg(method)
    variance<-match.arg(variance)
    plot<-match.arg(plot)
    
    if (method=="U-LLL") {
      
      if (variance=="121"){
        hats<-llsearch(x,y,n,jlo,jhi,klo,khi,plot)
        est<-p.est(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=c(est$sigma2,est$tau2),changepoints=c(est$xj,est$xk)))
      }
      
      
      if (variance == "Common"){
        hats<-llsearch.C(x,y,n,jlo,jhi,klo,khi,plot)
        est<-p.est.C(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=est$sigma2,changepoints=c(est$xj,est$xk)))
      }
      
      if (variance == "Differ"){
        hats<-llsearch.D(x,y,n,jlo,jhi,klo,khi,plot)
        est<-p.est.D(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=c(est$sigma2,est$tau2,est$u2),changepoints=c(est$xj,est$xk)))
      }
      
    }  
    
    if (method=="C-LLL") {
      
      if (variance=="121") {
        hats<-con.search(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con.vals(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma=c(sqrt(1/est$eta[1]),sqrt(1/est$eta[2])),coe=est$beta,changepoints=est$tau))
      }
      
      if (variance =="Common") {
        hats<-con.search.C(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con.vals.C(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],coe=est$beta,changepoints=est$tau))
      }
      
      if (variance =="Differ") {
        hats<-con.search.D(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con.vals.D(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2],1/est$eta[3]),coe=est$beta,changepoints=est$tau))
      }
      
    } 
    
    if (method=="C-LQL") {
      
      if (variance=="121") {
        hats<-con2.search(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con2.vals(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2]),coe=est$beta,changepoints=est$tau))
      }
      
      if (variance =="Common") {
        hats<-con2.search.C(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con2.vals.C(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma2=1/est$eta[1],coe=est$beta,changepoints=est$tau))
      }
      
      if (variance =="Differ") {
        hats<-con2.search.D(x,y,n,jlo,jhi,klo,khi,plot)
        est<-con2.vals.D(x,y,n,hats$jhat,hats$khat)
        return(list(maxloglik=hats$value,sigma2=c(1/est$eta[1],1/est$eta[2],1/est$eta[3]),coe=est$beta,changepoints=est$tau))
      }
      
    }
    
    if(method=="MMP") {
      if (variance=="121") {
        rr<-r1
        ss<-s1
        vv<-r2
        ww<-s2
        hats<- mmpiter(x,y,n,jlo,jhi,klo,khi,sigma21,sigma22,rr,ss,vv,ww)
        est<-p.est(x,y,n,hats$jhat,hats$khat)
        return(list(coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=c(est$sigma2,est$tau2),changepoints=c(est$xj,est$xk)))
      }
      
      if (variance =="Common") {
        rr<-r1
        ss<-s1
        hats<- mmpiter.C(x,y,n,jlo,jhi,klo,khi,sigma21,rr,ss)
        est<-p.est.C(x,y,n,hats$jhat,hats$khat)
        return(list(coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=est$sigma2,changepoints=c(est$xj,est$xk)))
      }
      
      if (variance =="Differ") {
        rr<-r1
        ss<-s1
        vv<-r2
        ww<-s2
        mm<-r3
        nn<-s3
        hats<- mmpiter.D(x,y,n,jlo,jhi,klo,khi,sigma21,sigma22,sigma23,rr,ss,vv,ww,mm,nn)
        est<-p.est.D(x,y,n,hats$jhat,hats$khat)
        return(list(coe=c(est$a0,est$a1,est$b0,est$b1,est$c0,est$c1),sigma2=c(est$sigma2,est$tau2,est$u2),changepoints=c(est$xj,est$xk)))
      }
      
    }
  }