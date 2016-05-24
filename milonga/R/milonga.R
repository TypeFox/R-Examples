milonga<-function(data, null_p=NULL, weight=10, burnin=NULL, rep=10){

  data<-as.matrix(data)
  n<-nrow(data)
  m<-ncol(data)
  pattern<-do.call(expand.grid, replicate(m, 0:1, simplify=FALSE))
  colnames(pattern)<-colnames(data)
  K<-nrow(pattern)
  cdata<-data[complete.cases(data), ]
  
  impi<-function(datai){
      subn<-NULL
      subp<-NULL
      wt.n<-NULL
      for(k in 1:K){
            a<-pattern[k,]
            srow<-grepl(paste(a, collapse="-"), 
                        apply(datai, 1, paste, collapse="-"))
            subn[k]<-sum(srow)
            if(is.null(null_p)){
                subp[k]<-1/K
            }else{
                subp[k]<-ifelse(sum(pattern[k,])==0, null_p, (1-null_p)/(K-1))
            }
            wt.n[k]<-(subn[k]+subp[k]*weight)/(n+weight)*n    
      } 
      
      ptheta<-rdirichlet(1,wt.n)

      datax<-matrix(NA, ncol=m, nrow=n)
      for (i in 1:n){
            nna<-sum(is.na(data[i, ]))
            if(nna==0){
              datax[i, ]<-data[i, ]
            } else if (nna==m){
                  sp<-rmultinom(1,1,ptheta)==1
                  datax[i,]<-as.matrix(pattern[sp,])
            } else {
                  b<-data[i,]
                  b[is.na(b)]<-"."
                  prow<-grepl(paste(b, collapse="-"), 
                              apply(pattern, 1, paste, collapse="-"))
                  sptheta<-ptheta[prow]
                  spattern<-pattern[prow,]
                  ss<-rmultinom(1,1,sptheta/sum(sptheta))==1
                  datax[i,]<-as.matrix(spattern[ss,])
            }
       }
    return(datax)
  }

  if(is.null(burnin)){
      datalist<-vector("list", rep)
      datalist[[1]]<-impi(cdata)

      for(i in 1:(rep-1)){
          colnames(datalist[[i]])<-colnames(data)
          datalist[[i+1]]<-impi(datalist[[i]])
      }

      return(datalist)

  }else{
      burninlist<-vector("list", burnin)
      burninlist[[1]]<-impi(cdata)

      for(i in 1:(burnin-1)){
          colnames(burninlist[[i]])<-colnames(data)
          burninlist[[i+1]]<-impi(burninlist[[i]])
      }

      datalist<-vector("list", rep)
      datalist[[1]]<-impi(burninlist[[burnin]])

      for(i in 1:(rep-1)){
          colnames(datalist[[i]])<-colnames(data)
          datalist[[i+1]]<-impi(datalist[[i]])
      }

      return(list(burnin=burninlist, imp=datalist))
  }

}


