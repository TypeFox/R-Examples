############# Functions from lsd.package!
sample.max.depth <-
  function(y,iter=100,eps=0.0001,p.length=10){
    # Calculates the maximum sample location-scale depth  
    # for the data set y 
    # Uses function sample.depth.cont.for.mu and sample.max.depth.for.mu
    # p.length is the maximum length of the precision step at the end
    y<-sort(y)
    N<-length(y)
    d.min<-floor(N/3)
    n.mid<-round(N/2)
    res<-sample.max.depth.for.mu(mu=y[n.mid],y=y,d.min=d.min,iter=iter,eps=eps)
    d<-res["d"]
    s<-res["sigma"]
    difb<-res["difbound"]
    i.mu<-res["i"]
    #cat("d: ", d, "s: ", s, "difb: ", difb, "\n")
    if(d<N/2){
      i<-1
      n.up<-ceiling(N*2/3)
      n.low<-floor(N/3)
      #  cat("i: ",i, "depth: ", d, "n.mid: ", n.mid, "s: ", s, "\n")
      #  cat( "    n.up: ", n.up, "n.low: ", n.low, "i.mu: ", res["i"], "difb: ", difb, "\n")
      dec<-1
      while(i<iter & (n.up-n.low>1 & dec>0)){
        # cat(sprintf("Nup i low %i %i %i ",d,n.up,n.low))
        i<-i+1
        n.mid.low<-ceiling((n.mid+n.low)/2)
        n.mid.up<-floor((n.mid+n.up)/2)
        #    cat("i: ", i, " n.low:", n.low, " n.mid.low", n.mid.low, " n.mid:", n.mid,
        #           " n.mid.up:" , n.mid.up, " n.up:", n.up, "\n")
        ##   cat(sprintf("|| %f %i ||", n.mid.low, n.mid.up))
        #    cat(sprintf("|| %f %i ||", y[n.mid.low], d.min))
        res.low<-sample.max.depth.for.mu(mu=y[n.mid.low],y=y,d.min=d.min,
                                         iter=iter,eps=eps)
        res.up<-sample.max.depth.for.mu(mu=y[n.mid.up],y=y,d.min=d.min,
                                        iter=iter,eps=eps)
        i.mu<-cbind(i.mu,res.low["i"])
        i.mu<-cbind(i.mu,res.up["i"])
        d.low<-res.low["d"]
        d.up<-res.up["d"]
        if(d.low> d ){
          
          
          
          d<-d.low
          s<-res.low["sigma"]
          difb<-res.low["difbound"]
          n.up<-n.mid
          n.mid<-n.mid.low
          #       cat("1.Fall (low): ", "depth: ", d,  "s: ", s,  "n.low: ", n.low,  
          #              "n.mid: ", n.mid," n.up: ", n.up,"\n")
        } # d.low>d
        else{
          if(d.up >d){
            
            d<-d.up
            s<-res.up["sigma"]
            difb<-res.up["difbound"]
            n.low<-n.mid
            n.mid<-n.mid.up
            #         cat("2. Fall (up) ", "depth: ", d, "s: ", s,  "n.low: ", n.low,  
            #                  "n.mid: ", n.mid," n.up: ", n.up,"\n")
          } # d.up
          else{
            if(d.low< d | d.up < d){
              
              if(d.low < d){
                
                n.low<-n.mid.low
                #               cat("3. Fall (d.low<d): ", "depth: ", d, "s: ", s,  
                #                   "n.low: ", n.low,"n.mid: ", n.mid, " n.up: ", n.up, "\n")
              } # d.low< d 
              if(d.up < d){
                
                n.up<-n.mid.up
                #               cat("4. Fall (d.up<d): ", "depth: ", d,"s: ", s, "n.low: ", 
                #                      n.low,  "n.mid: ", n.mid," n.up: ", n.up,  "\n")
              } # d.up < d
            } # d.low< d | d.up < d
            else{
              dec<-0
              #            cat("5. Fall: ", "depth: ", d, "s: ", s,  "n.low: ", n.low, 
              #                            "n.mid: ", n.mid, " n.up: ", n.up,"\n")
            } # d.low< d | d.up < d else
          } # d.up else
        }# d.low else
        #    cat(" KONIEC \n")
      } # while
    } # d< N/2
    # Precision step
    length<-max(c(p.length,2*(n.up-n.low+1)))
    #cat("Precision step", "length: ", length, "\n")
    mu<-seq(y[n.low],y[n.up],length=length)
    res<-sample.max.depth.for.mu(mu=mu[1],y=y,d.min=d.min,iter=iter,eps=eps)
    i.mu<-cbind(i.mu,res["i"])
    for(i in 2:length){
      res<-rbind(res,sample.max.depth.for.mu(mu=mu[i],y=y,d.min=d.min,iter=iter,eps=eps))
      #   cat("i: ", i, "res: ", res[i,"i"], "\n")
      i.mu<-cbind(i.mu, res[i,"i"])
    }
    d<-max(res[,"d"])
    dmax<-rep(1,length)
    for(i in 1:length){
      if(res[i,"d"]<d){
        dmax[i]<-0
      }
    }
    #cat("res[,d]: ", res[,"d"], "\n")
    #cat("dmax:      ", dmax, "\n")
    dmax<-as.logical(dmax)
    #cat("dmax:      ", dmax, "\n")
    res<-res[dmax,]
    mu<-mu[dmax]
    #cat("mu: ", mu, "\n")
    if(length(mu)>1){
      length.half<-round(length(mu)/2) 
      mu<-mu[length.half]   
      s<-res[length.half,"sigma"]
      difb<-res[length.half,"difbound"]
    }
    else{
      s<-res["sigma"]
      difb<-res["difbound"]
    }
    #cat(i.mu,"\n")
    cat("Sum of all iterations: ", sum(i.mu), " Maximum number of iterations in substeps: ", max(i.mu)," Number of substeps: ", length(i.mu), "\n")
    res<-c(d,mu,s)
    names(res)<-c("max.depth","mu","sigma")
    res
  }




sample.max.depth.for.mu <-
  function(mu,y,d.min=0,iter=100,eps=0.0000001){
    # Calculates the maximum sample location-scale depth for a given mu 
    # for the data set y 
    # Uses function sample.depth.cont.for.mu
    y<-sort(y)
    M<-0
    N<-length(y)
    
    for(i in 1:N){
      if(y[i]<mu){
        M<-M+1
      }
      
    }
    i = 0
    if(y[M+1]>mu){ 
      d<-min(M,N-M)
    }
    else{
      d<-min(M,N-M-1)
    }  
    cont<-sample.depth.cont.for.mu(d,mu,y=y)
    difbound<-cont["ubound"]-cont["lbound"]
    if(abs(difbound)>eps){
      
      i<-1
      d.up<-d
      d.low<-d.min
      # cat(paste(d.up, d.low, i<iter , abs(difbound)>eps , d.up-d.low>1))
      #  cat("\n")
      while(i<iter & (abs(difbound)>eps & d.up-d.low>1)){
        i<-i+1
        if(difbound < -eps){
          d.up<-d
          d<-round((d.up+d.low)/2) 
          cont<-sample.depth.cont.for.mu(d,mu,y=y)
          difbound<-cont["ubound"]-cont["lbound"]
        }
        else{
          d.low<-d
          d<-round((d.up+d.low)/2)
          cont<-sample.depth.cont.for.mu(d,mu,y=y)
          difbound<-cont["ubound"]-cont["lbound"]
        }
      }
    }
    if(difbound< -eps){
      
      d<-d-1
      cont<-sample.depth.cont.for.mu(d,mu,y=y)
      difbound<-cont["ubound"]-cont["lbound"]
    }
    res<-c(d,(cont["ubound"]+cont["lbound"])/2,i,difbound)
    names(res)<-c("d","sigma","i","difbound")
    res 
  }



sample.depth.cont.for.mu <-
  function(d,mu,y,length=100){
    # Calculates the depth contours with respect to sigma for a given mu
    # for the the data set y for depth d
    # The result is a vector, where lbound denoted the lower bound, ubound the upperbound and 
    # tbound is a logic variable which is true if lbound<=ubound 
    y<-sort(y)
    M<-0
    N<-length(y)
    for(i in 1:N){
      if(y[i]<mu){
        M<-M+1
      }
      i<-i+1
    }
    # cat(sprintf("M value %i \n", M))
    
    if(y[M+1]>mu){       
      if((d<=M & d<=N-M) & d>0){
        k<-seq(M-d+1,M,1)
        lk<-(mu-y[k])*(y[d+k]-mu)
        #print(k)
        lbound<-max(sqrt(abs(lk)))
        k<-seq(1,d,1)
        uk<-(mu-y[k])*(y[N-d+k]-mu) 
        ubound<-min(sqrt(abs(uk)))
        
        case<-0
        if(lbound <= ubound){tbound<-T}
        else{tbound<-F}
      }
      else{
        lbound<-1
        ubound<-0
        case<-0
        tbound<-F
      }
    }
    if(y[M+1]==mu){       
      if((d<=M & d<=N-M-1) & d>0){
        k<-seq(M-d+1,M,1)
        lk<-(mu-y[k])*(y[d+k]-mu)
        lbound<-max(sqrt(abs(lk)))
        k<-seq(1,d,1)
        uk<-(mu-y[k])*(y[N-d+k]-mu)
        
        ubound<-min(sqrt(abs(uk)))
        
        case<-1
        if(lbound <= ubound){tbound<-T}
        else{tbound<-F}
      }
      else{
        lbound<-1
        ubound<-0
        case<-1
        tbound<-F
      }
    }
    bound<-c(lbound,ubound,tbound, case,M)
    names(bound)<-c("lbound","ubound","tbound","case","M")
    bound 
  }
