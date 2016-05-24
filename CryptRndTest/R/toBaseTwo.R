toBaseTwo=function(x,m=128,prec=512,num.CPU=4){  
  M=length(x)  
  if (m>64){
    x=new('mpfr',unlist(x))
  }
  cl <- makePSOCKcluster(num.CPU)
  setDefaultCluster(cl)
  clusterExport(cl=cl, varlist=c("x"), envir=environment())
  bit=parLapply(cl,1:M,function(i){     
    if (m>64){
      bit=mpfrArray(0, prec, dim = c(m,M))
    }else{
      bit=array(0, dim = c(m,M))
    }
    j=1
    while (x[i]>0){
      if (j>m){      
        stop(paste0("Required bit length is greater than ",m,"!"))
      }    
      bit[j,i] = x[i]%%2 
      x[i] = floor(x[i]/2)    
      j=j+1
    }      
    return(bit)
  })
  if (m>64){
    r.bit=mpfrArray(0,prec,dim=c(m,M)) 
    r.bit=parLapply(cl,1:M,function(i){          
      r.bit[,i]=rev(new('mpfr',unlist(bit[[i]][,i])))    
    })
  }else{
    r.bit=array(0,dim=c(m,M)) 
    r.bit=parLapply(cl,1:M,function(i){          
      r.bit[,i]=rev(bit[[1]][,i])    
    })
  }
  stopCluster(cl) 
  return(r.bit)
}