sinc<-function(fh=NULL,fl=NULL,tw=0,sf,type,n){
  df<-tw/sf
  if(is.null(fl)!=TRUE){fcl<-fl/sf+df/2}
  if(is.null(fh)!=TRUE){fch<-fh/sf-df/2}
  
  if(type=="low"){pulse<-2*fcl*sin(n*2*pi*fcl)/(n*2*pi*fcl)
                  pulse[n==0]<-2*fcl
  }
  if(type=="high"){pulse<- -2*fch*sin(n*2*pi*fch)/(n*2*pi*fch)
                   pulse[n==0]<-1-2*fch
  }
  if(type=="band"){pulse<- (2*fcl*sin(n*2*pi*fcl)/(n*2*pi*fcl)) - (2*fch*sin(n*2*pi*fch)/(n*2*pi*fch))
                   pulse[n==0]<-2*(fcl-fch)
  }
  return(pulse)
}