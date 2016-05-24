ModelEff=function(nfac=3,mod=1,nproc=0,dir=1,ufunc=mod,dimensions = list(NULL), pvslice=c(1,1,1),
                  lc=c(0,0,0,0,0,0,0,0,0,0,0,0),uc=c(1,1,1,1,1,1,1,1,1,1,1,1))
{
  library(mixexp)
  # glimit is the resolution for the plot
  glimit<-25
  
  # Check for allowed number of process variables
  if (nproc>3)
    stop("No more than 3 process variables allowed")
  
  # check for valid mod
  if(mod < 1 | mod > 6) 
    stop("mod must be a number between 1 and 6")
  #  if(mod==4 | is.null(ufunc))
  #    stop("When mod=4, you must supply the user model as an lm function through ufunc=object")
  
  # checks for valid number of factors
  if (nfac<2 | nfac>12)
    stop("The number factors must be between 2 and 12")
  if (nfac>=7 & mod>=2)
    stop("Linear models only when the number of factors is greater than 6")
  if (nfac<3 & mod==3)
    stop("cubic or special cubic models require at least 3 factors")
  if (nfac<3 & mod==4)
    stop("cubic or special cubic models require at least 3 factors")
  if(mod>=5 & nproc==0)
    stop("for models 5 or 6 there must be at least one process variable")
  
  # get the beta coefficients from the user function supplied
  if (! is.null(ufunc)) {
    Beta <- ufunc$coeff
  }
  # Set Upper and Lower Bounds for the components
  
  x1=c(0,1)
  x2=c(0,1)
  x3=c(0,1)
  x4=c(0,1)
  x5=c(0,1)
  x6=c(0,1)
  x7=c(0,1)
  x8=c(0,1)
  x9=c(0,1)
  x10=c(0,1)
  x11=c(0,1)
  x12=c(0,1)
  if(! is.null(lc)){
    if (nfac>=1) {
      x1[1]<-lc[1]
      x1[2]<-uc[1]
    }
    if (nfac>=2) {
      x2[1]<-lc[2]
      x2[2]<-uc[2]
    }
    if (nfac>=3) {
      x3[1]<-lc[3]
      x3[2]<-uc[3]
    }
    if (nfac>=4) {   
      x4[1]<-lc[4]
      x4[2]<-uc[4]
    }
    if (nfac>=5) {
      x5[1]<-lc[5]
      x5[2]<-uc[5]
    }    
    if (nfac>=6) {
      x6[1]<-lc[6]
      x6[2]<-uc[6]
    }    
    if (nfac>=7) {
      x7[1]<-lc[7]
      x7[2]<-uc[7]
    }    
    if (nfac>=8) {
      x8[1]<-lc[8]
      x8[2]<-uc[8]
    }    
    if (nfac>=9) {
      x9[1]<-lc[9]
      x9[2]<-uc[9]
    }    
    if (nfac>=10) {
      x10[1]<-lc[10]
      x10[2]<-uc[10]
    }    
    if (nfac>=11) {
      x11[1]<-lc[11]
      x11[2]<-uc[11]
    }    
    if (nfac>=12) {
      x12[1]<-lc[12]
      x12[2]<-uc[12]
    }    
    
  }
  # get constraint matrix for crvtave
  
  ck<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
  v<-c(-x1[1],x1[2])
  for (i in 2:nfac){
    v<-c(v,-ck[1,i],ck[2,i])
  }
  
  
  
  Ip<-diag(nfac)
  In<--1*Ip
  conmx<-interleave(Ip,In)
  conmx<-cbind(conmx,v)
  
  
  # calls crvtave to create exteme vertices design plus centroid
  ndm<-0
  des<-crvtave(ndm,conmx)
  des<-as.matrix(des)
  
  # selects the centroid from the design
  cnr<-des[,nfac+1]>0
  cent<-des[cnr,1:nfac]
  
  # gets pseudo component transformation constants
  J<-matrix(rep(1,times=12),ncol=1)
  Bds<-ck%*%J
  Den<-1-Bds[1,1]
  
  # transforms centroid to pseudo-component space
  pcent<-c(rep(0,times=nfac))
  for (i in 1:nfac) {
    pcent[i]<-(cent[i]-ck[1,i])/Den
  }
  
  # create grid
  grid<-glimit:1
  
  # create top half of x1 grid
  ifac<-1
  Xgrid<-array(rep(0,times=glimit*nfac*nfac))
  dim(Xgrid)<-c(nfac,glimit,nfac)
  
  
  # gets direction for Piepel method
  
  if (dir==1) {
    
    
    # get xi max at pseudo component vertex
    Vi<-Den+x1[1]
    Delta<-(Vi-cent[ifac])/glimit
    w<-Delta*grid+cent[ifac]
    Xgrid[ifac,,ifac]<-w
    if (ifac>1) {
      for (i in 1:(ifac-1)) {
        Xgrid[ifac,,i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
      }
    }
    if (ifac<nfac){
      for (i in (ifac+1):nfac) {
        Xgrid[ifac,,i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
      }
    }
    
    
    
    # when ifac=1 augment grid below with centroid
    if (ifac==1) {
      Temp<-array(rep(0,times=(glimit+1)*nfac*nfac))
      dim(Temp)<-c(nfac,(glimit+1),nfac)
      Temp[ ,1:glimit, ]<-Xgrid[, 1:glimit, ]
      
      
      
      for (i in 1:nfac) {
        Temp[i,,]<-rbind(Xgrid[i,,],cent)
      }
      Xgrid<-Temp
    }
    
    
    
    # create bottom half of x1 grid
    Temp<-array(rep(0,times=(2*glimit+1)*nfac*nfac))
    dim(Temp)<-c(nfac,(2*glimit+1),nfac)
    Temp[ , 1:(glimit+1), ]<-Xgrid[ , 1:(glimit+1), ]
    Xgrid<-Temp
    
    
    
    # get reverse grid
    grid<--1*(1:glimit)
    Delta<-(cent[ifac]-x1[1])/glimit
    w<-Delta*grid+cent[ifac]
    Xgrid[ifac,((glimit+2):(2*glimit+1)),ifac]<-w
    if (ifac>1) {
      for (i in 1:(ifac-1)) {
        Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
      }
    }
    if (ifac<nfac) {
      for (i in (ifac+1):nfac) {
        Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
      }
    }
    
    
    
    for (ifac in 2:nfac) {
      grid<-glimit:1
      
      # create top half of xi grid
      
      
      # get xi max at pseudo component vertex
      Vi<-Den+ck[1,ifac]
      Delta<-(Vi-cent[ifac])/glimit
      w<-Delta*grid+cent[ifac]
      Xgrid[ifac,(1:glimit),ifac]<-w
      if (ifac>1) {
        for (i in 1:(ifac-1)) {
          Xgrid[ifac,(1:glimit),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
        }
      }
      if (ifac<nfac){
        for (i in (ifac+1):nfac) {
          Xgrid[ifac,(1:glimit),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
        }
      }
      
      
      
      # create bottom half of xi grid
      grid<--1*(1:glimit)
      Delta<-(cent[ifac]-ck[1,ifac])/glimit
      w<-Delta*grid+cent[ifac]
      Xgrid[ifac,((glimit+2):(2*glimit+1)),ifac]<-w
      if (ifac>1) {
        for (i in 1:(ifac-1)) {
          Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
        }
      }
      if (ifac<nfac) {
        for (i in (ifac+1):nfac) {
          Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-cent[i]-(Delta*grid)*(cent[i]/(1-cent[ifac]))
        }
      }
    }
    
    mTitle<-"Effect Plot (Piepel direction)"
  } 
  # end of code for Piepel direction coordinates
  
  
  
  # get direction for Cox method
  
  if (dir==2) {
    
    # get xi max at pseudo component vertex
    Vi<-1.0
    Delta<-(Vi-pcent[ifac])/glimit
    w<-Delta*grid+pcent[ifac]
    
    Xgrid[ifac,,ifac]<-w
    if (ifac>1) {
      for (i in 1:(ifac-1)) {
        Xgrid[ifac,,i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
      }
    }
    if (ifac<nfac){
      for (i in (ifac+1):nfac) {
        Xgrid[ifac,,i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
      }
    }
    
    
    
    # when ifac=1 augment grid below with centroid
    if (ifac==1) {
      Temp<-array(rep(0,times=(glimit+1)*nfac*nfac))
      dim(Temp)<-c(nfac,(glimit+1),nfac)
      Temp[ ,1:glimit, ]<-Xgrid[, 1:glimit, ]
      
      
      
      for (i in 1:nfac) {
        Temp[i,,]<-rbind(Xgrid[i,,],pcent)
      }
      Xgrid<-Temp
    }
    
    
    
    # create bottom half of x1 grid
    Temp<-array(rep(0,times=(2*glimit+1)*nfac*nfac))
    dim(Temp)<-c(nfac,(2*glimit+1),nfac)
    Temp[ , 1:(glimit+1), ]<-Xgrid[ , 1:(glimit+1), ]
    Xgrid<-Temp
    
    
    
    # get reverse grid
    grid<--1*(1:glimit)
    Delta<-(pcent[ifac]-0)/glimit
    w<-Delta*grid+pcent[ifac]
    Xgrid[ifac,((glimit+2):(2*glimit+1)),ifac]<-w
    if (ifac>1) {
      for (i in 1:(ifac-1)) {
        Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
      }
    }
    if (ifac<nfac) {
      for (i in (ifac+1):nfac) {
        Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
      }
    }
    
    
    
    for (ifac in 2:nfac) {
      grid<-glimit:1
      
      # create top half of xi grid
      
      
      # get xi max at pseudo component vertex
      Vi<-1.0
      Delta<-(Vi-pcent[ifac])/glimit
      w<-Delta*grid+pcent[ifac]
      Xgrid[ifac,(1:glimit),ifac]<-w
      if (ifac>1) {
        for (i in 1:(ifac-1)) {
          Xgrid[ifac,(1:glimit),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
        }
      }
      if (ifac<nfac){
        for (i in (ifac+1):nfac) {
          Xgrid[ifac,(1:glimit),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
        }
      }
      
      
      
      # create bottom half of xi grid
      grid<--1*(1:glimit)
      Delta<-(pcent[ifac]-0)/glimit
      
      w<-Delta*grid+pcent[ifac]
      Xgrid[ifac,((glimit+2):(2*glimit+1)),ifac]<-w
      if (ifac>1) {
        for (i in 1:(ifac-1)) {
          Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
        }
      }
      if (ifac<nfac) {
        for (i in (ifac+1):nfac) {
          Xgrid[ifac,((glimit+2):(2*glimit+1)),i]<-pcent[i]-(Delta*grid)*(pcent[i]/(1-pcent[ifac]))
        }
      }
    }
    
    # transform from psuedo components to regular components
    Temp<-Xgrid
    for  (i in 1:nfac) {
      for (j in 1:nfac)  {
        for (k in 1:(2*glimit+1)) {
          Xgrid[i,k,j]<-Den*Temp[i,k,j]+ck[1,j]
        }
      }
      
    }
    
    
    mTitle<-"Effect Plot (Cox direction)"
  }
  
  
  
  # end of code for Cox direction 
  
  # set up the plot matrix
  PX<-array(rep(0,times=nfac*2*(2*glimit+1)))
  dim(PX)<-c((2*glimit+1),2*nfac)
  
  
  
  # get predicted values
  for (ifac in 1:nfac) {
    Xtemp<-Xgrid[ifac, , ]
    if (mod==1 | mod==3| mod==4) {
      Xmat<-Xtemp
    }
    if (mod==2) {
      if (nfac==2) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2])            
      }
      if (nfac==3) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3])
      }
      if (nfac==4) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,4]*Xtemp[,5])
      }
      if (nfac==6) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,1]*Xtemp[,6],Xtemp[,2]*Xtemp[,3],
                    Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,2]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,3]*Xtemp[,6],
                    Xtemp[,4]*Xtemp[,5],Xtemp[,4]*Xtemp[,6],Xtemp[,5]*Xtemp[,6])
      }
    }
    
    
    if (mod==3) {
      if (nfac==3) {
        Xmat<-Xtemp
        a<-Xtemp[,1]
        b<-Xtemp[,2]
        c12<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,3]
        c13<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,3]
        c23<-cubic(a,b)
        Xmat<-cbind(Xmat,c12,c13,c23,
                    Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3], 
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3])
      }
      
      if (nfac==4) {
        Xmat<-Xtemp
        a<-Xtemp[,1]
        b<-Xtemp[,2]
        c12<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,3]
        c13<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,4]
        c14<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,3]
        c23<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,4]
        c24<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,4]
        c34<-cubic(a,b)
        Xmat<-cbind(Xmat,c12,c13,c14,c23,c24,c34,
                    Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,4])
      }
      
      if (nfac==5) {
        Xmat<-Xtemp
        a<-Xtemp[,1]
        b<-Xtemp[,2]
        c12<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,3]
        c13<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,4]
        c14<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,5]
        c15<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,3]
        c23<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,4]
        c24<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,5]
        c25<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,4]
        c34<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,5]
        c35<-cubic(a,b)
        a<-Xtemp[,4]
        b<-Xtemp[,5]
        c45<-cubic(a,b)
        Xmat<-cbind(Xmat,c12,c13,c14,c15,c23,c24,c25,c34,c35,c45,
                    Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],
                    Xtemp[,4]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,3],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,4],Xtemp[,1]*Xtemp[,2]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,4]*Xtemp[,5],Xtemp[,2]*Xtemp[,3]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,5],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,3]*Xtemp[,4]*Xtemp[,5])
      }
      if (nfac==6) {
        Xmat<-Xtemp
        a<-Xtemp[,1]
        b<-Xtemp[,2]
        c12<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,3]
        c13<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,4]
        c14<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,5]
        c15<-cubic(a,b)
        a<-Xtemp[,1]
        b<-Xtemp[,6]
        c16<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,3]
        c23<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,4]
        c24<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,5]
        c25<-cubic(a,b)
        a<-Xtemp[,2]
        b<-Xtemp[,6]
        c26<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,4]
        c34<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,5]
        c35<-cubic(a,b)
        a<-Xtemp[,3]
        b<-Xtemp[,6]
        c36<-cubic(a,b)
        a<-Xtemp[,4]
        b<-Xtemp[,5]
        c45<-cubic(a,b)
        a<-Xtemp[,4]
        b<-Xtemp[,6]
        c46<-cubic(a,b)
        a<-Xtemp[,5]
        b<-Xtemp[,6]
        c56<-cubic(a,b)
        Xmat<-cbind(Xmat,c12,c13,c14,c15,c16,c23,c24,c25,c26,c34,c35,c36,c45,c46,c56,
                    Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,1]*Xtemp[,6],Xtemp[,2]*Xtemp[,3],
                    Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,2]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,3]*Xtemp[,6],
                    Xtemp[,4]*Xtemp[,5],Xtemp[,4]*Xtemp[,6],Xtemp[,5]*Xtemp[,6],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,6],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,6],Xtemp[,1]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,4]*Xtemp[,6],Xtemp[,1]*Xtemp[,5]*Xtemp[,6],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,6],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,2]*Xtemp[,4]*Xtemp[,6],Xtemp[,2]*Xtemp[,5]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,4]*Xtemp[,5],Xtemp[,3]*Xtemp[,4]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,5]*Xtemp[,6],Xtemp[,4]*Xtemp[,5]*Xtemp[,6])
      }
    }
    if(mod==4) {
      if (nfac==3) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3])
      }
      if (nfac==4) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],
                    Xtemp[,4]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,3],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,4],Xtemp[,1]*Xtemp[,2]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,4]*Xtemp[,5],Xtemp[,2]*Xtemp[,3]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,5],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,3]*Xtemp[,4]*Xtemp[,5])
      }
      if (nfac==6) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,1]*Xtemp[,6],Xtemp[,2]*Xtemp[,3],
                    Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,2]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,3]*Xtemp[,6],
                    Xtemp[,4]*Xtemp[,5],Xtemp[,4]*Xtemp[,6],Xtemp[,5]*Xtemp[,6],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,2]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,6],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,3]*Xtemp[,6],Xtemp[,1]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,1]*Xtemp[,4]*Xtemp[,6],Xtemp[,1]*Xtemp[,5]*Xtemp[,6],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,5],
                    Xtemp[,2]*Xtemp[,3]*Xtemp[,6],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],
                    Xtemp[,2]*Xtemp[,4]*Xtemp[,6],Xtemp[,2]*Xtemp[,5]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,4]*Xtemp[,5],Xtemp[,3]*Xtemp[,4]*Xtemp[,6],
                    Xtemp[,3]*Xtemp[,5]*Xtemp[,6],Xtemp[,4]*Xtemp[,5]*Xtemp[,6])
      }
    }
    
    ##### Define Models 5 and 6 here according to the order in MixModel
    if (mod==5 & nfac >5)
      stop("No more than 5 mixture components allowed with model 5")
    if (mod==6 & nfac >5)
      stop("No more than 5 mixture components allowed with model 6")
    if (mod==5 | mod==6) {
      Xmat<-Xtemp
    }
    ########### Model 5  #######################################################
    if (mod==5) {
      if (nfac==2) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2])            
      }
      if (nfac==3) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3])
      }
      if (nfac==4) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,4]*Xtemp[,5])
      }
      # Creates new beta vector at constant values of process variables
      betaT<-Beta[1:(nfac+(nfac*(nfac-1)/2))]
      # here are the constant values of the process variables
      pvslice<-pvslice[1:nproc]
      # modifies the coeficients for the linear mixture terms by adding effect of linear process variables
      i1<-nfac
      i2<-nfac*(nfac-1)/2
      i3<-nfac*nproc
      i4<-(nfac*(nfac-1)/2)*nproc
      i5<-(nproc*(nproc-1)/2)
      i6<-i2*i5
      strt1<-i1+i2+1
      stp1<-i1+i2+i3
      betaL<-matrix(Beta[(strt1):(stp1)],nrow=nfac,ncol=nproc)
      betates<-as.matrix(betaL,nrow=1) 
      if(length(pvslice)>1) {vecL<-pvslice%*%betaL}
      else {vecL<-pvslice*betaL}
      betaT[1:nfac]<-betaT[1:nfac]+vecL
      betates<-as.matrix(betaT,nrow=1)
      # modifies the coefficients for quadratic mixture terms by adding effect of linear process variables
      strt2<-i1+i2+i3+1
      stp2<-i1+i2+i3+i4
      betaQL<-matrix(Beta[strt2:stp2],nrow=nfac*((nfac-1)/2),ncol=nproc) 
      vecQ<-pvslice%*%betaQL
      betaT[(i1+1):(i1+i2)]<-betaT[(i1+1):(i1+i2)]+vecQ
      betates<-as.matrix(betaT,nrow=1)
      # modifies the coeficients for the linear mixture terms by adding effect of quadratic process variables
      strt3<-i1+i2+i3+i4+1
      stp3<-i1+i2+i3+i4+i5*i1
      pvsq<-c(pvslice[1]*pvslice[2],pvslice[1]*pvslice[3],pvslice[2]*pvslice[3])
      betates<-as.matrix(Beta,nrow=1)
      betaLQ<-matrix(Beta[strt3:stp3],nrow=nfac,ncol=nproc*(nproc-1)/2)
      vecLQ<-pvsq%*%betaLQ
      betaT[1:i1]<-betaT[1:i1]+vecLQ
      # modifies the coefficients for quadratic mixture terms by adding effect of quadratic process variables
      strt4<-stp3+1
      stp4<-stp3+i6
      betates<-as.matrix(Beta,nrow=1)
      betaQQ<-matrix(Beta[strt4:stp4],nrow=nfac*((nfac-1)/2),ncol=nproc*(nproc-1)/2)
      vecQQ<-pvsq%*%betaQQ
      betaT[(i1+1):(i1+i2)]<-betaT[(i1+1):(i1+i2)]+vecQQ
      betates<-as.matrix(betaT,nrow=1)
      ##############################################################################
    }
    
    if(mod<=3 ) {
      yhX<-Xmat%*%Beta
    }
    if(mod<=4 ) {
      yhX<-Xmat%*%Beta
    }
    if(mod==5 ) {
      yhX<-Xmat%*%betaT
    }  
    ############ Model 6 ############################################################
    if (mod==6) {
      if (nfac==2) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2])            
      }
      if (nfac==3) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3])
      }
      if (nfac==4) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],
                    Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],
                    Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,4]*Xtemp[,5])
      }
      # Creates new beta vector at constant values of process variables
      betaT<-Beta[c(1:nfac,(nfac+nproc+1):(nfac+nproc+nfac*(nfac-1)/2))]
      # here are the constant values of the process variables
      pvslice<-pvslice[1:nproc]
      # modifies the coeficients for the linear mixture terms by adding effect of linear process variables
      i1<-nfac
      i2<-nfac*(nfac-1)/2
      i3<-nfac*nproc
      i4<-(nfac*(nfac-1)/2)*nproc
      i5<-(nproc*(nproc-1)/2)
      i6<-i2*i5
      strt1<-i1+nproc+i2+1
      stp1<-i1+nproc+i2+i3
      betaL<-matrix(Beta[(strt1):(stp1)],nrow=nfac,ncol=nproc)
      betates<-as.matrix(betaL,nrow=1) 
      if(length(pvslice)>1) {vecL<-pvslice%*%betaL}
      else {vecL<-pvslice*betaL}
      betaT[1:nfac]<-betaT[1:nfac]+vecL
      ## calculates constant to add to predicted values
      st<-i1+1
      sp<-i1+nproc
      betaQ<-(Beta[st:sp])
      if(length(pvslice)>1) {const<-(pvslice^2)%*%betaQ}
      else {const<-pvslice*betaQ}
    }
    if(mod==6 ) {
      yhX<-Xmat%*%betaT+const
      
    }  
    
    
    # get deviations from centroid
    DX<-c(rep(0,times=(2*glimit+1)))
    for (i in 1:(2*glimit+1)) {
      DX[i]<-Xgrid[ifac,i,ifac]-cent[ifac]
    }
    PX[,(2*ifac)]<-yhX
    PX[,(2*ifac-1)]<-DX
  }
  
  # select ranges for each the lines for each factor
  sel<-c(rep(0,times=nfac*(2*glimit+1)))
  dim(sel)<-c((2*glimit+1),nfac)
  for (i in 1:nfac) {
    sel[,i]<-(ck[1,i] <= Xgrid[i,,i]& Xgrid[i,,i]<=ck[2,i])
  }
  
  
  # get min and max x and y -coordinates for plot
  xmin<-c(rep(0,times=nfac))
  ymin<-c(rep(0,times=nfac))
  xmax<-c(rep(0,times=nfac))
  ymax<-c(rep(0,times=nfac))
  for (ifac in 1:nfac) {
    w1<-PX[(sel[,ifac]==1),(2*ifac-1)]
    w2<-PX[(sel[,ifac]==1),(2*ifac)]
    xmin[ifac]<-min(w1)
    xmax[ifac]<-max(w1)
    ymin[ifac]<-min(w2)
    ymax[ifac]<-max(w2)
  }
  
  xaxismin<-min(xmin)
  xaxismax<-max(xmax)
  yaxismax<-max(ymax)
  yaxismin<-min(ymin)
  
  
  # makes the plot
  plabs<-c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12")
  xvec<-c(xaxismin,xaxismax)
  yvec<-c(yaxismin,yaxismax)
  plot(xvec,yvec,type="n",ylab="Predicted Response",xlab="Deviation from centroid",main=mTitle)
  for (i in 1:nfac) {
    xvec<-PX[(sel[,i]==1),(2*i-1)]
    yvec<-PX[(sel[,i]==1),(2*i)]
    lines(xvec,yvec,lty=i,col=i)
    text(xvec[length(xvec)],yvec[length(yvec)],plabs[i])
    text(xvec[1],yvec[1],plabs[i])
    
  }
  
  
  #  return(PX)
}
##############################################################
