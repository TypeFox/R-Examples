EffPlot = function(des=NULL,nfac=3,mod=1,dir=1)
{
  # glimit is the resolution for the plot
  glimit<-25
  
  # check for valid mod
  if(mod < 1 | mod == 3| mod > 4) 
    stop("mod must be one of the following: 1 = linear, 2 = quadratic, or 4 = Special Cubic")
  
  # get the number of factors from the design
  if(! is.null(des)) {
    dd<-dim(des)
    nfac<-dd[2]-1
  }
  
  # checks for valid number of factors
  
  
  if (nfac<2 | nfac>12)
    stop("The number factors must be between 2 and 12")
  if (nfac>=7 & mod>=2)
    stop("Linear models only when the number of factors is greater than 6")
  if (nfac<3 & mod==4)
    stop("Special cubic model requires at least 3 factors")
  
  # initialize lower and upper bounds
  x1<-c(0,1)
  x2<-c(0,1)
  x3<-c(0,1)
  x4<-c(0,1)
  x5<-c(0,1)
  x6<-c(0,1)
  x7<-c(0,1)
  x8<-c(0,1)
  x9<-c(0,1)
  x10<-c(0,1)
  x11<-c(0,1)
  x12<-c(0,1)
  
  
  # set up lower and upper bounds for factors if design is in des
  if(! is.null(des)) {
    if (nfac>=1) {
      x1[1]<-min(des[,1])
      x1[2]<-max(des[,1]) 
    }
    if (nfac>=2) {
      x2[1]<-min(des[,2])
      x2[2]<-max(des[,2]) 
    }
    if (nfac>=3) {
      x3[1]<-min(des[,3])
      x3[2]<-max(des[,3]) 
    }
    if (nfac>=4) {
      x4[1]<-min(des[,4])
      x4[2]<-max(des[,4]) 
    }
    if (nfac>=5) {
      x5[1]<-min(des[,5])
      x5[2]<-max(des[,5]) 
    }
    if (nfac>=6) {
      x6[1]<-min(des[,6])
      x6[2]<-max(des[,6]) 
    }
    if (nfac>=7) {
      x7[1]<-min(des[,7])
      x7[2]<-max(des[,7]) 
    }
    if (nfac>=8) {
      x8[1]<-min(des[,8])
      x8[2]<-max(des[,8]) 
    }
    if (nfac>=9) {
      x9[1]<-min(des[,9])
      x9[2]<-max(des[,9]) 
    }
    if (nfac>=10) {
      x10[1]<-min(des[,10])
      x10[2]<-max(des[,10]) 
    }
    if (nfac>=11) {
      x11[1]<-min(des[,11])
      x11[2]<-max(des[,11]) 
    }
    if (nfac>=12) {
      x12[1]<-min(des[,12])
      x12[2]<-max(des[,12]) 
    }
    
    
    # get the beta coefficients
    if (mod==1) {
      modl<-lm(y~.-1,data=des)
      Beta<-modl$coeff
    }
    
    if (mod==2) {
      
      if (nfac==2) {
        modl<-lm(y~x1+x2+x1*x2-1,data=des)
        Beta<-modl$coeff
      }
      
      if (nfac==3)  {
        modl<-lm(y~x1+x2+x3+x1*x2+x1*x3+x2*x3-1,data=des)
        Beta<-modl$coeff
      }
      if (nfac==4)  {
        modl<-lm(y~x1+x2+x3+x4+x1*x2+x1*x3+x1*x4+x2*x3+x2*x4+x3*x4-1,data=des) 
        Beta<-modl$coeff
        
      }
      if (nfac==5)  {
        modl<-lm(y~x1+x2+x3+x4+x5+x1*x2+x1*x3+x1*x4+x1*x5+x2*x3+x2*x4+x2*x5+x3*x4+x3*x5+x4*x5-1,data=des)
        Beta<-modl$coeff 
      }
      if (nfac==6)  {
        modl<-lm(y~x1+x2+x3+x4+x5+x6+x1*x2+x1*x3+x1*x4+x1*x5+x1*x6+x2*x3+x2*x4+x2*x5+x2*x6+x3*x4+x3*x5+x3*x6+x4*x5+x4*x6+x5*x6-1,data=des)
        Beta<-modl$coeff 
      }
    }
    
    if (mod==4) {
      if (nfac==3)  {
        modl<-lm(y~x1+x2+x3+x1*x2+x1*x3+x2*x3+x1*x2*x3-1,data=des)
        Beta<-modl$coeff
      }
      if (nfac==4)  {
        modl<-lm(y~x1+x2+x3+x4+x1*x2+x1*x3+x1*x4+x2*x3+x2*x4+x3*x4+x1*x2*x3+x1*x2*x4+x1*x3*x4+x2*x3*x4-1,data=des) 
        Beta<-modl$coeff
      }
      if (nfac==5)  {
        modl<-lm(y~x1+x2+x3+x4+x5+x1*x2+x1*x3+x1*x4+x1*x5+x2*x3+x2*x4+x2*x5+x3*x4+x3*x5+x4*x5+x1*x2*x3+x1*x2*x4+x1*x2*x5+x1*x3*x4+x1*x3*x5+x1*x4*x5+x2*x3*x4+x2*x3*x5+x2*x4*x5+x3*x4*x5-1,data=des)
        Beta<-modl$coeff 
      }
      if (nfac==6)  {
        modl<-lm(y~x1+x2+x3+x4+x5+x6+x1*x2+x1*x3+x1*x4+x1*x5+x1*x6+x2*x3+x2*x4+x2*x5+x2*x6+x3*x4+x3*x5+x3*x6+x4*x5+x4*x6+x5*x6+x1*x2*x3+x1*x2*x4+x1*x2*x5+x1*x2*x6+x1*x3*x4+x1*x3*x5+x1*x3*x6+x1*x4*x5+x1*x4*x6+x1*x5*x6+x2*x3*x4+x2*x3*x5+x2*x3*x6+x2*x4*x5+x2*x4*x6+x2*x5*x6+x3*x4*x5+x3*x4*x6+x3*x5*x6+x4*x5*x6-1,data=des)
        Beta<-modl$coeff 
      }
    }   
    
    ## end of calculations involving the data frame des
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
    if (mod==1 | mod==4) {
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
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,4]*Xtemp[,5])
      }
      if (nfac==6) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,1]*Xtemp[,5],Xtemp[,1]*Xtemp[,6],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,2]*Xtemp[,6],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,3]*Xtemp[,6],Xtemp[,4]*Xtemp[,5],Xtemp[,4]*Xtemp[,6],Xtemp[,5]*Xtemp[,6])
      }
    }
    
    
    if (mod==4) {
      if (nfac==3) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,3])
      }
      if (nfac==4) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,4])
      }
      if (nfac==5) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,1]*Xtemp[,5],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,4]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],Xtemp[,1]*Xtemp[,2]*Xtemp[,5],Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],Xtemp[,1]*Xtemp[,4]*Xtemp[,5],Xtemp[,2]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,5],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],Xtemp[,3]*Xtemp[,4]*Xtemp[,5])
      }
      if (nfac==6) {
        Xmat<-cbind(Xtemp,Xtemp[,1]*Xtemp[,2],Xtemp[,1]*Xtemp[,3],Xtemp[,1]*Xtemp[,4],Xtemp[,1]*Xtemp[,5],Xtemp[,1]*Xtemp[,6],Xtemp[,2]*Xtemp[,3],Xtemp[,2]*Xtemp[,4],Xtemp[,2]*Xtemp[,5],Xtemp[,2]*Xtemp[,6],Xtemp[,3]*Xtemp[,4],Xtemp[,3]*Xtemp[,5],Xtemp[,3]*Xtemp[,6],Xtemp[,4]*Xtemp[,5],Xtemp[,4]*Xtemp[,6],Xtemp[,5]*Xtemp[,6],Xtemp[,1]*Xtemp[,2]*Xtemp[,3],Xtemp[,1]*Xtemp[,2]*Xtemp[,4],Xtemp[,1]*Xtemp[,2]*Xtemp[,5],Xtemp[,1]*Xtemp[,2]*Xtemp[,6],Xtemp[,1]*Xtemp[,3]*Xtemp[,4],Xtemp[,1]*Xtemp[,3]*Xtemp[,5],Xtemp[,1]*Xtemp[,3]*Xtemp[,6],Xtemp[,1]*Xtemp[,4]*Xtemp[,5],Xtemp[,1]*Xtemp[,4]*Xtemp[,6],Xtemp[,1]*Xtemp[,5]*Xtemp[,6],Xtemp[,2]*Xtemp[,3]*Xtemp[,4],Xtemp[,2]*Xtemp[,3]*Xtemp[,5],Xtemp[,2]*Xtemp[,3]*Xtemp[,6],Xtemp[,2]*Xtemp[,4]*Xtemp[,5],Xtemp[,2]*Xtemp[,4]*Xtemp[,6],Xtemp[,2]*Xtemp[,5]*Xtemp[,6],Xtemp[,3]*Xtemp[,4]*Xtemp[,5],Xtemp[,3]*Xtemp[,4]*Xtemp[,6],Xtemp[,3]*Xtemp[,5]*Xtemp[,6],Xtemp[,4]*Xtemp[,5]*Xtemp[,6])
      }
    }
    
    if(mod <=2 | mod == 4) {
      yhX<-Xmat%*%Beta
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
  
  
    return(PX)
}
##############################################################

