simesi3d<-function(n,p1,p2,p3,s1,s2,s3,siemen,noisedim=0){
#muodostaa n*2 havaintomatriisin
#
#p1+p2+p3+p4=1
#
#tiheysfunktio on paloittain vakio 4 palassa:
#I: [0,split1]*[0,1]*[0,1]
#II: [split1,1]*[0,split2]*[0,1] 
#III: [split1,1]*[split2,1]*[0,split3]
#IV: [split1,1]*[split2,1]*[split3,1]
#
#p1 on palan I tn. ja p2 on palan II tn.,...
#p1=f1*s1, p2=(1-s1)*s2*f2, p3=(1-s1)*(1-s2)*s3*f3, 
#p4=(1-s1)*(1-s2)*(1-s3)*f4
#kokeillaan esim. p1=0.1, p2=0.2, p3=0.3, p4=0.4 
#
set.seed(siemen)
d<-noisedim+3
x<-matrix(0,n,d)     #havaintomatriisi
i<-1
while (i<=n){         
  U<-runif(1)        #arpoo mihin palaan havainto tulee
  U1<-runif(1)
  U2<-runif(1)
  U3<-runif(1)
  if (U<p1){         #tn:lla p1 palaan I
    x[i,1]<-s1*U1    #x-koord skaalataan valiin [0,split1]
    x[i,2]<-U2       #y-koordinaatti valiin [0,1]
    x[i,3]<-U3       #z-koordinaatti valiin [0,1]
  }
  else{ 
    if ((U>=p1) && (U<p1+p2)){     #pala II
       x[i,1]<-s1+(1-s1)*U1  #x-koord skaalataan valiin [split1,1]
       x[i,2]<-s2*U2         #y-koord skaalataan valiin [0,split2]
       x[i,3]<-U3            #z-koordinaatti valiin [0,1] 
    }
    else{
       if ((U>=p1+p2) && (U<p1+p2+p3)){    #pala III
          x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
          x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]
          x[i,3]<-s3*U3          #z-koord skaalataan valiin [0,split3]
       }
       else{
          x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
          x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]
          x[i,3]<-s3+(1-s3)*U3   #z-koord skaalataan valiin [split3,1]
       }
    }
  }
  if (noisedim>0){
     nd<-3
     while (nd<=(3+noisedim)){
         x[i,nd]<-runif(1)
         nd<-nd+1
     }
  }
i<-i+1
}
return(x)
}




