simesi<-function(n,p1,p2,s1,s2,siemen,noisedim=0){
#muodostaa n*2 havaintomatriisin
#
#p1+p2+p3=1
#split1, split2 in (0,1)
#
#tiheysfunktio on paloittain vakio 3 palassa:
#I: [0,split1]*[0,1], II: [split1,1]*[0,split2], III: [split1,1]*[split2,1]
#p1 on palan I tn. ja p2 on palan II tn., 
#p1=f1*s1, p2=(1-s1)*s2*f2, p3=(1-s1)*(1-s2)*f3 
#kokeillaan esim. p1=0.1, p2=0.3, p3=0.6 
#
#tulos: 
#val: 0.5 NA 0.5 NA NA, 
#vec: 1   NA 2   NA NA,
#mean: havlkm/(n*tilavuus), likimain: 1  2*p1 2*(p2+p3) 4*p2 4*p3 
#nelem:  n  p1*n  (p1+p2)*n  p2*n  p3*n
#ssr: havlkm*log(mean) = havlkm*log(havlkm/(n*tilavuus)), 
#    likimain: n*log(1)=0      p1*n*log(2*p1)   (p1+p2)*n*log(2*(p2+p3))
#              p2*n*log(4*p2)  p3*n*log(4*p3) 
#
#p1<-0.1
#p2<-0.3
#p3<-0.6
#s1<-0.75
#s2<-0.5
#values<-c(p1/s1,p2/((1-s1)*s2),p3/((1-s1)*(1-s2)))
#recs<-matrix(0,3,4)
#recs[1,]<-c(0,s1,0,1)
#recs[2,]<-c(s1,1,0,s2)
#recs[3,]<-c(s1,1,s2,1)
#koe<-drawgene(values,recs,plkm=30)
#persp(koe$x,koe$y,koe$z,phi=30,theta=60) 
#
set.seed(siemen)
d<-noisedim+2
x<-matrix(0,n,d)     #havaintomatriisi
i<-1
while (i<=n){         
  U<-runif(1)        #arpoo mihin palaan havainto tulee
  U1<-runif(1)
  U2<-runif(1)
  if (U<p1){         #tn:lla p1 palaan I
    x[i,1]<-s1*U1   #x-koord skaalataan valiin [0,split1]
    x[i,2]<-U2       #y-koordinaatti valiin [0,1]
  }
  else{ if (U>p1+p2){          #tn:lla p3 palaan III
            x[i,1]<-s1+(1-s1)*U1   #x-koord skaalataan valiin [split1,1]  
            x[i,2]<-s2+(1-s2)*U2   #y-koord skaalataan valiin [split2,1]

          }
           else{                 #tn:lla p2 palaan II
             x[i,1]<-s1+(1-s1)*U1  #x-koord skaalataan valiin [split1,1]
             x[i,2]<-s2*U2      #y-koord skaalataan valiin [0,split2]
           }
  }
  if (noisedim>0){
     nd<-3
     while (nd<=(2+noisedim)){
         x[i,nd]<-runif(1)
         nd<-nd+1
     }
  }
i<-i+1
}
return(x)
}




