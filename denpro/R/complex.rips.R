complex.rips<-function(dendat,rho)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

lkm<-0
complex<-matrix(0,1,d+1)

nn<-nn.indit(dendat)
maxk<-n-1
nnr<-nn.radit(dendat,maxk)

for (i in 1:n){
   rcur<-nnr[i,1]
   j<-1
   while ( (rcur<=rho) && (j<=(n-1)) ){
        cind<-nn[i,j]

        if (cind>i){
        # find the connection
        rcur1<-nnr[i,j]
        rcur2<-nnr[cind,1]
        kk<-j+1
        found<-FALSE
        while ( (rcur1<=rho) && (kk<=(n-1)) && (!found) ){
           koe1<-nn[i,kk]
           ll<-1
           while ( (rcur2<=rho) && (ll<=(n-1)) && (!found) ){
                koe2<-nn[cind,ll]
                if (koe1==koe2){
                     found<-TRUE
                     addi<-matrix(c(i,cind,koe1),1,3)
                     if (lkm==0) complex<-matrix(c(i,cind,koe1),1,3)
                     else complex<-rbind(complex,addi)
                     lkm<-lkm+1 
                } 
                ll<-ll+1 
                if (ll<=(n-1)) rcur2<-nnr[cind,ll]     
           }
           kk<-kk+1
           if (kk<=(n-1)) rcur1<-nnr[i,kk]
        }        
        # connection search end
        }

        j<-j+1
        if (j<=(n-1)) rcur<-nnr[i,j]
   }
}

if (lkm==0) complex<-NULL
return(complex)
}

