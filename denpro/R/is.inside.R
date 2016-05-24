is.inside<-function(simp1,simp2)
{
# simp1, simp2 (d+1)*d matrices
# returns 1 if simp1 is inside simp2

d<-2  #dim(simp1)[2]

deet<-matrix(0,3,1)
lk<-1
for (ii in 1:(d+1)){
    v1<-simp2[ii,]
    jj<-ii+1
    while (jj<=(d+1)){
         v2<-simp2[jj,]
         deet[lk]<-sqrt( sum((v1-v2)^2) )
         jj<-jj+1
         lk<-lk+1
    }
}
rho<-max(deet)

tulos<-1
i<-1
while ( (i<=(d+1)) && (tulos==1) ){
    vertice1<-simp1[i,]
    j<-1
    while ( (j<=(d+1)) && (tulos==1) ){
        vertice2<-simp2[j,]
        eta<-sqrt( sum((vertice1-vertice2)^2) )
        if (eta>rho) tulos<-0
        j<-j+1
    }
    i<-i+1
}

return(tulos)
}


