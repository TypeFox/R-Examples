touchi<-function(rec1,rec2,rho=0)
{
#Checks whether rectangles rec1, rec2 touch.
#rec1,rec2 are 2*d vectors, discrete rectangles (grid)

#Returns 0 if intersection is empty

d<-length(rec1)/2
if (length(rho)==1) rho<-rep(rho,d)

tulos<-1
i<-1
while ((i<=d) && (tulos==1)){  
    ala<-max(rec1[2*i-1],rec2[2*i-1])
    yla<-min(rec1[2*i],rec2[2*i])
    if (yla+2*rho[i]<ala) tulos<-0
    i<-i+1
}
return(tulos)
}




