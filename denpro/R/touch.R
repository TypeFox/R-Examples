touch<-function(rec1,rec2,epsi=0.000001)
{
# Checks whether rectangles rec1, rec2 touch.
# rec1,rec2 are 2*d vectors

# Returns FALSE if intersection is empty

d<-length(rec1)/2
tulos<-TRUE
i<-1
while ((i<=d) && (tulos)){  
    ala<-max(rec1[2*i-1],rec2[2*i-1])
    yla<-min(rec1[2*i],rec2[2*i])
    if (yla+epsi<ala) tulos<-FALSE
    i<-i+1
}
return(tulos)
}




