leikkaa<-function(rec1,rec2){
#Makes an intersection of rectangles rec1, rec2
#rec1,rec2 are 2*d vectors
#
#Returns 2*d-vector or NA if intersection is empty
#
d<-length(rec1)/2
tulos<-matrix(0,2*d,1)
i<-1
while ((i<=d) && (!is.na(tulos))){  
    tulos[2*i-1]<-max(rec1[2*i-1],rec2[2*i-1])
    tulos[2*i]<-min(rec1[2*i],rec2[2*i])
    if (tulos[2*i]<=tulos[2*i-1]) tulos<-NA
    i<-i+1
}
return(tulos)
}

