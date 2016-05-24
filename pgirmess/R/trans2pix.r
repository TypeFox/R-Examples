trans2pix<-function(vect) {

# Giraudoux 24.1.2004
# Convert a transect coordinate file with some landmarks 
# and "NA" values in between into a matrix with 
# intermediate coordinates.
# The argument passed is a matrix or data.frame of two columns
# each row is a transect interval; each column must start (first row)
# and end (last row) with a landmark ; intermediate landmarks must have
# numbers in the two columns of the row. Other rows must be NA values
# The returned value is a matrix of the same dimensions

vect1<-vect[,1]
vect2<-vect[,2]

j<-0
deb1<-NULL;deb2<-NULL
vecfin1<-NULL;vecfin1<-NULL;vecfin2<-NULL
for (i in 1:length(vect1)){
    j<-j+1
    if (i==1){
        deb1<-vect1[i]
        deb2<-vect2[i]
        vecfin1<-deb1
        vecfin2<-deb2
        }
    if ((i!=1)&&(!is.na(vect1[i]))) {
        seg1<-seq(deb1,vect1[i],l=j)
        seg2<-seq(deb2,vect2[i],l=j)
        vecfin1<-c(vecfin1,seg1[2:length(seg1)])
        vecfin2<-c(vecfin2,seg2[2:length(seg2)])
        deb1<-vect1[i]
        deb2<-vect2[i]
        j<-1
        }
    }
cbind(long=vecfin1,lat=vecfin2)
}
