"shift" <-
function(v,places,dir="right"){

vnew<-NULL

d<-substring(dir,1,1)

if(d=="r" & places<0 ){
d<-"l"}
else {if (d=="l" & places<0){d<-"r"}}


n<-length(v)

p<-abs(places)

if (p==0){vnew<-v}
else{
if (d=="r"){
vnew<-c(v[(n-p+1):n],v[1:(n-p)])
}
else{
vnew<-c(v[(p+1):n],v[1:p])
}
}

vnew
}

