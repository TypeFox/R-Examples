xlogy <-
function(x,y){if(length(x)==1) x<-rep(x,length(y)); z<-x*log(y); z[y==0]<-x[y==0]*(-230);z}
