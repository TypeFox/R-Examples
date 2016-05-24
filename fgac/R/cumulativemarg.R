"cumulativemarg" <-
function(cumulative,x,a){
	if(missing(a)){resp<-cumulative(x)}
	else{
		d<-length(a)
		if(d==1){resp<-cumulative(x,a[1])}
		else{
		if(d==2){resp<-cumulative(x,a[1],a[2])}
		else{
		if(d==3){resp<-cumulative(x,a[1],a[2],a[3])}	
			}
		}
		}
		}

