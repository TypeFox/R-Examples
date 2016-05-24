circlePoint <- function(circle, T){
	if(is.vector(T) && length(T) > 1){
		r <- matrix(NA, nrow=length(T), ncol=length(circle$C))
		for(i in 1:length(T)) r[i, ] <- circle$C + circle$R*cos(T[i])*circle$U + circle$R*sin(T[i])*circle$V		
	}else{
		r <- circle$C + circle$R*cos(T)*circle$U + circle$R*sin(T)*circle$V
	}
	r
}