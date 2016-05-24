runs <- 8; size <- 40
time <- replicate(runs, cumsum(rexp(size)))
df <- data.frame(time = as.vector(time), run = rep(1:runs,each=size) )
stop <- min(apply(time,2,max)) 
stop <- 5 * trunc(stop/5)
df <- df[time <= stop,]
myplot <- stripplot(run~time, df, pch=1, cex=.7, col='black',
	panel=function(x,y,...){
		panel.abline(h=seq(1.5,7.5,by=1),col='gray60') 
		panel.abline(v=seq(0,stop,by=5), col='gray60') 
		panel.stripplot(x,y,...)
	})
