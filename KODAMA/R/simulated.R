

spirals = function (n=c(100,100,100),sd=c(0,0,0)){
   clusters=length(n)
   x=NULL
   y=NULL
   for(i in 1:clusters){
      t=seq(1/(4*pi),1,length.out=n[i])^0.5*2*pi
      a=rnorm(n[i],sd=sd[i])
      x=c(x,cos(t+(2*pi*i)/clusters)*(t+a))
      y=c(y,sin(t+(2*pi*i)/clusters)*(t+a))
   }
   cbind(x,y)   
}



dinisurface = function(N=1000){
	u=sort(runif(N)*4*pi)	
	v=runif(N)
	a=1
	b=0.2
	x=a*cos(u)*sin(v)
	y=a*sin(u)*sin(v)
	z=a*(cos(v)+log(tan(v/2)))+b*u
	data=cbind(x,y,z)
	return(data)
}



swissroll = function(N=1000){
	n <- 3
    m <- 2
    tt <- sort((3 * pi/2) * (1 + 2 * runif(N)))
    height <- 21 * runif(N)
    x <- tt * cos(tt)
    y <- height
    z <- tt * sin(tt)
	data=cbind(x,y,z)
	return(data)
}




helicoid = function(N=1000){
	a=1
	p=sample((seq(1,-1,length.out=N)))
	t=seq(-pi,pi,length.out=N)
	x=p*cos(a*t)
	y=p*sin(a*t)
	z=t
	data=cbind(x,y,z)
	return(data)
}







