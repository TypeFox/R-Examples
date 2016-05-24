generate.benchmark.data <-
function(typ #choice of dependence to generate, number between 1 and 8
		,noises #vector of noises
		,n #sample size
		#data for torus projection
		,project=FALSE #project onto torus?
		,windx = 1 #how many windings in x-direction
		,windy = 1 #how many windings in y-direction
){
	
	x= replicate(length(noises),runif(n))
	
	#lin+noise
	if(typ==1){
		y=x+ matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#parabolic+noise
	if(typ==2){
		y=4*(x-.5)^2+  matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#cubic+noise
	if(typ==3){
		y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10* matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#sin+noise
	if(typ==4){
		y=sin(4*pi*x) + 2*matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#their sine + noise
	if(typ==5){
		y=sin(16*pi*x) + matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#x^(1/4) + noise
	if(typ==6){
		y=x^(1/4) + matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#circle
	if(typ==7){
		y=(2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2)) + 3/4*matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	#step function
	if(typ==8){
		y = (x > 0.5) + 5*matrix(noises * rnorm(n*length(noises)),ncol=length(noises),byrow=T)
	}
	if(project){
		proj_torus = function(dat,wind){
			return((dat*wind) %% 1)
		}
		x = proj_torus(x,windx)
		y = proj_torus(y,windy)
	}
	return(list(x=x,y=y))
}
