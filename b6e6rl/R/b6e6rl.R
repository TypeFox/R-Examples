b6e6rl <-
function(fn_name, a, b, N, my_eps, max_evals, n0, delta){
h <- 12
#initialization
d <- length(a)
cni <- rep(0,h)			
ni <- rep(0,h) + n0		
nrst <- 0
func_evals <- N
success <- 0
func_near <-0
P <- matrix(0,N,d+1)		
for (i in 1:N){
	P[i,1:d] <- a+(b-a)*runif(d)
	P[i,d+1] <- do.call(fn_name,list(P[i,1:d]))
}#0-th generation initialized
fmax <- max(P[,d+1])
fmin <- min(P[,d+1])
indmin <- which.min(P[,d+1])
Q <- P
p2 <- 0.5*(1+1/d)
p1 <- 0.5*(p2+1/d)
p3 <- 0.5*(p2+1)
CR1 <- CRexp_set(p1,d)
CR2 <- CRexp_set(p2,d)
CR3 <- CRexp_set(p3,d)
while ((fmax-fmin > my_eps) && (func_evals < d*max_evals)){	#main loop
	for (i in 1:N){
		hh <- roulete(ni)[1]		#select hh-th set of parameters
		p_min <- roulete(ni)[2]		
		if (p_min < delta){
			cni <- cni + ni - n0	
			ni <- rep(0,h) + n0	
			nrst <- nrst + 1
		}#reset
		switch(toString(hh),		#number of selected heuristics
			'1'={
				F <- 0.5
				CR <- 0
				y <- derand_RL(P,F,CR,i)				
			},
			'2'={
				F <- 0.5
				CR <- 0.5
				y <- derand_RL(P,F,CR,i)				
			},
			'3'={
				F <- 0.5
				CR <- 1
				y <- derand_RL(P,F,CR,i)				
			},
			'4'={
				F <- 0.8
				CR <- 0
				y <- derand_RL(P,F,CR,i)				
			},
			'5'={
				F <- 0.8
				CR <- 0.5
				y <- derand_RL(P,F,CR,i)				
			},
			'6'={
				F <- 0.8
				CR <- 1
				y <- derand_RL(P,F,CR,i)				
			},
			'7'={
				F <- 0.5
				CR <- CR1
				y <- derandexp_RL(P,F,CR,i)				
			},
			'8'={
				F <- 0.5
				CR <- CR2
				y <- derandexp_RL(P,F,CR,i)				
			},
			'9'={
				F <- 0.5
				CR <- CR3
				y <- derandexp_RL(P,F,CR,i)				
			},
			'10'={
				F <- 0.8
				CR <- CR1
				y <- derandexp_RL(P,F,CR,i)				
			},
			'11'={
				F <- 0.8
				CR <- CR2
				y <- derandexp_RL(P,F,CR,i)				
			},
			'12'={
				F <- 0.8
				CR <- CR3
				y <- derandexp_RL(P,F,CR,i)				
			}
		)
		y <- zrcad(y,a,b)
		fy <- do.call(fn_name, list(y))
		func_evals <- func_evals+1
		if (fy < P[i,d+1]){#trial point y is good for renewing population
			Q[i,] <- c(y,fy)			
			success <- success+1
			hh <- as.integer(hh)
			ni[hh] <- ni[hh]+1	#zmena prsti qi	
		}
	}#end of generation
	P <- Q
	fmax <- max(P[,d+1])
	fmin <- min(P[,d+1])
	indmin <- which.min(P[,d+1])	
}# end of main loop
cni <- cni + ni - n0			
x_star <- P[indmin,1:d]
fn_star <- fmin
result <- list(x_star=x_star,fn_star=fn_star,fn_name=fn_name,func_evals=func_evals,success=success,nrst=nrst,cni=cni)

return(result)
}
