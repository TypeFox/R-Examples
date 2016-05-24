# Description: 	Implementation of hit.run algorithm, p. 377
# Usage:	hit.run(theta.mat,reps,I.mat)



hit.run <- function(theta.mat,reps,I.mat)  {
        		    for (i in 2:reps)  {
			        u.vec <- c(runif(1,0,pi/2),runif(1,pi/2,pi),
			                   runif(1,pi,3*pi/2),
			                   runif(1,3*pi/2,2*pi))
            			u.dr <- sample(u.vec,size=1,
			                   prob=c(1/3,1/6,1/3,1/6))
            			g.ds <- rgamma(1,1,1)  
            			xy.theta <- c(g.ds*cos(u.dr),g.ds*sin(u.dr)) 
			                    + theta.mat[(i-1),]
            			a <- dmultinorm(xy.theta[1],xy.theta[2],
			                        c(0,0),I.mat)/
                 		     dmultinorm(theta.mat[(i-1),1],
			                        theta.mat[(i-1),2],c(0,0),I.mat)
            			r.uniform <- runif(1)
            			if (a > r.uniform) theta.mat[i,] <- xy.theta
            			else theta.mat[i,] <- theta.mat[(i-1),]
        		    }
			    theta.mat
			}

