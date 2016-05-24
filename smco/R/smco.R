####################################################################################
#                                                                                  #   
# smco                                                                             # 
#                                                                                  #
# A simple Monte Carlo optimizer using adaptive coordinate sampling                #
# Prof. Juan David Velasquez H.                                                    #
#                                                                                  #
####################################################################################
smco <-
function(par = NULL, fn, gr = NULL, ..., N = length(par), LB, UB, maxiter = 1000, 
	Co = 0.01, Cmin = 0.0001, Cmax = 0.5, trc = FALSE, lambda = 20,
    useBFGS = FALSE, control = list(), hessian = FALSE)
{
    #-------------------------------------------------------------------------------
    #
    # Local improvement using gradient-based optimization
    # 
    runBFGS <- function( u )
    {
        z.curr[idim.curr] = u
        x.curr = LB + z.curr * (UB - LB)
		return (fn( x.curr ))	
    }
    #-------------------------------------------------------------------------------
    #
    # Optimization algorithm
    #
    if (trc) { print( match.call() ) }
    #
    #-------------------------------------------------------------------------------
    #
    if(is.null(par)) { z.curr = runif(N) } else	{ z.curr = (par - LB) / (UB - LB) }
    #
	z.min       = z.curr                   # optimo local  en [0, 1]^n 
	z.opt       = z.curr
	x.min       = LB + z.curr * (UB - LB)  # optimo local  en [L, U]^n 
    #    
    f.min       = rep(0, maxiter);  f.min[1]  = fn( x.min ) # f evaluada en el optimo local	    
	f.curr      = rep(0, maxiter);  f.curr[1] = f.min[1]    # f evaluada en el punto actual	
    f.opt       = rep(0, maxiter);  f.opt[1]  = f.min[1]    # f evaluada en el optimo global
	#
	idim.curr   = 1
	bad.counter = 0    
	C           = Co
	#
	#-------------------------------------------------------------------------------
    #
	# Cycle
	#
	for (iter in 2:maxiter) 
	{
		z.curr = z.min
        u = z.min[idim.curr]			
        C =  max(Cmin, min(C * exp( 0.0001 + rnorm(1) / sqrt(2 * sqrt(N)+ 2 * N) ), Cmax))
        z.curr[idim.curr] = u + C * qnorm(runif(1) * (pnorm((1 - u) / C) - pnorm((0 - u) / C)) + pnorm((0 - u) / C))        
		#
		x.curr = LB + z.curr * (UB - LB)
		f.curr[iter] = fn( x.curr )	
		#		
        #--------------------------------------------------------------------------
        #
        if (useBFGS) 
        {
            result = optim(par = z.curr[idim.curr], fn = runBFGS, 
                gr = NULL, method='L-BFGS-B', lower = 0, upper = 1, 
                control = control, hessian = hessian)
            #
            if (result$value < f.curr[iter])
            {
                z.curr[idim.curr] = result$par
                x.curr = LB + z.curr * (UB - LB)
                f.curr[iter] = fn( x.curr )                 
            }
        }
        #
        #--------------------------------------------------------------------------
        #
		if (f.curr[iter] < f.min[iter - 1]) 
		{						
			x.min       = x.curr
			f.min[iter] = f.curr[iter]
			z.min       = z.curr
            bad.counter = 0
		}		
		else 
		{
			f.min[iter] = f.min[iter - 1]			
            bad.counter = bad.counter + 1
		}	
        #
        #--------------------------------------------------------------------------
		#
		# Restart
		#
        if (bad.counter > lambda * N)      
        {
			for (idim.curr in 1:N) {
				u = z.min[idim.curr]		
				z.min[idim.curr] = u + Co * qnorm(runif(1) * (pnorm((1 - u) / Co) - pnorm((0 - u) / Cmax)) + pnorm((0 - u) / Co))
			}
            z.min = runif(N)             
            f.curr[iter] = fn( LB + z.min * (UB - LB) )	     
            f.min[iter] = f.curr[iter]
            bad.counter = 0            
            C = Co
            if (trc) cat("Restarting ...", "\n")
        }
        #
        #--------------------------------------------------------------------------
        #
        if (f.curr[iter] < f.opt[iter - 1]) 
		{						
			x.opt       = x.curr
			f.opt[iter] = f.curr[iter]
			z.opt       = z.curr            
		}		
		else 
		{
			f.opt[iter] = f.opt[iter - 1]			            
		}	
        #
        idim.curr = idim.curr + 1 
        if (idim.curr > N) 
        {
            idim.curr = 1    			            		
            if (trc) 
            {
                cat(iter, f.curr[iter], f.min[iter], f.opt[iter], "\n")
            }
        }
        #
    }
	#
	return (list(par = x.opt, value = f.opt[maxiter], 
            f.opt = f.opt, f.min = f.min, f.curr=f.curr,
            call = match.call()))
	#
}



#
###########################################################
#                                                         #
#                 Load library function                   #
#                                                         #
###########################################################
#
.onAttach <- 
function(...)
{
	version = library(help = smco)$info[[1]] 
	version = version[pmatch("Version",version)]
	um = strsplit(version, " ")[[1]]
	version = um[nchar(um) > 0][2]
	#
	cat(paste("This is smco package", version, "\n"))
}