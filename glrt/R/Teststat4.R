Teststat4 <-
function(data0, n0, est0, cens0, data1, n1, est1, cens1, rho = 0, gamma = 0, c0 = 1)
{
    k0 = ncol(est0$intmap)
    k1 = ncol(est1$intmap)
    F0 = est0$sigma
    F1 = est1$sigma
    FF0 = c(0, F0)
    FF1 = c(0, F1)
    indices0 = est0$ppairs
    indices1 = est1$ppairs
    alpha0cross = AlphaCross(data0[,1], data0[,2], n0, est1, k1)
    alpha1cross = AlphaCross(data1[,1], data1[,2], n1, est0, k0)

    indices0cross = StartEnd(alpha0cross) #for sample 0 based estimated F1
    indices1cross = StartEnd(alpha1cross) #for sample 1 based estimated F0

    u = rep(0, 2)
    tiny = .Machine$double.eps
    for (i in 1:n0) 
    {
        if (cens0[i] == 1) { #LC
	    if(!is.na(indices0cross[i,1])) # only LC obs that is at risk at > 0 interval of intmap 1; o/w contribution to u[1] is zero
	    {
                F0u = F0[indices0[i, 2]]
    	        F1u = F1[indices0cross[i, 2]]
                if (F0u > 0 + tiny) 
                    u[1] = u[1] + (Linkfunc(F1u, rho, gamma) - c0)/F0u
	    }
        }
        else if (cens0[i] == 2) { #IC
	    if(!is.na(indices0cross[i,1])) # only IC obs that is at risk at > 0interval of intmap 1	    
	    {
                F0u = FF0[indices0[i, 1]]
                F0v = F0[indices0[i, 2]]
	        	F1u = FF1[indices0cross[i, 1]]
                F1v =  F1[indices0cross[i, 2]]
                if (F0v - F0u > 0 + tiny) 
                    u[1] = u[1] + (Linkfunc(F1v, rho, gamma) - Linkfunc(F1u, rho, gamma))/(F0v - F0u)
	    }
        }
        else { #RC
	    if(!is.na(indices0cross[i,1])) # only RC obs that is at risk at > 0interval of intmap 1
	    {	
                F0v = FF0[indices0[i, 1]]
	        	F1v = FF1[indices0cross[i, 1]]
                if (F0v < 1 - tiny) 
                    u[1] = u[1] + (c0 - Linkfunc(F1v, rho, gamma))/(1 - F0v)
	    }
        }
    }
    u[1] = u[1] / sqrt(n0)

    for (i in 1:n1) 
    {
        if (cens1[i] == 1) { #LC
	    if(!is.na(indices1cross[i,1])) # only LC obs that is at risk at > 0 interval of intmap 0; o/w contribution to u[1] is zero
	    {
                F1u = F1[indices1[i, 2]]
	        F0u = F0[indices1cross[i, 2]]
                if (F1u > 0 + tiny) 
                    u[2] = u[2] + (Linkfunc(F0u, rho, gamma) - c0)/F1u
	    }
        }
        else if (cens1[i] == 2) { #IC
	    if(!is.na(indices1cross[i,1])) # only IC obs that is at risk at > 0interval of intmap 0
            {
				F1u = FF1[indices1[i, 1]]
				F1v = F1[indices1[i, 2]]
	    		F0u = FF0[indices1cross[i, 1]]
            	F0v =  F0[indices1cross[i, 2]]
            	if (F1v - F1u > 0 + tiny) 
                    u[2] = u[2] + (Linkfunc(F0v, rho, gamma) - Linkfunc(F0u, rho, gamma))/(F1v - F1u)
	    }
        }
        else { #RC
	    if(!is.na(indices1cross[i,1])) # only RC obs that is at risk at > 0interval of intmap 1
	    {
                F1v = FF1[indices1[i, 1]]
	        F0v = FF0[indices1cross[i, 1]]
                if (F1v < 1 - tiny) 
                    u[2] = u[2] + (c0 - Linkfunc(F0v, rho, gamma))/(1 - F1v)
	    }
        }
    }
    u[2] = u[2] / sqrt(n1)
    u
}