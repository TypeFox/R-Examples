
###############
# Basic Brownian motion tree likelihood for any number of independent dimensions

bm_loglik_tree = function(tree, values, par, dimen) {
  
  A = par[1:dimen]  	
  Sig = par[(dimen + 1):(dimen*2)] 	
  if (any(Sig < 0)) return(-Inf)
  
  n = length(tree$tip.label)
  V = vcv(tree)  	    	
  
  res = mapply(function(values, A, Sig) dmvnorm(t(as.matrix(values)), mean = rep(A, n), sigma = V*Sig, 1), values, A, Sig)
  
  return(sum(res))
}


#--------------
# Optimize
# 'values' has to be a list, with each element being one-dimensional values

point.like.bm = function(tree, values, start_values = NA, dimen = NA) {
	
	if (is.na(dimen)) dimen = length(values)
	if (is.na(start_values)) start_values = unlist(c(lapply(values, mean), lapply(values, sd)))	
	
	opt.res = nlm(function(p) -1*bm_loglik_tree(tree, values, p, dimen), p = start_values)	
	
	return(list(mrcas = opt.res$estimate[1:dimen],
		 rates = opt.res$estimate[(dimen + 1):(dimen*2)],
		 nlm.details = opt.res))
}


###############
# Rectangle constrained Brownian motion tree likelihood for any number of independent dimensions

bm_loglik_tree_constrained = function(tree, lower_bounds, upper_bounds, par, dimen) { 
		
  	A = par[1:dimen]  	
  	Sig = par[(dimen + 1):(dimen*2)]  
  	if (any(Sig < 0)) {
      return(-1e30)   		  	
    }
  	
  	n = length(tree$tip.label)
  	V = vcv(tree)	
	
	  res = mapply(function(lower_bounds, upper_bounds, A, Sig) pmvnorm(lower_bounds, upper_bounds, mean = rep(A, n), sigma = V*Sig), lower_bounds, upper_bounds, A, Sig)  	  	
  	
  	area_cor = mapply(function(upper_bounds,lower_bounds) sum(log(upper_bounds - lower_bounds)), upper_bounds,lower_bounds)	

    v = sum(log(as.numeric(res))) - sum(area_cor)
    if(is.nan(v) || abs(v)==Inf) v = -1e30
  	
  	return(v)
}

#--------------
# Optimize
# 'lower/upper_bounds' have to be 2 separate lists, with each list element being one-dimensional values

ranges.like.bm = function(tree, lower_bounds, upper_bounds, start_values = NA, dimen = NA) {
	
	if (is.na(dimen)) dimen = length(lower_bounds)
	if (is.na(start_values)) { start_values = c((unlist(lapply(lower_bounds, mean)) + 
		unlist(lapply(upper_bounds, mean)))/2,
		(unlist(lapply(lower_bounds, sd)) + 
		unlist(lapply(upper_bounds, sd)))/2) }
	
	opt.res = nlm(function(p) -1*bm_loglik_tree_constrained(tree,lower_bounds, 
				upper_bounds, p, dimen), p = start_values)	
	
	return(list(mrcas = opt.res$estimate[1:dimen],
		          rates = opt.res$estimate[(dimen + 1):(dimen*2)],
		          nlm.details = opt.res))
}
  

###############
# Functions for internal use
###############
# Basic BM for ancestors

bm_loglik_ancestors = function(tree, values, params) {	 				

	ntaxa = length(tree$tip.label)
	nnode = tree$Nnode

	ax = params[1:nnode]
	ay = params[(nnode+1):(2*nnode)]
  	sigma2x = params[2*nnode+1]
  	sigma2y = params[2*nnode+2]
  	if (sigma2x<0 || sigma2y<0) return(-Inf)

	tips = tree$edge[,2] <= ntaxa

	ax1 = ax[tree$edge[,1] - ntaxa]
  	ax2 = rep(NA, length(ax1))
  	ax2[tips]  = values$x[tree$edge[tips,2]]
  	ax2[!tips] = ax[tree$edge[!tips,2]-ntaxa]

  	ay1 = ay[tree$edge[,1] - ntaxa]
  	ay2 = rep(NA, length(ay1))
  	ay2[tips]  = values$y[tree$edge[tips,2]]
  	ay2[!tips] = ay[tree$edge[!tips,2]-ntaxa]

  	lx = -1*sum((ax1 - ax2)^2/tree$edge.length)/(2 * sigma2x) - nnode * log(sigma2x)
  	ly = -1*sum((ay1 - ay2)^2/tree$edge.length)/(2 * sigma2y) - nnode * log(sigma2y)
  	
  	return(lx + ly)
}

######################################
# utility function to estimate sigma2x and sigma2y by maximum likelihood and compute the 
# standard errors of these estimates.  This is only used to generate a proposal distribution for 
# sigma2 in rase.  

bm_est_sigma2 = function(tree, values, params) {
    
    ntaxa = length(tree$tip.label)
    nnode = tree$Nnode
    
    ax = params[1:nnode]
    ay = params[(nnode+1):(2*nnode)]
    
    tips = tree$edge[,2] <= ntaxa
    
    ax1 = ax[tree$edge[,1] - ntaxa]
    ax2 = rep(NA, length(ax1))
    ax2[tips]  = values$x[tree$edge[tips,2]]
    ax2[!tips] = ax[tree$edge[!tips,2]-ntaxa]
    
    ay1 = ay[tree$edge[,1] - ntaxa]
    ay2 = rep(NA, length(ay1))
    ay2[tips]  = values$y[tree$edge[tips,2]]
    ay2[!tips] = ay[tree$edge[!tips,2]-ntaxa]
    
    sigma2xhat = mean( (ax1 - ax2)^2 / tree$edge.length )
    sigma2xhat_sd = sqrt( -1/sum(1/(2*sigma2xhat^2) - (ax1-ax2)^2/(sigma2xhat^3 * tree$edge.length)) ) 
    
    sigma2yhat = mean( (ay1 - ay2)^2 / tree$edge.length )
    sigma2yhat_sd = sqrt( -1/sum(1/(2*sigma2yhat^2) - (ay1-ay2)^2/(sigma2yhat^3 * tree$edge.length)) ) 
    
    return(list(sigma2xhat=sigma2xhat, sigma2xhat_sd=sigma2xhat_sd, sigma2yhat=sigma2yhat, sigma2yhat_sd=sigma2yhat_sd))
}


#########################
# Loglik for ancestors using polygons

bm_loglik_ancestors_poly = function(tree, polygons, params, nGQ) {   				

  ntaxa = length(tree$tip.label)
  nnode = tree$Nnode
  
  ax = params[1:nnode]
  ay = params[(nnode+1):(2*nnode)]
  sigma2x = params[2*nnode+1]
  sigma2y = params[2*nnode+2]
  
  if(sigma2x<0 || sigma2y<0) return(-Inf)
  
  dat = cbind(tree$edge, tree$edge.length)
  logliks = apply(dat, 1, function(r) {
    ax1 = ax[r[1]-ntaxa]
    ay1 = ay[r[1]-ntaxa]
    if(r[2]<= ntaxa) {
      v = polyCub.SV(polygons[[r[2]]], bigauss_pdf, mx = ax1, my = ay1, sx = sqrt(sigma2x*r[3]), sy = sqrt(sigma2y*r[3]), rho = 0, nGQ=nGQ)    
      logv = ifelse(is.nan(v), -1e30, log(v)) 
      loglik = logv - log(area(polygons[[r[2]]]))
      if(is.nan(loglik)) loglik = -1e30
      #cat("a =",c(ax1,ay1), "\n")
      #cat("r[3] =",r[3], "\n")
      #cat("sigmas =", c(sigma2x,sigma2y), "\n")
      #cat("polycub loglik =", loglik, "\n")
    } else {
      loglik = sum(dnorm(c(ax[r[1]-ntaxa], ay[r[1]-ntaxa]), c(ax1, ay1), 
                         sd=sqrt(r[3]*c(sigma2x,sigma2y)), log=TRUE))
    }
    return(loglik)
  })
  
  return(sum(logliks))
}

#########################
# propose ancestral values (x,y) at an internal node 
# the ancestral values is a=(ax,ay)
# the daughter values are d1=(d1x,d1y) and d2=(d2x,d2y) 
# s is the ancestral branch length
# t1,t2 are the daughter branch lengths

bm_propose_trio = function(a, d1, d2, s,t1, t2, sigma2x, sigma2y) {

	if (is(d1, 'owin')) { 
    	d1 = poly_center(d1)
  	}

  	if (is(d2, 'owin')) { 
    	d2 = poly_center(d2)
  	}

	ax = a[1]; ay = a[2]
  	d1x = d1[1]; d1y = d1[2]
  	d2x = d2[1]; d2y = d2[2]

  	xm1 = (d1x*t2 + d2x*t1)/(t1+t2)
  	xv1 = t1*t2*sigma2x/(t1+t2)

  	ym1 = (d1y*t2 + d2y*t1)/(t1+t2)
  	yv1 = t1*t2*sigma2y/(t1+t2)

  	if(is.na(ax) || is.na(ay)) { 
   		x = rnorm(1, mean=xm1, sd=sqrt(xv1))
    	y = rnorm(1, mean=ym1, sd=sqrt(yv1))
     	logfwdprob = dnorm(x,xm1, sd=sqrt(xv1), log=TRUE)
      	logbwdprob = dnorm(y,ym1, sd=sqrt(yv1), log=TRUE)
  	} else {
   		xm2 = (xm1*s*sigma2x + ax*xv1)/(xv1 + s*sigma2x)
   		xv2 = xv1*s*sigma2x/(xv1 + s*sigma2x)
    	x = rnorm(1,mean=xm2,sd=sqrt(xv2))
    	ym2 = (ym1*s*sigma2y + ay*yv1)/(yv1 + s*sigma2y)
    	yv2 = yv1*s*sigma2y/(yv1 + s*sigma2y)
    	y = rnorm(1,mean=ym2,sd=sqrt(yv2))
     	logfwdprob = dnorm(x,xm1, sd=sqrt(xv1), log=TRUE)
      	logbwdprob = dnorm(y,ym1, sd=sqrt(yv1), log=TRUE)
  	}

    # return value, logfwdprob, logbwdprob
  	return(list(value=c(x,y), logfwdprob=logfwdprob, logbwdprob=logbwdprob))
}

#########################
# compute the log-likelihood of a trio
# using fast integration
# d1 and/or d2 can be shapes OR points

bm_loglik_trio = function(a, v, d1, d2, s, t1, t2, sigma2x, sigma2y, areas, daughter_ids, nGQ) {
  
  if (is.na(a[1])) {
    la = 0
  } else {
    la = sum(dnorm(v, a, sd=sqrt(s*c(sigma2x, sigma2y)), log=TRUE))
  }
  
  if (!is(d1, 'owin')) { # d1 is an internal node
    l1 = sum(dnorm(d1, v, sd = sqrt(t1*c(sigma2x, sigma2y)), log=TRUE))
  } else { # d1 is a tip
    l1 = log(as.numeric(polyCub.SV(d1, bigauss_pdf, mx = v[1], my = v[2], sx = sqrt(sigma2x*t1), sy = sqrt(sigma2y*t1), rho = 0, nGQ = nGQ))) - 
      log(areas[daughter_ids[1]])
    if (is.nan(l1)) l1 = -1e30
    #if (l1==0) l1 = -1e30
  }
  
  if (!is(d2, 'owin')) { # d2 is an internal node
    l2 = sum(dnorm(d2,v,sd=sqrt(t2*c(sigma2x,sigma2y)), log=TRUE))
  } else { # d2 is a tip
    l2 = log(as.numeric(polyCub.SV(d2, bigauss_pdf, mx = v[1], my = v[2], sx = sqrt(sigma2x*t2), sy = sqrt(sigma2y*t2), rho = 0, nGQ=nGQ))) - log(areas[daughter_ids[2]])
    if(is.nan(l2)) l2 = -1e30
    #if (l2==0) l2 = -1e30
  }
  
  return(la + l1 + l2)
}




#########################
# gibbs sampling for ancestors in 2-dimensional point bm

bm_ase = function(tree, values, niter=1e3, logevery=10, sigma2_scale=0.05, screenlog=TRUE, params0 = NA) {

    if (!is(tree, "phylo")) {
        stop('tree should be of class phylo')
    }
    
    ntaxa = length(tree$tip.label)
    nnode = tree$Nnode   
       
    ax = array(NA, dim=c(niter,nnode))
    sigma2x = rep(NA, niter)
    ay = array(NA, dim=c(niter,nnode))
    sigma2y = rep(NA, niter)
     
    if (any(is.na(params0))) {
        ace_resx = ace(values$x, tree, method = 'ML')
        ace_resy = ace(values$y, tree, method = 'ML')
        ax[1, (1:nnode)] = ace_resx$ace[1:nnode]
        sigma2x[1] = ace_resx$sigma[1]
        ay[1, (1:nnode)] = ace_resy$ace[1:nnode]
        sigma2y[1] = ace_resy$sigma[1]
    } else {
        if (length(params0) != 2*nnode+2) stop("starting values not of correct length")
        ax[1,] = params0[1:nnode]
        sigma2x[1] = params0[2*nnode+1]
        ay[1,] = params0[(nnode+1):(2*nnode)]
        sigma2y[1] = params0[2*nnode+2]
    }
    
    for (iter in 2:niter) {
                
        # random permutation of internal nodes
    	nodelist = ntaxa + sample(1:nnode, nnode, replace=FALSE)
        
        # populate current node values
    	ax[iter,] = ax[iter-1,]
    	ay[iter,] = ay[iter-1,]
    	sigma2x[iter] = sigma2x[iter-1]
    	sigma2y[iter] = sigma2y[iter-1]
        
        for (node in nodelist) {
            
            # daughters
            daughter_ids = tree$edge[tree$edge[,1]==node,2]
            t = c(tree$edge.length[tree$edge[,2]==daughter_ids[1]],
                tree$edge.length[tree$edge[,2]==daughter_ids[2]])
            
            if (daughter_ids[1]<= ntaxa) {
                d1_value = c(values$x[daughter_ids[1]], values$y[daughter_ids[1]])
            } else {
                d1_value = c(ax[iter,daughter_ids[1]-ntaxa], ay[iter,daughter_ids[1]-ntaxa])
            }
            
            if (daughter_ids[2]<=ntaxa) {
    	    	d2_value = c(values$x[daughter_ids[2]], values$y[daughter_ids[2]])
            } else {
    	    	d2_value = c(ax[iter,daughter_ids[2]-ntaxa], ay[iter,daughter_ids[2]-ntaxa])
            }
            
            # ancestor
            a_id = tree$edge[tree$edge[,2]==node,1]
            if (length(a_id)==0) {
                a_value = c(NA,NA)
                s = NA
            } else {
                a_value = c(ax[iter,a_id-ntaxa], ay[iter,a_id-ntaxa])
    	    	s = tree$edge.length[a_id-ntaxa]
            }	
            
            xy = bm_propose_trio(a_value, d1_value, d2_value, s, t[1],t[2], sigma2x, sigma2y)
            ax[iter,node-ntaxa] = xy$value[1]
            ay[iter,node-ntaxa] = xy$value[2]
            
    	}
                
        obj = bm_est_sigma2(tree, list(x=values$x, y=values$y), c(ax[iter,], ay[iter,]))
      	s2x_cur = obj$sigma2xhat
      	s2x_sd = sigma2_scale*obj$sigma2xhat_sd
      	s2y_cur = obj$sigma2yhat
      	s2y_sd = sigma2_scale*obj$sigma2yhat_sd

        repeat {
            sigma2x_prop = rnorm(1, mean=s2x_cur, sd=s2x_sd)
            if(sigma2x_prop>0) break
        }
    	  logprobratiox = dnorm(s2x_cur, sigma2x_prop, sd=s2x_sd, log=TRUE) - dnorm(sigma2x_prop, s2x_cur, sd=s2x_sd, log=TRUE)
        
        repeat {
            sigma2y_prop = rnorm(1, mean=s2y_cur, sd=s2y_sd)
            if(sigma2y_prop>0) break
        }
    	  logprobratioy = dnorm(s2y_cur, sigma2y_prop, sd=s2y_sd, log=TRUE) - dnorm(sigma2y_prop, s2y_cur, sd=s2y_sd, log=TRUE)
        
      
       # flat prior
        loglik_cur = bm_loglik_ancestors(tree, values, c(ax[iter,], ay[iter,], sigma2x[iter], sigma2y[iter]))
        loglik_prop = bm_loglik_ancestors(tree, values, c(ax[iter,], ay[iter,], sigma2x_prop, sigma2y_prop))
	
        logratio = loglik_prop - loglik_cur + logprobratiox + logprobratioy
        
         if (log(runif(1)) < logratio) {
            sigma2x[iter] = sigma2x_prop
            sigma2y[iter] = sigma2y_prop
        }
  	
  	    if (screenlog && iter%%logevery==0) {
                
                if (iter == logevery) {
                    cat('Iteration sigma2x sigma2y\n')
                    cat(iter, sigma2x[iter], sigma2y[iter],'\n')
                } else {
                    cat(iter, sigma2x[iter], sigma2y[iter],'\n')
                }

            } 
    }
    
    colnames(ax) = paste('n', names(branching.times(tree)), '_x', sep = '')
    colnames(ay) = paste('n', names(branching.times(tree)), '_y', sep = '')
    
    return(cbind(ax, ay, sigma2x, sigma2y))  
}

   

############################
# Almost-Gibbs sampling of ancestors with polygons at tips
# using polyCub for the terminal branch likelihoods

# Bivariate gaussian pdf
bigauss_pdf = function(s, mx, my, sx, sy, rho) {
  x = s[,1]
  y = s[,2]
  tr1 = 1 / (2 * pi * sx * sy * sqrt(1 - rho^2))
  tr2 = (x - mx)^2 / sx^2
  tr3 = -(2 * rho * (x - mx)*(y - my))/(sx * sy)
  tr4 = (y - my)^2 / sy^2
  z = tr2 + tr3 + tr4
  tr5 = tr1 * exp((-z / (2 *(1 - rho^2))))
  return (tr5)
}


#rase = range ancestral state estimation

rase = function(tree, polygons, niter=1e3, logevery=10, sigma2_scale=0.05, screenlog=TRUE, params0 = NA, nGQ = 20) {
  
  if (!is(tree, "phylo")) {
    stop('tree should be of class phylo')
  }

  if (!any(sapply(polygons, is.owin))) {
    stop('one or more polygons are not in \'owin\' format')
  }
  
  if (any(sapply(polygons, is.empty))) {
    stop('one or more polygons are empty')
  }  

  if (any(is.na(match(tree$tip.label, names(polygons))))) {
    stop('tip labels and polygon names do not match')
  }
  
  ntaxa = length(tree$tip.label)
  nnode = tree$Nnode
  
  # initialize 
  areas = mapply(area, polygons)    
  xy.tips = t(mapply(poly_center, polygons))
  
  ax = array(NA, dim=c(niter,nnode))
  sigma2x = rep(NA, niter)
  ay = array(NA, dim=c(niter,nnode))
  sigma2y = rep(NA, niter)
  
  if (any(is.na(params0))) {
    ace_resx = ace(xy.tips[,1], tree, method = 'ML')
    ace_resy = ace(xy.tips[,2], tree, method = 'ML')
    ax[1, (1:nnode)] = ace_resx$ace[1:nnode]
    sigma2x[1] = ace_resx$sigma[1]
    ay[1, (1:nnode)] = ace_resy$ace[1:nnode]
    sigma2y[1] = ace_resy$sigma[1]
  } else {
    if (length(params0) != 2*nnode+2) stop("starting values not of correct length")
    ax[1,] = params0[1:nnode]
    sigma2x[1] = params0[2*nnode+1]
    ay[1,] = params0[(nnode+1):(2*nnode)]
    sigma2y[1] = params0[2*nnode+2]
  }
  
  for (iter in 2:niter) {
    # random permutation of internal nodes
    nodelist = ntaxa+sample(1:nnode, nnode, replace=FALSE)
    
    # populate current node values
    ax[iter,] = ax[iter-1,]
    ay[iter,] = ay[iter-1,]
    sigma2x[iter] = sigma2x[iter-1]
    sigma2y[iter] = sigma2y[iter-1]
    
    for (node in nodelist) {
      
      approx = 0
      # daughters
      daughter_ids = tree$edge[tree$edge[,1]==node,2]
      t = c(tree$edge.length[tree$edge[,2]==daughter_ids[1]],
            tree$edge.length[tree$edge[,2]==daughter_ids[2]])
      
      if (daughter_ids[1]<= ntaxa) { # tips
        d1_value = polygons[[daughter_ids[1]]]
        approx = 1
      } else {
        d1_value = c(ax[iter,daughter_ids[1]-ntaxa], ay[iter,daughter_ids[1]-ntaxa])
      }
      
      if (daughter_ids[2]<=ntaxa) { # tips
        d2_value = polygons[[daughter_ids[2]]]
        approx = 1
      } else {
        d2_value = c(ax[iter,daughter_ids[2]-ntaxa], ay[iter,daughter_ids[2]-ntaxa])
      }
      
      # ancestor
      a_id = tree$edge[tree$edge[,2]==node,1]
      if (length(a_id)==0) {
        a_value = c(NA,NA)
        s = NA
      } else {
        a_value = c(ax[iter,a_id-ntaxa], ay[iter,a_id-ntaxa])
        s = tree$edge.length[a_id-ntaxa]
      }	
      
      # proposal 
      xy_prop = bm_propose_trio(a_value, d1_value, d2_value, s, t[1], t[2], sigma2x[iter], sigma2y[iter])
      
      if (approx) {
        loglik_prop = bm_loglik_trio(a_value, xy_prop$value, d1_value, d2_value, 
                                     s, t[1], t[2], sigma2x[iter], sigma2y[iter], areas, daughter_ids, nGQ)
        
        loglik_cur = bm_loglik_trio(a_value, c(ax[iter,node-ntaxa], ay[iter,node-ntaxa]), 
                                    d1_value, d2_value, s, t[1], t[2], sigma2x[iter], sigma2y[iter], areas, daughter_ids, nGQ)
        
        logratio = loglik_prop - loglik_cur + xy_prop$logbwdprob - xy_prop$logfwdprob
        
        if (log(runif(1)) < logratio) {
          ax[iter,node-ntaxa] = xy_prop$value[1]
          ay[iter,node-ntaxa] = xy_prop$value[2]
        }
      } else {
        ax[iter,node-ntaxa] = xy_prop$value[1]
        ay[iter,node-ntaxa] = xy_prop$value[2]
      }
    }
    
    # sample sigma2x and sigma2y
    
    # new function estimates sigma2x and sigma2y, and their standard errors.
    # these are used in a normal approximation proposal below
    obj = bm_est_sigma2(tree, list(x=xy.tips[,1], y=xy.tips[,2]), c(ax[iter,], ay[iter,]))
    s2x_cur = obj$sigma2xhat
    s2x_sd = sigma2_scale*obj$sigma2xhat_sd
    s2y_cur = obj$sigma2yhat
    s2y_sd = sigma2_scale*obj$sigma2yhat_sd
    
    
    repeat {
      sigma2x_prop = rnorm(1, mean=s2x_cur, sd=s2x_sd)
      if(sigma2x_prop>0) break
    }
    logprobratiox = dnorm(s2x_cur, sigma2x_prop, sd=s2x_sd, log=TRUE) - dnorm(sigma2x_prop, s2x_cur, sd=s2x_sd, log=TRUE)
    
    repeat {
      sigma2y_prop = rnorm(1, mean=s2y_cur, sd=s2y_sd)
      if(sigma2y_prop>0) break
    }
    logprobratioy = dnorm(s2y_cur, sigma2y_prop, sd=s2y_sd, log=TRUE) - dnorm(sigma2y_prop, s2y_cur, sd=s2y_sd, log=TRUE)
    
    
    # 1/sigma2 prior
    loglik_cur = bm_loglik_ancestors_poly(tree,polygons,c(ax[iter,],ay[iter,],sigma2x[iter],sigma2y[iter]), nGQ) - (log(sigma2x[iter])+log(sigma2y[iter]))
    loglik_prop = bm_loglik_ancestors_poly(tree,polygons,c(ax[iter,],ay[iter,],sigma2x_prop,sigma2y_prop), nGQ) - (log(sigma2x_prop)+log(sigma2y_prop))
    #cat("loglik_cur =", loglik_cur, "\n")
    #cat("loglik_prop =", loglik_prop, "\n")
    
    logratio = loglik_prop - loglik_cur + logprobratiox + logprobratioy
    
    if (log(runif(1)) < logratio) {
      sigma2x[iter] = sigma2x_prop
      sigma2y[iter] = sigma2y_prop
    }
    
    if (screenlog && iter%%logevery==0) {
      
      if (iter == logevery) {
        cat('Iteration sigma2x sigma2y\n')
        cat(iter, sigma2x[iter], sigma2y[iter],'\n')
      } else {
        cat(iter, sigma2x[iter], sigma2y[iter],'\n')
      }
      
    } 
  }
  
  
  colnames(ax) = paste('n', names(branching.times(tree)), '_x', sep = '')
  colnames(ay) = paste('n', names(branching.times(tree)), '_y', sep = '')
  
  return(cbind(ax, ay, sigma2x, sigma2y))  	
}




#########################
# polygon center

poly_center = function(poly) {
	xy = poly$bdry
	centr = c()
	for (j in 1:length(xy)) {
		x = c(xy[[j]]$x, xy[[j]]$x[1])
		y = c(xy[[j]]$y, xy[[j]]$y[1])
		
	  	n = length(x)
	  	if (length(y) != n) stop("length of x and y must be equal")

	  	Cx = 0; Cy = 0; A = 0
	  	for(i in 1:(n-1)) {
	    	Cx = Cx + (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i]) 
	    	Cy = Cy + (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i]) 
	    	A = A + x[i]*y[i+1] - x[i+1]*y[i]
	  	}
	  	A = A/2
	  	Cx = Cx/(6*A)
	  	Cy = Cy/(6*A)
		
  		centr = rbind(centr,c(Cx,Cy, area(poly)))
	}
	return(c(weighted.mean(centr[,1], centr[,3]),
		weighted.mean(centr[,2], centr[,3])))
}



###################
# Locations at time slices
###################

###################
# takes a phylogenetic tree and a time slice
# and returns the pair of nodes whose branch 
# exists at that time

tree.slice = function(tree, slice) {

	if (!is(tree, 'phylo')) {
        stop('tree should be of class phylo')
    }
    if (!is.binary.tree(tree)) {
    	stop("tree should be fully dichotomous")
    }
	nnode = tree$Nnode
	ntips = min(tree$edge[,1]) - 1
	b.t = branching.times(tree)

	rage = cbind(tree$edge, NA, NA)

	rage[match(1:ntips, rage[,2]), 4] = 0
	rage = rage[order(rage[,1]),]
	rage[,3] = rep(b.t, each = 2)
	rage = rage[order(rage[,2]),]
	rage[which(rage[,2]>ntips),4] = b.t[2:(ntips - 1)]


	c1 = rage[which(rage[, 3] > slice),]
	c2 = c1[which(c1[, 4] < slice),]

	if (nrow(c2) == 0) cat('no branches at that time slice') else return(c2)
		
}


################
# MCMC of species ranges at time slice t


#########################
# compute the log-likelihood of the location of v
# at time slice u, given ancestor a, descendant d, 
# separated by time t

bm_loglik_duo = function(a, v, d, u, t, sx, sy, nGQ) {

   	la = sum(dnorm(v, a, sd=sqrt(u*c(sx, sy)), log=TRUE))

  	if (!is(d, 'owin')) { 
    	ld = sum(dnorm(d, v, sd = sqrt((t-u)*c(sx, sy)), log=TRUE))
  	} else { # d1 is a tip
    	ld =  log(as.numeric(polyCub.SV(d, bigauss_pdf, mx = v[1], my = v[2], sx = sqrt(sx*(t-u)), sy = sqrt(sy*(t-u)), rho = 0, nGQ = nGQ))) - 
        log(area(d))
    	if (is.nan(ld)) ld = -1e30
  	}

	return(la + ld)
}




#########################
# propose v values (x,y) at time slice u 
# the ancestral values is a=(ax,ay)
# the daughter values are d=(dx,dy) 
# u is the ancestral branch length
# t is the branch length

bm_propose_duo = function(a, d, u, t, sx, sy) {

	if (is(d, 'owin')) { 
    	d = poly_center(d)
  	}

	ax = a[1]; ay = a[2]
  	dx = d[1]; dy = d[2]

	xm1 = (dx*u + ax*(t-u))/t 
	xv1 = (t-u)*u*sx/t
	
	ym1 = (dy*u + ay*(t-u))/t 
	yv1 = (t-u)*u*sy/t
	
    x = rnorm(1, mean=xm1, sd=sqrt(xv1))
	y = rnorm(1, mean=ym1, sd=sqrt(yv1))

    logfwdprob = dnorm(x, mean=xm1, sd=sqrt(xv1), log=TRUE)
    logbwdprob = dnorm(y, mean=ym1, sd=sqrt(yv1), log=TRUE)

    # return value, logfwdprob, logbwdprob
  	return(list(value=c(x,y), logfwdprob=logfwdprob, logbwdprob=logbwdprob))
}



#################
# MCMC for locations at a given 
# slice of time, with the results
# from a rase run 

rase.slice = function(tree, slice, res, polygons, params0 = NA, niter=1e3, logevery=10, nGQ = 20) {

    if (!is(tree, "phylo")) {
        stop('tree should be of class phylo')
    }
    
    if (any(is.na(match(tree$tip.label, names(polygons))))) {
        stop('tip labels and polygon names do not match')
    }
    
                                        # Slice the tree with separate function
    a_d = tree.slice(tree, slice)
    
    ntaxa = length(tree$tip.label) # number of taxa
    nbranch = nrow(a_d) # number of branchs
    
    anc_x = res[, a_d[,1]-ntaxa] 				# all the ancestors posterior distributions from rase results in _x 
    anc_y = res[, tree$Nnode + a_d[,1]-ntaxa] 	# and _y
    
                                        # if no starting parameter values are provided,	
    if (is.na(params0)) { 	# use time weighted mean as starting values
	
        for (b in 1:nbranch) { # for each branch		
            if (a_d[b,4] == 0) { # if daughter is a tip, use the polygon centroid to calculate the mean
                pxy = poly_center(polygons[[b]])
                
                params0[b] = ((mean(anc_x[,b])*(slice - a_d[b, 4])) + 
                              (pxy[1]*(a_d[b, 3] - slice)))/(a_d[b, 3] - a_d[b, 4])
                
                params0[b+nbranch] = ((mean(anc_y[,b])*(slice - a_d[b, 4])) + 
                                      (pxy[2]*(a_d[b, 3] - slice)))/(a_d[b, 3] - a_d[b, 4])
                
            } else { # if daughter is not a tip
                params0[b] = ((mean(anc_x[,b])*(slice - a_d[b, 4])) + 
                              (mean(res[, a_d[b,2]-ntaxa])*(a_d[b, 3] - slice)))/(a_d[b, 3] - a_d[b, 4])
                
                              
                params0[b+nbranch] = ((mean(anc_y[,b])*(slice - a_d[b, 4])) + 
                                      (mean(res[, tree$Nnode + a_d[b,2] - ntaxa])*(a_d[b, 3] - slice)))/(a_d[b, 3] - a_d[b, 4])

            }		
        }
        
		cat("Using time weighted mean as starting parameter values \n")
	
    } else {		
                                        # Check for length of input initial parameter values
        if (length(params0) != (nrow(a_d)*2)) stop("starting values not of correct length")
    }
    
	
	bx = array(NA, dim=c(niter,nbranch))		# array of branch slice x's and y's  	
  	bx[1,] = params0[1:nbranch]					# populate starting nodes with starting values
  	by = array(NA, dim=c(niter,nbranch))
  	by[1,] = params0[(nbranch+1):(2*nbranch)]
  	
  	# use the posterior distribution given by rase of sigma x and y 
	sigma2x = res[,ncol(res)-1]
	sigma2y = res[,ncol(res)]

	for (iter in 2:niter) {
	
		if (logevery && iter%%logevery==0) cat("iter =", iter, "\n")
    	
    	# sample one single joint sampler
    	samp = sample(1:length(sigma2x),1)
    	
    	# sigma2x and sigma2y from this sample (sx & sy)
    	    	
    	sx = sigma2x[samp]
    	sy = sigma2y[samp]
    	
    	# populate current node values
    	bx[iter,] = bx[iter-1,]
    	by[iter,] = by[iter-1,]

		for (i in 1:nbranch) { # each i is a branch in the tree, a row in the output of tree.slice
			
			approx = 0
			# sample the ancestor of the branch from rase posterior distribution (a_value)
			a_value = c(anc_x[samp,i], anc_y[samp,i])
			
			# get time t - time from the ancestor to daughter
			t = a_d[i, 3] - a_d[i, 4]
			# get time u (slice since the ancestor)
    		u = a_d[i, 3] - slice
			
			daughter_id = a_d[i, 2]
    		       
			if (daughter_id <= ntaxa) { # if daughter is a tip
				d_value = polygons[[daughter_id]] # get the polygon as d_value
      			approx = 1
    		} else { # if not a tip, sample from the posterior distribution given by rase
    			d_value = c(res[samp, a_d[i,2]-ntaxa], res[samp, tree$Nnode + a_d[i,2]-ntaxa])			
			}

		    # proposal duo (x & y value)
 			xy_prop = bm_propose_duo(a_value, d_value, u, t, sx, sy)

			if (approx) {
				
				# likelihood of proposal
				loglik_prop = bm_loglik_duo(a_value, xy_prop$value, d_value, 
                                       u, t, sx, sy, nGQ)
          		
          		# likelihood of current node
				loglik_cur = bm_loglik_duo(a_value, c(bx[iter,i], by[iter,i]), 
                                      d_value, u, t, sx, sy, nGQ)

				#logratio
				logratio = loglik_prop - loglik_cur + xy_prop$logbwdprob - xy_prop$logfwdprob

				if (log(runif(1)) < logratio) {
    				bx[iter, i] = xy_prop$value[1]
    	    		by[iter, i] = xy_prop$value[2]
    			}
			} else {
  				bx[iter, i] = xy_prop$value[1]
        		by[iter, i] = xy_prop$value[2]
    		}
		
		}


  	}

	colnames(bx) = paste('b', 1:nbranch, '_x', sep = '')
	colnames(by) = paste('b', 1:nbranch, '_y', sep = '')

	return(cbind(bx, by))

}




################
# Data Handling
################


################
# Transform shapefiles from 'sp' package for use in rase

shape.to.rase = function(shape_poly) {
	if (!is(shape_poly, 'SpatialPolygonsDataFrame')) {
		stop('Error: object is not of class SpatialPolygonsDataFrame')
	}
	
	pols = list()
	for (i in 1:length(shape_poly)) {		
    	fp1 = shape_poly[i,]
    	fp2 = fp1@polygons[[1]]
		fp2 = lapply(fp2@Polygons, function(lst){owin(poly = list(x=rev(lst@coords[,1]),y=rev(lst@coords[,2]),hole=lst@hole), check = TRUE)})
		pols = c(pols, fp2)
	}	
	return(pols)	
}




################
# Name polygons, order them as the tree tips

name.poly = function(polygons, tree, poly.names = NA) {
	
	if (!is(tree, "phylo")) {
        stop("Error: tree is not of class phylo")
    }

	tips = tree$tip.label
	
	if (any(is.na(poly.names))) {
		nam = names(polygons)
	} else {
		nam = poly.names
		names(polygons) = nam
	}
	if (any(is.na(match(tips, nam)))) stop('tip labels and polygon names do not match')
			
	if (all(match(tips, nam) == 1:length(polygons))) {
		cat('tip labels and polygon names match and are in the same order \n') 
	} else {	
		polygons = polygons[match(tips, nam)]	
	}
		
	return(polygons)
}



################
# Handling the post Gibbs sampling
	
post.mcmc = function(res, burnin = 1000, thin = 10, as.ggmcmc = TRUE) {
	
	if (thin == 0) warning('Note that if thin = 0, then no iteration is saved at all.')
	
	mc = res[-(1:burnin),]
	mc = mc[seq(1, dim(mc)[1], thin),]

	if (as.ggmcmc == TRUE) {
		fgg = data.frame(Iteration = rep(1:dim(mc)[1], dim(mc)[2]),
				Chain = rep(1, dim(mc)[1]*dim(mc)[2]),
				Parameter = rep(colnames(mc), each = dim(mc)[1]),
				value = as.vector(mc))

		attr(fgg, 'nChains') = as.integer(1)
		attr(fgg, 'nParameters') = dim(mc)[2]
		attr(fgg, 'nIterations') = dim(mc)[1]
		attr(fgg, 'nBurnin') = burnin
		attr(fgg, 'nThin') = thin
		attr(fgg, 'description') = 'S'
		attr(fgg, 'parallel') = FALSE
	
		return(fgg)
	} else {
		return (mc)		
	}
}


#########################
# Visualization functions with rgl
#########################

#########################
# Function to transform data structure for 3D plotting

data.for.3d = function(res, tree, polygons) {  
		
	if (is.null(tree$node.label)) {  	
       	tree$node.label <- paste("n", names(branching.times(tree)), sep="")
    }
      
    xy.tips = t(mapply(poly_center, polygons)) 
    xy.nodes = matrix(colMeans(res)[1:(2*tree$Nnode)], nrow = tree$Nnode, ncol = 2)
    xy.all = data.frame(rbind(xy.tips, as.data.frame(xy.nodes)))

    z = c(rep(0, times=length(tree$tip.label)), branching.times(tree))

    names(xy.all) = c("x", "y")
    label = c(tree$tip.label, tree$node.label)
    return(list(xyz = data.frame(xy.all, z, label), edge = tree$edge, pol = polygons))
}


#########################
# Function that opens the 'rgl' 3d device and plots the phylogenetic tree 
# to the 3d setting
# Note: each node is placed to its estimated geographic position
  
phylo.3d = function(df3, z.scale = 1, pts = TRUE, ...) {
    
    edg = df3$edge
 
    if (pts == TRUE) points3d(df3$xyz$x, df3$xyz$y, z.scale*df3$xyz$z)
        
	for (i in 1:nrow(edg)) {
		x = df3$xyz$x[edg[i,]]
      	y = df3$xyz$y[edg[i,]]
      	z = df3$xyz$z[edg[i,]]
  
      	lines3d(x, y, z*z.scale, ...)
    }
}


#########################
# Function that adds the terminal polygons to the opened 'rgl' device

add.polygons = function(df3, axes = 2, ...) {
    
	for (j in 1:length(df3$pol)) {
    	poli = df3$pol[[j]]

      	polygon3d(x=c(poli$bdry[[1]]$x, poli$bdry[[1]]$x[1]), 
      			y=c(poli$bdry[[1]]$y, poli$bdry[[1]]$y[1]), 
                z=rep(0, times=length(poli$bdry[[1]]$x) + 1), 
                ...)
    }
	if (axes == 0) print('No axes displayed')
	if (axes == 1) {
		axes3d(edges = c('x'), xlab = 'Longitude')
		title3d('','','Longitude')
	} 
	if (axes == 2) {
		axes3d(edges = c('x', 'y'), xlab = 'Longitude', ylab = 'Latitude')
		title3d('','','Longitude','Latitude')
	}
	if (axes == 3) {
		axes3d(edges = c('x', 'y', 'z'), xlab = 'Longitude', ylab = 'Latitude', zlab = 'Time')
		title3d('','','Longitude','Latitude','Time')
	}
}


#########################
# Function that adds the ancestral posterior densities from the Gibbs sampling 
# to the opened 'rgl' device

add.dens = function(df3, res, nlevels = 20, z.scale = 1, col = c(1:nnode), ...) {
    node.ages = df3$xyz[df3$xyz[, 'z'] != 0, 3] 
    nnode = length(node.ages)
   
    for (i in 1:nnode)
    {
    	densy = sm.density(res[,c(i, i + nnode)], display = 'none')		
	    cont = contourLines(densy$eval.points[,1], densy$eval.points[,2], densy$estimate, nlevels = nlevels)
		polygon3d(x = cont[[1]]$x, y = cont[[1]]$y, 
		z = 0.02+(z.scale*rep(node.ages[i], times = length(cont[[1]]$x))), col = col[i], ...)
    }
}



