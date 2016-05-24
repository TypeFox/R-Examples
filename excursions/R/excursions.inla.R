## excursions.inla.R
##
##   Copyright (C) 2012, 2013, 2014, David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.


excursions.inla <- function(result.inla,
                            stack,
                            name=NULL,
                            tag=NULL,
                            ind=NULL,
                            method,
                            alpha=1,
                            F.limit,
                            u,
                            u.link = FALSE,
                            type,
                            n.iter=10000,
                            verbose=0,
                            max.threads=0,
                            seed=NULL)
{
  if (!requireNamespace("INLA", quietly=TRUE))
    stop('This function requires the INLA package (see www.r-inla.org/download)')
  if(missing(result.inla))
	  stop('Must supply INLA result object')
  if(missing(method)){
	  cat('No method selected, using QC\n')
	  method = 'QC'
  }
  if(missing(u))
	  stop('Must specify level u')
  if(missing(type))
	  stop('Must specify type of excursion')

  if(result.inla$.args$control.compute$config==FALSE)
	  stop('INLA result must be calculated using control.compute$config=TRUE')

 if(missing(F.limit)) {
    F.limit = alpha
  } else {
    F.limit = max(alpha,F.limit)
  }

  n = length(result.inla$misc$configs$config[[1]]$mean)

  #Get indices for the component of interest in the configs
  ind.stack <- inla.output.indices(result.inla, name=name, stack=stack, tag=tag)
  n.out = length(ind.stack)
  #Index vector for the nodes in the component of interest
  ind.int <- seq_len(n.out)

  #ind is assumed to contain indices within the component of interest
  if(!missing(ind) && !is.null(ind)){
    ind <- private.as.vector(ind)
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack

  # If u.link is TRUE, the limit is given in linear scale
  # then transform to the scale of the linear predictor
  u.t = rho = rep(0,n)
  if(u.link == TRUE){
    links = result.inla$misc$linkfunctions$names[
                                            result.inla$misc$linkfunctions$link]
    u.tmp = sapply(ind, function(i) private.link.function(u,links[i]))
    u.t[ind] = u.tmp
  } else {
    u.t = u
  }

  #Calculate marginal probabilities
  #If stack and tag are provided, we are interested in the predictor
  #Otherwise, we are interested in some of the effects

  if(verbose)
	  cat('Calculating marginal probabilities\n')

  #Are we interested in a random effect?
  random.effect = FALSE
  if(!missing(name)&&!is.null(name)&&name!="APredictor"&&name!="Predictor"){
    random.effect = TRUE
    if(is.null(result.inla$marginals.random))
    stop('INLA result must be calculated using return.marginals.random=TRUE if excursion sets to be calculated for a random effect of the model')
  }

  if(!random.effect && is.null(result.inla$marginals.linear.predictor))
    stop('INLA result must be calculated using return.marginals.linear.predictor=TRUE if excursion sets are to be calculated for the linear predictor.')

  if(random.effect) {
    rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal(ind.int[j],
                           u=u,result = result.inla,
                           effect.name=name, u.link = u.link, type = type))
  } else {
    rho.ind <- sapply(1:length(ind), function(j) inla.get.marginal(ind[j],
                           u=u,result = result.inla, u.link = u.link, type = type))
  }
  rho[ind] = rho.ind

  n.theta = result.inla$misc$configs$nconfig
  for(i in 1:n.theta){
    config = private.get.config(result.inla,i)
    if(config$lp == 0)
      break
  }

  if(verbose)
    cat('Calculating excursion function using the ', method, ' method\n')

  if(method == 'EB' || method == 'QC' ) {
	  res = excursions(alpha=alpha, u=0, mu=config$mu-u.t, Q=config$Q,
	                   type=type, method=method, F.limit = F.limit,
	                   vars=config$vars, rho=rho,
		  	  		       ind=ind, n.iter=n.iter, max.threads=max.threads,seed=seed)
  	F = res$F[ind]
  } else if (method =='NI' || method == 'NIQC') {
    qc = 'QC'
  	if(method == 'NI')
	  	qc = 'EB'

  	res = lw = NULL

	  for(i in 1:n.theta){
		  if(verbose)
			  cat('Configuration ', i, ' of ', n.theta, '\n')

  		conf.i = private.get.config(result.inla,i)
	  	lw[i] = conf.i$lp
		  res[[i]] = excursions(alpha=alpha,u=0,mu=conf.i$mu-u.t,Q=conf.i$Q,
		                        type=type,method=qc,rho=rho,vars=conf.i$vars,
							              ind=ind,n.iter=n.iter,
							              F.limit = F.limit,
							              max.threads=max.threads,seed=seed)
  	}

  	w = exp(lw)/sum(exp(lw))

  	F = w[1]*res[[1]]$F[ind]
	  for(i in 2:n.theta){
 		  F = F + w[i]*res[[i]]$F[ind]
	  }
  } else if (method == 'iNIQC') {

  	pfam.i = rep(-0.1,n)
  	pfam.i[ind] = rho.ind
  	reo = sort(pfam.i,index.return=TRUE)$ix
  	pfam.i[!ind] = 0

	  res = lw =  NULL

	  for(i in 1:n.theta){
	  	if(verbose)
	  		cat('Configuration ', i, ' of ', n.theta, '\n')

  		conf.i = private.get.config(result.inla,i)
  		lw[i] = conf.i$lp
  		r.i <- INLA::inla(result.inla$.args$formula,
                        family = result.inla$.args$family,
                        data=result.inla$.args$data,
                        control.compute = list(config = TRUE),
                        control.predictor=result.inla$.args$control.predictor,
                        control.mode = list(theta=
                         as.vector(result.inla$misc$configs$config[[i]]$theta),
                                  fixed=TRUE))

      if(random.effect) {
        p1.i <- sapply(1:length(ind), function(j) inla.get.marginal(ind.int[j],
                           u=u,result = r.i, effect.name=name,
                           u.link = u.link, type = type))

      } else {
        p1.i <- sapply(1:length(ind), function(j) inla.get.marginal(ind[j],
                           u=u,result = result.inla,
                           u.link = u.link, type = type))
      }
	  	pfam.i[ind] = p1.i

	  	res[[i]] = excursions(alpha=alpha,u=0,mu=conf.i$mu-u.t,
		  					            Q=conf.i$Q, type=type, method='QC',
	             						  rho=pfam.i,vars=conf.i$vars,
		  				          	  max.size=length(ind),reo=reo,
		  				          	  F.limit = F.limit,
			  				            n.iter=n.iter,max.threads=max.threads,seed=seed)
	  }

  	w = exp(lw)/sum(exp(lw))
  	F = w[1]*res[[1]]$F[ind]
  	for(i in 2:n.theta){
  		F = F + w[i]*res[[i]]$F[ind]
  	}
  } else {
  	stop('Method must be one of EB, QC, NI, NIQC, iNIQC')
  }

  if(type == "="){
    rho.ind = 1-pmax(rho.ind,1-rho.ind)
  }
  if(type == "!="){
    rho.ind = pmax(rho.ind,1-rho.ind)
  }

  F.out = mu.out = rho.out = vars.out = E.out = M.out = G.out = rep(NA,n.out)

  F.out[ind.int] = F
  vars.out[ind.int] = config$vars[ind]
  mu.out[ind.int] = config$mu[ind]
  rho.out[ind.int] = rho.ind


  G = rep(0,length(config$mu[ind]))
  if(type == "<") {
    G[config$mu[ind] > u.t[ind]] = 1
  } else {
    G[config$mu[ind] >= u.t[ind]] = 1
  }
  G.out[ind.int] = G

  E = rep(0,length(F))
  E[F>1-alpha] = 1
  E.out[ind.int] = E

  M = rep(-1,length(F))
  if (type=="<") {
    M[E==1] = 0
  } else if (type == ">") {
    M[E==1] = 1
  } else if (type == "!=" || type == "=") {
    M[E==1 & config$mu[ind]>u] = 1
    M[E==1 & config$mu[ind]<u] = 0
  }

  M.out[ind.int] = M

  output <- list(E = E.out,
                 F = F.out,
                 G = G.out,
                 M = M.out,
                 mean = mu.out,
                 vars = vars.out,
                 rho=rho.out,
                 meta=list(calculation="excursions",
                           type=type,
                           level=u,
                           level.link=u.link,
                           alpha=alpha,
                           F.limit = F.limit,
                           n.iter=n.iter,
                           method=method,
                           ind=ind.int,
                           call = match.call()))
  class(output) <- "excurobj"
  output
}
