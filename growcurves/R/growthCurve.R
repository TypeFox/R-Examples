######################################################
## growthCurve.R function,
## which produces within subject predicted
## growth curves.
######################################################

#' Within subject model-predicted growth curve
#'
#' Produces a set of predicted response values, by subject, at \code{T} time points.  The response values
#' are predicted by employing the posterior samples of model parameters where the resultant response
#' values for each subject are composed by averaging over all posterior samples in a Rao-Blackwellizing fashion.
#'
#' @param y.case The \code{N x 1} (subject-time case) vector of data response values. 
#' @param B The \code{M x P*q} matrix of subject random effect posterior samples.  \code{M} = number of MCMC samples,
#'	\code{P} = number of subjects, \code{q} = number of random effect parameters, per subject.
#' @param Alpha The \code{M x 1} vector for the model intercept parameter.
#' @param Beta The \code{M x F} matrix of model fixed effects parameters, where \code{F} = number of fixed effects
#' @param U The \code{M x S} matrix of univariate multiple membership random effects, where \code{S} = number of random effects.
#'	\code{U} is multivariate, then the input is of dimension \code{M x Nmv*S}, where \code{Nmv} is the multivariate dimension.
#'	Leave \code{NULL} is don't require the multiple membership effects. 
#'	Input as list of \code{M x S} matrices if have more than one mutiple membership term.
#' @param aff.clients Vector of length \code{P.aff} that identifies subjects affected by \code{U}.  Identical to \code{subj.aff} from \code{\link{dpgrowmm}}.
#'			Input as list of vectors, each comprised of affected subjects attached to the equivalent multiple membership term if have more than one term.
#' @param W.subj A \code{P x S} multiple membership weight matrix for \code{U} that expands \code{W.subj.aff} of \code{\link{dpgrowmm}} from affected subjects, \code{Paff} to all subjects, \code{P}.
#'			Input as list of \code{P[i] x S[i]} matrices, where i indexes an MM term, if have more than one multiple membership term.
#' @param X.n A design matrix with \code{N} rows (for subject-measure) cases providing nuisance fixed effects.  Will be expanded to
#' 	the \code{T} within sample predictions, but held constant between successive observed values (for generating expanded predictions).
#' @param Z.n A design matrix with \code{N} rows providing nuisance random effects.  Grouping is assumed to be by-subject.
#' @param trt.case The treatment group membership vector of length \code{N} (subject-time cases).  Assumed numeric with lowest group level == 0; 
#'		.e.g. \code{(0,0,0,1,1,2,2,2,2,)}.
#' @param trt.lab Associated labels for the numeric treatment groups.  Each distinct treatment group assumed to have a unique label. 
#' @param subject.case Vector of length \code{N} providing subject-measure cases.  Must be in numerical format with unique subjects sequential starting at 1.
#' @param subject.lab \code{N x 1} case length vector with user desired labels that map 1:1 to \code{subject.case}.
#' @param T Number of time points to build each subject curve.  \code{T = 10} is typically sufficient.
#' @param min.T The minimum time value that \code{T} will take.
#' @param max.T The maximum time value that \code{T} will take.
#' @param n.thin The gap between each MCMC sample used for the growth curve.
#' @param n.waves The maximum number of observed measurement waves, per subject.
#' @param time.case A vector of length \code{N} providing times for associated subject-measure observations.  Identical to \code{time} from \code{\link{dpgrowmm}}.
#' @param n.fix_degree The highest polynomial degree to employ for constructing time-based fixed effects covariates.  
#' @param Nrandom A scalar input providing the number of by-subject time-based random effect parameters. Only need to input if employ nuisance random effects.
#' @return A list object containing the following \code{data.frames} and plots:
#'     \item{plot.dat}{A \code{data.frame} object containing the within-subjects predicted growth curves.
#'			Fields are titled, \code{c("fit","time","subject","trt")}.}
#'     \item{dat.data}{A \code{data.frame} object containing actual by-subject data for each measurement occasion.
#'			Fields are titled, \code{c("fit","time","subject","trt")}.}
#'     \item{p.gcall}{A \code{ggplot2} object that aggregates growth curves by treatment type.}
#'     \item{p.gcsel}{A \code{ggplot2} object that plots growth curves with associated data points for 10 randomly selected subjects.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}, \code{\link{dpgrowmult}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrow}}, \code{\link{dpgrowmm}}, and \code{\link{dpgrowmult}}
#' @aliases growthCurve
#' @export
growthCurve	= function(y.case, B, Alpha, Beta, U = NULL, aff.clients = NULL, W.subj = NULL, X.n = NULL, Z.n = NULL, trt.case, 
				trt.lab, subject.case, subject.lab, T, min.T, max.T, n.thin, n.waves = NULL, time.case, n.fix_degree, Nrandom = NULL)
{
  ## function to produce set of growth curves for P clients over T time points using the MCMC output from the growth curve mixed
  ## effects models, with or without a multiple membership component.
  ## note: "trt.case" is of length N (number of subject-time cases), not number of subjects, P, and allows for multiple treatment groups
	
  subject		= unique(subject.case)
  P			= length(subject)

  ## create treatment vector, trt, of length subject
  tmp			= data.frame(trt.case,subject.case,stringsAsFactors = FALSE) 
  names(tmp)		= c("trt","subject")
  tmp			= unique(tmp) ## produces trt of length Nsubject
  trt			= tmp$trt
  rm(tmp)
  
  ## capture clients receiving additive session effect term
  if ( !is.null(W.subj) )
  {
  	if( is.list(W.subj) ) ## multiple MM terms, so cull unique clients
  	{
		nty = length(aff.clients) ## aff.clients is also a list object
		## aff.clients = do.call("rbind",aff.clients) ## assuming aff.clients is a list of column vectors
		aff.clients = unlist(aff.clients)
		aff.clients = unique(aff.clients) ## clients affected by ANY mm term
  	}else{ ## single MM term, either univariate or multivariate
		S <- ncol(W.subj) ## number of MM (session) effects
	}
  } ## end conditional statement on !is.null(W.subj)

  mm.clients		= matrix(0,P,1)
  if(!is.null(aff.clients))
  {
  	mm.clients[aff.clients]	= 1
  } ## else mm.clients == 0, which bypasses session effects for lgm and dp models

  ## capture number of random effects - both growth and nuisance			
  Q			= ncol(B)/P ## equivalent to Nrandom only if Z.n = NULL; otherwise, Q > Nrandom

  ## fix num.ran.time to be number of time-indexed random effects under either presence or absence of nuisance covariates
  if(is.null(Z.n)) 
  {
	num.ran.time	= Q
  }else{
	if( is.null(Nrandom) )
        {
		stop("Must input 'Nrandom', number of time-indexed random effects, because have additional nuisance random effects.")
	}else{ ## user inputs Nrandom
		num.ran.time 	= Nrandom
	}
  }

  ## determine number of fixed effect growth variables and nuisance variables
  len.gcurve.fe		= n.fix_degree + (n.fix_degree+1)*(max(trt))	
  F			= ncol(Beta)
  if (is.null(X.n) )
  {
	len.n.fe	= 0 ## no nuisance fixed covariates
  }else{ ## user input fixed covariates, X.n
  	len.n.fe		= F - len.gcurve.fe
  	stopifnot( len.n.fe == ncol(X.n) ) ## computed length of nuisance covariates must equal number of columns of X.n
  }

  ## subsample MCMC iterations
  M			= nrow(Beta) ## number of MCMC iterations
  K			= floor(M/n.thin) ## number of MCMC iterations to sample (after thinning)
  Nlevel		= length(unique(trt))
  samps			= as.integer(seq(1,M, length.out = K)) ## set of thinned MCMC iterations to sample

  ## define P*T growth curve variables for plotting
  trt.g			= vector(mode="numeric",length = P*T)
  time.g		= vector(mode="numeric",length = P*T)
  subject.g		= vector(mode="numeric",length = P*T)
  Y			= matrix(0,P*T,K) ## to capture fit values for each client*time 
  y			= vector(mode="numeric",length = P*T) ## capture fit values averaged over MCMC iterations (Rao-Blackwellized)
  time			= seq(min.T,max.T,length.out = T)

  ## define throw-away / written-over by-subject variables
  bmati			= matrix(0,M,Q)
  z.c			= as.matrix(0,num.ran.time,1)
  xc.trt_gen		= as.matrix(0,(n.fix_degree+1),1)

  if(!is.null(X.n) | !is.null(Z.n) )  ## there exists extra non-time-based columns of X that must be expanded
  { 
    #############################################################
    ## expand X from P*n.waves to P*T - allows different values for X.e, Z.e by wave, but holds those values constants for added time points
    #############################################################
    expander		= matrix(0,P*T,1) ## vector to expand nrows of X.n from P*n.waves to P*T
    markers		= matrix(0,n.waves,1)
    wave		= as.integer(as.factor(time.case)) 
    wave.u		= 1:n.waves
    if( !is.null(X.n) )
    {
    		X.look		= as.data.frame(cbind(subject.case,wave,1:nrow(X.n)))
    }else{ ## !is.null(Z.n)
		X.look		= as.data.frame(cbind(subject.case,wave,1:nrow(Z.n)))
    }
    names(X.look) 	= c("subject","wave","row")

    ## Fill in missing rows of X.look (for missing waves)
    X.look	= split(X.look,X.look$subject)
    X.look 	= lapply(X.look,function(x){
					num.waves 	= nrow(x)
					num.add  	= n.waves - num.waves 
                              		if(num.add > 0)
					{
					  	add.waves 			= setdiff(wave.u,x$wave)
						x[num.waves + 1:num.add,] 	= 0
						x[num.waves + 1:num.add,2] 	= add.waves
						x				= x[order(x$wave),]	
						for(i in 1:length(add.waves))
						{
							dist = x$wave - add.waves[i]
							loc  = which(dist == min(dist))
							x[add.waves[i],c(1,3)] = x[loc[1],c(1,3)]
						}
					}
					x})
    X.look			= do.call("rbind",X.look)
    row.names(X.look) 		= 1:nrow(X.look)

    ## compute wave markers in time
    for(j in 1:n.waves)
    {
	  markers[j] = min.T + (j-1)*(max.T/(n.waves-1))
    }

    ## Define X row 'expander' to pad non-time based covariates with duplicate rows
    for (i in 1:P)
    {
		for(t in 1:T)
		{	
			if( time[t] < markers[2] ) wave.t = 1
			for(j in 2:n.waves)
			{
				if( (time[t] < markers[j]) & (time[t] > markers[j-1]) )
				{
					wave.t = j - 1
				}else wave.t = j

			} ## end loop j over n.waves
			expander[(i-1)*T+t] = subset( X.look,(subject == i & wave == wave.t) )$row
      	} ## end loop t over time grid
    } ## end loop i over subjects

    if( !is.null(X.n) ) X.e		= as.matrix(X.n[expander,]) ## P*T version of X
    if ( !is.null(Z.n) ) Z.e	 	= as.matrix(Z.n[expander,]) ## P*T version of Z; Xlook applies here because Z is not time-based and nrow(Z) = nrow(X) = Ncase
  } ## end conditional statement on whether there are non-time-based (extra) X predictors that need to be expanded

  ## capture common fixed effects coefficients, beta, for all subjects (non-treatment covariates)
  pos.common	= c(1:n.fix_degree)
  for(i in 1:P) ## over subjects
  {
	for(j in 1:Q)
	{
		bmati[,j]	= B[,((j-1)*P + i)] ## M x 1
	}
	for(t in 1:T) ## over time grid points
	{
		## random effects
		for( j in 1:num.ran.time)
		{
			z.c[j]		= time[t]^(j-1)
		}
		if ( num.ran.time < Q )
		{
			z.c[(num.ran.time+1):Q]		= Z.e[(i-1)*T+t,]
		}

		## fixed effects
		## construct generic polynomial against which treatment effect is multiplied
    		for( k in 1:(n.fix_degree+1) )
    		{
    			xc.trt_gen[k]			= time[t]^(k-1)
    		}
    		x.c					= xc.trt_gen[-1] ## non-treatment effect fixed covariates (excludes intercept)
		this.trt				= as.integer(trt[i] != 0)
		if( this.trt == 1 )  x.c		= c(x.c,this.trt*xc.trt_gen)
		l.x.c					= length(x.c)
		if( len.n.fe > 0 ) ## there are nuisance effefcts
		{
			start				= l.x.c + 1
			end				= start + len.n.fe - 1
			x.c[start:end]			= X.e[(i-1)*T+t,]
		}

		x	= as.matrix(x.c)  ## note: x.c doesn't capture full columns of X, only those growth curve columns for subject i + nuisance
		z 	= as.matrix(z.c)

		## identifiers (for plotting in ggplot2)
		trt.g[(i-1)*T+t]		= trt[i]
		time.g[(i-1)*T+t]		= time[t]
		subject.g[(i-1)*T+t]		= subject[i]

		##
		## generate curves for over K thinned MCMC iterations
		##

  		## In the case there ARE nuisance FE: starting and end columns of Beta for nuisance fixed effects
		start.e			= len.gcurve.fe + 1
		end.e			= start.e + len.n.fe - 1

		if( mm.clients[i] > 0 ) ## client i receives multiple membership effect
		{
			if(trt[i] == 0) ## control group - no treatment*time covariates
			{
				if (len.n.fe > 0 ) ## there are nuisance covariates
				{
					if( is.list(W.subj) ) ## there is more than 1 MM term
					{
						common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.e:end.e)] %*% x + bmati[samps,] %*% z 
						for(block in 1:nty)
						{
							common.term = common.term + U[[block]][samps,] %*% W.subj[[block]][i,]
						}
						Y[((i-1)*T+t),]		= common.term
					}else{ ## single MM term
						if( ncol(U) > ncol(W.subj) ) ## U is multivariate
						{
							common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.e:end.e)] %*% x + bmati[samps,] %*% z
							for( q in 1:num.ran.time ) ## fix Nmv = Nrandom, number of time-indexed random effect terms
							{
								indx.q 		<- (q-1)*S + 1:S
								common.term 	<- common.term + (time[t]^(q-1))*U[samps,indx.q] %*% W.subj[i,]
							}
							Y[((i-1)*T+t),]		= common.term
						}else{ ## U is univariate
							Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.e:end.e)] %*% x +
								 		bmati[samps,] %*% z + U[samps,] %*% W.subj[i,]
						}
					} ## end conditional statement on whether multiple MM effects

				}else{ ## no nuisance covariates
					if( is.list(W.subj) )
					{
						common.term	=  Alpha[samps] + Beta[samps,pos.common] %*% x + bmati[samps,] %*% z
						for(block in 1:nty)
						{
							common.term = common.term + U[[block]][samps,] %*% W.subj[[block]][i,]
						}
						Y[((i-1)*T+t),]		= common.term
					}else{ ## single MM term
						if( ncol(U) > ncol(W.subj) ) ## U is multivariate
						{
							common.term	=  Alpha[samps] + Beta[samps,pos.common] %*% x + bmati[samps,] %*% z
							for( q in 1:num.ran.time )
							{
								indx.q 		<- (q-1)*S + 1:S
								common.term 	<- common.term + (time[t]^(q-1))*U[samps,indx.q] %*% W.subj[i,]
							}
							Y[((i-1)*T+t),]		= common.term
						}else{ ## U is univariate
							Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,pos.common] %*% x +
								 		bmati[samps,] %*% z + U[samps,] %*% W.subj[i,]
						}
					} ## end conditional statement on whether multiple MM effects
				}
		
			}else{ ## find the specific set of num.ran.time fixed effects associated to treatment trt[i] for subject i
				## growth curve fixed effects
				start.c			= n.fix_degree + (n.fix_degree+1)*(trt[i] - 1)	+ 1
				end.c			= start.c + (n.fix_degree+1) - 1
				## nuisance fixed effects
				if ( len.n.fe > 0 )
				{
					if( is.list(W.subj) )
					{
						common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c,start.e:end.e)] %*% x + bmati[samps,] %*% z
						for(block in 1:nty)
						{
							common.term = common.term + U[[block]][samps,] %*% W.subj[[block]][i,]
						}
						Y[((i-1)*T+t),]		= common.term
					}else{ ## single MM term
						if( ncol(U) > ncol(W.subj) ) ## U is multivariate
						{
							common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c,start.e:end.e)] %*% x + bmati[samps,] %*% z
							for( q in 1:num.ran.time )
							{
								indx.q 		<- (q-1)*S + 1:S
								common.term 	<- common.term + (time[t]^(q-1))*U[samps,indx.q] %*% W.subj[i,]
							}
							Y[((i-1)*T+t),]		= common.term
						}else{ ## U is univariate
							Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c,start.e:end.e)] %*% x +
								 		bmati[samps,] %*% z + U[samps,] %*% W.subj[i,]
						}
					} ## end conditional statement on whether multiple MM effects
				}else{ ## no nuisance covariates
					if( is.list(W.subj) )
					{
						common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c)] %*% x + bmati[samps,] %*% z 
						for(block in 1:nty)
						{
							common.term = common.term + U[[block]][samps,] %*% W.subj[[block]][i,]
						}
						Y[((i-1)*T+t),]		= common.term
					}else{ ## single MM term
						if( ncol(U) > ncol(W.subj) ) ## U is multivariate
						{
							common.term	=  Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c)] %*% x + bmati[samps,] %*% z 
							for( q in 1:num.ran.time )
							{
								indx.q 		<- (q-1)*S + 1:S
								common.term 	<- common.term + (time[t]^(q-1))*U[samps,indx.q] %*% W.subj[i,]
							}
							Y[((i-1)*T+t),]		= common.term
						}else{ ## U is univariate
							Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c)] %*% x +
								 			bmati[samps,] %*% z + U[samps,] %*% W.subj[i,]
						}
					} ## end conditional statement on whether multiple MM effects
				} ## end conditional statement on whether nuisance covariates
			} ## end conditional statement on whether control group or one of treatment groups
		}else{ ## no multiple membership effects

			if(trt[i] == 0) ## control group - no treatment*time covariates
			{
				if (len.n.fe > 0 )
				{
					Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.e:end.e)] %*% x + bmati[samps,] %*% z 
				}else{ ## no nuisance covariates
					Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,pos.common] %*% x + bmati[samps,] %*% z 
				}
		
			}else{ ## find the specific set of Nrandom fixed effects associated to treatment trt[i] for subject i
				## growth curve fixed effects
				start.c			= n.fix_degree + (n.fix_degree+1)*(trt[i] - 1)	+ 1
				end.c			= start.c + (n.fix_degree+1) - 1
				## nuisance fixed effects
				if ( len.n.fe > 0 )
				{
					Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c,start.e:end.e)] %*% x  + 
									bmati[samps,] %*% z 
				}else{ ## no nuisance covariates
					Y[((i-1)*T+t),]		= Alpha[samps] + Beta[samps,c(pos.common,start.c:end.c)] %*% x +
								  	bmati[samps,] %*% z 
				} ## end conditional statement on whether nuisance covariates
			} ## end conditional statement on whether control group or one of treatment groups
		} ## end conditional statement on whether client[i] receives an additive multiple membership term

		## compute Rao-Blackwellized average response, y, for subject i at time t
		y[(i-1)*T+t]		= mean(Y[(i-1)*T+t,])
	} ## end loop over t
  } ## end loop over i

  plot.dat			= as.data.frame(cbind(y,time.g,subject.g,trt.g))
  names(plot.dat) 		= c("fit","time","subject","trt")

  ##
  ## plot growth curves for ALL clients, grouped by CBT (treatment or no) - mm_car_gcurves.eps
  ## for MM(CAR) method (or any other single method).
  ##

  ## construct trt labels for plotting
  dat.gc		= plot.dat
  trt.lab		= unique(trt.lab)
  Nlevel		= length(trt.lab) ## number of unique treatment groups
  labs			= c()
  if(Nlevel > 1 )
  {
  	for(k in 1:Nlevel)
   	{
		labs	= c( labs, paste("TRT",trt.lab[k], sep = "_") )
	}
  }else{ ## only a single group
	labs	= "baseline"
  }

  slabs			= unique(subject.lab) ## subject labels for plotting

  ##
  ## for ALL clients - MEAN growth curve, by TREATMENT category
  ## 					
  dat.gc$trt		= factor(dat.gc$trt,labels = labs)
  dat.gc$subject	= factor(dat.gc$subject,labels = slabs) ## NOTE: the user-input labels, not the modeling labels, are returned in dat.gc and dat.data
  p			= ggplot(data=dat.gc,aes(x=time, y = fit, group = subject)) ## no individual subject visibility
  ## l			= geom_line(colour = alpha("black",1/5), linetype=5)
  l			  = geom_line(colour = "black", alpha = 1/5, linetype=5)
  l.2			= geom_smooth(aes(group=1),method = "loess", size = 1.1, colour = "black")
  axis	 	= labs(x = "Time", y = "Model Fit", type = "Cluster")
  f			  = facet_wrap(~trt, scales="fixed") 
  p.gcall		= p + l + l.2 + f + axis

  ##
  ## for SELECTED clients - WITH DATA - maybe looking at 2 methods - gcurves_subject.eps
  ##
  num.subjects			= max(subject.case)
  num.samp			= min(5,num.subjects) ## number of subjects to sample
  dat.data			= as.data.frame(cbind(y.case,time.case,subject.case,trt.case))  ## use subject.lab, instead of subject to return user input labels
  names(dat.data)		= c("fit","time","subject","trt")
  dat.data$subject		= factor(dat.data$subject, labels = slabs)
  dat.data$trt			= factor(dat.data$trt,labels = labs)
  dat.by.trt			= split(dat.gc,dat.gc$trt) ## sample num.samp subjects from each of the 5 clusters to plot 
  per.samp			= as.vector(sapply(dat.by.trt,function(x){
					choices = sample(unique(x$subject),num.samp)
					return(choices)}))
  dat.p				= subset(dat.gc,subject %in% per.samp)
  datreal.p			= subset(dat.data,subject %in% per.samp)
  p				= ggplot(data=dat.p,aes(x=time, y = fit) )
  l				= geom_line(size = 1.2)
  l.2				= geom_point(data=datreal.p,size=3,shape=1) ## Add the data
  axis	 			= labs(x = "Time", y = "Fit")
  f				= facet_wrap(trt~subject, scales="fixed") 
  p.gcsel			= p + l + f + axis + l.2 

  return(list(plot.dat = dat.gc, dat.data = dat.data, trt = trt.g, subject = subject.g, time = time.g, p.gcall = p.gcall, p.gcsel = p.gcsel))

  fit <- time <- subject <- trt <- NULL; rm(fit); rm(time); rm(subject); rm(trt);

  gc()

} ## end of function growthCurve