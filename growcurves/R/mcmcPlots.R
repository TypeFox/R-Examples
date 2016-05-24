#' generate plots of model(s) posterior results
#'
#' Constructs plots of subject and multiple membership effects, as well as traceplots
#' for model precision and clustering parameters.  Returns a list of objects of class \code{ggplot}.
#' 
#' @param subjecti.u A vector of length \code{P}, number of unique subjects, containing unique set of user input values for \code{subject}.
#' @param subj.aff A vector of length \code{P.aff} identifying the unique subjects (which are a subset of variable, subject) receiving multiple membership random effects. 
#'			Applies only to case of a single set of multiple membership random effects. 
#' @param subjaff.input User input version of \code{subj.aff} that may be character or numeric format.  (Again, this is a strict subset of subjecti.u). 
#'			Applies only to case of a single set of multiple membership random effects. 
#' @param bmat.summary A list object of \code{q} elements, each containing an \code{P x 3} matrix of c(2.5\%,50\%,97.5\%) quantile summaries
#'		for each subject of the applicable subject random effect parameter.  \code{P} = number of subjects, \code{q} = number of random effect parameters, per subject.
#' @param group  An \code{S x 1} vector of group identifiers for the multiple membership random effects,
#'			where \code{S} is the number of multiple membership random effects.  The format is sequential numeric, starting at 1.
#'			Applies only to case of a single set of multiple membership effects.
#' @param groupi.u A vector of user input unique values for the multiple membership effect group identifiers where employ 1 multiple membership term.
#'			Input as a list of \code{S x 1} vectors in the case of more than one set of multiple membership effects.
#' @param u.summary An \code{S x 3} matrix of of quantile summaries for each multiple membership session effect where employ 1 multiple membership term.
#'			Input as list of \code{S x 3} quantiles in the case of more than one set of multiple membership effects.
#' @param Nmv The order for the multiple membership effects.  Defaults to \code{Nmv = 1} for univariate effects.  Otherwise, \code{Nmv > 1}
#'		indicates that \code{u.summary} is dimensioned as \code{Nmv*S x 3}.
#' @param ulabs An \code{nty} vector of labels for each term (block) in the case of more than one set of multiple membership effects.
#' @param mm.summary A \code{P.aff x 3} matrix of quantile summaries.  \code{mm} was created by multiple the set of \code{S} multiple membership
#'		effects, \code{u}, on each MCMC iteration by the multiple membership design matrix, \code{W.subj.aff}.
#' @param M The \code{iter.keep x 1} matrix of posterior samples for the parameter capturing the number of clusters formed under the DP prior on the client effects.
#' @param Tauu \code{iter.keep x 1} matrix of posterior samples capturing the precision parameter for \code{"mmcar", "mmi" and "mmigrp"}.
#'		Input as \code{iter.keep x nty} matrix in the case of \code{nty} multiple membership effect terms.
#' @param Taub \code{iter.keep x Nrandom} matrix of posterior samples capturing the precision parameter for each of the sets of subject random effects.
#' @param Taue \code{iter.keep x 1} matrix of posterior samples capturing the precision parameter for the model error term.
#' @param Deviance \code{iter.keep x 1} matrix of posterior samples for the model deviance.
#' @return A list of plot objects of class \code{ggplot2} including:
#'     \item{p.U}{by group plot of session effects, u[1:Nsession]. Plot is faceted for more than one set of effect terms.}
#'     \item{p.Umm}{plot of "mm = W.subj.aff \%*\% u" for those clients attending assessions.}
#'     \item{p.Ub0}{plot of " mm + b0", the total random intercept, for those clients attending sessions.} 
#'     \item{p.Ub}{plot of "mm + b" for multivariate MM effects with order equal to "Nrandom".}
#'     \item{p.b}{stacked plots of b0,...,b(q-1) - vertical lines for each client span 2.5\% - 97.5\% values with mean noted.} 
#'     \item{p.M}{MCMC trace plot of M, number of clusters.} 
#'     \item{p.tauu}{MCMC trace plots of tau.u. Plot is faceted for more than one set of effect terms.} 
#'     \item{p.taue}{MCMC trace plots of tau.e.} 
#'     \item{p.taub}{MCMC faceted trace plot for each of the q components of tau.b.} 
#'     \item{p.dev}{MCMC trace plots of deviance.} 
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrow}}, \code{\link{dpgrowmm}}, and \code{\link{dpgrowmult}}
#' @aliases mcmcPlots plotmcmc
#' @export
mcmcPlots = function(subjecti.u, subj.aff = NULL, subjaff.input = NULL, bmat.summary, group = NULL, groupi.u = NULL, u.summary = NULL, Nmv = 1,
			ulabs = NULL, mm.summary = NULL, M = NULL, Tauu = NULL, Taub, Taue, Deviance)
{  
   ##
   ## some dimensions
   ##
   if( !is.null(u.summary) ) 
   {
	if( is.null(group) & is.null(groupi.u) ) stop("\nMust enter group, either through 'group' or 'groupi.u' to plot MM term effects.\n")
	if( is.null(groupi.u) ) groupi.u = group
	if( !is.list(u.summary) )
	{
		if( is.null(Nmv) )
		{
			Nsession 	= nrow(u.summary)
		}else{
			Nsession	= nrow(u.summary)/Nmv
		}
		G		= length(unique(group))
	}else{ ## multiple MM terms
		nty		= ncol(Tauu)
		numt		= vector("numeric",nty)
		for(i in 1:nty)
		{
			numt[i] = nrow(u.summary[[i]])
		}
		Nsession	= sum(numt)

		if( !is.null(ulabs) ) ## session block labels
		{
			stopifnot(length(ulabs) == nty)
		}else{
			ulabs = as.character(1:nty)
		}	
	} ## end conditional statment oon !is.list(u.summary)
   }

   iter.keep	= length(Taue)
   Nrandom	= length(bmat.summary)  ## bmat.summary is a list object of length Nrandom

   ##
   ## produce plots
   ##
   if( !is.null(u.summary) ) ## no use of subj.aff b/c not mapping back to an mm term for plotting
   {
	if( is.list(u.summary) ) ## assume group is a list of length nty, each containing some vector of group labels of length numt[i]
	{
		block		= factor(ulabs) ## must be of length nty
		dat		= vector(mode = "list", length = nty)
		for(i in 1:nty)
		{
			## note: using unique values of subject labels of length P, not case formatted of length N.
                        ## note: groupi.u is here of length numt[i] (number of sessions), so it contains duplicate values
			dat[[i]]			= data.frame(1:numt[i],block[i],groupi.u[[i]],u.summary[[i]]) ## using treat to facet plot
			names(dat[[i]])			= c("session","block","group","low","mean","high")
			## dat.long[[i]]		= melt(dati,id=c("session","block","group"), measure = c("low","mean","high"))
   		} ## end loop i to create long data.frame with label identifier for each bi
      			datU			= as.data.frame(do.call("rbind",dat))
		## plot
   		p.U		 		= ggplot(data=datU,aes(x=session, y=mean, colour=factor(group)))
   		l 				= geom_point()
		l.2				= geom_smooth(aes(group=block),method = "loess", span = 1.0, size = 1.1, se = FALSE, colour = "pink")
   		f				= facet_wrap(~block,scales="free",ncol=2)
   		axis	 			= labs(x = "Effect ID", y = "Effect Size", colour = "Group")
   		options		 		= labs(title="Random Multiple Membership Effects under Blocks")
   		p.U		 		= p.U + l + l.2 + f + axis + options

		## Tau.u - MCMC trace
   		dattauu			= as.data.frame(cbind(1:iter.keep,Tauu))
  		names(dattauu)		= c("iteration",ulabs)
   		dattauu.long		= melt(dattauu,id="iteration")
   		names(dattauu.long) 	= c("iteration","label","value")
   		p.tauu			= ggplot(data=dattauu.long,aes(x=iteration,y=value))
   		l				= geom_line()
   		f 				= facet_wrap(~label, scales="free_y",ncol=1)
   		axis				= labs(x = "Iterations", y = "Sampled Value")
   		options			= labs(title="MCMC Trace plots for Tauu")
   		p.tauu			= p.tauu + l + f + axis + options

		rm(datU); rm(dat); rm(dattauu); rm(dattauu.long)
	
	}else{
		if( Nmv > 1 ) ## multivariate MM effects
		{
			## construct trt labels for plotting
   			labs = c()
			groupi.u = unique(groupi.u) ## in the case entered with length equal to number of sessions
   			for(g in 1:G)
   			{
				labs	= c( labs, paste("Group",groupi.u[g], sep = "_") )
   			}

   			## session effects - by session number - without credible ranges
   			## data.frame
   			group				= factor(group,labels = labs)
			group				= rep(group, times = Nmv)
   			session 			= rep(1:Nsession, times = Nmv)
			order				= rep(1:Nmv,each = Nsession)
   			datU				= data.frame(session,group,order,u.summary)  ##u.summary is Nmv*Nsession x 3 (3 statistics for each component of u)
   			names(datU)			= c("session","cluster","order","low","mean","high")
 
   			## plot
   			p.U		 		= ggplot(data=datU,aes(x=session, y=mean, colour=factor(cluster)))
   			l 				= geom_point(aes(shape=factor(order)))
   			f				= facet_wrap(~cluster,scales="free",ncol=2)
   			axis	 			= labs(x = "Session", y = "Effect Size", shape = "Poly Order", colour = "Group")
   			options		 		= labs(title="Session Effects (u)")
   			p.U		 		= p.U + l + f + axis + options


   			## mm(u) - by client - with credible ranges
   			## long data.frame
   			client				= rep(subjaff.input,times = Nmv)  ## direct correspondence between subj.aff and subjaff.input
			order.mm			= rep(1:Nmv,each = length(subjaff.input))
   			datUmm				= data.frame(client,order.mm,mm.summary)
   			names(datUmm)			= c("client","order","low","mean","high")
   			datUmm.long			= melt(datUmm,id.vars=c("client","order"),measure = c("low","mean","high"))
 
   			## plot
   			p.Umm 		= ggplot(data=datUmm.long,aes(x=client, y=value , group=factor(client)))
   			l.1 		= geom_line(colour = "steelblue4")
   			l.2		= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
			f		= facet_wrap(~order,scales="free",ncol=2)
   			axis 		= labs(x = "Subject", y = "Effect Size")
   			options 	= labs(title="Session-induced client effects - mm = H (W * U)")
   			p.Umm 		= p.Umm + l.1 + l.2 + f + axis + options 


   			## mm(u) + b - by client - with credible ranges
   			## long data.frame
   			## Note: used subj.aff instead of subjaff.input to select subjects from bmat.summary 
			num.subject			= length(unique(subjecti.u))
			sel.subject			= rep(1:num.subject,times=Nmv)
			b.summary			= do.call("rbind",bmat.summary)
			b.summary			= b.summary[sel.subject,] ## just in case there are nuisance random effects such that Nrandom > Nmv
			pos.sel				= which(sel.subject %in% subj.aff)
   			Ub.summary			= b.summary[pos.sel,] + mm.summary  ## only using subset of clients that attended sessions
   			datUb				= data.frame(client,order.mm,Ub.summary)
   			names(datUb)			= c("client","order","low","mean","high")
   			datUb.long			= melt(datUb,id.vars=c("client","order"),measure = c("low","mean","high"))

   			p.Ub 	= ggplot(data=datUb.long,aes(x=client, y=value , group=factor(client)))
   			l.1 	= geom_line(colour = "steelblue4")
   			l.2	= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
			f	= facet_wrap(~order,scales="free",ncol=2)
   			axis 	= labs(x = "Subject", y = "Effect Size")
   			options = labs(title="Session-induced client effects - mm + b")
   			p.Ub 	= p.Ub + l.1 + l.2 + f + axis + options 

	
   			## tau.u - MCMC trace
   			dattauu			= data.frame(1:iter.keep,Tauu)
   			names(dattauu) 		= c("iteration",paste("order",1:Nmv,sep="_"))
			dattauu			= melt(dattauu,id.vars = "iteration")
			dattauu$order		= dattauu$variable; dattauu$variable <- NULL
   			p.tauu			= ggplot(data=dattauu,aes(x=iteration,y=value))
   			l			= geom_line()
			f			= facet_wrap(~order,scales="free",ncol=2)
   			axis			= labs(x = "Iterations", y = "Sampled Value")
   			options			= labs(title="MCMC Trace plot for tau.u")
   			p.tauu			= p.tauu + l + f + axis + options

			rm(datU); rm(datUmm); rm(datUmm.long); rm(datUb); rm(datUb.long); 

		}else{ ## univariate MM effects
   			## construct trt labels for plotting
   			labs = c()
			groupi.u = unique(groupi.u) ## in the case entered with length equal to number of sessions
   			for(g in 1:G)
   			{
				labs	= c( labs, paste("Group",groupi.u[g], sep = "_") )
   			}

   			## session effects - by session number - without credible ranges
   			## data.frame
   			group				= factor(group,labels = labs)
   			session 			= 1:Nsession
   			datU				= as.data.frame(cbind(session,group,u.summary))  ##u.summary is Nsession x 3 (3 statistics for each component of u)
   			names(datU)			= c("session","group","low","mean","high")
 
   			## plot
   			p.U		 		= ggplot(data=datU,aes(x=session, y=mean, colour=factor(group)))
   			l 				= geom_point()
   			f				= facet_wrap(~group,scales="free",ncol=2)
   			axis	 			= labs(x = "Session", y = "Effect Size", colour = "Group")
   			options		 		= labs(title="Session Effects (u)")
   			p.U		 		= p.U + l + f + axis + options


   			## mm(u) - by client - with credible ranges
   			## long data.frame
   			client				= subjaff.input  ## direct correspondence between subj.aff and subjaff.input
   			datUmm				= as.data.frame(cbind(client,mm.summary))
   			names(datUmm)			= c("client","low","mean","high")
   			datUmm.long			= melt(datUmm,id="client",measure = c("low","mean","high"))
 
   			## plot
   			p.Umm 		= ggplot(data=datUmm.long,aes(x=client, y=value , group=factor(client)))
   			l.1 		= geom_line(colour = "steelblue4")
   			l.2		= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
   			axis 		= labs(x = "Subject", y = "Effect Size")
   			options 	= labs(title="Session-induced client effects - mm = (W * u)")
   			p.Umm 		= p.Umm + l.1 + l.2 + axis + options 


   			## mm(u) + b0 - by client - with credible ranges
   			## long data.frame
   			## Note: used subj.aff instead of subjaff.input to select subjects from bmat.summary 
   			Ub0.summary			= bmat.summary[[1]][subj.aff,] + mm.summary  ## only using subset of clients that attended sessions
   			datUb0				= as.data.frame(cbind(client,Ub0.summary))
   			names(datUb0)			= c("client","low","mean","high")
   			datUb0.long			= melt(datUb0,id="client",measure = c("low","mean","high"))

   			p.Ub0 	= ggplot(data=datUb0.long,aes(x=client, y=value , group=factor(client)))
   			l.1 	= geom_line(colour = "steelblue4")
   			l.2	= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
   			axis 	= labs(x = "Subject", y = "Effect Size")
   			options = labs(title="Session-induced client effects - mm0 + b0")
   			p.Ub0 	= p.Ub0 + l.1 + l.2 + axis + options 
	
   			## tau.u - MCMC trace
   			dattauu			= as.data.frame(cbind(1:iter.keep,Tauu))
   			names(dattauu) 		= c("iteration","value")
   			p.tauu			= ggplot(data=dattauu,aes(x=iteration,y=value))
   			l			= geom_line()
   			axis			= labs(x = "Iterations", y = "Sampled Value")
   			options			= labs(title="MCMC Trace plot for tau.u")
   			p.tauu			= p.tauu + l + axis + options

			rm(datU); rm(datUmm); rm(datUmm.long); rm(datUb0); rm(datUb0.long);
		}

	} ## end conditional statment on is.list(u.summary) - are there multiple MM terms

   } ## end conditional statement on !is.null(u.summary)

   ## b_0, .., b_nr - by client - with credible ranges
   ## data.frame
   dat.long		  	= vector(mode = "list", length = Nrandom)
   for(i in 1:Nrandom)
   {
	## note: using unique values of subject labels of length P, not case formatted of length N.
	dat			= as.data.frame(cbind(subjecti.u,bmat.summary[[i]])) ## using unique values of subject.input for plotting
	names(dat)		= c("client","low","mean","high")
	dat.long[[i]]		= melt(dat,id="client",measure = c("low","mean","high"))
	records			= nrow(dat.long[[i]])
	labs			= rep(paste("b",i-1,sep=""),records)
	dat.long[[i]]		= as.data.frame(cbind(dat.long[[i]],labs))
	names(dat.long[[i]])[4] = "type"
   } ## end loop i to create long data.frame with label identifier for each bi
      	datB			= as.data.frame(do.call("rbind",dat.long))
	datB$type		= factor(datB$type)
	rm(dat); rm(dat.long)

   ## plot
   p.b 		= ggplot(data=datB,aes(x=client, y=value , group=factor(client)))
   l.1 		= geom_line(colour = "steelblue4")
   ## l.2	= stat_summary(fun.y = mean, geom = "point", size = 3, shape = 16, colour = "pink")
   l.2		= geom_smooth(aes(group=1),method = "loess", span = .2, size = 1.2, colour = "pink")
   f 		= facet_wrap(~type, scales="free_y",ncol=1)
   axis		= labs(x = "Subject", y = "Effect Size")
   options 	= labs(title="Client Effects, b")
   p.b 		= p.b + l.1 + l.2 + f + axis + options 

   ## Tau.b - MCMC trace
   dattaub		= as.data.frame(cbind(1:iter.keep,Taub))
   names(dattaub)	= c("iteration",as.character(0:(Nrandom-1)))
   dattaub.long		= melt(dattaub,measure=as.character(0:(Nrandom-1)))
   names(dattaub.long) 	= c("iteration","label","value")
   p.taub		= ggplot(data=dattaub.long,aes(x=iteration,y=value))
   l			= geom_line()
   f 			= facet_wrap(~label, scales="free_y",ncol=1)
   axis			= labs(x = "Iterations", y = "Sampled Value")
   options		= labs(title="MCMC Trace plot for Taub")
   p.taub		= p.taub + l + f + axis + options

   rm(datB); rm(dattaub); rm(dattaub.long);

   ## tau.e - MCMC trace
   dattaue		= as.data.frame(cbind(1:iter.keep,Taue))
   names(dattaue) 	= c("iteration","value")
   p.taue		= ggplot(data=dattaue,aes(x=iteration,y=value))
   l			= geom_line()
   axis			= labs(x = "Iterations", y = "Sampled Value")
   options		= labs(title="MCMC Trace plot for tau.e")
   p.taue		= p.taue + l + axis + options

   ## deviance - MCMC trace
   datdev		= as.data.frame(cbind(1:iter.keep,Deviance))
   names(datdev) 	= c("iteration","value")
   p.dev		= ggplot(data=datdev,aes(x=iteration,y=value))
   l			= geom_line()
   axis			= labs(x = "Iterations", y = "Sampled Value")
   options		= labs(title="MCMC Trace plot for Deviance")
   p.dev		= p.dev + l + axis + options

   if( !is.null(M) )
   {
	
   	## M - MCMC trace
   	datM		= as.data.frame(cbind(1:iter.keep,M))
   	names(datM) 	= c("iteration","value")
   	p.M		= ggplot(data=datM,aes(x=iteration,y=value))
   	l		= geom_line()
   	axis		= labs(x = "Iterations", y = "Sampled Value")
   	options		= labs(title="MCMC Trace plot for M")
   	p.M		= p.M + l + axis + options
   }

   ##
   ## return plot objects of class ggplot
   ##

   if( !is.null(u.summary) )  ## DP + MM models
   {
	if( is.list(u.summary) ) ## multiple MM terms
	{
		return(list(p.U = p.U, p.b = p.b, p.tauu = p.tauu, p.taue = p.taue, p.taub = p.taub, p.dev = p.dev, p.M = p.M))
	}else{ ## univariate or multivariate MM term
		if( Nmv > 1 ) ## multivariate MM term
		{
			return(list(p.U = p.U, p.Umm = p.Umm, p.Ub = p.Ub, p.b = p.b, p.tauu = p.tauu, p.taue = p.taue, p.taub = p.taub, p.dev = p.dev, p.M = p.M))
		}else{ ## univariate MM term
			return(list(p.U = p.U, p.Umm = p.Umm, p.Ub0 = p.Ub0, p.b = p.b, p.tauu = p.tauu, p.taue = p.taue, p.taub = p.taub, p.dev = p.dev, p.M = p.M))
		}
	}
   }else{ ## either DP or LGM
	if( !is.null(M) ) ## DP
        {
		return(list(p.b = p.b, p.taue = p.taue, p.taub = p.taub, p.dev = p.dev, p.M = p.M))
        }else{ ## LGM
		return(list(p.b = p.b, p.taue = p.taue, p.taub = p.taub, p.dev = p.dev))
	}
   } ## end conditional fork on returning plot objects

   value <- iteration <- variable <- order <- cluster <- NULL; rm(value); rm(iteration);


} ## end function mcmcPlots