#' Plot comparison of Effect parameters of a Multiple membership (MM) term under varied prior formulations
#'
#' Produces an overlaid plot of multiple membership effect values under chosen prior formulations.  The effect values
#' for a focus (or all) multiple membership term(s) or group(s) may be compared under different model formulations (one vs. multiple terms, or
#' multivariate vs. univariate) for differing prior formulations (e.g. "mmcar" vs. "mmi").
#'
#' @param objects A list object where each entry is an object with class within any of  either \code{dpgrowmm} or \code{dgrowmult}.
#'	returned from a model run.  The list length is equal to the number of models on which a comparison of multiple membership term effects is desired.
#' @param mm.terms An optional character vector of MM term names or numeric vector of MM term positions targeted for plotting.
#'		   The entries a term name for each \code{dpgrowmult} object and group name for each \code{dpgrowmm} object.  
#'		   If a character vector is used, names must correspond to those used in package engine 
#'		   functions that generated each object.  Length is equal to number of \code{objects}.   
#' @param prior.labs An optional character vector providing names for the prior formulations under comparison.  Length is equal to number of objects.
#' @param center A scalar boolean input to indicate whether effect values in all objects should be centered.  Defaults to \code{TRUE}.
#' @param axis.labs An optional vector input of length \code{2} corresponding to \code{x.axis} and \code{y.axis} labels, respectively. 
#' @param group.plot An optional numeric vector that groups or re-groups the MM effects in each term selected in \code{mm.terms}.
#' @param smoother A scalar boolean input indicating whether to co-plot a smoother line with point values.
#' @param order A scalar boolean input indicating whether to sort the plots of effects in order of increasing (mean) value.
#' @param orderto A list object of length equal to number of MM terms (for \code{dpgrowmult} objects), 
#'	or groups (for \code{dpgrowmm} objects), each holding a vector of length \code{S[n]}, the number of MM effects for object \code{n}
#'	that conveys a desired order for the MM effects under option \code{order = TRUE} as an alternative to ordering by increasing mean value.
#' @param label.mm An optional list object, each containing a vector of labels for the MM random effects in each term (\code{dpgrowmult}) or each group (\code{dpgrowmm}).
#'			Used for plotting all terms or groups for input objects.   The length of \code{label.mm} is equal to the number of terms or groups per object.
#' @return A list object containing a plot objects and associated plot data.frames
#'     \item{dat.term}{A \code{data.frame} object used to generate a plot comparing a specific term or group from an MM object
#'		with field titles, \code{c("session","group","low","mean","high","prior")}, where \code{session} denotes the
#'		MM effect labels, \code{group} the sub-group labels within the focus MM terms, and \code{prior}, the model identifiers.}
#'     \item{p.term}{A \code{ggplot2} object of effect plots for chosen object terms.}
#'	\item{p.all}{A \code{ggplot2} plots effects, by term and/or compatible group, for each object.  Only renders if the number of effects in terms and/or groups are compatible between objects.}
#'	 \item{dat.all}{A \code{data.frame} object used to generate \code{p.all} 
#'		with field titles, \code{c("prior","session","block","low","mean","high")},
#'		where \code{block} denotes labels for MM terms or groups included.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrowmult}}, \code{\link{ddpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases effectsplot
#' @export
effectsplot <- function(objects, mm.terms = NULL, prior.labs = NULL, center = TRUE, axis.labs = NULL, group.plot = NULL, smoother = FALSE, order = FALSE, orderto = NULL, label.mm = NULL)
{
  ## check inputs
  if( is.null(objects) ) stop("\n Must input 'objects' to plot. \n")
  if( !is.list(objects) ) stop("\n 'objects' must be input as a list.\n")
  nty 			= length(objects)
  if( !is.null(prior.labs) )
  {
  	if( length(prior.labs) != nty ) stop("\n Number of objects and elements of prior.labs must be equal.\n")
  }else{
	prior.labs = 1:nty ## prior formulations are sequentially numbered if omitted by the user
  }

  if( !is.null(mm.terms) )
  {
	if( length(unique(class(mm.terms))) > 1 ) ## mixture of character and numeric
	{
		stop("\n All entries of 'mm.terms' must be of class 'character' or 'numeric', not both.\n")
	}

  }

  if( !is.null(group.plot) & is.null(mm.terms) ) warning("\n The 'group.plot' input will be ignored when not employing 'mm.terms'.")

  if( is.null(axis.labs) )
  {
      	axis.labs 	= as.character(axis.labs)
	axis.labs[1] 	= "Effect Number"
	axis.labs[2] 	= "Effect Size"
  }


  ##
  ##  plot specific term for "dpgrowmult" objects or a specific group for "dpgrowmm" objects
  ## 

  term.flag <- 0
  if( !(is.null(mm.terms))  )
  {
	term.flag <- 1

	## memo: mm.terms is either a numeric vector of term positions or a character vector of term names
  	if( !(is.null(mm.terms)) )
  	{
		if( length(mm.terms) != nty ) stop("Length of 'mm.terms' must equal number of objects.")

  		## convert term label to position (among sampled terms) where one or more objects is of class "dpgrowmult"
  		if( is.character(mm.terms) )
  		{
			mm.nterms 	= vector("numeric", nty)
			for( i in 1:nty )
			{
				if( class(objects[[i]]) == "dpgrowmult" )
				{
					ulabs		= objects[[i]]$summary.results$ulabs
				}else{ ## class(object) = "dpgrowmm", univariate or multivariate
					ulabs		= unique(objects[[i]]$summary.results$map.grp$label.input) ## unique group labels for dpgrowmm objects
				}
				mm.nterms[i]	= which(ulabs == mm.terms[i])
			}
  		}else{  ## is.numeric(mm.terms) == TRUE
			mm.nterms = mm.terms
		}
  	}

  	## build plot data.frame
  	groups 			= vector("list",nty)
  	us.summary 		= vector("list",nty)
	Nmv			= vector("numeric",nty)
  	term			= 1
  	for( i in 1:nty)
  	{
		if( class(objects[[i]]) == "dpgrowmult" )
		{
			groups[[i]]		= objects[[i]]$summary.results$group.i[[mm.nterms[i]]] ## of length equal to number of sessions
			us.summary[[i]]		= objects[[i]]$summary.results$u.summary[[mm.nterms[i]]]
			Nmv[i]			= 1
		}else{ ## object is of class = "dpgrowmm" 
			all.groups.i		= as.matrix(objects[[i]]$summary.results$map.grp$label.new)  ## numeric, converted from user-input label. "as.matrix" matches groups[[i]].
			groups[[i]]		= all.groups.i[all.groups.i == mm.nterms[i]] 
			all.us.i		= objects[[i]]$summary.results$u.summary
			Nmv[i]			= objects[[i]]$summary.results$Nmv
			if( is.null(Nmv[i]) ) stop("\n Object produced by growcurves version 0.1.0.  Re-run under current version to used 'effectsplot'.")
			all.groups.i		= rep(all.groups.i,times=Nmv[i])
			pos.sel.i		= which(all.groups.i %in% mm.nterms[i])
			us.summary[[i]]		= all.us.i[pos.sel.i,]
		}
  	}

	## check consistency of inputs
  	## test to ensure group assignments are identical for all sets of session effects since comparison is restricted to the same term
  	for( i in 1:nty )
  	{
		j = i + 1
		while( j <= nty )
		{
			if( nrow(us.summary[[i]])/Nmv[i] != nrow(us.summary[[j]])/Nmv[j] ) stop("\nCan only compare session effects values for the same multiple membership term.\n")
			if( !identical(groups[[i]],groups[[j]]) ) ## one of the dpgrowmm objects left with default grouping - a vector of single values.
			{
				if(length(unique(groups[[i]])) > 1) ## use the object with multiple group identifier for plotting.
				{
					groups[[j]] <- groups[[i]]
				}else{
					groups[[i]] <- groups[[j]]
				}
			}
			j = j + 1
		}
  	}

  	sessions		= 1:(nrow(us.summary[[1]])/Nmv[1])
  	if( !is.null(group.plot) ) 
  	{
		if( length(group.plot) != nrow(us.summary[[1]])/Nmv[1] ) stop("\nLength of 'group.plot' must be equal to number of effects (of any order) for each term.")
  		group 			= group.plot
  	}else{
		group			= groups[[1]]
  	}
  	dat				= vector("list", nty)
  	for( i in 1:nty )
  	{
		sessions.i			= rep(sessions,times=Nmv[i])
		group.i				= rep(group,times=Nmv[i])
		order.i				= rep(1:Nmv[i],each = length(sessions))
  		dat[[i]] 			= data.frame(sessions.i,group.i,order.i,us.summary[[i]],prior.labs[i])
  		names(dat[[i]])			= c("session","groups","order","low","mean","high","prior") ## use name 'groups' to not confuse ggplot group factor
		if( center == TRUE )
		{
			dat[[i]]			= split(dat[[i]],list(dat[[i]]$groups,dat[[i]]$order)) ## group mean center
  			dat[[i]]			= lapply(dat[[i]],function(x){
								x$low = x$low - mean(x$low)
								x$mean = x$mean - mean(x$mean)
								x$high = x$high - mean(x$high)
								x})
			dat[[i]]			= do.call("rbind", dat[[i]])
		}
  	}
  	dat.term		= do.call("rbind", dat)

  	## output plot

  	## shows session effects for each prior on same within group plot
	if( order == TRUE )
	{
		## first create variable, session_groups, that confines reordering to within group.
		dat.term			<- transform( dat.term, session_groups = factor(paste(session,groups,sep=".")) )
  		p.term				= ggplot(data=dat.term,aes(x=reorder(session_groups,mean), y = mean, colour = prior))
	}else{ ## don't order
		p.term				= ggplot(data=dat.term,aes(x=session, y = mean, colour = prior))
	}
	if( all(Nmv == 1) )
	{
		l 			= geom_point(aes(shape = prior))
		l.2			= geom_smooth(aes(group=prior, colour = prior), method = "loess", span = 1.0, se = FALSE)##, colour = "pink")
  		axis	 		= labs(x = axis.labs[1], y = axis.labs[2], shape = "Prior", colour = "Prior")
	}else{ ## at least one object is multivariate
		l			= geom_point(aes(shape = interaction(prior,order), colour = interaction(prior,order)))
		l.2			= geom_smooth(aes(group=interaction(prior,order), colour = interaction(prior,order)), method = "loess", span = 1.0, se = FALSE)##, colour = "pink")
  		axis	 		= labs(x = axis.labs[1], y = axis.labs[2], shape = "Prior-Order", colour = "Prior-Order")
	}
  	f				= facet_wrap(~groups, scales="free")
	
	if( length(label.mm) == 1 ) ## entered for a single term
	{
		if( length(label.mm[[1]]) != length(sessions) ) ## not entered for this term
		{
			labels.x = sessions
		}else{
			labels.x = label.mm[[1]]
		}
	}else{ ## no supplemental labels entered for just this term
		labels.x = sessions
	}
	if( order == FALSE )
	{
		options 			<- scale_x_discrete(labels=labels.x)
	}else{ ## reorder sessions within group - although plotting session_block, set x-axis tick labels to session (within block)
		options 			<- scale_x_discrete(labels=labels.x, breaks=dat.term$session_groups)
	}
  	if( smoother == TRUE)
  	{
  		p.term			= p.term + l + l.2 + f + axis + options
  	}else{
		p.term			= p.term + l + f + axis + options
  	}
  	dev.new()
  	print(p.term)

 } ## end conditional statement on !(is.null(mm.terms)) | !(is.null(group.plot))



  ##
  ## if (dpgrowmult or dpgrowmm) objects have same number of treatments and dosages per treatment, then plot effects from all of the treatments across objects.
  ## convert "group" structure for dpgrowmm objects to "block" structure of dpgrowmm objects.
  ##

  all.flag		<- 1 ## start with assumption that all objects have identical treatments and dosages

  ## check for identical number of treatments and dosages/treatment across objects
  ulabs.block		<- vector("list",nty) 
  numt.block		<- vector("list",nty)
  Nmv			<- vector("numeric",nty)
  for( i in 1:nty)
  {
	if( class(objects[[i]]) == "dpgrowmult" )
	{
		ulabs.block[[i]]	<- summary(objects[[i]])$summary.results$ulabs
		numt.block[[i]]		<- summary(objects[[i]])$summary.results$numt
		Nmv[i]			<- 1
	}else{ ## object of class "dpgrowmm" 
		ulabs.block[[i]]	<- unique(objects[[i]]$summary.results$map.grp$label.input) ## label.input maps labels to each effect, so take unique values
		numt.block[[i]]		<- as.numeric(table(objects[[i]]$summary.results$map.grp$label.input)) ## number of session identifiers per group.  we assign them to a block so they may be plotted.
		Nmv[i] 			<- objects[[i]]$summary.results$Nmv  ## object MM effects are univariate if Nmv = 1, else multivariate
	}
  }

  for( i in 1:nty )
  {
	j = i + 1
	while( j <= nty )
	{
		if( !identical(numt.block[[i]],numt.block[[j]]) ) all.flag = 0 ## not the same dosages within treatment (MM matrix)
		j = j + 1
	}
  }

  if( all.flag == 1 ) ## identical length treatments and dosages within treatments
  {
	## since treatment length and numbers are identical across objects, use the first object values
	ulabs.char		<- sapply(ulabs.block,class)
	pos.char		<- which(ulabs.char == "character")
	if( length(pos.char) > 0 )
	{
	 	ulabs			= ulabs.block[[pos.char[1]]]
	}else{ ## none of the labels is of type character
		ulabs			= ulabs.block[[1]]
	}
	numt			= numt.block[[1]]
	block			= factor(ulabs) 
	nty.n			= length(numt) ## length of treatments within each object (versus number of objects, nty).

	## input dummy value for orderto if not input by user
	if( is.null(orderto) ) 
  	{
		orderto		<- vector("list",nty.n)
		for( i in 1:nty.n){ orderto[[i]] <- rep(0,numt[i]) }
  	}
	
	## check whether orderto input conforms to term or group structure
	for( i in 1:nty.n )
	{
		if( length(orderto[[i]]) != numt[i] ) stop("Length of each vector element in 'orderto' must be equal to number of MM effects for the MM term or group.")
	}

	## check input label.mm to relabel x.axis tick labels on the faceted plots
	if( !is.null(label.mm)  )
  	{
		if( !is.list(label.mm) ) stop("\n Input 'label.mm' as a list object with length equal to number of terms per object.")
		if( length(label.mm) != length(numt) ) stop("\n Length of list object, label.mm, must be equal to number of terms or groups per object.")
  	}

	dat.all			= vector("list",nty)
	for( n in 1:nty )
	{
		if( class(objects[[n]]) == "dpgrowmm" )
		{
			## Convert dpgrowmm objects to format of dpgrowmult.  Use for plotting of all groups where the groups of a "dpgrowmm" object correspond to the terms of a "dpgrowmult" object.
			u.summary.n		= vector("list",nty.n)
			## split u.summary into list objects corresponding to blocks - to be consistent with "dpgrowmult" objects.
			grouping.n		<- rep(objects[[n]]$summary.results$map.grp$label.input, times = Nmv[n]) ## will duplicate effect sets by order for multivariate objects (where Nmv[n] > 1)
			for( t in 1:nty.n ) 
			{
				posn.t			<- which(grouping.n %in% ulabs.block[[n]][t])
				u.summary.n[[t]]	<- objects[[n]]$summary.results$u.summary[posn.t,]
			}
		}else{ ## objects are of class "dpgrowmult"
			u.summary.n		= objects[[n]]$summary.results$u.summary
		}
	
		if( center == TRUE )
		{
			if( Nmv[n] == 1 ) ## u.summary is univariate
			{
				for(i in 1:nty.n) ## number of "treatments "(groups) per block
				{
					mn			<- colMeans(u.summary.n[[i]])
					u.summary.n[[i]] 	= apply(u.summary.n[[i]],1,function(x){
									x <- x - mn
									x})
					u.summary.n[[i]]	<- t(u.summary.n[[i]])
				}
			}else{ ## u.summary is multivariate - so center by order
				for(i in 1:nty.n) ## number of "treatments "(groups) per block
				{
					Nsessionn.i		<- nrow(u.summary.n[[i]])/Nmv[n]
					for( m in 1:Nmv[n])
					{
						indx.m				<- (m-1)*Nsessionn.i + 1:Nsessionn.i
						u.summaryni.m			<- u.summary.n[[i]][indx.m,]
						mn.m				<- colMeans(u.summaryni.m)
						u.summaryni.m 			<- apply(u.summaryni.m,1,function(x){
											x <- x - mn.m
											x})
						u.summaryni.m			<- t(u.summaryni.m)
						u.summary.n[[i]][indx.m,]	<- u.summaryni.m
					} ## end loop m over Nmv[n] orders - separately centering for each order
				} ## end loop over nty.n dosages within block n
			} ## end conditional statement on whether u.summary is univariate or multivariate
		} ## end conditional statement on whether to center effects

		dat			= vector(mode = "list", length = nty.n)
		for(i in 1:nty.n)
		{
			Nsessionn.i		<- nrow(u.summary.n[[i]])/Nmv[n]
			ordern.i		<- rep(1:Nmv[n],each=Nsessionn.i)
			if( !is.null(label.mm)  )
			{
				if(length(label.mm[[i]]) != Nsessionn.i){stop("Each vector entry in 'label.mm'  must have same length as the number of MM effects in that term or group.")}
				session.i		<- rep(label.mm[[i]],times=Nmv[n])
			}else{
				session.i		<- rep(1:Nsessionn.i,times=Nmv[n])
			}
			orderto.i		<- rep(orderto[[i]],times=Nmv[n])
			dat[[i]]		= data.frame(prior.labs[n],session.i,ordern.i,block[i],orderto.i,u.summary.n[[i]]) ## using treat to facet plot
			names(dat[[i]])		= c("prior","session","order","block","ordering","low","mean","high")
   		} ## end loop i to create long data.frame with label identifier for each bi

		
      	dat.all[[n]]					= as.data.frame(do.call("rbind",dat))
	} ## end loop n over nty treatments or blocks

	dat.all						= as.data.frame(do.call("rbind",dat.all))
	## plot
	if( order == TRUE)
	{
		if( sum(orderto[[i]]) == 0 )
		{
			## first create variable, session_block, that confines reordering to within block.
			dat.all				<- transform( dat.all, session_block = factor(paste(session,block,sep=".")) )
			## although plotting session_block, set x-axis tick labels to session (within block)
			options 			<- scale_x_discrete(labels=dat.all$session, breaks=dat.all$session_block)
   			p.all		 		= ggplot(data=dat.all,aes(x=reorder(session_block,mean), y=mean)) + options
		}else{
			p.all		 		= ggplot(data=dat.all,aes(x=reorder(session,ordering), y=mean))
		}
	}else{ ## don't order
		p.all		 		= ggplot(data=dat.all,aes(x=session, y=mean))
	}
	if( all(Nmv[n] == 1) ) ## no multivariate objects
	{
   		l 				= geom_point(aes(shape=factor(prior)), size = 2.3)
		l.2				= geom_smooth(aes(group=prior), method = "loess", span = 1.0, se = FALSE)
		axis	 			= labs(x = axis.labs[1], y = axis.labs[2], shape = "Prior")
	}else{ ## at least one multivariate object
		l				= geom_point(aes(shape=interaction(prior,order), colour = interaction(prior,order)),size = 2.3)
		l.2				= geom_smooth(aes(group=interaction(prior,order),colour = interaction(prior,order)), method = "loess", span = 1.0, se = FALSE)
		axis	 			= labs(x = axis.labs[1], y = axis.labs[2], shape = "Prior-Order", colour = "Prior-Order")
	}
   	f				= facet_wrap(~block,scales="free",ncol=2)
	if( smoother == TRUE )
	{
   		p.all		 		= p.all + l + l.2 + f + axis ## + scale_x_continuous( breaks = seq(1:max(dat.all$session)) )
	}else{
		p.all		 		= p.all + l + f + axis ## + scale_x_continuous( breaks = seq(1:max(dat.all$session)) )
	}
	dev.new()
  	print(p.all)

 	if( term.flag == 1 )
	{
		return(invisible(list(p.term = p.term, dat.term = dat.term, p.all = p.all, dat.all = dat.all)))
	}else{
		return(invisible(list(p.all = p.all, dat.all = dat.all)))
	}

  }else{ ## all.flag == 0
	if( term.flag == 1 )
	{
		return(invisible(list(p.term = p.term, dat.term = dat.term)))
	}
  }

  session <- session_groups <- session_block <- group <- block <- ordern.i <- ordering <- session.i <- Nsessionn.i <- sessions.i <- group.i <- order.i <- low <- mean <- high <- prior <- NULL; rm(session); rm(group); rm(low); rm(mean); rm(high); rm(prior); rm(block)

  ## gc()

} ## end of function effectsplot