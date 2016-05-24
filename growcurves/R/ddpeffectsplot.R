#' Plot by-subject and by-treatment posterior mean values for dossage random effects
#'
#' Each \code{ddpgrow} object contains posterior mean estimates for the  \code{q} x \code{T} 
#' matrix of by-subject random effects that is extracted from the \code{ddpgrow}
#' object that is input to \code{ddpEffectsplot}.
#' This function produces a \code{q} x \code{T.m} heat map plot of posterior mean effect values 
#' for the dosages in treatment \code{m} faceted on a set of chosen subjects.  The resulting plot
#' produces a heatmap for each trt-subject combination.  Both a ggplot2 plot object
#' and a \code{data.frame} object are returned.
#'
#' @param object A \code{ddpgrow} object. 
#' @param subjects.plot A vector of subjects for performing plots that is composed of some subset of the \code{subject} vector input for modeling.
#'			If left blank, a random subset is chosen from \code{subject}.
#' @param n.plot An optional scalar input for number of randomly generated subjects to plot (if \code{subjects.plot} is left blank).
#' @param trts.plot A vector of focus treatments to use for plotting.  
#' @param x.axis.label An optional scalar character entry to label the treatment(s) dosages 
#' @param smoother  A scalar boolean input indicating whether to co-plot a smoother line with point values.
#' @param re.order A scalar boolean input indicating whether to sort the plots of effects in order of increasing value.
#' @param cred.intervals A boolean scalar indicating whether the by-subject effects plots should include credible intervals.
#' @param map.group A \code{matrix} or \code{data.frame} object containing a grouping of subjects that will be used
#'			to produce an additional set of effect plots that aggregate subjects by the grouping structure.
#'			The first column containing subject identifiers for all subjects modeled in \code{object}.
#'			The second column contains the desired desired group identifiers that may be of type character or numeric.
#' @param n.dose.plot Optional numeric input for number of randomly chosen doses for which to plot effects growth curves.
#' @param orderto A numeric vector of length equal to the total number of dosages across all treatments that conveys an order to be used
#'			by-dosage growth curve plots within cluster and treatment.
#' @return A list object containing a faceted set of heat maps (one per subject),
#'		a faceted set of effect point plots, and the associated \code{data.frame} objct.
#'     \item{dat.se}{A \code{data.frame} object used to generate the trt-subject faceted plots for effect means
#'			Fields are titled, \code{c("order","dose","trt","subject","effects")}.}
#'     \item{dat.ci}{A \code{data.frame} object used to generate the trt-subject faceted plots for effect credible intervals
#'			Fields are titled, \code{c("order","dose","trt","subject","quantile","effects")}.}
#'     \item{dat.clust}{A \code{data.frame} object used to generate the trt-group faceted plots for effect means
#'			Fields are titled, \code{c("order","dose","trt","cluster","effects")}.}
#'     \item{dat.clust.ci}{A \code{data.frame} object used to generate the trt-group faceted plots for effect credible intervals
#'			Fields are titled, \code{c("order","dose","trt","cluster","quantile","effects")}.}
#'     \item{dat.gc}{A \code{data.frame} object used to generate the by-dose growth curves (for multivariate polynomial effects)
#'			Fields are titled, \code{c("fit","time","cluster","trt","dose")}.}
#'     \item{p.hm}{A \code{ggplot2} object of heat maps for mean random effect values, faceted by trt and subject combinations.}
#'     \item{p.pp}{A \code{ggplot2} object of point plots for mean random effect values or credible intervals, faceted by trt and subject combinations.}
#'     \item{pc.m}{A \code{ggplot2} object of point plots for mean random effect values, faceted by trt and group combinations.}
#'     \item{pc.ci}{A \code{ggplot2} object of point plots for credible intervals of random effect values, faceted by trt and cluster combinations.}
#'     \item{pc.gc}{A \code{ggplot2} growth curve plots for each dose where the does effects are multivariate polynomial.}
#' @seealso \code{\link{ddpgrow}}, \code{\link{dpgrow}}, \code{\link{dpgrowmm}}, \code{\link{dpgrowmult}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases ddpEffectsplot
#' @export
ddpEffectsplot <- function(object, subjects.plot = NULL, n.plot = 3, trts.plot = NULL, x.axis.label = NULL, smoother = TRUE, re.order = TRUE, cred.intervals = TRUE, map.group = NULL, n.dose.plot = 5, orderto = NULL)
{
  ## check if object is of class ddpgrow
  if( class(object) != "ddpgrow" ) stop("\nThe input 'object' must be of class 'ddpgrow' to use this function.\n")

  ## check if subjects.plot are in subjects modeled
  map.subject 		= summary(object)$summary.results$map.subject ## data.frame
  map.subject		= unique(map.subject)
  subjects.input 	= map.subject$label.input
  if( is.null(subjects.plot) ) ## is.null(subjects.plot) = TRUE, randomly sample subject set
  {
	  subjects.num	= sample(1:nrow(map.subject), n.plot, replace = FALSE)
	  subjects.plot	= subjects.num ## labels are same as numerical count of subjects used to model
  }else{ ## check that subjects.plot is a strict subset of subjects.input
  	if( length(setdiff(subjects.plot,subjects.input)) > 0 ) stop("\nSubjects to plot must be a subset of those used in model.\n")
	    tmp		= subset(map.subject, label.input %in% subjects.plot)
	  subjects.num	= tmp$label.new ## user inputs names, we extract numerical equivalent used in modeling
	  rm(tmp)
  }
  ## if( !is.null(map.group) ), check that subjects mapped are in subjects modeled
  if( !is.null(map.group) )
  {
 	  map.group		<- as.data.frame(map.group)
	  names(map.group)	<- c("label.input","cluster")
	  if( length(setdiff(map.group$label.input,subjects.input)) > 0 ) stop("\nSubjects to group for plotting must be those used to model.  Please check that first column of 'map.group' contains
										subjects and second contains group membership.\n")
	  ## replace label.input in map.group with numerical label.new
	  tmp			<- merge(map.group, map.subject, by="label.input", all.x = TRUE, sort = FALSE)
	  map.newgroup		<- data.frame(tmp$label.new,tmp$cluster)
	  names(map.newgroup)	<- c("subjects.num","cluster")
	  ## capture number of clusters
	  clusters		<- unique(map.newgroup$cluster)
	  num.clust		<- length(clusters)
  }
  ## look for conflicts on re-ordering instructions
  if( re.order == TRUE & !is.null(orderto) ) warning("\nSince re.order is set to TRUE, the orderto input will be ignored.  If want to use orderto, then set re.order to FALSE.\n")
  if( is.null(orderto) ) orderto = 1:sum(summary(object)$summary.results$numt) ## each entry in numt is the number of dosages associated to each treatment.  length(numt) is the number of treatments.

 ## extract numerical trt as user inputs label
 map.trtcov		= summary(object)$summary.results$map.trtcov ## data.frame
 map.trtcov		= unique(map.trtcov)
 if( !is.null(trts.plot) )
 {
	if( length(setdiff(map.trtcov$label.input,trts.plot)) > 0) stop("\nTreatments selected for plotting must be those used to model.\n")
  	tmp			    = subset(map.trtcov, label.input %in% trts.plot)
  	trt.num			= tmp$label.new
  	rm(tmp)
  }else{ ## is.null(trts.plot), so randomly select a set of treatments
	  typet			  = summary(object)$summary.results$typet
	  samp.treat	= min(length(typet),2)
	  trt.num			= sample(1:length(typet),samp.treat,replace=FALSE)
	  trts.plot		= unique(map.trtcov$label.input)[trt.num]
  }

  ## create plot data.frame
  Nrandom			= ncol(summary(object)$summary.results$Z)
  Nsubject			= nrow(unique(summary(object)$summary.results$map.subject))
  numt				= summary(object)$summary.results$numt
  nt				= length(trt.num)
  np				= length(subjects.num)
  dat				= vector("list",nt)
  effect.mat			= vector("list",nt)

  for(m in 1:nt)
  {
	effect.mat[[m]]	= summary(object)$summary.results$dosetrt.mean[[trt.num[m]]] ## select the trt.num[m] list element of an P x q*numt[m] matrix
	effect.sum	= summary(object)$summary.results$dosetrt.summary[[trt.num[m]]] 
	leff.m		= t(matrix(effect.sum[,1], (Nrandom*numt[trt.num[m]]), Nsubject, byrow = FALSE))
	heff.m		= t(matrix(effect.sum[,3], (Nrandom*numt[trt.num[m]]), Nsubject, byrow = FALSE))
	dose.m		= colnames(effect.mat[[m]])[1:numt[trt.num[m]]]
  	dose.m	 	= sapply(strsplit(dose.m,""),function(x){
					x <- x[-c(1,2)] ## get rid of the "q." indicator of random effects in the names
					x <- paste(x,collapse="")
					x}) ## now we just have the names of the numt[trt.num[m]] doses associated to treatment trt.num[m]
	dose.m		= rep(dose.m, Nrandom)
	order.m		= rep(1:Nrandom, each = numt[trt.num[m]])
	dat.m			= vector("list",np)
	for(i in 1:np)
	{
		effectsm.i		= effect.mat[[m]][subjects.num[i],] ## extract subject (row), subjects.num[i], from the P x q*numt[m] effects.mat
		leffm.i		= leff.m[subjects.num[i],]
		heffm.i		= heff.m[subjects.num[i],]
		dat.m[[i]]		= data.frame(order.m, dose.m, trts.plot[m], subjects.plot[i], effectsm.i, leffm.i, heffm.i, stringsAsFactors = FALSE)
 	}
	dat[[m]]	= do.call("rbind",dat.m)
  }
  dat			= do.call("rbind", dat)
  names(dat)	= c("order","dose","trt","subject","mean.effects","low.effects","high.effects")	

  ## output plots
  ## tile geom 
  if( re.order == TRUE )
  {
	## first create variable, dose_trtsubject, that confines reordering to within trt-subject.
	dat			<- transform( dat, dose_trtsubject = factor(paste(dose,trt,subject,sep=".")) )
	## although plotting dose_trtsubject, set x-axis tick labels to dose (within trt-subject)
	options 		<- scale_x_discrete(labels=dat$dose, breaks=dat$dose_trtsubject)
  	p.t			= ggplot(data=dat,aes(x=reorder(dose_trtsubject,mean.effects), y = factor(order), fill = mean.effects)) + options
  }else{
	p.t			= ggplot(data=dat,aes(x=dose, y = factor(order), fill = mean.effects))
  }

  l			= geom_tile()
  yaxis 		= ylab("Polynomial Order")
  if( is.null(x.axis.label) ) 
  {
	xaxis = xlab("Dose")
  }else xaxis = xlab(eval(x.axis.label))
  f			= facet_wrap(trt~subject, scales="free", nrow=min(nt,5) )
  dev.new()
  p.t			= p.t + l + xaxis + yaxis + f + labs(fill = "effects") + scale_fill_gradient(low = "white", high = "steelblue")
			## + scale_fill_gradient( limits = c(-3002,-2980) )
  print(p.t)

  ## point geom 
  if( cred.intervals == FALSE )
  {
	  if( re.order == TRUE )
	  {
		  ## first create variable, dose_trtsubject, that confines reordering to within trt-subject.
		  dat			<- transform( dat, dose_trtsubject = factor(paste(dose,trt,subject,sep=".")) )
		  ## although plotting dose_trtsubject, set x-axis tick labels to dose (within trt-subject)
		  options 		<- scale_x_discrete(labels=dat$dose, breaks=dat$dose_trtsubject)
  		p.p			= ggplot(data=dat,aes(x=reorder(dose_trtsubject,mean.effects), y = mean.effects, colour = factor(order), shape = factor(order))) + options
	  }else{
		  p.p			= ggplot(data=dat,aes(x=dose, y = mean.effects, colour = factor(order), shape = factor(order)))
	  }
  	l			      = geom_point()
  	l.2			    = geom_smooth(aes(group=factor(order)), method = "loess", span = 1.0, alpha = 0.1, se = FALSE)##, colour = "pink")
  	yaxis 			= ylab("Effect Values")
  	options			= labs(colour = "Order", shape="Order")
  	f			= facet_wrap(trt~subject, scales="free", nrow=min(nt,n.plot) )
  }else{
	  dat.ci			      = melt(dat, id = c("order","dose","trt","subject","dose_trtsubject"))
    dat.ci$quantile		= dat.ci$variable
	  dat.ci$effects		= as.numeric(dat.ci$value)
	  dat.ci$variable   <- dat.ci$value <- NULL

	if( re.order == TRUE )
	{
		## first create variable, dose_trtsubject, that confines reordering to within trt-subject.
		dat.ci			<- transform( dat.ci, dose_trtsubject = factor(paste(dose,trt,subject,sep=".")) )
		## although plotting dose_trtsubject, set x-axis tick labels to dose (within trt-subject)
		options 		<- scale_x_discrete(labels=dat.ci$dose, breaks=dat.ci$dose_trtsubject)
		p.p			= ggplot(data=dat.ci,aes(x=reorder(dose_trtsubject,effects), y = effects, colour = factor(order), shape = factor(order))) + options
	}else{
		p.p			= ggplot(data=dat.ci,aes(x=dose, y = effects, colour = factor(order), shape = factor(order)))
	}
  	l			= geom_point(na.rm=TRUE)
  	## l = geom_line
  	l.2			= geom_smooth(aes(group=factor(order)), method = "loess", span = 1.0, alpha = 0.1, 
  	                    se = FALSE, na.rm=TRUE)##, colour = "pink")
  	yaxis 			= ylab("Effect Values")
  	options			= labs(colour = "Order", shape="Order")
  	f			= facet_wrap(trt~subject, scales="free", nrow=min(nt,n.plot) )
  }	
  if( smoother == TRUE)
  {
  	p.p			= p.p + l + l.2 + xaxis + yaxis + f + options + theme(axis.text.x=element_text(size=6, angle = 90))
  }else{
	p.p			= p.p + l + xaxis + yaxis + f + options + theme(axis.text.x=element_text(size=6, angle = 90))
  }
  dev.new()
  print(p.p)	

  if(!is.null(map.group)) ## plot subject effects aggregated to group
  {
	dat.clust				= vector("list",nt)
	dose.lab				= vector("list",nt)
	dose.ord				= vector("list",nt)
	md					<- sapply(1:nt,function(m){md = numt[trt.num[m]]}) ## vector of dosages for each plotted treatment.
	start.m					<- 1 ## to produce user-input ordering of dosages for dosage growth curve plots (within cluster and treatment).
	for( m in 1:nt )
	{
		## develop posterior draws for cluster effects means
		DoseEffects.m			<- object$DoseEffects[[trt.num[m]]] ## B x Nsubject*(Nrandom*numt[m])
		tDE.m					<- t(DoseEffects.m)
		## some dimensions
		step.m				<- Nrandom*numt[trt.num[m]] ## effects per subject (and cluster once we average over subjects_
		B					<- nrow(DoseEffects.m)
		## some structures
		tDC.m				<- matrix(0,(num.clust*step.m),B)
		for( c in 1:num.clust )
		{
			subj.c		<- subset(map.newgroup, cluster == clusters[c])$subjects.num
			nums.c		<- length(subj.c)
			cols.c		<- vector("numeric",(nums.c*step.m))
			for( s in 1:nums.c )
			{
				indx.s		<- ((s-1)*step.m + 1):(s*step.m)
				cols.c[indx.s]	<- c( ((subj.c[s]-1)*step.m + 1):(subj.c[s]*step.m) )
			}
			indx.c			<- ((c-1)*step.m + 1):(c*step.m)
			## we use 'rowsum' to chunk sum for each member in step.m and then transpose the result
			chunk.c			<- rep(1:step.m, times = nums.c) ## recall subjects are the slow moving index
			## composing average of the step.m block of effects by subjects in cluster c for each posterior sample in B
			tDC.m[indx.c,]		<- rowsum(tDE.m[cols.c,],group = chunk.c) * (1/nums.c) ## num.clust*step.m x B
		}
		DoseCluster.m			<- t(tDC.m) ## B x num.clust*step.m
		## compose quantile summaries
		mcmc.data			= DoseCluster.m
		mcmc.mean			= colMeans(mcmc.data)
  		mcmc.low			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.025)})
 		mcmc.high			= apply(mcmc.data,2,function(x){quantile(x,probs = 0.975)})
        	## for doseclust.summary.m, cluster is slowest moving index, followed by polynomial order, followed by numt[m]
		## <--> num.clust*(Nrandom*numt[m]) x 1, index speed order: cluster, order, dose
		## compose dose effect quantiles for the clusters - has num.clust rows (which is probably pretty small)
  		doseclust.summary.m		= t(rbind(mcmc.low,mcmc.mean,mcmc.high)) ## transpose to 3 columns for c("mcmc.low","mcmc.mean","mcmc.high") 
		## define indices
		order.m				= rep(1:Nrandom, each = numt[trt.num[m]], times = num.clust)  ## of polynomial order
		dose.m				= colnames(effect.mat[[m]])[1:numt[trt.num[m]]] ## use effect.mat from above b/c contains dose names
  		dose.mlab	 			= sapply(strsplit(dose.m,""),function(x){
								x <- x[-c(1,2)] ## get rid of the "q." indicator of random effects in the names
								x <- paste(x,collapse="")
							x}) ## now we just have the names of the numt[trt.num[m]] doses associated to treatment trt.num[m]
		end.m					<- start.m + md[m] - 1
		dose.mord				= orderto[start.m:end.m] ## ordering dose names across treatments
		start.m				<- end.m + 1
		dose.m				= rep(dose.mlab, times = (Nrandom*num.clust))
		clust.m				= rep(clusters,each=(Nrandom*numt[trt.num[m]]))
		dat.clust[[m]]			<- data.frame(order.m, dose.m, trts.plot[m], clust.m, doseclust.summary.m, stringsAsFactors = FALSE)	
		dose.lab[[m]]			<- dose.mlab
		dose.ord[[m]]			<- dose.mord
	} ## end loop over nt treatments
  	dat.clust		= do.call("rbind", dat.clust)
  	names(dat.clust)	= c("order","dose","trt","cluster","low.effects","mean.effects","high.effects")	

	## perform point plotting
	## of means
	if( re.order == TRUE )
	{
		## first create variable, dose_trtcluster, that confines reordering to within trt-cluster.
		dat.clust		<- transform( dat.clust, dose_trtclust = factor(paste(dose,trt,cluster,sep=".")) )
		## although plotting dose_trtclust, set x-axis tick labels to dose (within trt-cluster)
		options 		<- scale_x_discrete(labels=dat.clust$dose, breaks=dat.clust$dose_trtclust)
		pc.m			= ggplot(data=dat.clust,aes(x=reorder(dose_trtclust,mean.effects), y = mean.effects, colour = factor(order), shape = factor(order))) + options
	}else{
		pc.m			= ggplot(data=dat.clust,aes(x = dose, y = mean.effects, colour = factor(order), shape = factor(order)))
	}
  	l				  = geom_point()
  	l.2				= geom_smooth(aes(linetype = factor(order), group=factor(order)), method = "loess", span = 1.0, alpha = 0.1, se = FALSE)##, colour = "pink")
  	yaxis 			= ylab("Effect Values")
  	options			= labs(colour = "Order", shape="Order",linetype = "Order")
  	f			= facet_wrap( trt~cluster, scales = "free_x", nrow=nt )
	if( smoother == TRUE)
  	{
  		pc.m			= pc.m + l + l.2 + xaxis + yaxis + f + options + theme_bw() + theme(axis.text.x=element_text(size=6, angle = 90)) 
  	}else{
		pc.m			= pc.m + l + xaxis + yaxis + f + options + theme_bw() + theme(axis.text.x=element_text(size=6, angle = 90)) 
  	}
  	dev.new()
  	print(pc.m)


	if( re.order == TRUE )
	{
		## of credible intervals, by cluster
		datclust.ci					= melt(dat.clust, id = c("order","dose","trt","cluster","dose_trtclust"))
      		datclust.ci$quantile				= datclust.ci$variable
		datclust.ci$effects				= datclust.ci$value
		datclust.ci$variable <- datclust.ci$value <- NULL

		## although plotting dose_trtclust, set x-axis tick labels to dose (within trt-cluster)
		options 		<- scale_x_discrete(labels=datclust.ci$dose, breaks=datclust.ci$dose_trtclust)
		pc.ci			= ggplot(data=datclust.ci,aes(x=reorder(dose_trtclust,effects), y = effects, colour = factor(order), shape = factor(order))) + options
	}else{
		## of credible intervals, by cluster
		datclust.ci					= melt(dat.clust, id = c("order","dose","trt","cluster"))
    datclust.ci$quantile				= datclust.ci$variable
		datclust.ci$effects				= datclust.ci$value
		datclust.ci$variable <- datclust.ci$value <- NULL

		pc.ci			= ggplot(data=datclust.ci,aes(x = dose, y = effects, colour = factor(order), shape = factor(order)))
	}
  	l			= geom_line()
  	l.2			= geom_smooth(aes(group=factor(order)), method = "loess", span = 1.0, alpha = 0.1, se = FALSE)##, colour = "pink")
  	yaxis 			= ylab("Effects")
  	options.2		= labs(colour = "Order", shape="Order")
  	f			= facet_wrap( trt~cluster, scales= "free_x", nrow=nt )
	if( smoother == TRUE)
  	{
  		pc.ci			= pc.ci + l + l.2 + xaxis + yaxis + f + options.2 + theme(axis.text.x=element_text(size=6, angle = 90))
  	}else{
		pc.ci			= pc.ci + l + xaxis + yaxis + f + options.2 + theme(axis.text.x=element_text(size=6, angle = 90))
  	}
  	dev.new()
  	print(pc.ci)	
	
	if( Nrandom > 1 ) ## more than one effect per dose
	{
		##
		## plot within sample predicted growth curves for each dose, within cluster
		##

		## set-up structures
		C				<- length(clusters)
		M				<- length(trts.plot)
		MD				<- sum( md )
		md.cum				<- cumsum( md )
		T				<- 10 ## length of time intervals
		cluster.gc			<- trt.gc <- dose.gc <- vector("character", C*MD*T)
		y.gc				<- time.gc <- dose.order <- vector("numeric", C*MD*T)
		X				<- object$summary.results$X
		measure.times			<- unique( X[,(colnames(X) == "time^1")] )
		min.T				<- min(measure.times); max.T <- max(measure.times)
		time				<- seq(min.T,max.T,length.out = T)
	
		## compute "response" from effect orders for each dose within treatment within cluster at T timempoints
		for( c in 1:C )
		{
			for( m in 1:M )
			{
				for( d in 1:md[m] )
				{
					dat.cmd				<- subset( dat.clust, (cluster %in% clusters[c]) & (trt %in% trts.plot[m]) & (dose %in% dose.lab[[m]][d]) )
					effect.cmd			<- vector("numeric", Nrandom)
					for( r in 1:Nrandom )
					{
						effect.cmd[r]		<- dat.cmd$mean.effects[dat.cmd$order == r]
					}
					for( t in 1:T )
					{
						## fill return objects
						indx.cmdt			<- (c-1)*(MD*T) + ((md.cum[m]-md[m])*T) + (d-1)*T + t
						cluster.gc[indx.cmdt]		<- as.character(clusters[c])
						dose.gc[indx.cmdt]		<- dose.lab[[m]][d]
						dose.order[indx.cmdt]		<- dose.ord[[m]][d]
						trt.gc[indx.cmdt]		<- trts.plot[m]
						time.gc[indx.cmdt]		<- time[t]

						## generate response
						y.cmdt				<- 0
						for( r in 1:Nrandom )
						{
							y.cmdt		<- y.cmdt + effect.cmd[r]*(time[t]^(r-1))
						}
						y.gc[indx.cmdt]		<- y.cmdt
					}  ## end loop over T time points
				} ## end loop over D doses 
			} ## end loop over M treatments
		} ## end loop over C clusters

  		dat.gc				= data.frame(y.gc,time.gc,cluster.gc,trt.gc,dose.gc,dose.order, stringsAsFactors = FALSE)
  		names(dat.gc) 			= c("fit","time","cluster","trt","dose","dose.order")

		## select doses to plot
		dose.plot						<- vector("list",nt)
		for(m in 1:nt) {dose.plot[[m]]		<- sample(dose.lab[[m]],min(length(dose.lab[[m]]),n.dose.plot),replace = FALSE)}
		dose.plot						<- unlist(dose.plot)
		## dose.plot					<- c(4,5,6,7,9,10,11,12,16,18,22,23,25,26,29,35,42,45,47,52)
		dat.gcplot			<- subset(dat.gc, dose %in% dose.plot)
		n.clust			<- length(unique(dat.gcplot$cluster))

		##
  		## for ALL clients - MEAN growth curve, by TREATMENT category
  		## 					
  		p.gc			= ggplot(data=dat.gcplot,aes(x=time, y = fit, colour = factor(reorder(dose,dose.order)))) ## here, looking for individual subject visibility
  		## l			= geom_line(aes(linetype = factor(dose)))
		l			= geom_line()
		l.2			= geom_point()
  		## l.2			= geom_smooth(aes(group=1),method = "loess", size = 1.1, colour = "black")
  		axis	 		= labs(x = "Time", y = "Response Units", colour = eval(x.axis.label), linetype = eval(x.axis.label), shape = eval(x.axis.label)) ## x.axis.label captures "dose" label
  		f			= facet_wrap(trt~cluster, ncol = n.clust, scales="fixed") 
  		p.gc			= p.gc + l + l.2 + f + axis + scale_colour_grey(start = 0.1, end = 0.5) + theme_bw()
		dev.new()
  		print(p.gc)
	}else{ ## Nrandom = 1, so don't plot a growth curve
		p.gc		<- NULL
		dat.gc	<- NULL
	} ## end conditional statement on Nrandom > 1 to plot by-dose effect growth curves	

  } ## end conditional statement on whether is.null(map.group)

  ## return objects
  if( is.null(map.group) )
  {
	if( cred.intervals == FALSE )
	{
		return(invisible(list(dat.se = dat, p.hm = p.t, p.pp = p.p)))
	}else{
		return(invisible(list(dat.se = dat, dat.ci = dat.ci, p.hm = p.t, p.pp = p.p)))
	}
  }else{
	if( cred.intervals == FALSE )
	{
		return(invisible(list(dat.se = dat, dat.clust = dat.clust, dat.clust.ci = datclust.ci, dat.gc = dat.gc, p.hm = p.t, p.pp = p.p, pc.m = pc.m, pc.ci = pc.ci, p.gc = p.gc)))
	}else{
		return(invisible(list(dat.se = dat, dat.ci = dat.ci, dat.clust = dat.clust, dat.clust.ci = datclust.ci, dat.gc = dat.gc, p.hm = p.t, p.pp = p.p, pc.m = pc.m, pc.ci = pc.ci, p.gc = p.gc)))
	}
  }

  subject <- subjects.num <- label.input <- mean.effects <- dose <- cluster <- fit <- trt  <- dose_trtsubject <- dose_trtclust <- effects <- quantile <- NULL  

} ## end function ddpEffectsplot