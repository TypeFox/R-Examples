#' Plot by-subject and by-group growth curves
#'
#' Produces a by-subject plot of predicted growth curves with associated data values and an aggregated by-group
#' plot of growth curves along with a smoother line for each group based on user input.  Facilates inference
#' for different growth curve patterns based on subsets of subjects and various subject groupings.
#'
#' @param object A \code{dpgrow}, \code{dgrowmm} or \code{dpgrowmult} object obtained from running the appropriate modeling function. 
#' @param compare.objects An optional list of \code{dpgrow}, \code{dgrowmm} or \code{dpgrowmm} objects to employ as comparison models for by-subject plotting
#' @param subjects.plot A vector of subjects to use for by-group plot that is composed of some subset of \code{subject}.  Defaults to all subjects modeled in \code{object}.
#' @param groups.plot A vector with associated group identifiers for \code{subjects.plot} if the grouping is different from \code{trt}
#'		used to run the \code{dpgrow} or \code{dpgrowmm} model.  The entered grouping does not have to relate to that used for modeling.
#'		Defaults to use of treatment groups modeled in \code{object} (if not entered). 
#' @param subjects.subset A vector of a subset of \code{subjects.plot} to use for by-subject plotting for readability.  The full \code{subjects.plot}
#'	set is used for by group plotting. If left blank, the full \code{subjects.plot} vector is used for by-subject plotting.
#' @param subjects.random A boolean scalar.  If \code{TRUE} a random subset (of 10) is selected from \code{subjects.plot} for by-subject plotting.
#'		Leave blank if enter \code{subjects.subset}.  Defaults to \code{TRUE}.
#' @param x.lab Optional title for x-axis.  Default = "Time"
#' @param y.lab Optional title for y-axis.  Default = "Fit"
#' @param main.label An optional character model label for \code{object} to use in by-subject plots for comparison with models from \code{compare.objects}
#' @param title.lab option plot title. A vector of 2 character entries is allowed.  The first entry is the title for the group aggregated plot.
#'	The second entry is the title for the plot of selected subject growth curves.  If \code{title.lab} contains a single entry, then it is used
#'	as the title for both plots.  Otherwise, defaults to \code{NULL}.
#' @return A list object containing a plot objects, data.frame object from which it is constructed, and data.frame with actual data values co-plotted.
#'     \item{dat.gc}{A \code{data.frame} object used to generate the within subject predicted growth curves for \code{object}.
#'			Used for aggregated plot of growth curves by group returned in \code{p.gctrt}
#'			Fields are titled, \code{c("fit","time","subject","trt")}.}
#'     \item{dat.igc}{A \code{data.frame} object containing within subject predicted growth curves under models in \code{object}
#'			and \code{compare.objects}.  Used for by-subject growth curves plot returned in \code{p.gcsub}.}
#'     \item{dat.data}{A \code{data.frame} object containing the actual data observations for plotted subjects.
#'			Field titles are the same as for \code{dat.gc}.}
#'     \item{p.gctrt}{A \code{ggplot2} object of subjects aggregated by group.}
#'     \item{p.gcsub}{A \code{ggplot2} object of subjects.}
#' @seealso \code{\link{dpgrowmm}}, \code{\link{dpgrow}}, \code{\link{dpgrowmult}}, \code{\link{growthCurve}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @aliases growplot
#' @export

growplot <- function(object, compare.objects = NULL, subjects.plot = NULL, groups.plot = NULL, subjects.subset = NULL, subjects.random = TRUE, x.lab = "Time", y.lab = "Fit", 
			main.label = "Main_Model", title.lab = NULL)
{
  ## check inputs
  if( !is.null(compare.objects) ) ## ensure class of compare.objects within allowed set
  {
       L = length(compare.objects)
       for( l in 1:L )
       {
            if(class(compare.objects[[l]]) != "dpgrow" & class(compare.objects[[l]]) != "dpgrowmm" & class(compare.objects[[l]]) != "dpgrowmult" & class(compare.objects[[l]]) != "ddpgrow") 
            {
                 stop("Class of each element of 'compare.objects' must be either 'dpgrow', 'dpgrowmm', 'dpgrowmult' or 'ddpgrow'.")
            }
       } ## end loop l over elements in compare.objects[[l]]
  } ## end conditional statement on whether is.null(compare.objects)
 
  if(class(object) != "dpgrow" & class(object) != "dpgrowmm" & class(object) != "dpgrowmult" & class(object) != "ddpgrow") stop("Class of 'object' must be either 'dpgrow', 'dpgrowmm', 'dpgrowmult' or 'ddpgrow'.")
  if( !is.null(groups.plot) & !is.null(subjects.plot) )
  {
  	if( length(groups.plot) != length(subjects.plot) ) stop("Variable 'groups.plot' must be of same length as 'subjects.plot'.")
  }
  if( !is.null(subjects.subset) )
  {
	subjects.random = FALSE
  }
  ## check if subjects.plot are in subjects modeled
  map.subject 		= summary(object)$summary.results$map.subject ## data.frame
  map.subject		= unique(map.subject)
  subjects.input 	= map.subject$label.input
  if( is.null(subjects.plot) ) ## is.null(subjects.plot) = TRUE, so set equal to all subjects (as entered by user for modeling)
  {
	subjects.plot = subjects.input
  }else{ ## check that subjects.plot is a strict subset of subjects.input
  	if( length(setdiff(subjects.plot,subjects.input)) > 0 ) stop("\nSubjects to plot must be a subset of those used in model.\n")
  }

  if( !is.null(compare.objects) )
  {
       ## check that same set of subjects modeled in 'object' are modeled in 'compare.objects'
       num.objects = length(compare.objects)
       for( m in 1:num.objects )
       {
            map.subject 		= summary(compare.objects[[m]])$summary.results$map.subject ## data.frame
            map.subject		= unique(map.subject)
            subjects.input 		= map.subject$label.input
            if( length(setdiff(subjects.plot,subjects.input)) > 0 ) stop("\nObjects in 'compare.objects' must employ the same subjects as the main 'object'.\n")
       }  
  } ## end statements on to ensure same subjects in both object and compare.objects
  

  ## check if subjects.subset are in subjects.plot
  if( !is.null(subjects.subset) )
  {
	if( length(setdiff(subjects.subset,subjects.plot)) > 0 ) 
          stop("\nsubjects.subset must be a subset of subjects.plot.\n")
  }

  ## check if labels input properly
  if(length(title.lab) < 2)
  {
	title.lab[2] = title.lab
  }


  ## extract plot data.frame - names(dat.growcurve = c("fit","time","subject","trt")
  ## dat.gc will contain 'subject' and 'trt' in the form entered by the user
  ## both will be factor variables.  So no need to translate from user to internal format.

  ## pull-in predicted and actual growth curve data for all subjects
  plot.object 	<- plot(object,plot.out = FALSE)
  dat.pred	<- plot.object$dat.growcurve
  dat.gcdata	<- plot.object$dat.gcdata

  ## subset entered data to those subjects in subjects.plot
  dat.gc 			<- subset(dat.pred, subject %in% subjects.plot)  ## predicted growth curve data
  dat.data			<- subset(dat.gcdata, subject %in% subjects.plot) ## real data

  ## create new "trt" variable for grouping if !is.null(groups.plot)
  if( !is.null(groups.plot) )
  {
	## create map of subjects.plot to groups.plot
	map.group			<- data.frame(subjects.plot,groups.plot)
      	names(map.group)		<- c("subject","trt")
	map.group			<- unique(map.group) ## in the case duplicate values were entered
  	## merge map.group with dat.pred to get new "trt" field to aggregate subjects in plots	
     dat.gc$trt 		<- NULL
	dat.data$trt		<- NULL
    	dat.gc			<- merge(dat.gc, map.group, by = "subject", all.x = T, sort = FALSE)
	dat.data			<- merge(dat.data, map.group, by = "subject", all.x = T, sort = FALSE)
  }
  
  ##
  ## construct plots of individual growth curves, faceted by "trt" group.
  ##
  p				= ggplot(data=dat.gc,aes(x=time, y = fit, group = subject)) ## no individual subject visibility
  l				= geom_line(linetype=5,colour = "slategray")
  ## l			= geom_line(colour = alpha("black",1/5),aes(linetype=5))
  l.2			= geom_smooth(aes(group=1),method = "loess", size = 1.1, colour = "black")
  axis	 		= labs(x = x.lab, y = y.lab)
  f				= facet_wrap(~trt, scales="fixed")
  if( !is.null(title.lab) )
  { 
  	options		 		= labs(title = title.lab[1])
  }else{
	options				= NULL
  }
  dev.new()
  p.gctrt		     = p + l + l.2 + f + axis + options
  print(p.gctrt)


  ##
  ## for SELECTED clients - WITH DATA - including comparison models
  ##
  
  ## prepare data.frame that combines in comparison models
  dat.igcdata     				<- dat.data
  if( !is.null(compare.objects) )
  {
       ## set first entry to main model from object
       dat.igc					<- vector("list",(L+1))
       dat.igc[[1]]				<- dat.gc
       dat.igc[[1]]$model		<- main.label
       dat.igcdata$model			<- main.label
       
       for( l in 2:(L+1)  )
       {
            plot.cobject 		<- plot(compare.objects[[l-1]],plot.out = FALSE)
            dat.igc[[l]]		<- subset(plot.cobject$dat.growcurve, subject %in% subjects.plot) 
            
            if( !is.null(groups.plot) )
            {
                 ## reset treatment variable to main model
                 dat.igc[[l]]$trt 		<- NULL
                 dat.igc[[l]]			<- merge(dat.igc[[l]], map.group, by = "subject", all.x = T, sort = FALSE)
            }
            
            ## adding model id variable
            if( !is.null(names(compare.objects)) )
            {
                 dat.igc[[l]]$model		<- names(compare.objects)[[l-1]]
            }else{
                 dat.igc[[l]]$model		<- paste("Comp_Model_",l-1,sep="")
            }
       }
       dat.igc		<- do.call("rbind", dat.igc)
  }else{ ## no comparison models
       dat.igc		     <- dat.gc
       dat.igc$model     <- eval( main.label )
  } ## end conditional statement on !is.null(compare.objects)
  
  ## reorder model such that main.label is first
  dat.igc$model		<- factor(dat.igc$model)
  dat.igc$model		<- relevel(dat.igc$model, ref = main.label)
  
  if( !is.null(subjects.subset) | subjects.random == TRUE )
  {
       if( is.null(subjects.subset) )
       {
            ## sample from each treatment group roughly equally, but make total sampled add to 12
            trt.u		= unique(dat.igc$trt)
            num.trt		= length(trt.u)
            if( num.trt >= 12)
            {
                 subj.per.trt = 1
            }else{
                 subj.per.trt = floor(12/num.trt)
            }
            subjects.subset	= vector("list", num.trt)
            for( i in 1:num.trt )
            {
                 tmp			= subset(dat.igc, trt %in% trt.u[i])
                 subj.trt		= tmp$subject
                 subjects.subset[[i]]	= sample(subj.trt, subj.per.trt, replace = FALSE)
            }
            subjects.subset		= unlist(subjects.subset)
            subj.sampled		= subj.per.trt*num.trt
            if(subj.sampled < 12) add.subj = sample(setdiff(subjects.plot, subjects.subset),(12-subj.sampled), replace = FALSE)
            else add.subj = NULL
            subjects.subset		= c(subjects.subset, add.subj)
       }
       dat.igc 			<- subset(dat.igc, subject %in% subjects.subset)  ## predicted growth curve data
       dat.igcdata			<- subset(dat.igcdata, subject %in% subjects.subset) ## real data
  }
  
  p				= ggplot(data=dat.igc,aes(x=time, y = fit) )
  l				= geom_line(aes(linetype = model), size = 1.2)
  l.2				= geom_point(data=dat.igcdata,size=3,shape=1) ## Add the data
  axis	 			= labs(x = x.lab, y = y.lab, type = "Model")
  f				= facet_wrap(trt~subject, scales="fixed")
  if( !is.null(title.lab) )
  { 
       options		 		= labs(title = title.lab[2])
  }else{
       options				= NULL
  }
  dev.new()
  p.gcsub			= p + l + f + axis + l.2 + options + scale_linetype_discrete(name = "Model")
  print(p.gcsub)
  
  return(invisible(list(dat.gc = dat.gc, dat.data = dat.data, dat.igc = dat.igc, p.gctrt = p.gctrt, p.gcsub = p.gcsub)))
  
  time <- fit <- trt <- subject <- model <- NULL; rm(time); rm(fit); rm(subject); rm(tmp)
  
  gc()

} ## end of function growplot