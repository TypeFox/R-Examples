

tfrPredTab <- function(tfr.pred, main.win, wpp.year) {
	eTFRp <- new.env()

	eTFRp$wpp.year <- wpp.year
	
	create.sim.dir.widget(env=eTFRp, parent=tfr.pred, type=c('tfr', 'tfr3'), 
				main.win=main.win,
				default=eval(formals(run.tfr.mcmc)$output.dir))
	
	nb <- bDem.gnotebook(container=tfr.pred, expand=TRUE)

	# Run MCMC group
	# color chosen from http://html-color-codes.info
	mcmc.g <- ggroup(label="  <span color='#B40404'>Run MCMC</span>  ", markup=TRUE, 
							horizontal=FALSE, container=nb, spacing=10)
	mcmc.g.env <- TFRrunMCMCgroup(mcmc.g, main.win, parent=eTFRp)	
	# Continue MCMC group
	cont.mcmc.g <- ggroup(label="  <span color='#B40404'>Continue MCMC</span>  ", markup=TRUE, 
							horizontal=FALSE, container=nb, spacing=10)
	cont.mcmc.g.env <- TFRcontinueMCMCgroup(cont.mcmc.g, main.win, parent=eTFRp)

	# Predictions group
	pred.g <- ggroup(label="  <span color='#B40404'>Make predictions</span>  ", markup=TRUE,
							horizontal=FALSE, container=nb, spacing=10)
	pred.g.env <- TFRnewPred.group(pred.g, main.win, parent=eTFRp)
	
	# Result group
	result.g <- ggroup(label="  <span color='#B40404'>Explore results</span>  ", markup=TRUE,
							horizontal=FALSE, container=nb)
	result.g.env <- TFRresults.group(result.g, main.win, parent=eTFRp)
	
	svalue(nb) <- 1
	label <- glabel(paste('Dependency in use: bayesTFR  v.', packageVersion("bayesTFR")), container=tfr.pred)
	font(label) <- c(style='italic', family='serif')
}






