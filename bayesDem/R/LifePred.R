e0PredTab <- function(e0w, main.win, wpp.year) {
	eLifeExp <- new.env()
	eLifeExp$wpp.year <- wpp.year
	
	create.sim.dir.widget(env=eLifeExp, parent=e0w, type='e0', 
				main.win=main.win,
				default=eval(formals(run.e0.mcmc)$output.dir))
				
	nb <- bDem.gnotebook(container=e0w, expand=TRUE)
	# Run MCMC group
	# color chosen from http://html-color-codes.info
	mcmc.g <- ggroup(label="  <span color='#B40404'>Run MCMC</span>  ", markup=TRUE, 
							horizontal=FALSE, container=nb, spacing=10)
	mcmc.g.env <- e0RunMCMCgroup(mcmc.g, main.win, parent=eLifeExp)	
	# Continue MCMC group
	cont.mcmc.g <- ggroup(label="  <span color='#B40404'>Continue MCMC</span>  ", markup=TRUE, 
							horizontal=FALSE, container=nb, spacing=10)
	cont.mcmc.g.env <- e0ContinueMCMCgroup(cont.mcmc.g, main.win, parent=eLifeExp)

	# Predictions group
	pred.g <- ggroup(label="  <span color='#B40404'>Make predictions</span>  ", markup=TRUE,
							horizontal=FALSE, container=nb, spacing=10)
	pred.g.env <- e0NewPred.group(pred.g, main.win, parent=eLifeExp)
	
	# Result group
	result.g <- ggroup(label="  <span color='#B40404'>Explore results</span>  ", markup=TRUE,
							horizontal=FALSE, container=nb)
	result.g.env <- e0Results.group(result.g, main.win, parent=eLifeExp)
	
	svalue(nb) <- 1
	label <- glabel(paste('Dependency in use: bayesLife  v.', packageVersion("bayesLife")), container=e0w)
	font(label) <- c(style='italic', family='serif')
	
}