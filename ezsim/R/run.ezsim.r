#' Run the Simulation of an ezsim object. The simulation result is store into the ezsim object in the argument directly, reassignment is not needed.
#' @name run.ezsim
#' @aliases run.ezsim
#' @title Run the Simulation
#' @method run ezsim
#' @param x An ezsim object
#' @param \dots not used
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @export 
#' @examples              
#' \dontrun{
#' ezsim_basic<-ezsim(
#'     m             = 100,
#'     run           = FALSE,
#'     run_test      = TRUE,
#'     display_name  = c(mean_hat="hat(mu)",sd_mean_hat="hat(sigma[hat(mu)])"),
#'     parameter_def = createParDef(list(n=seq(20,80,20),mu=c(0,2),sigma=c(1,3,5))),
#'     dgp           = function() rnorm(n,mu,sigma),
#'     estimator     = function(x) c(mean_hat = mean(x), 
#'                                  sd_mean_hat=sd(x)/sqrt(length(x)-1)),
#'     true_value    = function() c(mu, sigma / sqrt(n-1))
#' )
#' run(ezsim_basic)
#' }

run.ezsim <-function(x,...){
	x$parameter_list <- generate(x$parameter_def)
	## TODO tryCatch
	## TODO flexible way to store time used
	
	create_cluster_flag <- FALSE
	
	## Decide whether we need to stop the cluster after the run
	create_cluster_flag <- FALSE
	
	if (x$use_core > 1 & is.null(x$cluster) ){
		x$cluster <- makeCluster(x$use_core)
		if (!is.null(x$cluster_packages)){
			for (i in 1:length(x$cluster_packages))
				eval(substitute( 
					clusterEvalQ(x$cluster , require(w, character.only=TRUE)   )  ,
					list( w=x$cluster_packages[i] )
				))
		}
		clusterSetRNGStream(x$cluster,x$use_seed)
		create_cluster_flag<-TRUE
	}

	if (x$use_core == 1) 
		set.seed(x$use_seed)

	# if (x$parallel & is.null(x$cluster)){ 
	# 	x$cluster<-makeCluster(x$number_of_workers)
	# 	create_cluster_flag<-TRUE
	# }
	
	##  A local function to conduct simulation
	## fix parameter, repeat for m times
	i=j=NULL
	compute_simulation<-function(par,ezsim_object,m){
		parser<-function(i,par,ezsim_object){
			ezsim_object$estimator(Jmisc::evalFunctionOnList(ezsim_object$dgp,par))
		}
		
		out<-
		if (x$use_core == 1){
			lapply(1:m,parser,par=par,ezsim_object=ezsim_object)
		}
		else{
			parLapply(ezsim_object$cluster,1:m,parser,par=par,ezsim_object=ezsim_object)
		}
		return(out)
	}
	
	tryCatch({
		# if no autosave
		if (x$auto_save==0){
			time_used<-system.time({
				x$simulation_result <- lapply(x$parameter_list, function(i) compute_simulation(i,x,x$m)) 
			})
		} else{ 
			# compute each autosave interval
			m<-rep(trunc(x$m/x$auto_save),x$auto_save-1)
			m<-c(m,x$m-sum(m))
			
			time_stamp <- toString(Sys.time())
			time_stamp <- gsub(' ','_',time_stamp)
			time_stamp <- gsub('-|:','',time_stamp)
			# run the simulation one by one

			ezsim_list<-
			foreach (i = 1:length(m), j=m) %do% {
				obj_name<-paste('ezsim_part_',i,sep='')
				file_name<-paste('part_',i,'_',time_stamp,'.rData',sep='')
				
				# clone the ezsim object and change the number of simulation
				temp_ezsim<-x
				temp_ezsim$m<-j
				temp_ezsim$auto_save<-0
				
				# run the simulation
				cat('Part ',i,'. Number of replication:',j,'\n',sep='')
				temp_ezsim <- run(temp_ezsim)

				# rename the name of the ezsim object and save it to a file.
				assign(obj_name,temp_ezsim)
				cat('Saving part ',i,' to ',file_name,'\n\n',sep='')
				save(file=file_name,list=obj_name)
				temp_ezsim
			}
			## merge them back and return
			temp_x<-foreach ( i=ezsim_list , .combine=merge.ezsim) %do% i
			temp_x$auto_save <- x$auto_save
			x<-temp_x
		}
		## create simulation table
		tryCatch({
			x<-createSimulationTable(x)
		}, error = function(e){
			cat("Error in creating simulation table\n")
			stop(e)
		})
	}, finally = {
		if (create_cluster_flag){
			tryCatch({
				stopCluster(x$cluster)
			}, finally = {
				x$cluster<-NULL
			})
		}
	})
	return(x)
}
