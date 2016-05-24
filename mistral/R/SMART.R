## -----------------------------------------------------------------------------
## Fonction SMART
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

SMART = function(dimension,
			lsf,
			N1 = 10000,                     # Number of samples for the (L)ocalisation step
			N2 = 50000,                     # Number of samples for the (S)tabilisation step
			N3 = 200000,                    # Number of samples for the (C)onvergence step
			Nu = 50,                        # Size of the first Design of Experiments
			lambda1 = 7,                    # Relaxing parameter for MH algorithm at step L
			lambda2 = 3.5,                  # Relaxing parameter for MH algorithm at step S
			lambda3 = 1,                    # Relaxing parameter for MH algorithm at step C
			tune_cost = c(1,10,100,1000),   # Input for tuning cost paramter of the SVM
			tune_gamma = c(0.5,0.2,0.1,0.05,0.02,0.01), # Input for tuning gamma parameter of the SVM
			clusterInMargin = TRUE,         # Enforce selected clusterised points to be in margin
			alpha_margin = 1,
			#Localization, Convergence and Stabilization bounds
			k1 = round(6*(dimension/2)^(0.2)), # Rank of the first iteration of step S
			k2 = round(12*(dimension/2)^(0.2)), # Rank of the first iteration of step C
			k3 = k2 + 16,                  # Rank of the last iteration of step C
			#Arguments for SMART used in Subset simulation
      learn_db  = NULL,              # Coordinates of alredy known points
      lsf_value = NULL,              # Value of the LSF on these points
      failure   = 0,                 # Failure threshold
      limit_fun_MH = NULL,           # Define an area of exclusion with a limit function, eg in metaSS
      sampling_strategy = "MH",      # Either MH for Metropolis-Hastings of AR for accept-reject
      seeds = NULL,                  # If some points are already known to be in the appropriate subdomain, eg in metaSS
      seeds_eval = NULL,             # Value of the metamodel on these points
      burnin = 30,                   # Burnin parameter for MH
      thinning = 4,                  # Thinning parameter for MH
      plot = FALSE,                  # Set to TRUE for a full plot, ie refresh at each iteration
      limited_plot = FALSE,          # Set to TRUE for a final plot with final DOE, metamodel and LSF
      add = FALSE,                   # If plots are to be added to a current device
      output_dir = NULL,             # If plots are to be saved in jpeg in a given directory
      z_MH = NULL,                   # For plots, if metamodel has already been evaluated on the grid
      z_lsf = NULL,                  # For plots, if LSF has already been evaluated on the grid
      verbose = 0) {                 # Either 0 for almost no output, 1 for medium size output and 2 for all outputs
	
	cat("==========================================================================================\n")
	cat("                              Beginning of SMART algorithm\n")
	cat("==========================================================================================\n\n")
	
	# Fix NOTE for R CMD check
	x <- y <- z <- ..level.. <- NULL
	
	## STEP 0 : INITIALISATION
	
	Ncall = 0;
	z_meta = NA;
	
	#Define the proportion of margin, switching and closest points at the beginning and at the end of each stage
#current proportion at each iteration is given by linear regression between these two extremities
	proportion = list(L_stage=list(Nmargin=c(100,100),Nswitch=c(0,0),Nclose=c(0,0)),
			S_stage=list(Nmargin=c(80,40),Nswitch=c(20,50),Nclose=c(0,10)),
			C_stage=list(Nmargin=c(0,0),Nswitch=c(90,90),Nclose=c(10,10)))
	
	#Define the list variables used in the core loop
	#U stands for the points coordinates, one column per points, one row per dimension
	U = list(N1=matrix(nrow=dimension,ncol=N1),
		N2=matrix(nrow=dimension,ncol=N2),
		N3=matrix(nrow=dimension,ncol=N3),
		Nu=matrix(nrow=dimension,ncol=Nu),
		Nmargin=NA,
		Nswitch=NA,
		Nclose=NA)
  
	#G stands for the value of the limit state function on these points
	G = list(g=NA,#value on learn_db points, ie all the value already calculated at a given iteration
		Nmargin=NA,
		Nswitch=NA,
		Nclose=NA,
		Nu=NA*c(1:Nu))
	#G_meta stands for the value of the surrogate on these points
	G_meta = list(N1=NA*c(1:N1),
		N2=NA*c(1:N2),
		N3=NA*c(1:N3))
	
	cat(" ============================================= \n")
	cat(" STEP 1 : EVALUATION OF A FIRST METAMODEL \n")
	cat(" ============================================= \n\n")
  
	if(is.null(seeds) || sampling_strategy=="AR"){
	  if(verbose>0){cat("* Get Ru maximum radius of N3 standard Gaussian samples \n")}
	  tmp = rnorm(dimension*N3,mean=0,sd=1)
	  dim(tmp) = c(dimension,N3)
	  radius = apply(tmp,2,function(x) {sqrt(sum(x^2))})
	  Ru = max(radius)
	  if(verbose>1){cat("   - Ru maximum radius of",N3,"standard samples =",Ru,"\n")}
	}
  
	if(is.null(limit_fun_MH)) {
	  if(verbose>0){cat(" * Generate Nu =",Nu,"uniformly distributed samples in a sphere of radius Ru =",Ru,"\n")}
		U$Nu = t(runifSphere(dimension=dimension,N=Nu,radius=Ru))
	}
	else{
		if(sampling_strategy=="MH"){
		  if(verbose>0){cat(" * Generate N1 =",N1,"points from seeds with l1r-mM algorithm\n")}
			gen_pop = generateWithlrmM(seeds=seeds$N1,seeds_eval=seeds_eval$N1,N=N1,lambda=lambda1,limit_f=limit_fun_MH,burnin=0,thinning=0)
			U$N1 = gen_pop$points
			G_meta$N1 = gen_pop$eval
		  if(verbose>0){cat(" * Get Nu =",Nu,"points by clustering of the N1 points\n")}
			U$Nu = t(kmeans(t(U$N1), centers=Nu,iter.max=30)$centers)
		}
		else{
		  if(verbose>0){cat(" * Generate Nu =",Nu,"uniformly distributed samples with an accept-reject strategy\n")}
			rand = function(dimension,N) {t(runifSphere(dimension,N,radius=Ru))}
			U$Nu = generateWithAR(dimension=dimension,N=Nu,limit_f=limit_fun_MH,rand=rand)
		}
	}
	rownames(U$Nu) <- rep(c('x', 'y'), length.out = dimension)

	if(verbose>0){cat(" * Assessment of g on these points\n")}
	G$Nu = lsf(U$Nu);Ncall = Ncall + Nu

	#plotting part
	if(plot==TRUE){
	  
	  xplot <- yplot <- c(-80:80)/10
	  df_plot = data.frame(expand.grid(x=xplot, y=yplot), z = lsf(t(expand.grid(x=xplot, y=yplot))))
	  
	  if(add==TRUE){
	    p <- last_plot()
	  } else {
	    p <- ggplot(data = df_plot, aes(x,y)) +
	      geom_contour(aes(z=z, color=..level..), breaks = failure) +
	      theme(legend.position = "none") +
	      xlim(-8, 8) + ylim(-8, 8)
	  }
	  if(verbose>0){cat("2D PLOT : FIRST DoE \n")}
		if(is.null(limit_fun_MH)) {
		  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
		    r = diameter / 2
		    tt <- seq(0,2*pi,length.out = npoints)
		    xx <- center[1] + r * cos(tt)
		    yy <- center[2] + r * sin(tt)
		    return(data.frame(x = xx, y = yy))
		  }
		  
		  dat <- circleFun(diameter = Ru*2,npoints = 100)
		  #geom_path will do open circles, geom_polygon will do filled circles
		  p <- p + geom_path(data = dat, aes(alpha = 0.3))
		  print(p)
		}
		else{
		  z_MH = limit_fun_MH(t(df_plot[,1:2]))
		  df_plot_MH <- data.frame(df_plot[,1:2], z = z_MH$mean)
		  p <- p + geom_contour(data = df_plot_MH, aes(z=z, color=..level..), breaks = 0)
		}
	  p <- p + geom_point(data = data.frame(t(U$Nu), z = lsf(U$Nu)), aes(color = z))
	  print(p)
	}

	#add points to the learning database
	if(verbose>0){cat(" * Add points to the learning database\n")}
	if(is.null(learn_db)){
		learn_db = cbind(seq(0,0,l=dimension),U$Nu)
		g0 = lsf(as.matrix(seq(0,0,l=dimension)));Ncall = Ncall + 1;#this is to insure consistency between model sign and lsf sign
		G$g = c(g0,G$Nu)
	}
	else{
		learn_db = cbind(learn_db,U$Nu)
		G$g = c(lsf_value,G$Nu)
	}

	#Train the model
	if(verbose>0){cat(" * Train the model\n")}
	if(verbose<1){capture.output(meta <- trainModel(design=learn_db,response=(G$g-failure),type="SVM",cost=tune_cost,gamma=tune_gamma))}
	else{meta = trainModel(design=learn_db,response=(G$g-failure),type="SVM",cost=tune_cost,gamma=tune_gamma)}

	#Update meta_fun & meta_model
	if(verbose>0){cat(" * UPDATE quantities based on surrogate model \n")}
	meta_model = meta$model
	meta_fun = meta$fun

	#plotting part
	if(plot==TRUE){
	  if(verbose>0){cat("2D PLOT : FIRST METAMODEL\n")}
		z_meta = meta_fun(t(df_plot[,1:2]))
		df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean)
		print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = c(0, -1, 1)))
	}

	cat(" ============================= \n")
	cat(" STEP 2 : REFINEMENT PROCEDURE \n")
	cat(" ============================= \n\n")
  
	for (k in c(1:k3)) {
	  cat("\n * ITERATION ",k,"in [",1,k3,"]\n")
		stage = 1*(k<k1) + 2*(k>=k1)*(k<k2) + 3*(k>=k2)
    stage_str = switch(stage, "1" = {"A- Localisation"}, "2" = {"B- Stabilisation"}, "3" = {"C- Convergence"})

		cat("\n",stage_str,"STAGE \n")
		cat("    ======================================================== \n\n")
    
		#Get the parameters corresponding to the current stage
	  if(verbose>0){cat(" * Get the parameters corresponding to the current stage\n")}
		k_start = 1*(stage==1) + k1*(stage==2) + k2*(stage==3)
		k_end = (k1-1)*(stage==1) + (k2-1)*(stage==2) + k3*(stage==3)
		N = N1*(stage==1) + N2*(stage==2) + N3*(stage==3)
	  if(verbose>1){cat("   - Number of points to sample =",N,"\n")}
		lambda = lambda1*(stage==1) + lambda2*(stage==2) + 1*(stage==3)

		#Get the parameters corresponding to the current iteration
	  if(verbose>0){cat(" * Get the parameters corresponding to the current iteration\n")}
		Nsup = round((3+(k-1)/(k1-1))*sqrt(dimension))*(stage==1) +
			round((4+(k-k1)/(k2-k1))*sqrt(dimension))*(stage==2) +
			round(5*sqrt(dimension))*(stage==3)
		aMargin = proportion[[stage]]$Nmargin[1]; bMargin = proportion[[stage]]$Nmargin[2]
		aSwitch = proportion[[stage]]$Nswitch[1]; bSwitch = proportion[[stage]]$Nswitch[2]
		aClose = proportion[[stage]]$Nclose[1]; bClose = proportion[[stage]]$Nclose[2]
		Nmargin = round(Nsup/100*(aMargin+(k-k_start)/(k_end-k_start)*(bMargin-aMargin)))
		Nswitch = round(Nsup/100*(aSwitch+(k-k_start)/(k_end-k_start)*(bSwitch-aSwitch)))
		Nclose = round(Nsup/100*(aClose+(k-k_start)/(k_end-k_start)*(bClose-aClose)))
	  if(verbose>1){cat("   - Nsup =",Nsup,"\n")
		cat("   - Nmargin =",Nmargin,"\n")
		cat("   - Nswitch =",Nswitch,"\n")
		cat("   - Nclose =",Nclose,"\n\n")}

		if(is.null(limit_fun_MH)){
			if(stage==3) {
				#Generate standard gaussian samples
			  if(verbose>0){cat(" * Generate standard gaussian samples\n")}
				U[[stage]] = matrix(rnorm(dimension*N, mean=0, sd=1),dimension,N)
			}
			else {
				#Generate samples with a uniform distribution on a sphere of radius Ru
			  if(verbose>0){cat(" * Generate samples with a uniform distribution on a sphere of radius Ru =",Ru,"\n")}
				U[[stage]] = t(runifSphere(dimension=dimension,N=N,radius=Ru))
			}
		}
		else{
			if(sampling_strategy=="MH"){
			  if(verbose>0){cat(" * Generate N =",N,"points with lr-mM algorithm\n")}
				U[[stage]] = generateWithlrmM(seeds=seeds[[stage]],seeds_eval=seeds_eval[[stage]],N=N,lambda=lambda,limit_f=limit_fun_MH,burnin=0,thinning=0)$points
			}
			else{
			  if(verbose>0){cat(" * Generate N =",N,"points with an accept-reject strategy\n")}
				if(stage==3) {
					rand = function(dimension,N) {matrix(rnorm(dimension*N),dimension,N)}
					U[[stage]] = generateWithAR(dimension=dimension,N=N,limit_f=limit_fun_MH,rand=rand)
				}
				else {
					rand = function(dimension,N) {t(runifSphere(dimension=dimension,N=N,radius=Ru))}
					U[[stage]] = generateWithAR(dimension=dimension,N=N,limit_f=limit_fun_MH,rand=rand)
				}
			}
		}

		#Evaluate g on U[[stage]] using the meta model
	  if(verbose>0){cat(" * Evaluate the metamodel on these points \n")}
		meta_pred = meta_fun(U[[stage]])
		G_meta[[stage]] = meta_pred$mean

		#Selection of Nsup = Nmargin + Nswitch + Nclose new points where g is to be assessed
		if(Nmargin>0) {
			#Selection of Nmargin new points where g is to be assessed
		  if(verbose>0){cat(" * Selection of Nmargin =",Nmargin,"new points where g is to be assessed\n")}
			isMargin = inMargin(meta_pred,type="SVM",alpha=alpha_margin)
		  if(verbose>1){cat("   - Number of margin points =",sum(isMargin*1),"\n")}
			if(!sum(isMargin)==0){
				U$Nmargin = tryCatch(
					clusterize(data=U[[stage]][,isMargin],Ncluster=Nmargin,inMargin=clusterInMargin),
					error = function(cond) {
						message(cond)
						message("\nall margin points are kept")
						return(U[[stage]][,isMargin])
					})
				U$Nmargin = as.matrix(U$Nmargin)
				rownames(U$Nmargin) <- rep(c('x', 'y'), length.out = dimension)
        
				#Assessment of g
				if(verbose>0){cat(" * Assessment of g on these points \n")}
				G$Nmargin <- lsf(U$Nmargin);Ncall = Ncall + Nmargin
				#Add points U$Nmargin to the learning database
				if(verbose>0){cat(" * Add points to the learning database\n")}
				learn_db <- cbind(learn_db,U$Nmargin)
				G$g = c(G$g,G$Nmargin)
	
				#plotting part
				if(plot==TRUE){
				  if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
				  print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nmargin), z = G$Nmargin), aes(color = z)))
				  p <- p + geom_point(data = data.frame(t(U$Nmargin), z = G$Nmargin), aes(color = z))
				}
				
				#Train the model
				if(verbose>0){cat(" * Train the model\n")}
				meta = trainModel(meta_model,design=learn_db,response=(G$g-failure),updesign=U$Nmargin,upresponse=(G$Nmargin-failure),type="SVM")
				#Update meta_fun & meta_model
				meta_model.prev = meta_model
				meta_model = meta$model
				meta_fun.prev = meta_fun
				meta_fun = meta$fun
			}
		}

		if(Nswitch>0) {
			#Selection of Nswitch new points where g is to be assessed
		  if(verbose>0){cat(" * Selection of Nswitch =",Nswitch,"new points where g is to be assessed\n")}
			G_meta_prev = meta_fun.prev(U[[stage]])$mean
			isSwitched = as.logical(sign(G_meta_prev*G_meta[[stage]])<0)
		  if(verbose>1){cat("   - Number of switching points =",sum(isSwitched*1),"\n")}
			if(!sum(isSwitched)==0){
			  U$Nswitch <- tryCatch(t(kmeans(t(U[[stage]][,isSwitched]), centers=Nswitch, iter.max=30)$centers),
			           error = function(cond) {
			             message(cond);
			             message("\nall switching points are kept");
			             return(U[[stage]][,isSwitched]);
			           })
				U$Nswitch = as.matrix(U$Nswitch)
				rownames(U$Nswitch) <- rep(c('x', 'y'), length.out = dimension)
				#Assessment of g
				if(verbose>0){cat(" * Assessment of g\n")}
				G$Nswitch = lsf(U$Nswitch);Ncall = Ncall + Nswitch
				#Add points U$Nmargin + U$Nswitch to the learning database
				if(verbose>0){cat(" * Add points U$Nswitch  to the learning database\n")}
				learn_db = cbind(learn_db,U$Nswitch)
				G$g = c(G$g,G$Nswitch)
	
				if(plot==TRUE){
				  if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
				  print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nswitch), z = G$Nswitch), aes(color = z)))
				  p <- p + geom_point(data = data.frame(t(U$Nswitch), z = G$Nswitch), aes(color = z))
				}
				
				#Train the model
				if(verbose>0){cat(" * Train the model\n")}
				meta = trainModel(meta_model,design=learn_db,response=(G$g-failure),updesign=U$Nswitch,upresponse=(G$Nswitch-failure),type="SVM")
				#Update meta_fun & meta_model
				meta_model.prev = meta_model
				meta_model = meta$model
				meta_fun.prev = meta_fun
				meta_fun = meta$fun
			}
		}
		
		if(Nclose>0) {
			#Selection of Nclose new points where g is to be assessed
		  if(verbose>0){cat(" * Selection of Nclose =",Nclose,"new points where g is to be assessed\n")}
			isClose = sort(abs(G_meta[[stage]]),index=TRUE)$ix[1:Nclose]
			U$Nclose = as.matrix(U[[stage]][,isClose])
      rownames(U$Nclose) <- rep(c('x', 'y'), length.out = dimension)
      
			#Assessment of g
		  if(verbose>0){cat(" * Assessment of g\n")}
			G$Nclose = lsf(U$Nclose);Ncall = Ncall + Nclose
			#Add points U$Nclose to the learning database
		  if(verbose>0){cat(" * Add points U$Nclose to the learning database\n")}
			learn_db = cbind(learn_db,U$Nclose)
			G$g = c(G$g,G$Nclose)

			if(plot==TRUE){
			  if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
			  print(p_meta <- p_meta + geom_point(data = data.frame(t(U$Nclose), z = G$Nclose), aes(color = z)))
			  p <- p + geom_point(data = data.frame(t(U$Nclose), z = G$Nclose), aes(color = z))
			}
			
			#Train the model
		  if(verbose>0){cat(" * Train the model\n")}
			meta = trainModel(meta_model,design=learn_db,response=(G$g-failure),updesign=U$Nclose,upresponse=(G$Nclose-failure),type="SVM")
			#Update meta_fun & meta_model
			meta_model.prev = meta_model
			meta_model = meta$model
			meta_fun.prev = meta_fun
			meta_fun = meta$fun
		}
    
	  #plotting part
	  if(plot==TRUE){
	    if(verbose>0){cat(" * 2D PLOT : UPDATE\n")}
	    z_meta = meta_fun(t(df_plot[,1:2]))
	    df_plot_meta <- data.frame(df_plot[,1:2], z = z_meta$mean)
	    print(p_meta <- p + geom_contour(data = df_plot_meta, aes(z=z, color=..level.., alpha = 0.5), breaks = c(0, -1, 1)))
	  }
	}

	if( (plot || limited_plot) & add==FALSE & !is.null(output_dir)) {dev.off()}

	cat("\n ======================================= \n")
	cat(" STEP 3 : FAILURE PROBABILITY ESTIMATION \n")
	cat(" ======================================= \n\n")
  
	#estimation of probability P using metamodel classifier values
	if(verbose>0){cat(" * estimation of probability P using metamodel classifier values\n")}
	if(!is.null(limit_fun_MH)) {
		capture.output(gen_pop <- generateWithlrmM(seeds=seeds$N3,seeds_eval=seeds_eval$N3,
                                               N=N,lambda=lambda,limit_f=limit_fun_MH,
                                               burnin=burnin,thinning=thinning,
                                               VA_function=function(x){1*(meta_fun(x)$mean<0)})
                   )
		fail_points = gen_pop$VA_values
		points=gen_pop$points[,fail_points]
		meta_eval= meta_fun(points)$mean
		MC_est = gen_pop$est
		MC_var = gen_pop$var
		MC_delta = gen_pop$delta
		MC_gamma = gen_pop$gamma
	}
	else{
		G_meta$N3 = meta_fun(U$N3)$mean
		fail_points = (G_meta$N3<0)
		points=U$N3[,fail_points]
		meta_eval=G_meta$N3[fail_points]
		MC_est = sum(1*fail_points)/N3
		MC_var = MC_est*(1-MC_est)/N3
		MC_delta = sqrt(MC_var)/MC_est
		MC_gamma = 0
	}

	G_meta$N1 = meta_fun(U$N1)$mean
	indN1 = (G_meta$N1<0)
	G_meta$N2 = meta_fun(U$N2)$mean
	indN2 = (G_meta$N2<0)
	points = list(N1=U$N1[,indN1],N2=U$N2[,indN2],N3=points)
	meta_eval = list(N1=G_meta$N1[indN1],N2=G_meta$N2[indN2],N3=meta_eval)

	cat("==========================================================================================\n")
	cat("                              End of SMART algorithm\n")
	cat("==========================================================================================\n\n")
  
	cat(" * proba =",MC_est,"\n")
  cat(" * variance =",MC_var,"\n")
  cat(" * cov =",MC_delta,"\n")
	cat(" * Markov chain gamma =",MC_gamma,"\n")

	if(plot+limited_plot) {
		res = list(proba=MC_est,
               cov=MC_delta,
               gamma=MC_gamma,
               Ncall=Ncall,
               learn_db=learn_db,
               lsf_value=G$g,
               meta_fun=meta_fun,
               meta_model=meta_model,
               points=points,
               meta_eval=meta_eval,
               z_meta=z_meta$mean)
	}
	else {res = list(proba=MC_est,
	                 cov=MC_delta,
	                 gamma=MC_gamma,
	                 Ncall=Ncall,
	                 learn_db=learn_db,
	                 lsf_value=G$g,
	                 meta_fun=meta_fun,
	                 meta_model=meta_model,
	                 points=points,
	                 meta_eval=meta_eval)
	}
	
	return(res)
}
