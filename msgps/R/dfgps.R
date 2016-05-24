dfgps <-
function(X,y,penalty="enet", ex_para=c(0),  STEP=10000, STEP.max=100000, DFtype="MODIFIED", p.max=300){
		#####check if the algorithm works
		candidate <- c("enet","genet","alasso")
	 	if(sum(candidate == penalty) != 1)	stop('penalty must equal "enet", "genet" or "alasso".')

	   if(mode(X)!="numeric") stop(" X must be numeric.")
	   if(!is.matrix(X)) stop(" X must be a matrix.")
	   if(mode(y)!="numeric") stop(" y must be numeric.")
		if (!is.vector(y)) stop("y must be a vector.")
		if (nrow(X)!=length(y)) stop("The number of sample must not be differenet between X and y.")
		if (sum(complete.cases(X)==FALSE)>0)  stop("X must be complete data.")
		if (sum(complete.cases(y)==FALSE)>0)  stop("y must be complete data.")

	    penalty_int <- which(candidate == penalty)
		penalty_int <- as.integer(penalty_int - 1)
		if((penalty_int==0 ||  penalty_int==1)  &&   length(ex_para)!=1)	stop('"ex_para" must be a scalar.')
		if(!is.numeric(ex_para))	stop('"ex_para" must be a numeric.')
		if(penalty_int==0  &&   (ex_para < 0 || ex_para >=1))	stop('"ex_para" must be in [0,1).')
		if(penalty_int==1  &&   (ex_para <= 0 || ex_para >=1))	stop('"ex_para" must be in (0,1).')
		if(penalty_int==2  &&   length(ex_para)!=2) stop("ex_para must be a 2-dimensional vector.")
		if(penalty_int==2  &&   length(ex_para)!=2) stop("ex_para must be a 2-dimensional vector.")
		if(penalty_int==2  &&   ex_para[1]<0 ) stop("ex_para[1] must be non-negative.")
		if(penalty_int==2  &&   ex_para[2]<0 ) stop("ex_para[2] must be non-negative.")

		if(mode(STEP)!="numeric")	stop('"STEP" must be numeric.')
		if(STEP < 500)	stop('"STEP" must be greater than or equal to 500.')
		if(STEP >= 1e+8)	stop('"STEP" must be less than 1e+8.')

		candidate_DFtype <- c("NAIVE","MODIFIED")
	 	if(sum(candidate_DFtype == DFtype) != 1)	stop('DFtype must be "MODIFIED"  or "NAIVE".')

		if(mode(STEP)!="numeric")	stop('"STEP" must be numeric.')
		if(p.max < 1)	stop('"p.max" must bea  positive integer.')
		if(p.max >= 10000)	stop('"p.max" must be less than 10000.')

		if(mode(STEP.max)!="numeric")	stop('"STEP.max" must be numeric.')
		if(STEP.max< 500)	stop('"STEP.max" must be greater than or equal to 500.')
		if(STEP.max >= 1e+8)	stop('"STEP.max" must be less than 1e+8.')

		#weight for adaptive lasso
		if(penalty_int==2){ 
			if(ncol(X) < nrow(X)) PLS <- solve(t(X)%*%X + ex_para[2]*diag(ncol(X)))%*%t(X)%*%y
			if(ncol(X) >= nrow(X)) PLS <- t(X)%*%solve((X)%*%t(X) + ex_para[2]*diag(nrow(X)))%*%y
			weight_vec <- (abs(PLS))^(-ex_para[1])
		}else{
			weight_vec <- rep(1,ncol(X))
		}

		#standardize
		#mean, varを計算
		
		meanX <- apply(X,2,mean)
		meanX_mat <- sweep(X, 2, meanX)
		standardize_vec <- 1 / sqrt(apply(meanX_mat^2,2,sum))
		standardize_vec[standardize_vec==Inf] <- 0
		standardize_mat <- matrix(rep(standardize_vec,nrow(X)),nrow(X),ncol(X),byrow=T)
		meany <- mean(y)

		X0 <- X
		X <- meanX_mat * standardize_mat
#		X <- scale(X) / sqrt(nrow(X)-1)
		y0 <- y
		y <- y-meany

		

		#####calculating delta_t
		if( nrow(X) <= ncol(X)){
			beta_OLS <- t(X)%*%solve(X%*%t(X) + 0.001*diag(nrow(X)))%*%y
		}else if(det(t(X)%*%X) < 1e-3){
			beta_OLS <- solve(t(X)%*%X + 0.001*diag(ncol(X))) %*%t(X)%*%y
		}else{
			beta_OLS <-  solve(t(X)%*%X)%*%t(X)%*%y 
			}
		 delta_t <- sum(abs(beta_OLS)) / STEP
		#if(penalty=="enet") delta_t <- sum(ex_para[1]*beta_OLS^2/2 + (1-ex_para[1]) * abs(beta_OLS)) / STEP
		#if(penalty=="genet") delta_t <- (sum(log(ex_para[1]  + (1-ex_para[1]) * abs(beta_OLS))) - ncol(X)*log(ex_para[1])) / STEP
		#if(penalty=="alasso") delta_t <- sum(weight_vec * abs(beta_OLS)) / STEP

		#main
		STEP.max=as.integer(STEP.max)
         p.max=as.integer(p.max)
		 N <- length(y)
        #library.dynam("msgps")
        gps_C=.Call("gps", X,y,delta_t,penalty_int,ex_para, STEP.max,p.max,weight_vec,standardize_vec)
        #library.dynam.unload("msgps")
		betagps_matrix <- gps_C[1][[1]]
        STEP_adj  <- gps_C[2][[1]]
        RSS  <- gps_C[3][[1]]
        increment_covpenalty_vec_power <- gps_C[4][[1]]
        selected_variable_index_vec <- gps_C[5][[1]]
        betahat_index_vec_adj <- gps_C[6][[1]]
        tuning <- gps_C[7][[1]]
        tuning_stand <- gps_C[8][[1]]
        selected_variable_index_vec <- selected_variable_index_vec[selected_variable_index_vec>0]
        betagps_matrix <- betagps_matrix[1:STEP_adj]
        RSS <- RSS[1:STEP_adj]
        increment_covpenalty_vec_power <- increment_covpenalty_vec_power[1:STEP_adj]
        betahat_index_vec_adj <- betahat_index_vec_adj[1:STEP_adj]
        tuning <- tuning[1:STEP_adj]
		tuning_stand <- tuning_stand[1:STEP_adj]
        ex_X_selected = X[,selected_variable_index_vec]
        betagps_matrix=as.integer(betagps_matrix)
        STEP_adj2=as.integer(STEP_adj-1)
        STEP_adj=as.integer(STEP_adj)
        betahat_index_vec_adj = as.integer(betahat_index_vec_adj)
        selected_variable_index_vec = as.integer(selected_variable_index_vec)
		increment_covpenalty_vec0 <- (1-2*delta_t/N)^ increment_covpenalty_vec_power
		increment_covpenalty_vec1 <- 1-(1-2*delta_t/N)^ increment_covpenalty_vec_power
		t0 <- log( 1-2*delta_t/N*increment_covpenalty_vec_power ) / log(1-2*delta_t/N)
		increment_covpenalty_vec2 <- 1-(1-2*delta_t/N)^ t0
		increment_covpenalty_vec3 <- 2*delta_t/N * increment_covpenalty_vec_power
        
if(DFtype=="NAIVE"){
         #DFNAIVE
        dfgps_vec=.Call("DFNAIVE2", ex_X_selected,y,betahat_index_vec_adj,STEP_adj2,increment_covpenalty_vec3)
		if(sum( abs(dfgps_vec)>1e+10 | is.na(dfgps_vec) | dfgps_vec<0 ) > 1 ) stop("DF is not correct")
         
         
}else if(DFtype=="MODIFIED"){
        #DFMODIFIED
        #QR decomposition
        qr_X <- qr(ex_X_selected)
        qr_X_R <- qr.R(qr_X)
		qr_X_R[is.nan(qr_X_R)] <- 0
        #dfgps_vec=.Call("DFMODIFIED", qr_X_R, y, betahat_index_vec_adj, STEP_adj2, increment_covpenalty_vec, selected_variable_index_vec)
        dfgps_vec=.Call("DFMODIFIED2", qr_X_R, y, betahat_index_vec_adj, STEP_adj2, increment_covpenalty_vec3, selected_variable_index_vec)
		if(sum( abs(dfgps_vec)>1e+10 | is.na(dfgps_vec) | dfgps_vec<0 ) > 1 ) stop("DF is not correct")
		}

        p<-ncol(X)
        N<-nrow(X)
        if(DFtype=="NAIVE" || DFtype=="MODIFIED" || DFtype=="MODIFIED2"){
                ans <- list(N=N,p=p,delta_t=delta_t,coefficient_index=betagps_matrix,df=dfgps_vec,STEP_adj=STEP_adj2,RSS=RSS,tuning=tuning,tuning_stand=tuning_stand,X=X0,y=y0,ex_para=ex_para,beta_OLS=beta_OLS,Xstand=X,ystand=y)
                }else{
                ans <- list(N=N,p=p,delta_t=delta_t,coefficient_index=betagps_matrix,STEP_adj=STEP_adj2,RSS=RSS,tuning=tuning,tuning_stand=tuning_stand,X=X0,y=y0,ex_para=ex_para,beta_OLS=beta_OLS,Xstand=X,ystand=y)
                }
		class(ans) <- "dfgps"
        return(invisible(ans))
}

