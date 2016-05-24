#' guess_stnderr
#'
#' @title Bootstrapped standard errors of effect size estimates
#' 
#' @param pre_test data.frame carrying pre_test items
#' @param pst_test data.frame carrying pst_test items
#' @param nsamps  number of resamples, default is 100
#' @param seed    random seed, default is 31415
#' @param force9       Optional. There are cases where DK data doesn't have DK. But we need the entire matrix. By default it is FALSE.
#' @return  list with standard error of parameters, estimates of learning, standard error of learning by item
#' @export
#' @examples
#' pre_test <- data.frame(pre_item1=c(1,0,0,1,0), pre_item2=c(1,NA,0,1,0)) 
#' pst_test <- data.frame(pst_item1=pre_test[,1] + c(0,1,1,0,0), 
#'						  pst_item2 = pre_test[,2] + c(0,1,0,0,1))
#' \dontrun{guess_stnderr(pre_test, pst_test, nsamps=10, seed = 31415)}

guess_stnderr <- function(pre_test=NULL, pst_test=NULL, nsamps=100, seed = 31415, force9=FALSE) 
{
	
	# pre_test <- alldat[,t1]; pst_test <-  alldat[,t2]; nsamps=10; seed = 31415
	# pre_test <- alldat_dk[,t1]; pst_test <-  alldat_dk[,t2]; nsamps=10; seed = 31415

	# build a df
	df 		<- data.frame(cbind(pre_test, pst_test))

	# set nitems and nparams based on df		
	nitems 	    <- ncol(df)/2	
	transmatrix <- multi_transmat(pre_test, pst_test, force9=force9)
	nparams     <- ifelse(ncol(transmatrix)==4, 4, 8)	
			
	#define matrices to store samples and sample results		
	resamps.t1 		 <- resamps.t3 <- resamps.results <- list()
	resamps.lca.eff  <- matrix(ncol=nitems+1, nrow=nsamps)
	resamps.agg      <- matrix(ncol=ncol(df), nrow=nsamps)
			
	stnderrs.lca.params <- matrix(ncol=nitems, nrow=nparams)
			
	stnderrs.effects <- avg.effects <- matrix(ncol=nitems+1, nrow=1)
				
	resamps.lca.params <- rep(list(matrix(nrow = nsamps, ncol = nparams)), nitems)
				
	#extracting samples from the data	
	set.seed(seed)
	resamples <- lapply(1:nsamps, function(i) df[sample(1:nrow(df), replace = T),])

	# Looping through the samples; estimating based one each
	for(i in 1:length(resamples)) {
		print(i)
		transmatrix_i           <- multi_transmat(resamples[[i]][,1:nitems], resamples[[i]][,(nitems+1):(2*nitems)], force9=force9)
		resamps.results[[i]] 	<- guesstimate(transmatrix_i)
		resamps.lca.eff[i,] 	<- resamps.results[[i]]$est.learning
		resamps.agg[i,] 		<- transmatrix_i[nitems,]

	    for(j in 1:nitems) {
			resamps.lca.params[[j]][i,]    <- resamps.results[[i]]$param.lca[,j]
		}	
	}
			
	# Now getting standard error and means of different effects
	stnderrs.effects[1,]	<- sapply(as.data.frame(resamps.lca.eff), sd, na.rm=T)			
	avg.effects[1,]			<-	sapply(as.data.frame(resamps.lca.eff), mean, na.rm=T)
				
	for (j in 1:nitems) {
		stnderrs.lca.params[,j]	<- sapply(as.data.frame(resamps.lca.params[[j]]), sd, na.rm=T)						
	}
					
	# Assigning row names
	row.names(avg.effects) <- row.names(stnderrs.effects) <- c("lca")
				    			   
	if(nrow(stnderrs.lca.params)==8) {
		row.names(stnderrs.lca.params) <- c("lgg", "lgk", "lgc", "lkk", "lcg", "lck", "lcc", "gamma")
	} else { 
		row.names(stnderrs.lca.params) <- c("lgg", "lgk",  "lkk", "gamma")
	}
			
	# Get results out			    
	res 		<- list(stnderrs.lca.params, avg.effects, stnderrs.effects)
	names(res)	<- c("stnderrs.lca.params", "avg.effects", "stnderrs.effects") 

	res
}	