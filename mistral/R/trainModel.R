## -----------------------------------------------------------------------------
## Fonction trainModel
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------
#' @import e1071
#' @import DiceKriging
#' 
trainModel = function(model=NULL,
                      kernel,
                      design,
                      response,
                      updesign=design,
                      upresponse=response,
                      type="Kriging",
                      optim.method = "BFGS",
                      cost=NA,
                      gamma=NA,
                      coef.trend) {

    ## SVM
    
	if(type=="SVM") {
		Y = sign(response)
		threshold_ind = (Y==0)
		Y = factor(Y[!threshold_ind])
		design = t(design[,!threshold_ind])
		row.names(design) <- NULL
		if(is.null(model)) {
			model.new = tune.svm(design,
                                 Y,
                                 scale=FALSE,
                                 type="C-classification",
                                 cost=cost,
                                 gamma=gamma)$best.model
			cat("SVM parameters tunned :\n  gamma =",model.new$gamma,"\n  C =",model.new$cost,"\n")
		}
		else {
			gamma = model$gamma
			cost = model$cost
			model.new = svm(design,Y,scale=FALSE,type="C-classification",cost=cost,gamma=gamma)
		}

	}

    ## KRIGEAGE
    
	if(type=="Kriging") {
		if(is.null(model)) {
			design = t(design)
			dimnames(design) = list(NULL,NULL)
			if(missing(kernel)) kernel="matern5_2";
            ## KRIGING : CONSTANT TREND !!!
      if(missing(coef.trend)) {
        capture.output(model.new <- km(design   = design,
                                       response = response,
                                       covtype  = kernel,
                                       nugget.estim=TRUE,
                                       optim.method = optim.method #, #control = list(max.generations = 30, print.level = 0),
                                       # estim.method = "LOO",
                                       # lower   = rep(0.01, dim(design)[2]),
                                       # upper   = rep(100.0, dim(design)[2])#,
                                       # parinit = rep(0.5,   dim(design)[2])
        ))
      }
      else{
  			capture.output(model.new <- km(design   = design,
                             response = response,
                             covtype  = kernel,
                             coef.trend = coef.trend,
                             nugget.estim=TRUE,
                             optim.method = "gen"#, #control = list(max.generations = 30, print.level = 0),
#                              estim.method = "LOO",
  #                            lower   = rep(0.01, dim(design)[2]),
  #                            upper   = rep(100.0, dim(design)[2])#,
  #                            parinit = rep(0.5,   dim(design)[2])
                             ))
      }
			covariance = model.new@covariance

			cat("Kriging model parameters tunned\n")
            cat("-------------------------------\n")
            cat("  - cov_type =",covariance@name,"\n")
            cat("  - theta    =",covariance@range.val,"\n")
            cat("  - sd2      =",covariance@sd2,"\n")
            cat("  - trend    =",model.new@trend.coef,"\n")
		}
		else {
			updesign = t(updesign)
			dimnames(updesign) = list(NULL,NULL)
			model.new = update(model, 
                               newX=updesign, 
                               newy=upresponse, 
                               newX.alreadyExist = FALSE, 
                               cov.reestim = FALSE, 
                               trend.reestim = FALSE, 
                               nugget.reestim = FALSE)
		}
	}

	limit.meta = limitFunction(model.new)
	
	
	contour = function(x,y){
	  grid = expand.grid(x=x, y=y)
	  z_meta = limit.meta(t(grid))
	  z_crit = abs(z_meta$mean)/z_meta$sd
	  return(data.frame(grid, z = z_meta$mean, margin = z_crit))
	}
	
	res = list(model=model.new,fun=limit.meta, contour=contour)
	return(res)
}
