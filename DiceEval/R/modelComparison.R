modelComparison <- function (X,Y,type="all",K=10,test = NULL,...) {

    if (length(type) == 1 && type == "all") {
        type <- c("Linear", "StepLinear","Additive", "PolyMARS", "MARS","Kriging")
    }
	
    data <- data.frame(X, Y)
    X <- data.frame(X)
    f <- dim(X)[2]
    n <- dim(data)[1]
	
    # Learning criteria
	learningCrit 			<- data.frame(matrix(0, nrow = 2, ncol = length(type)))
    rownames(learningCrit) <- c("R2", "RMSE")
	# Cross Validation criteria
	CVCrit <- data.frame(matrix(0, nrow = 2, ncol = length(type)))
    rownames(CVCrit) <- c("Q2", "RMSE")
	# Test criteria
	testCrit <- data.frame(matrix(0, nrow = 2, ncol = length(type)))
	rownames(testCrit) <- c("R2", "RMSE")
	
    # Vérification des données (apprentissage et test de même dim.)
	if (!is.null(test)){
		if (dim(test)[2] != dim(data)[2]) {
	        warning("The dimensions of the test set are not correct.")
			test <- NULL
	    }
	}

	# Pour les labels besoin de connaître le nombre de modèles à comparer
	Nb_model <- data.frame(matrix(0,nrow=2,ncol=6))
	colnames(Nb_model) <- c("Linear", "StepLinear","Additive","MARS","PolyMARS","Kriging")
	
	argList <- list(...)
	id_fmla <- 0; id_deg <- 0; p <- 0; id_gcv <- 0
	for (i in 1:length(type)){
	  type_ <- type[i]
	  switch(type_, Linear   = {Nb_model[1,1] = Nb_model[1,1]+1}, 
					Additive = {Nb_model[1,3] = Nb_model[1,3]+1},
					StepLinear = {Nb_model[1,2] = Nb_model[1,2]+1},
					MARS ={Nb_model[1,4] = Nb_model[1,4]+1},
					PolyMARS = {Nb_model[1,5] = Nb_model[1,5]+1},
					Kriging = {Nb_model[1,6] = Nb_model[1,6]+1})
					
	  if(type_ == "Linear" | type_ == "Additive"){
		# Récupère la formule + vérification cvompatibilité des entrées
		# (autant de modèles que de formules -- idem pour le stepwise et
		# le krigeage)
		if(!is.null(argList$formula)){
			# Récupère le nb de formules (l_fmla) en entrée
			if (class(argList$formula)=="formula"){
				l_fmla <- 1
			} else l_fmla <- length(argList$formula)
			
			# Récupère la formule (fmla) et vérifie que tout est ok.
			if(l_fmla==1){
				fmla <- argList$formula
			} else if(l_fmla>=2 & l_fmla<=length(type)){
				id_fmla <- id_fmla+1
				if (id_fmla > length(type) | id_fmla > l_fmla){
					stop("Argument \'formula'\ has not a valid length.")
				}
				fmla <- argList$formula[[id_fmla]]
			} else if (l_fmla > length(type)){
					stop("Argument \'formula'\ has not a valid length.")
			}
		} else { # dans le cas où rien n'est indiqué ==> formules par défaut
			switch(type_, Linear   = {warning("[Linear model] Argument \'formula\' not found, set at \'Y~.\'")
										fmla <- formulaLm(X,Y)}, 
					  Additive = {  warning("[Additive model] Argument \'formula\' not found, set at \'Y~s(X1)+...+s(Xp)\'")
									fmla <- formulaAm(X,Y)})
		}
		modTmp <- modelFit(X, Y, type = type_, formula = fmla)
	  } else if (type_ == "StepLinear"){
		if(!is.null(argList$formula)){
			if (class(argList$formula)=="formula"){
				l_fmla <- 1
			} else l_fmla <- length(argList$formula)
			
			if(l_fmla==1){
				fmla <- argList$formula
			} else if(l_fmla>=2 & l_fmla<=length(type)){
				id_fmla <- id_fmla+1
				if (id_fmla > length(type) | id_fmla > l_fmla){
					stop("Argument \'formula'\ has not a valid length.")
				}
				fmla <- argList$formula[[id_fmla]]
			} else if (l_fmla > length(type)){
				stop("Argument \'formula'\ has not a valid length.")
			}
		} else {
			warning("[StepLinear model] Argument \'formula\' not found, set at \'Y~.\'")
			fmla <- formulaLm(X,Y)
		}
		init <- modelFit(X,Y,type="Linear",formula=fmla)
			if (length(init$coefficients)>n) {
				stop("There are too many terms into the full length model")
			}
			if (!is.null(argList$penalty)){
				if(length(argList$penalty)==1){
					penalty <- argList$penalty
				} else if(length(argList$penalty)>=2 & length(argList$penalty)<= length(type)){
					p <- p+1
					penalty <- argList$penalty[p]
				} else if(length(argList$penalty)> length(type)){
					stop("Argument \'penalty'\ has not a valid length.")
				}
			} else {
				warning("[StepLinear model] Argument \'penalty\' not found, set at \'2\' (AIC criteria)")
				penalty <- 2
			}
			modTmp <- modelFit(X, Y, type = type_, formula = fmla,penalty=penalty)
		} else if (type_ == "MARS"){
			if(!is.null(argList$degree)){
				l_deg <- length(argList$degree)
				if(l_deg == 1){ id_deg <- 1}
				else if(l_deg >= 2 & l_deg <= length(type) ){
					id_deg <- id_deg+1
					if (id_deg > length(type) | id_deg > l_deg){
						stop("Argument \'formula'\ has not a valid length.")
					}
				}
				degree <- argList$degree[[id_deg]]
		    } else {
			  warning("[MARS model] Argument \'degree\' not found, set at \'2\'")
			  degree <- 2
			}
			modTmp <- modelFit(X, Y, type = type_, degree=degree)
		} else if (type_ == "PolyMARS"){
			if(!is.null(argList$gcv)){
				l_gcv <- length(argList$gcv)
				if(l_gcv == 1){ 
					id_gcv <- 1
				} else if (l_gcv >= 2 & l_gcv <= length(type)){
					id_gcv <- id_gcv+1
					if (id_gcv > length(type) | id_gcv > l_gcv){
						stop("Argument \'gcv'\ has not a valid length.")
					}
				} else if (l_gcv > length(type)){
					stop("Argument \'gcv'\ has not a valid length.")
				}
				gcv <- argList$gcv[[id_gcv]]
			} else {
				warning("[PolyMARS model] Argument \'gcv\' not found, set at \'4\'")
				gcv <- 4
			}
			modTmp <- modelFit(X, Y, type = type_, gcv=gcv)
		} else if (type_ == "Kriging"){
		  if(!is.null(argList$formula)){
			if (class(argList$formula)=="formula"){ l_fmla <- 1
			} else l_fmla <- length(argList$formula)
			
			if(l_fmla==1){	fmla <- argList$formula
			} else if(l_fmla>=2 & l_fmla<=length(type)){
				id_fmla <- id_fmla+1
				if (id_fmla > length(type) | id_fmla > l_fmla){
					stop("Argument \'formula'\ has not a valid length.")
				}
				fmla <- argList$formula[[id_fmla]]
			} else if (l_fmla > length(type)){
					stop("Argument \'formula'\ has not a valid length.")
			}
		  } else {
			warning("[Kriging model] Argument \'formula\'  set at \'~1\'")
			fmla <- ~1
		  }
		  id_covtype <- 0
			if(!is.null(argList$covtype)){
				if (class(argList$formula)=="vector"){ 
					l_covtype <- length(argList$covtype)
				} else l_covtype <- 1
				if(l_covtype==1){
					covtype <- argList$covtype
				} else if(l_covtype>=2 & l_covtype<=length(type)){
					id_covtype <- id_covtype+1
					if (id_covtype > length(type) | id_covtype > l_covtype){
						stop("Argument \'covtype'\ has not a valid length.")
					}
					covtype <- argList$covtype[[id_covtype]]
					} else if (l_fmla > length(type)){
						stop("Argument \'formula'\ has not a valid length.")
					}
				} else {
					warning("[Kriging model] Argument \'covtype\'  set at \"matern5_2\"")
							covtype="matern5_2"
				}	
			modTmp <- modelFit(X, Y, type = type_,formula=fmla,covtype=covtype)
		}
		
		if(type_ == "Linear" | type_ == "StepLinear" | type_ == "Additive" | type_ == "MARS"){
			learningCrit[1, i] <- R2(Y, modTmp$model$fitted.values)
			learningCrit[2, i] <- RMSE(Y, modTmp$model$fitted.values)
		} else if(type_ == "PolyMARS"){
			learningCrit[1, i] <- R2(Y, modTmp$model$fitted)
			learningCrit[2, i] <- RMSE(Y, modTmp$model$fitted)
		} else if(type_=="Kriging"){
			# Vérification que les données prédites aux points du plan interpolent
			# les données d'apprentissage
			tmp_ <- modelPredict(modTmp,newdata=X)
			if(sum(abs(Y-tmp_)) >= 10e-10){
				warning("Warning: Kriging model doesn't interpolate the data.")
			}
			learningCrit[1, i] <- R2(Y, tmp_)
			learningCrit[2, i] <- RMSE(Y, tmp_)
		}
		cvTmp <- crossValidation(modTmp,K=10)
		CVCrit[1,i] <- cvTmp$Q2
		CVCrit[2,i] <- cvTmp$RMSE_CV
		if (!is.null(test)){
			Ytest 	<- modelPredict(model = modTmp, newdata = test[, 1:f])
			testCrit[1, i] <- R2(test[,f+1], Ytest)
			testCrit[2, i] <- RMSE(test[,f+1], Ytest)
		}
	}
		
    if (min(testCrit[1, ]) < 0) {
        warning("R2 calculated for the test data is negative.")
    }
    if (min(CVCrit[1, ]) < 0) {
        warning("R2 estimated by cross validation is negative.")
    }
	
    if (!is.null(test)){
		R2Crit <- data.frame(t(learningCrit[1, ]), t(CVCrit[1,]), t(testCrit[1, ]))
	    colnames(R2Crit) <- c("learning", "cross-validation","test")
	    RMSECrit <- data.frame(t(learningCrit[2, ]), t(CVCrit[2,]), t(testCrit[2, ]))
	    colnames(RMSECrit) <- c("learning", "cross-validation","test")
	} else {
		R2Crit <- data.frame(t(learningCrit[1, ]), t(CVCrit[1,]))
	    colnames(R2Crit) <- c("learning", "cross-validation")
	    RMSECrit <- data.frame(t(learningCrit[2, ]), t(CVCrit[2,]))
	    colnames(RMSECrit) <- c("learning", "cross-validation")
	}
	
	
	c_names <- NULL
	for (i in 1:length(type)){
		type_ <- type[i]
		
			switch(type_, Linear 	= {val <- 1}, 
					Additive 		= {val <- 3},
					StepLinear		= {val <- 2},
					MARS 			= {val <- 4},
					PolyMARS 		= {val <- 5},
					Kriging 		= {val <- 6})
		Nb_model[2,val] <- Nb_model[2,val]+1
		if(Nb_model[1,val]==1){
			c_names <- c(c_names,colnames(Nb_model)[val])
		} else c_names <- c(c_names,paste(colnames(Nb_model)[val],Nb_model[2,val],sep=" "))
	}
	
	colnames(learningCrit) 	<- c_names
	colnames(CVCrit) 		<- c_names
	colnames(testCrit) 		<- c_names
	rownames(R2Crit) <- c_names
	rownames(RMSECrit) <- c_names
	
#    op <- par(mfrow = c(1, 2), cex = 0.7)
    op <- par(ask=TRUE, cex = 0.7)
    dotchart(as.matrix(R2Crit), pch = 19, xlim = c(0, 1),xlab = list("Comparison of R2", cex=1.2, font=2))
	axis(3)
	abline(v = axTicks(1), lty = "dotted", col = "gray60")
    dotchart(as.matrix(RMSECrit), pch = 19, xlab = list("Comparison of RMSE", cex=1.2, font=2), xlim = range(RMSECrit))
	abline(v = axTicks(1), lty = "dotted", col = "gray60")
	axis(3)
    par(op)

	if (is.null(test)){
        return(list(Learning = learningCrit, CV = CVCrit))
    } else return(list(Learning = learningCrit, CV = CVCrit, Test = testCrit))
}