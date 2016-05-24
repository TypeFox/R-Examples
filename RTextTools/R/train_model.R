train_model <- function(container,algorithm=c("SVM","SLDA","BOOSTING","BAGGING","RF","GLMNET","TREE","NNET","MAXENT"),
						method="C-classification", cross=0, cost=100, kernel="radial",  # SVM PARAMETERS
						maxitboost=100, # BOOSTING PARAMETERS
						maxitglm=10^5, # GLMNET PARAMETERS
						size=1,maxitnnet=1000,MaxNWts=10000,rang=0.1,decay=5e-4,trace=FALSE, # NNET PARAMETERS
						ntree=200, # RF PARAMETERS
						l1_regularizer=0.0,l2_regularizer=0.0,use_sgd=FALSE,set_heldout=0,verbose=FALSE, # MAXENT PARAMETERS
						...) {
        
        # CLEAN UP FROM PREVIOUS MODEL TRAINED
        gc()
        
        # CONDITIONAL TRAINING OF MODEL
        if (algorithm=="SVM") {
			model <- svm(x=container@training_matrix, y=container@training_codes, method=method, cross=cross, cost=cost, probability=TRUE, kernel=kernel)
		} else if (algorithm=="SLDA") {
           model <- slda(container.training_codes ~ ., data=data.frame(as.matrix(container@training_matrix),container@training_codes))
        } else if (algorithm=="BOOSTING") {
            model <- LogitBoost(xlearn=as.matrix(container@training_matrix), ylearn=container@training_codes, nIter=maxitboost)
        } else if (algorithm=="BAGGING") {
            model <- bagging(container.training_codes ~ ., data=data.frame(as.matrix(container@training_matrix),container@training_codes))
        } else if (algorithm=="RF") {
            model <- randomForest(x=as.matrix(container@training_matrix), y=container@training_codes, ntree=ntree)
        } else if (algorithm=="GLMNET") {
			training_matrix <- as(container@training_matrix,"sparseMatrix")
            model <- glmnet(x=training_matrix, y=container@training_codes, family="multinomial", maxit=maxitglm)
        } else if (algorithm=="TREE") {
            model <- tree(container.training_codes ~ ., data=data.frame(as.matrix(container@training_matrix),container@training_codes))
        } else if (algorithm=="NNET") {
            model <- nnet(container.training_codes ~ ., data=data.frame(as.matrix(container@training_matrix),container@training_codes), size=size, maxit=maxitnnet, MaxNWts=MaxNWts, rang=rang, decay=decay, trace=trace)
        } else if (algorithm=="MAXENT") {
			model <- maxent(container@training_matrix,as.vector(container@training_codes),l1_regularizer,l2_regularizer,use_sgd,set_heldout,verbose)
		} else {
			stop("ERROR: Invalid algorithm specified. Type print_algorithms() for a list of available algorithms.")
		}
		
		# RETURN TRAINED MODEL
		gc() # CLEAN UP AFTER MODEL
		return(model)
}