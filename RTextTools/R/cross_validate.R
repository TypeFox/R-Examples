cross_validate <- function(container,nfold,algorithm=c("SVM","SLDA","BOOSTING","BAGGING","RF","GLMNET","TREE","NNET","MAXENT"),seed=NA,
							method="C-classification", cross=0, cost=100, kernel="radial",  # SVM PARAMETERS
							maxitboost=100, # BOOSTING PARAMETERS
							maxitglm=10^5, # GLMNET PARAMETERS
							size=1,maxitnnet=1000,MaxNWts=10000,rang=0.1,decay=5e-4, # NNET PARAMETERS
							ntree=200, # RF PARAMETERS
							l1_regularizer=0.0,l2_regularizer=0.0,use_sgd=FALSE,set_heldout=0,verbose=FALSE # MAXENT PARAMETERS
							) {

    options(warn=-1) #Supress warnings
    if (!is.na(seed)) #Set seed for replicability
        set.seed(seed)
    extract_label_from_prob_names <- function(x) return(rownames(as.matrix(which.max(x))))
    #Bring in info from the container function    
	alldata <- rbind(container@training_matrix,container@classification_matrix) #put all data together
	allcodes <- as.factor(c(container@training_codes,container@testing_codes))
    rand <- sample(nfold,dim(alldata)[1], replace=T) #replace
    
    cv_accuracy <- NULL
    for (i in sort(unique(rand))) {
        if (algorithm=="SVM") {
            model <- svm(x=alldata[rand!=i,], y=allcodes[rand!=i],method=method,cross=cross,cost=cost,kernel=kernel) #put function here
            pred <- predict(model,alldata[rand==i,])
        } else
		if (algorithm=="SLDA") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			colnames(alldata) <- container@column_names
			data_and_codes <- cbind(as.matrix(alldata),allcodes)
			model <- slda(as.factor(allcodes)~.,data=data.frame(data_and_codes[rand!=i,]))
            pred <- predict(model,data.frame(alldata[rand==i,]))
			pred <- as.numeric(pred$class)
		} else
        
        if (algorithm=="RF") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			#colnames(alldata) <- container@column_names
			data_and_codes <-cbind(alldata,allcodes)
            model <- randomForest(as.factor(allcodes)~.,data=data_and_codes[rand!=i,], ntree=ntree)
            pred <- predict(model,newdata=alldata[rand==i,])
        } else
        if (algorithm=="GLMNET") {
			sparsedata <- as(as.matrix.csc(alldata[rand!=i,]),"dgCMatrix")
            model <- glmnet(x=sparsedata, y=as.vector(allcodes[rand!=i]),family="multinomial", maxit=maxitglm)
			
            prob <- predict(model,sparsedata,s=0.01,type="response")            
            pred <- apply(prob[,,1],1,extract_label_from_prob_names)
            pred <- as.numeric(pred)
 
        } else
        if (algorithm=="BOOSTING") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			colnames(alldata) <- container@column_names
			data_and_codes <- cbind(alldata,allcodes)
            model <- LogitBoost(xlearn=alldata[rand!=i,], ylearn=allcodes[rand!=i],maxitboost)
            pred <- predict(model,data.frame(alldata[rand==i,]))
        } else
		if (algorithm=="BAGGING") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			#colnames(alldata) <- container@column_names
			data_and_codes <-cbind(alldata,allcodes)
            model <- bagging(as.factor(allcodes)~.,data=data.frame(data_and_codes[rand!=i,]))
            pred <- predict(model,newdata=alldata[rand==i,])
		} else
        if (algorithm=="TREE") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			colnames(alldata) <- container@column_names
            data_and_codes <- cbind(alldata,allcodes)
            model <- tree(as.factor(allcodes)~ ., data = data.frame(data_and_codes[rand!=i,]))
            prob <- predict(model,newdata=data.frame(alldata[rand==i,]), type="vector")
            pred <- apply(prob,1,which.max)

        } else
        if(algorithm=="NNET") {
			alldata <- rbind(as.matrix(container@training_matrix),as.matrix(container@classification_matrix))
			colnames(alldata) <- container@column_names
			data_and_codes <- cbind(alldata,allcodes)
            model <- nnet(as.factor(allcodes)~ ., data = data.frame(data_and_codes[rand!=i,]),size=size,maxit=maxitnnet,MaxNWts=MaxNWts,rang=rang,decay=decay,trace=FALSE)
            prob <- predict(model,newdata=data.frame(alldata[rand==i,]))
            pred <- apply(prob,1,which.max)
        } else
		if (algorithm=="MAXENT") {
			model <- maxent(container@training_matrix,as.vector(container@training_codes),l1_regularizer,l2_regularizer,use_sgd,set_heldout,verbose)
			pred <- predict(model,alldata[rand==i,])
			pred <- pred[,1]
		}

        cv_accuracy[i] <- recall_accuracy(allcodes[rand==i],pred)
		
        cat("Fold ",i," Out of Sample Accuracy"," = ",cv_accuracy[i],"\n",sep="")
    }

	return(list(cv_accuracy,meanAccuracy=mean(cv_accuracy)))
}
