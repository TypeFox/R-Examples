`.run.cv` <-
function(x,y, fs.method, cross.outer,lambda1.set=NULL, class.weights, parms.coding="none", seed=123,maxIter=700, 
							verbose=TRUE	){

# Stratified cross-validation
  #  Same as cross-validation, except that the folds are stratified so that
  # they contain approximately the same proportions of labels as the original dataset.  

	if (verbose) print(paste("lambda1.set=", lambda1.set))
	
	nn<- nrow(x)
	nlevels.class<- nlevels(as.factor(y)) 
	levels.class <- levels(as.factor(y))
	
	if (cross.outer > nn | cross.outer < 2)  stop(paste("You have to specify at	least two different groups and not more than", nn))
	the.cut <- cut(1:nn,	cross.outer, 1:cross.outer) 
	# set seed 
	if (!is.null(seed)) set.seed(seed)
	
	# vote matrix 		
	votal.matrix <- matrix(0, ncol = nlevels.class, nrow = nn)
	the.vector.of.all.parameters <- vector(length = 0, mode = "list")
	model.info.list <- list()
	
	# maybe will implement in the future 		# several starts, like in 	MCRestimate # 
	#for (l in 1:cross.repeat) { 
	
	the.votes.per.cv <- matrix(NA, 	ncol = nlevels.class, nrow = nn) 
	rownames(the.votes.per.cv) <- y
	#rownames(the.votes.per.cv) <- names(y) 
	colnames(the.votes.per.cv) <- levels.class
	
	
	# create cut.list  for cross.outer - fold cv, which is not a leave one out cv!!!!
	if (cross.outer < nn){
		the.cut.list<- list() 	
		for (i in levels.class){ 
			the.cut.list[[i]] <- 	cut(1:sum( y == i ), cross.outer, 1:cross.outer)
		} 
		if (verbose){ 
			print("the.cut.list")
			print(the.cut.list)
		}
		
		# permute in each class separatly 	--> statifed cross validation
		perm.list<- lapply(the.cut.list, function (class.elems) sample (class.elems, length(class.elems))   )
		
		permutated.cut<- rep (NA, length(y))
		names(permutated.cut)<- names(y) 
		for (i in 1:length(perm.list)){ 
			permutated.cut[which(y == names(perm.list)[i] )]<- perm.list[[i]] 
		}
	}
	
	# if we have loo cv, just take one samle to test after the other, nothing to permute!  
	if (cross.outer == nn) {
	    permutated.cut<- c(1:nn)
	}
	
	#### end of create cut.list 
	
	# for Elastic NET
	set.opt.lam1<-NULL
	
	#################################### 	
	# do cv (change sample --> sample.i)
	#################################### 
	for (sample.i in 1:cross.outer) {
		if (verbose){ 
			print(paste("cv step ",sample.i, "of", cross.outer )) 
			print("")
		}
		block <- 	permutated.cut == sample.i
		train.matrix <- x[!block,, drop = FALSE]
		test.matrix <- x[ block, , drop = FALSE] 
		train.factor <- y[!block]
		test.factor<- y[block]
		table(train.factor)
		table(test.factor)
		
		
		### scad ###########################################################
		if (fs.method=="scad"){ 
			
			# scad for cv step # 1.do scv scad  + gacv --> 		optimal lambda # 2. prediction
			#1. construct model
			tmp.list<-.calc.scad(parms=lambda1.set, x.svm=train.matrix, y.svm=train.factor,  class.weights=class.weights,
									verbose= verbose, inner.val.method="none", parms.coding=parms.coding, maxIter=maxIter)
			
			fit<-tmp.list$model
			# 2. prediction (still in cv) : use test.matrix, test.factor
			# check if fit is not an empty model!
			if (fit[[1]][1] !="No variable selected."){
			
				predict.list<-predict(object=tmp.list, newdata=test.matrix, newdata.labels=as.factor(test.factor), labels.universe= levels(as.factor(y)))
					
				model.info<-list(class=test.factor, 
											pred.class= predict.list$pred.class,		
											tab=predict.list$tab, 
											accurancy =1 - predict.list$error, 
											w =fit$w,
											b = fit$b,
											lambda1 = fit$lambda1,
											lambda2 = NULL )
			} else # end of check if fit is not an empty model!	
			{
				if (verbose) print(paste("we have an empty model in the cv step ",sample.i, "of", cross.outer ))
			}
										
		}		# end of scad svm for cv step
		
		### 1norm ###########################################################
		if (fs.method=="1norm"){# 1norm svm nu = EstNuShort(train.matrix,		train.factor) 
			
			# 1. construct model 
			#		epsi - tuning parameter !
			#		find the best epsi via k-fold cv,  get the finla model with optimal epsi
			#fit<-run.1norm(x=train.matrix,y=train.factor,k=5,nu=0, output=1, seed=seed)
			fit<-.calc.1norm(parms=lambda1.set, x.svm=train.matrix, y.svm=train.factor, maxIter=maxIter,
											 class.weights=class.weights, verbose= verbose, 
											 inner.val.method="none", cross.inner=0,
											  parms.coding=parms.coding, seed=seed )
		
		
			# 2. prediction (still in cv) : use test.matrix, test.factor
			if (fit[[1]][1] !="No variable selected."){
				predict.list<-predict(object=tmp.list, newdata=test.matrix, newdata.labels=as.factor(test.factor), labels.universe= levels(as.factor(y)))
		
				model.info<-list(class=test.factor,
										pred.class= predict.list$pred.class,
										tab=predict.list$tab, 
										accurancy =1- predict.list$error, 
										nu =fit$nu,
										w =fit$w, 
										b= fit$b,
										lambda1 = fit$epsi,
										lambda2 = NULL )
			} else # end of check if fit is not an empty model!	
			{
				if (verbose) print(paste("we have an empty model in the cv step ",sample.i, "of", cross.outer ))
			}							
		}	#end of 1norm svm for cv step
			

		### drHSVM TODO   ###########################################################
		if (fs.method=="DrHSVM"){ 
			# scad for cv step # 1.do scv scad  + gacv --> 		optimal lambda # 2. prediction
			
					
			#1. construct model
			tmp.list<-.calc.DrHSVM(parms=lambda1.set, x.svm=train.matrix, y.svm=train.factor,  class.weights=class.weights,
									inner.val.method="none", verbose= verbose,  parms.coding=parms.coding, maxIter=maxIter)
			
			# 2. prediction (still in cv) : use test.matrix, test.factor
			predict.list<-.DrHSVM.predict(object=tmp.list$model, newx=test.matrix, newy=as.factor(test.factor), eps = 1e-10)
			predict.list$pred.class<- sign(predict.list$fit)
				
			# find optimal lambda1 with min error
			position.min.err<- which(predict.list$err == min(predict.list$err) )
			opt.lam1s<- tmp.list$model$lambda1[position.min.err ]
			
			if (verbose) print("# FS: ")
			len.w<-  apply(tmp.list$model$beta,1, function(betas) sum(betas !=0))
			sel.w<-  len.w[position.min.err]
			names(sel.w)<-opt.lam1s
			if (verbose)  print(sel.w)
			
			if (verbose) print("chose the model with min num of FS ")
			opt.lam1<- opt.lam1s[which.min(sel.w)]
			if (verbose) print(opt.lam1)
			
			# end model for opt.lam1
			pos.opt.lam1<- which(tmp.list$model$lambda1 == opt.lam1	)
			
			# coefs w,b
			w<- tmp.list$model$beta[pos.opt.lam1, ]
			names(w)<-colnames(train.matrix)
			w<-w[w!=0]
			b<- tmp.list$model$beta0[pos.opt.lam1 ]
			
			# table
			pred.class<- predict.list$pred.class[,pos.opt.lam1 ]
			tab<-table(predict.list$pred.class[, pos.opt.lam1], test.factor)
		
			# misclass error
			err<-	min(predict.list$err)
						
			model.info<-list(class=test.factor, 
											pred.class= pred.class,		
											tab=tab, 
											accurancy =1 - err, 
											w =w,
											b = b,
											lambda1 = opt.lam1,
											lambda2 = lambda1.set )
											
			set.opt.lam1<-c(set.opt.lam1, opt.lam1) 								
											
		}		# end of drHSVM svm for cv step
		
		### scad + L2   ###########################################################
		if (fs.method=="scad+L2"){ 
			# scad for cv step # 1.do scv scad  + gacv --> 		optimal lambda # 2. prediction
			
			
			if (verbose) print("lambda1.set")
			if (verbose) print(lambda1.set)
			
			# 1.do scv scad  + gacv --> optimal lambda 			ff.list<-list()
			tmp.list<- .calc.scad_L2(parms=lambda1.set, x.svm=train.matrix, y.svm=train.factor,  maxIter=maxIter, 
										class.weights=class.weights,verbose= verbose, inner.val.method="none", cross.inner=0,
											parms.coding=parms.coding, seed=seed )
			
			fit<-tmp.list$model
			# 2. prediction (still in cv) : use test.matrix, test.factor and fixed destingiged lambda.path.drHSVM from the whole data set
			if (fit[[1]][1] !="No variable selected."){
				predict.list<-predict(object=tmp.list, newdata=test.matrix, newdata.labels=as.factor(test.factor), labels.universe= levels(as.factor(y)) )
			
				model.info<-list(class=test.factor, 
											pred.class= predict.list$pred.class,		
											tab=predict.list$tab, 
											accurancy =1 - predict.list$error, 
											w =fit$w,
											b = fit$b,
											lambda1 = fit$lambda1,
											lambda2 = fit$lambda2 )
			} else # end of check if fit is not an empty model!	
			{
				if (verbose) print(paste("we have an empty model in the cv step ",sample.i, "of", cross.outer ))
			}									
		}		# end of scad svm for cv step
		### end of scad + L2 ###########################################################				
						
		if (exists("model.info")){		
			model.info.list[[ sample.i]] <- model.info 
							
			pred.vector <- 	model.info$pred.class 
			vote.matrix <- t(sapply(1:length(pred.vector), 	function(j) as.numeric(levels.class == pred.vector[j])))
			colnames(vote.matrix) <- levels.class 
			the.votes.per.cv[block, ] <- 		vote.matrix
			
			rm( model.info)
		}	
	} # end of cv
			
	votal.matrix <- votal.matrix + the.votes.per.cv 		# creating  the confusion table 
	
	# count NA as a non-correct prediction
	res <- .whatiscorrect(votal.matrix, count.na=TRUE)
		
	vote.table <- table(rownames(votal.matrix), res$best.vote) 	
	
	
	new.table <- 	matrix(0, 
											ncol=nrow(vote.table), 
											nrow=nrow(vote.table),
											dimnames=list(rownames(vote.table),rownames(vote.table)))
	new.table[,colnames(vote.table)] <- vote.table
	
	
	normed.table <- new.table/rowSums(new.table) 		
	confusion <- cbind(new.table, 1-diag(normed.table)) 
	colnames(confusion) <- 	c(levels(as.factor(y)), "class error")
	
	
	# missclassification error from the vote table 
	cv.error<- 1- sum(diag(vote.table))/sum(vote.table) 

	res.cv<- list(vote.table=vote.table,
								cv.error = cv.error,
								res= res,
								normed.table =normed.table,
								confusion =confusion,
								seed=seed,
								set.opt.lam1 =set.opt.lam1 
								, model.info.list = model.info.list )
	
	return(res.cv)
}

`.whatiscorrect` <-
function(votematrix, count.na=TRUE) {
# count.na=TRUE  count NA as a non correct prediction
#					= FALSE skip NAs
  correct.class.vote <- numeric(nrow(votematrix))
  correct.prediction <- logical(nrow(votematrix))
  best.vote          <- character(nrow(votematrix))
  
  votematrix.orig<-votematrix # for debugging
  
  for(i in 1:nrow(votematrix)) {
   
   # if count NA --> set the predicted value to the wrong one!
   if (count.na){
	   if (any(is.na(votematrix[i,]))){
	    wrong.pred<-setdiff( levels(as.factor(rownames(votematrix))) , rownames(votematrix)[i] )
	   	
	   	for (col in 1:ncol(votematrix)) 
	   	 votematrix[i,col]<- ifelse (wrong.pred == colnames(votematrix)[col], 1,0)  
	   } 
   }
    correct.class.vote[i] <- votematrix[i, rownames(votematrix)[i]==colnames(votematrix) ]
    correct.prediction[i] <- rownames(votematrix)[i] %in% colnames(votematrix)[votematrix[i,]==max(votematrix[i,])]
    if (!count.na) correct.prediction[i][is.na( correct.class.vote[i])] <- NA
   
   best.vote[i] <- ifelse(correct.prediction[i],
                          rownames(votematrix)[i],
                          colnames(votematrix)[which.max(votematrix[i,])])
  }
  
    if (count.na)  	correct.class.vote[is.na( correct.class.vote)]<- 0
 
  
  return(list(correct.class.vote=correct.class.vote,
              correct.prediction=correct.prediction,
              best.vote=best.vote,
              count.na=count.na))
}
