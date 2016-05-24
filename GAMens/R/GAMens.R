GAMens <- 
function(formula, data, rsm_size=2, autoform=FALSE, iter=10, df=4, bagging=TRUE, rsm=TRUE, fusion="avgagg") 
{

	formula<- as.formula(formula)
	setdiff.data.frame <- function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]
	if (!(as.character(fusion) %in% c("avgagg","majvote","w.avgagg","w.majvote"))){
		stop("fusion must be either 'avgagg', 'majvote', 'w.avgagg' or 'w.majvote' ")
	}
	if ((rsm == TRUE) & !(rsm_size > 0)){
		stop("Enter a valid random feature subspace size (rsm_size)")
	}
	if (!(rsm == TRUE) & !(bagging == TRUE)){
		stop("Either bagging, rsm or both must be set to TRUE")
	}
	if ((autoform == TRUE) & !(formula[3] == ".()")){
		warning("Autoform option will be overridden by the formula specification. In order to use autoform, set formula to '[outcome variable name]~.'")
		autoform <- FALSE
	}
	if ((autoform == FALSE) & (formula[3] == ".()")){
		warning("Autoform option is FALSE and a generic formula is used. The ensemble classifier will consist out of GAMs without nonparametric terms.")
		autoform <- FALSE
	}
	if ((rsm==TRUE) & (autoform == FALSE) & (!(formula[3] == ".()") & length(all.vars(formula[[3]]))<rsm_size)){
		stop("rsm_size should be smaller than the number of explanatory variables specified in formula")
	}

	targets <- data[,as.character(formula[[2]])]  
	nclasses <- nlevels(targets)
	depvarname <- as.character(formula[[2]])

	target_classes <- unique(targets)
	target_classes_s <- target_classes[order(target_classes)]
	newtargets <- as.numeric(data[,ncol(data)] == target_classes_s[2])
	newdata = cbind(data[,1:ncol(data)-1],newtargets)
	names(newdata)[ncol(newdata)] <- depvarname

	
	n <- length(data[,1])
	if (formula[3] == ".()"){
		p <- (length(newdata[1,])-1)
 		varnames <- as.matrix(setdiff.data.frame(as.matrix(names(data)),as.matrix(depvarname)))
	}else {p <- nrow(as.matrix(all.vars(formula)))-1
		varnames <- as.matrix(setdiff.data.frame(as.matrix(cbind(all.vars(formula))),as.matrix(depvarname)))
		formula_terms <- attr(terms(formula),"term.labels")
	}
	
     	gam_models <- list() 
	oob_predictions_p <- data.frame(cbind(1:nrow(newdata)))
	oob_predictions_c <- data.frame(cbind(1:nrow(newdata)))
	names(oob_predictions_p) <- "ID"
	names(oob_predictions_c) <- "ID"
	bootstraps <- array(0, c(n,iter))
	errors <- array(0,iter)
	
	id_added <- cbind(newdata,1:nrow(newdata))
	names(id_added)[ncol(id_added)] <- "ID"
	cutoff <- 0.5
	treshold <- 2	
	for (m in 1:iter) {
		flag <- FALSE
		gam_model <- NA
		while (is.na(gam_model[1])) {
			if (rsm==TRUE) {rfs <- sample(1:p,rsm_size,replace=FALSE)} else {rfs <- 1:p}
			if (bagging==TRUE) {bootstrap<- sample(1:n,replace=TRUE)} else {bootstrap <- 1:n}
			if (formula[3] == ".()") {
				selectedvars <- names(newdata[bootstrap,rfs])
				if (autoform==TRUE) {
					rfs_linvars <- character(0)
					rfs_nparvars <- character(0)
					for (vn in 1:length(selectedvars)) {
						varname <- selectedvars[vn]
						uniq <- dim(unique(as.data.frame(data[bootstrap, rfs[vn]])))[1]
						if (uniq<=treshold) {rfs_linvars <- rbind(rfs_linvars,varname) } 
						if (uniq>treshold) {rfs_nparvars <- rbind(rfs_nparvars,varname) }
					}
					if (length(rfs_nparvars > 0)) {
						npar_form1 <- paste("s(",rfs_nparvars,",",df,")")
						npar_form2 <- paste(npar_form1,collapse = "+")
					} else {npar_form2 <- character(0)}
					if (length(rfs_linvars > 0)) {
						lin_form <- paste(rfs_linvars,collapse= "+")
						if (length(rfs_nparvars > 0)) {finalstring = paste(depvarname,"~",npar_form2,"+",lin_form)} 
						else {finalstring <- paste(depvarname,"~",lin_form)}
					} else {finalstring = paste(depvarname,"~",npar_form2)}
					fmla <- as.formula(finalstring)
				}
				if (autoform==FALSE) {
					lin_form <- paste(selectedvars,collapse= "+")
					finalstring = paste(depvarname,"~",lin_form)
					fmla <- as.formula(finalstring)
				}
			} else {
				selected_terms <- formula_terms[rfs]
				lin_form <- paste(selected_terms,collapse= "+")
				finalstring = paste(depvarname,"~",lin_form)
				fmla <- as.formula(finalstring)
			}

		try(gam_model <- gam(fmla,data=newdata[bootstrap,], family=binomial(link="logit")),silent=TRUE)	
		}
		if (bagging == TRUE) {oob_data <- setdiff.data.frame(id_added,id_added[bootstrap,])} else {oob_data <- id_added[bootstrap,]}
		oob_cnt <- nrow(oob_data)
		oob_predict <- predict.gam(gam_model,oob_data,type="response")
		oob_predict_r <- as.data.frame(cbind(as.numeric(oob_predict),oob_data[,"ID"]))
		names(oob_predict_r)[2] <- "ID"
		coln <- paste("pred",m,sep="")
		names(oob_predict_r)[1] <- coln
		oob_predict_c <- as.data.frame(cbind((oob_predict_r[,1] > cutoff),oob_predict_r[,2]))
		names(oob_predict_c)[2] <- "ID"
		names(oob_predict_c)[1] <- coln
		ind<-as.numeric(oob_data[depvarname] != oob_predict_c[,1]) 	
		err<- sum(ind)/oob_cnt
		errors[m] <- err
		                     
		oob_predictions_p <- merge(oob_predictions_p,oob_predict_r, by = "ID", all= "TRUE")
		oob_predictions_c <- merge(oob_predictions_c,oob_predict_c, by = "ID", all= "TRUE")
		gam_models[[m]] <- gam_model			
		bootstraps[,m]<-bootstrap

	}
	temp_oob_pred_p <- oob_predictions_p
	temp_oob_pred_p[is.na(temp_oob_pred_p)] <- 0
	temp_oob_pred_c <- oob_predictions_c
	temp_oob_pred_c[is.na(temp_oob_pred_c)] <- 0
	temp_sums <- rowSums(temp_oob_pred_p[,2:ncol(temp_oob_pred_p)])
	temp_sums_cl <- rowSums(temp_oob_pred_c[,2:ncol(temp_oob_pred_p)])
	temp_n <- rowSums(!is.na(oob_predictions_p[,2:ncol(temp_oob_pred_p)]))
	
	if (fusion == "avgagg") {
		pred <- cbind(temp_sums / temp_n)
		class <- as.numeric(pred > cutoff)
	}else if (fusion == "w.avgagg") {
		temp_sums_weighted <- as.matrix(temp_oob_pred_p[,2:ncol(temp_oob_pred_p)]) %*% (1 - errors)
		temp_n_weighted <- as.matrix(!is.na(oob_predictions_p[,2:ncol(temp_oob_pred_p)])*1) %*% (1 - errors)
		pred <- temp_sums_weighted / temp_n_weighted
		class <- as.numeric(pred > cutoff)
	}else if (fusion == "majvote") {
		pred <- cbind(temp_sums_cl / temp_n)
		class <- as.numeric(pred > cutoff)
	}else if (fusion == "w.majvote") {
		temp_sums_weighted <- as.matrix(temp_oob_pred_c[,2:ncol(temp_oob_pred_c)]) %*% (1 - errors)
		temp_n_weighted <- as.matrix(!is.na(oob_predictions_c[,2:ncol(temp_oob_pred_c)])*1) %*% (1 - errors)
		pred <- temp_sums_weighted / temp_n_weighted
		class <- as.numeric(pred > cutoff)}
	class <-  cbind(gsub(1,target_classes_s[[2]],class,fixed=FALSE))
	class <-  cbind(gsub(0,target_classes_s[[1]],class,fixed=FALSE))
	
	ans<- list(GAMs=gam_models, formula=fmla, iter=iter, df=df, rsm=rsm, bagging=bagging, rsm_size=rsm_size, fusion_method=fusion, probs=pred, class=class, samples=bootstraps, weights=1-errors)
	class(ans) <- "GAMens"
	ans	
}
