#'@keywords internal
#'
#'@import reshape2
#'@importFrom stringr str_split
#'@importFrom splines ns
#'@importFrom stats quantile

TcGSA.dataLME<-
	function(expr, design, subject_name="Patient_ID", time_name="TimePoint", 
					 covariates_fixed="", time_covariates="",
					 group_name="", time_func="linear"){
		
		data_temp <- cbind.data.frame(design, expr)
		
		ID_vars <- c(subject_name, time_name, covariates_fixed, time_covariates, group_name)
		if(!(time_func %in% c("linear", "cubic", "splines")) && !(time_func %in% ID_vars)){
			time_func_vars <- gsub(" ", "", unlist(lapply(unlist(lapply(unlist(str_split(time_func, "\\+")), FUN=str_split, pattern="\\/")), FUN=str_split, pattern="\\*")))
			ID_vars <- c(ID_vars, time_func_vars)
		}
		if(length(which(ID_vars==""))>0){ID_vars <- ID_vars[-which(ID_vars=="")]}
		
		data_lm <- reshape2::melt(data_temp, id.vars=ID_vars, variable.name ="probe", value.name="expression", measure.vars=colnames(expr))
		data_lm$t1 <- data_lm[, time_name]
		data_lm$t1 <- data_lm$t1/10 # fixed effect estimations reduction in order to better estimate the variances (that are otherwise too small in regards of fixed effects)
		data_lm$t2 <- (data_lm$t1)^2
		data_lm$t3 <- (data_lm$t1)^3
		
		nk <- ceiling(length(unique(design[,time_name]))/4)
		if(length(unique(design[,time_name]))==2){
			noeuds <- min(design[,time_name])
			cat("Only 2 time-points here:\n longitudinal analysis is probably not the best idea...\n")
		}else{
			noeuds <- stats::quantile(design[,time_name], probs=c(1:(nk))/(nk+1))
		}
		NCsplines <- as.data.frame(ns(design[,time_name], knots = noeuds, Boundary.knots = range(design[,time_name]), intercept = FALSE))
		colnames(NCsplines) <- paste("spline_t",colnames(NCsplines) , sep="")
		NCsplines <- NCsplines*10
		data_lm <- cbind.data.frame(data_lm, NCsplines)
		return(data_lm)
	}
