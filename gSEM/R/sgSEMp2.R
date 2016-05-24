##' This function builds an gSEM model using gSEM principle 2. Principle 2 resembles the multiple regression principle in the way multiple predictors are considered simultaneously. Specifically, the first-level predictors to the system level variable, such as, Time and unit level variables, acted on the system level variable collectively by an additive model. This collective additive model can be found with a generalized stepwise variable selection (using the step() function in R, which performs variable selection on the basis of AIC) and this proceeds iteratively.
##' 
##' Data is analysed first using Principle 1 to find the best models. If needed, transformations based on the best models are applied to the predictors. Starting from the system response variable, each variable is regressed on all other variables except for the system response in an additive multiple regression model, which is reduced by a stepwise selection using stepAIC(). Then, for each selected variable, fitted regression for 6 selected functional forms and pick the best.
##'
##' @title Semi-supervised Generalized Structural Equation Modelling (gSEM) - Principle 2

##' @param x A dataframe, requiring at least 2 columns. By default its first column stores the main or primary influencing predictor, or exogenous variable e.g.., time, or a main predictor, the second column stores the response variable, and other columns store intermediate variables.
##' @param predictor A character string of the column name of the system predictor OR a numeric number indexing the column of the main predictor.
##' @param response A character string of the column name of the main response OR a numeric number indexing the column of the system response.
##' @return A list of the following items:
##'
##' \itemize{
##' \item "Graph": A network graph that contains the group and individual relationships between response and predictors determined by principle 2.
##' \item "res.print": A matrix. For each row, first column is the response variable, second column is the predictor, the other columns show corresponding summary information.
##'}

##' @seealso sgSEMp1() and plot.sgSEMp2()
##'
##' @export
##' @importFrom 'stats' 'coef'
##' @examples
##' # Using built-in dataset
##' data(acrylic)
##' ans <- sgSEMp2(acrylic)
##' ans$res.print
##' plot(ans)
##' 
##' \dontrun{
##' # Using simulated data
##' x4=runif(100,0,2)
##' x3=1+2.5*x4+rnorm(100,0,0.5)
##' x1=runif(100,1,4)
##' x2=-1-x1+x3+rnorm(100,0,0.3)
##' y=2+2*exp(x1/3)+(x2-1)^2-x3+rnorm(100,0,0.5)
##' # Check the pairwise plot 
##' sim=cbind(x4,y,x1,x2,x3)
##' pairs(sim)
##' ans <- sgSEMp2(as.data.frame(sim))
##' plot(ans)
##' }

sgSEMp2 <- function(x, predictor = NULL, response = NULL){
  
	if(!missing(predictor)){
		if(!is.character(predictor)){
			if(is.wholenumber(predictor)){
				if(predictor < 1 | predictor > length(colnames(x)))
					stop("Predictor location out of range!")
				predictor.loc <- predictor
			}
		}else{ 
			if(!(predictor %in% colnames(x))){
				stop(paste0("Predictor '", predictor, "' does not exist!"))
			}
			predictor.loc <- which(colnames(x) == predictor)
		}
		neworder <- 1:length(colnames(x))
		neworder[1] <- predictor.loc
		neworder[predictor.loc] <- 1
		x <- x[neworder]
	}
	if(!missing(response)){
		response.loc <- which(colnames(x) == response)
		if(!is.character(response)){
			if(is.wholenumber(response)){
				if(response < 1 | response > length(colnames(x)))
					stop("Response location out of range!")
				response.loc <- response
			}
		}else{
			if(!(response %in% colnames(x))){
				stop(paste0("Response '", response, "' does not exist!"))
			}
			response.loc <- which(colnames(x) == response)
		}
		neworder <- 1:length(colnames(x))
		neworder[2] <- response.loc
		neworder[response.loc] <- 2
		x <- x[neworder]
	}

	###############################
	## The following function does multiple selection
	############################### 
	Multiple.relation <- function(Mx){

		status_ind <- status_matrix[i_status,]
		Resp <- colnames(Mx)[which(status_ind==1)]   
		Var <- colnames(Mx)[which(status_ind==0)]
		## get rid of NA's
		x1 <- Mx[c(Resp,Var)]
		x1 <- x1[apply(x1, 1, FUN = function(x) sum(is.na(x)) == 0),]
		
		################## START: APPLY TRANSFORMATION TO PREDICTORS ###################
		trans <- rep(NA,length=length(Var))
		trans1 <- rep(NA,length=length(Var))
		names(trans) <- Var
		for ( i in 1:length(Var) ) {
			# If "SL", x ---> x
			if ( p1.bestModels[Var[i],Resp] == "SL" ) {
				trans1[i] <- NA
				names(trans1)[i] <- Var[i]
				Quad_value <- coef(p1.allModels[Var[i],Resp,'SL'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*",Var[i],sep="")
				}
			}
			# If "SQuad", x ---> x^2
			if ( p1.bestModels[Var[i],Resp] == "SQuad" ) {
				x1[,Var[i]] <- (x1[,Var[i]])^2
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'SQuad'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*",Var[i],"^2",sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*",Var[i],"^2",sep="")
				}
			}
			# If "Exp", x ---> exp(x)
			if ( p1.bestModels[Var[i],Resp] == "Exp" ) {
				x1[,Var[i]] <- exp(x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=e^", Var[i],sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Exp'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*e^",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*e^",Var[i],sep="")
				}
			}
			# If "Log", x ---> log(x)
			if ( p1.bestModels[Var[i],Resp] == "Log" ) {
				x1[,Var[i]] <- log(x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=log$", Var[i],"$",sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Log'][[1]])
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*log_",Var[i],sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*log_",Var[i],sep="")
				}
			}
			# If "Quad", x ---> (x+b/2c)^2
			if ( p1.bestModels[Var[i],Resp] == "Quad" ) {
			    temp_column <- x1[,Var[i]]^2
				x1 <- cbind(x1,temp_column)
				names(x1)[ncol(x1)] <- paste(Var[i],"__2",sep="")
			
				Quad_value <- coef(p1.allModels[Var[i],Resp,'Quad'][[1]])
				names(trans1)[i] <- Var[i]
				if ( (Quad_value[2]>0)&(Quad_value[3]>0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),'*',Var[i],'+',format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
				}
				if ( (Quad_value[2]>0)&(Quad_value[3]<0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),'*',Var[i],format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
				}
				if ( (Quad_value[2]<0)&(Quad_value[3]>0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),'*',Var[i],'+',format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
				}
				if ( (Quad_value[2]<0)&(Quad_value[3]<0) ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),'*',Var[i],format(Quad_value[3],scitific=T,digits=2),'*',Var[i],"^2",sep="")
				}
			}
			# If "nls", x ---> exp(cx)
			if ( p1.bestModels[Var[i],Resp] == "nls" ) {
				Quad_value <- coef(p1.allModels[Var[i],Resp,'nls'][[1]])
				x1[,Var[i]] <- exp(Quad_value[3]*x1[,Var[i]])
				trans1[i] <- paste(Var[i],"_t=exp$", format(Quad_value[3],scitific=T,digits=2),"*",Var[i],"$",sep="")
				names(trans1)[i] <- Var[i]
				names(x1)[i+1] <- paste(names(x1)[i+1],"_t",sep="")
				if ( Quad_value[2]>0 ) {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),"+",format(Quad_value[2],scitific=T,digits=2),"*e^",format(Quad_value[3],scitific=T,digits=2),Var[i],sep="")
					trans1[i] <- paste(Var[i],"_t=exp$", format(Quad_value[3],scitific=T,digits=2),"*",Var[i],"$",sep="")
				} else {
					trans[i] <- paste(Resp,"=",format(Quad_value[1],scitific=T,digits=2),format(Quad_value[2],scitific=T,digits=2),"*e^",format(Quad_value[3],scitific=T,digits=2),Var[i],sep="")
				}
			}
		}
		################## END: APPLY TRANSFORMATION TO PREDICTORS ###################
		Var_char <- paste(names(x1)[-1], collapse = "+")
		Rel <- paste0(Resp, "~", Var_char)
		lm.full.model <- do.call("lm", list(Rel, data=as.name("x1")))
		lm.best <- stepAIC(lm.full.model, direction = "backward", trace =FALSE)
		R2 <- summary(lm.best)$r.squared
		adj.R2 <- summary(lm.best)$adj.r.squared
		best.vars <- names(coef(lm.best))[-1] 
		
		coeff_number <- format(lm.best$coeff,digits=2)
		names(coeff_number) <- NULL
		coeff_names <- names(lm.best$coeff)
		lm.best.text <- paste(Resp,"=",sep="")
		for ( ii in 1:length(coeff_names) ) {
			if (as.numeric(gsub(" ","",coeff_number[ii]))<0) {
			    if (coeff_names[ii]=='(Intercept)') {
					lm.best.text <- paste(lm.best.text,coeff_number[ii],sep="")
				} else {
					if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))==0) {
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',coeff_names[ii],sep="")
					}
					if (length(grep("_t",coeff_names[ii]))!=0) {
						names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1)
						names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp2,sep="")
					}
					if (length(grep("__2",coeff_names[ii]))!=0) {
						names_temp <- gsub("__","^",coeff_names[ii])
						lm.best.text <- paste(lm.best.text,coeff_number[ii],'*',names_temp,sep="")
					}
				}
			} else {
				if (coeff_names[ii]=='(Intercept)') {
					lm.best.text <- paste(lm.best.text,gsub(" ","",coeff_number[ii]),sep="")
				} else {
					if (length(grep("_t",coeff_names[ii]))+length(grep("__2",coeff_names[ii]))==0) {
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',coeff_names[ii],sep="")
					}
					if (length(grep("_t",coeff_names[ii]))!=0) {
						names_temp <- substr(coeff_names[ii],1,regexpr("_t",coeff_names[ii])[1]-1)
						names_temp2 <- substr(trans1[names_temp],regexpr("=",trans1[names_temp])[1]+1,nchar(trans1[names_temp]))
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp2,sep="")
					}
					if (length(grep("__2",coeff_names[ii]))!=0) {
						names_temp <- gsub("__","^",coeff_names[ii])
						lm.best.text <- paste(lm.best.text,'+',gsub(" ","",coeff_number[ii]),'*',names_temp,sep="")
					}
				}
			}
		}
		
		# Make a new row in the Res.print table for each of the predictor variables in the additive model
		Res.print.local <- matrix(NA, nrow = 1, ncol = nRes)
		colnames(Res.print.local) <- c("resp",          "var",          "ar2",      "Tran",     "GModel", 
	                                   "Gpvalue",       "GR2aR2",       "IModel",   "Ipvalue",  "IR2aR2",
							           "Rank",          "Gcoeff",       "GModelB" )
		for (n in 1:length(best.vars)){
			Res.print.newrow <- matrix(NA, nrow = 1, ncol = nRes)
			Res.print.newrow <- as.data.frame(Res.print.newrow)
			colnames(Res.print.newrow) <- c("resp",          "var",          "ar2",      "Tran",     "GModel", 
	                                        "Gpvalue",       "GR2aR2",       "IModel",   "Ipvalue",  "IR2aR2",
							                "Rank",          "Gcoeff",       "GModelB" )
			Res.print.newrow[1,1] <- Resp
			Res.print.newrow[1,3] <- format(adj.R2,scitific=T,digits=2)		
			Res.print.newrow[1,6] <- paste(round(summary(lm.best)$coeff[,4],digits=3),collapse="  ")		
			Res.print.newrow[1,7] <- paste(round(R2,digits=3),round(adj.R2,digits=3),sep="  ")			
			Res.print.newrow[1,12] <- paste(format(lm.best$coeff,digits=2),collapse=" ")
			best.vars.rep <- best.vars
			for (iitemp in 1:length(best.vars)) {
				if (regexpr("_t",best.vars[iitemp])!=-1) {
					temp_str <- trans1[substr(best.vars[iitemp],1,regexpr("_t",best.vars[iitemp])-1)]
					best.vars.rep[iitemp] <- substr(temp_str,regexpr("=",temp_str)+1,nchar(temp_str))
				}
			}
			Res.print.newrow[1,13] <- paste(Resp,"=",gsub("__2","^2",paste(best.vars.rep,collapse="+")),sep="")
			if ((length(grep("__2",best.vars[n]))==0)&(length(grep("_t",best.vars[n]))==0)) {
				Res.print.newrow[1,2]<- best.vars[n]
				Res.print.newrow[1,5] <- lm.best.text
				Res.print.newrow[1,8] <- trans[best.vars[n]]
				if (p1.bestModels[best.vars[n],Resp]!='-1') {
					Individual_summary <- summary(p1.allModels[best.vars[n],Resp,p1.bestModels[best.vars[n],Resp]][[1]])
					Res.print.newrow[1,9] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
					if (p1.bestModels[best.vars[n],Resp]!="nls") {
						Res.print.newrow[1,10] <- paste(round(Individual_summary$r.squared,digits=3),round(Individual_summary$adj.r.squared,digits=3),sep="  ")
					}
				}
				Res.print.newrow[1,11] <- which(order(summary(lm.best)$coeff[-1,4])==n)
				Res.print <<- rbind(Res.print, Res.print.newrow)
				Res.print.local <- rbind(Res.print.local,Res.print.newrow)
			} else {
				if (length(grep("_t",best.vars[n]))!=0) {
					new_string <- substr(best.vars[n],1,regexpr("_t",best.vars[n])[1]-1)
				} else {
					new_string <- substr(best.vars[n],1,regexpr("__2",best.vars[n])[1]-1)
				}
				Res.print.newrow[1,2]<- new_string
				Res.print.newrow[1,5] <- lm.best.text
				Res.print.newrow[1,8] <- trans[new_string]
				Individual_summary <- summary(p1.allModels[new_string,Resp,p1.bestModels[new_string,Resp]][[1]])
				Res.print.newrow[1,9] <- paste(round(Individual_summary$coeff[,4],digits=3),collapse="  ")
				if (p1.bestModels[new_string,Resp]!="nls") {
					Res.print.newrow[1,10] <- paste(round(Individual_summary$r.squared,digits=3),round(Individual_summary$adj.r.squared,digits=3),sep="  ")
				}
				Res.print.newrow[1,11] <- which(order(summary(lm.best)$coeff[-1,4])==n)
				insert_ind <- 0
				Res.temp <- Res.print[-1,]
				if (nrow(Res.temp)>0) {
					for ( mm in 1:nrow(Res.temp) ) {
						if ( Res.temp[mm,1]==Resp ) {
							if ( Res.temp[mm,2]==new_string ) {
								insert_ind <- 1
							}
						}
					}
				}
				if (insert_ind==0) {
					Res.print <<- rbind(Res.print, Res.print.newrow)
					Res.print.local <- rbind(Res.print.local,Res.print.newrow)
				}
			}	
		}
		sig_var <- Res.print.local[-1,2]
		status_matrix[i_status,Resp] <<- 2
		status_matrix0 <- matrix(NA,nrow=1,ncol=dim(status_matrix)[2])
		colnames(status_matrix0) <- colnames(status_matrix)
		for ( sig_ind in 1:length(sig_var) ) {
			if ( (sig_var[sig_ind]!=colnames(status_matrix)[1])&(max(status_matrix[,sig_var[sig_ind]])==0) ) {
				status_temp <- matrix(NA,nrow=1,ncol=dim(status_matrix)[2])
				colnames(status_temp) <- colnames(status_matrix)
				status_temp[1,] <- status_matrix[i_status,]
				status_temp[1,sig_var[sig_ind]] <- 1
				new_sig <- sig_var[-sig_ind]
				if (length(new_sig)!=0) {
					for ( i_new in 1:length(new_sig) ) {
						if (max(status_matrix[,new_sig[i_new]])!=0) {
							status_temp[1,new_sig[i_new]] <- 2
						}
					}
				}
				status_temp[1,1] <- 0
				status_matrix0 <- rbind(status_matrix0,status_temp)
			}
		}
		status_matrix <<- rbind(status_matrix,status_matrix0[-1,])
	}
	
	###############################
	## Main scripts; Above two functions are called
	############################### 
  
	# Apply principle 1 on the data to find the best models
	if ( !is.null(predictor) ) {
		predictor.p1 <- predictor
	} else {
		predictor.p1 <- 1
	}
	if ( !is.null(response) ) {
		response.p1 <- response
	} else {
		response.p1 <- 2
	}
	p1.result <- sgSEMp1(x, predictor = predictor.p1, response = response.p1)
	p1.bestModels <- rbind(p1.result$bestModels[1,],rep(NA,length(colnames(p1.result$bestModels))),p1.result$bestModels[-1,])
	p1.bestModels <- cbind(rep(NA,length(rownames(p1.bestModels))),p1.bestModels)
	rownames(p1.bestModels) <- c(rownames(p1.result$bestModels)[1],colnames(p1.result$bestModels)[1],rownames(p1.result$bestModels)[-1])
	colnames(p1.bestModels) <- rownames(p1.bestModels)
	
	p1.allModels <- vector( "list", length(colnames(p1.bestModels)) * length(colnames(p1.bestModels)) * length(names(p1.result$allModels[1,1,])) )
	dim(p1.allModels) <- c(length(colnames(p1.bestModels)),length(colnames(p1.bestModels)), length(names(p1.result$allModels[1,1,])))
    p1.allModels[,,] <- NA
    dimnames(p1.allModels) <- list(colnames(p1.bestModels), colnames(p1.bestModels), names(p1.result$allModels[1,1,]))
	
	temp <- rbind(p1.result$allModels[1,,"SL"],rep(NA,length(colnames(p1.result$allModels[,,"SL"]))),p1.result$allModels[-1,,"SL"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"SL"] <- temp

	temp <- rbind(p1.result$allModels[1,,"Quad"],rep(NA,length(colnames(p1.result$allModels[,,"Quad"]))),p1.result$allModels[-1,,"Quad"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Quad"] <- temp
	
	temp <- rbind(p1.result$allModels[1,,"SQuad"],rep(NA,length(colnames(p1.result$allModels[,,"SQuad"]))),p1.result$allModels[-1,,"SQuad"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"SQuad"] <- temp
	
	temp <- rbind(p1.result$allModels[1,,"Exp"],rep(NA,length(colnames(p1.result$allModels[,,"Exp"]))),p1.result$allModels[-1,,"Exp"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Exp"] <- temp

	temp <- rbind(p1.result$allModels[1,,"Log"],rep(NA,length(colnames(p1.result$allModels[,,"Log"]))),p1.result$allModels[-1,,"Log"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"Log"] <- temp

	temp <- rbind(p1.result$allModels[1,,"nls"],rep(NA,length(colnames(p1.result$allModels[,,"nls"]))),p1.result$allModels[-1,,"nls"])
	temp <- cbind(rep(NA,length(rownames(temp))),temp)
	p1.allModels[,,"nls"] <- temp

	rownames(p1.allModels) <- colnames(p1.bestModels)
	# End apply principle 1
  
	nVar <- ncol(x)    #total number of variables
	nRVar <- nVar - 1  #total number of possible response variables
  
	nRes <- 13  # Number of cells in print variable; see below
	Res.print <- matrix(NA, nrow = 1, ncol = nRes)
	colnames(Res.print) <- c("resp",          "var",          "ar2",      "Tran",     "GModel", 
	                         "Gpvalue",       "GR2aR2",       "IModel",   "Ipvalue",  "IR2aR2",
							 "Rank",          "Gcoeff",       "GModelB" )
	status_matrix <- matrix(0, nrow=1, ncol=ncol(x))
	colnames(status_matrix) <- colnames(x)
	status_matrix[1,2] <- 1
	i_status <- 1
	while (i_status<=nrow(status_matrix)) {
		Multiple.relation(x)
		i_status <- i_status + 1
	}

	Res.print <- Res.print[-1,]
	res <- list(res.print = Res.print) #outputs are three tables
	class(res) <- c("sgSEMp2","list")
	invisible(res)
}