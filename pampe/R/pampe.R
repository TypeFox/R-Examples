pampe <-
function(time.pretr, time.tr, treated, controls = "All", data, nbest=1, nvmax=length(controls), select="AICc", placebos=FALSE){
    
    intercept=TRUE
    
    #Roadmap 2 to relax this restriction.   
    if (length(treated) != 1){
      stop("You can only have one treated unit.")
    } 
    
    #Avoid the overidentification problem
    if (length(c(time.pretr, time.tr)) > nrow(data)){
      stop("Treatment period has exceeded the data.")
    }  
    
    if (length(which(is.na(data)))>0){
      stop("There are missing values in the data. Pampe does not allow for unbalanced data at this point.")
    }
    
    #Choose the treated unit. This if is used to deal with the treated with numeric class.
    if (class(treated) == "numeric"){
      treated <- colnames(data)[treated]
    }
    
    if(is.null(rownames(data))){
      rownames(data) <- 1:nrow(data)
    }
    
    if(is.null(colnames(data))){
      colnames(data) <- 1:ncol(data)
    }
    
    #When user chooses not to use all the controls in the table, choose the specific one from the table.

      if(controls[1]!="All"){
        controls <- which(colnames(data) %in% controls)
      }
    

    #When user uses all the controls, then submit all the variables in the table, exluding the treated one.
    if (controls[1]=="All"){
      controls <- -which(colnames(data)==treated)
      #When the pretreatment duration is not enough, we just use as many control variables as possible so as not to violate the degree of freedom.
      if (nvmax==length(controls)){
        if (ncol(data[,controls])+3>=length(time.pretr)){
          nvmax <- length(time.pretr)-4
        }
        #As more control variable always means more significant, when the pretreatment duration is enough, we just set the more control variables the better.
        if (ncol(data[,controls])+3<length(time.pretr)){
          nvmax <- ncol(data[,controls])
        }
      }
    }
    
    
    
    if (length(controls)!=1){
      if (class(controls)=="character"){
        controls <- which(colnames(data) %in% controls)
      }
    }
    
    
    
    if (class(time.pretr)=="character"){
      #   While it could be boring to let user type in all of the time pretr, especially when the durition is long.
      #   How about we creat a function that let them type in like c(1993Q1,2003Q4) and transplant them into numerical 
      #   class? The user could use time.pretr <- c("1993Q1","2003Q4");time.tr <- c("2004Q1","2008Q1"), and for example:
      
      #   TimeCharToNum<-function(time,data){
      #     time<-which(rownames(data)==time[1]):which(rownames(data)==time[2]);
      #     return(time);
      #   };
      #   time.pretr<-TimeCharToNum(time=time.pretr,data=growth)
    }
    
    
    
    
    
    if (class(time.tr)=="character"){
      #   TimeCharToNum<-function(time,data){
      #     time<-which(rownames(data)==time[1]):which(rownames(data)==time[2]);
      #     return(time);
      #   };
      #   time.tr<-TimeCharToNum(time=time.tr,data=growth)
    }
    
    
    if (nvmax+3>=length(time.pretr)){
      stop("You have selected too many controls. Please change the 'nvmax' setting.")
    }
    
    treated.original <- treated
    possible.ctrls.all <- colnames(data)[controls]
    
    select.models   <- invisible(regsubsets(x=data[time.pretr, controls], y=data[time.pretr, c(treated)], names=colnames(data[,controls]), intercept=intercept, nbest=nbest, nvmax=nvmax, matrix.logical=TRUE, really.big=TRUE))
    #Calculate AIC, AICc, BIC, ... for the models computed in the previous step. Note that x and y have to match.
    select.criteria <- invisible(regsubsets2aic(x=data[time.pretr, controls], y=data[time.pretr, c(treated)], z=select.models))
    #From now on we keep the model selected according to the AICc criteria, but it's easy to modify it to keep the AIC one
    #Retrieve list of controls as a vector from the results obtained above
    controls <- unlist(strsplit(select.criteria[which.min(select.criteria[,select]),'model'], split=" + ", fixed=TRUE))
    #Set model, compute OLS
    
    fmla <- paste(paste("`", treated, "`", sep=""), " ~ ", paste(paste("`", controls, "`", sep=""), collapse= "+"))
    ols <- invisible(lm(fmla, data=data[time.pretr,]))
    
    ## TR EFFECT RESULTS
    results <- (cbind(matrix(data=1,ncol=1,nrow=nrow(data)), as.matrix(subset(data, select=controls))) %*% as.matrix(ols$coefficients))[c(time.pretr, time.tr),]
    tr.effect <- matrix(NA, nrow=length(c(time.pretr, time.tr)), ncol=1)
    tr.effect[,1] <- as.matrix(data[,treated])[c(time.pretr, time.tr),]-results
    results <- cbind(as.matrix(data[,treated])[c(time.pretr, time.tr),], results)
    colnames(results) <- c("Actual", "Counterfactual")
    
    
    
    ## MSPE-Mead Square Policy Effect?
    mspe <- matrix(NA, ncol=1, nrow=1)
    mspe[1,1] <- sum((tr.effect[(time.pretr-min(time.pretr)+1),1])^2)  #mspe for treated unit
    controls.estim <- controls
    
    tr.effect1 <- tr.effect
    mspe1 <- mspe
    
    
    #Check if user wants placebos, if not, show estim results only
    
    ###PLACEBOS
    if (("All" %in% placebos)==TRUE | ("Both" %in% placebos)==TRUE){
      placebos <- c("controls", "time")}
    
    #Control placebo study
    if (("controls" %in% placebos)==TRUE) {
      #First we orginaize all treatment effect into one data frame.
      tr.effect <- cbind(tr.effect1, matrix(NA, nrow=length(c(time.pretr, time.tr)), ncol=(length(possible.ctrls.all))))
      colnames(tr.effect) <- c(treated, possible.ctrls.all)
      rownames(tr.effect) <- rownames(data)[c(time.pretr, time.tr)]
      
      ## MSPE
      mspe <- cbind(mspe1, matrix(NA, ncol=(length(possible.ctrls.all)), nrow=1))
      colnames(mspe) <- c(treated, possible.ctrls.all)
      rownames(mspe) <- 'mspe'
      
      #ols.placebo <- list()
      placebo.r2 <- vector()
      
      for(i in 1:(length(possible.ctrls.all))) {
        treated <- possible.ctrls.all[i]
        possible.ctrls <- possible.ctrls.all[-i]
        ####### Same program as for Israel but in the loop. We'll only store the data necessary for the plot and the OLS results of each model
        select.models   <- regsubsets(x=data[time.pretr, possible.ctrls], y=data[time.pretr, treated], names=colnames(data[,possible.ctrls]), intercept=intercept, nbest=nbest, nvmax=nvmax, matrix.logical=TRUE, really.big=TRUE)
        
        #Then we calculate those models' AIC and AICc and choose the best one based on our criterion
        select.criteria <- regsubsets2aic(x=data[time.pretr, possible.ctrls], y=data[time.pretr, treated], z=select.models) 
        
        controls <- unlist(strsplit(select.criteria[which.min(select.criteria[,select]),'model'], split=" + ", fixed=TRUE))
        #Set model, compute OLS
        fmla <- paste(paste("`", treated, "`", sep=""), " ~ ", paste(paste("`", controls, "`", sep=""), collapse= "+"))
        #ols.placebo[[i]] <- invisible(lm(fmla, data=data[time.pretr,]))
        temp.placebo <- lm(fmla, data=data[time.pretr,])
        placebo.r2[i] <- summary(temp.placebo)$r.squared
        tr.effect[,(i+1)] <- as.matrix(data[,treated])[c(time.pretr, time.tr),]-(cbind(matrix(data=1,ncol=1,nrow=nrow(data)), as.matrix(subset(data, select=controls))) %*% as.matrix(temp.placebo$coefficients))[c(time.pretr, time.tr),]
      }#End control placebo loop
      
      
      for(i in 1:length(possible.ctrls.all)) { #mspe for ctrls
        mspe[1,(i+1)] <- sum(tr.effect[(time.pretr-min(time.pretr)+1),i+1]^2)
      }
      
      tr.effect.ctrl <- tr.effect
      mspe.ctrl <- mspe
      
    } #End control placebo study
    
    #Time-treatment reassign placebo study
    if (("time" %in% placebos)==TRUE) {
      time.pretr.original <- time.pretr
      time.tr.original <- time.tr
      nvmax.user <- nvmax
      placebo.r2 <- vector()
      time.reassign <- ceiling(mean(time.pretr)):(max(time.pretr)-1)
      possible.ctrls <- possible.ctrls.all
      treated <- treated.original
      
      tr.effect <- cbind(tr.effect1, matrix(NA, nrow=length(c(time.pretr, time.tr)), ncol=(length(time.reassign))))
      colnames(tr.effect) <- c(rownames(data)[time.tr[1]], rownames(data)[time.reassign+1])
      rownames(tr.effect) <- rownames(data)[c(time.pretr, time.tr)]
      
      ## MSPE
      mspe <- cbind(mspe1, matrix(NA, ncol=(length(time.reassign)), nrow=1))
      colnames(mspe) <- c(rownames(data)[time.tr[1]], rownames(data)[time.reassign+1])
      rownames(mspe) <- 'mspe'
      
      for(i in 1:length(time.reassign)) {
        time.pretr <- time.pretr[1]:time.reassign[i]
        time.tr <- (time.reassign[i]+1):time.tr[length(time.tr)]
        
        if (length(possible.ctrls)+3>=length(time.pretr)){
          nvmax = length(time.pretr) - 3
        }
        ####### Same program as for treated but in the loop. We'll only store the data necessary for the plot and the OLS results of each model
        select.models   <- regsubsets(x=data[time.pretr, possible.ctrls], y=data[time.pretr, treated], names=colnames(data[,possible.ctrls]), intercept=intercept, nbest=nbest, nvmax=nvmax, matrix.logical=TRUE, really.big=TRUE)
        
        #Then we calculate those models' AIC and AICc and choose the best one based on our criterion
        select.criteria <- regsubsets2aic(x=data[time.pretr, possible.ctrls], y=data[time.pretr, treated], z=select.models) 
        
        controls <- unlist(strsplit(select.criteria[which.min(select.criteria[,select]),'model'], split=" + ", fixed=TRUE))
        #Set model, compute OLS
        fmla <- paste(paste("`", treated, "`", sep=""), " ~ ", paste(paste("`", controls, "`", sep=""), collapse= "+"))
        #ols.placebo[[i]] <- invisible(lm(fmla, data=data[time.pretr,]))
        temp.placebo <- lm(fmla, data=data[time.pretr,])
        placebo.r2[i] <- summary(temp.placebo)$r.squared
        tr.effect[,(i+1)] <- as.matrix(data[,treated])[c(time.pretr, time.tr),]-(cbind(matrix(data=1,ncol=1,nrow=nrow(data)), as.matrix(subset(data, select=controls))) %*% as.matrix(temp.placebo$coefficients))[c(time.pretr, time.tr),]
        #Reset nvmax
        nvmax = nvmax.user
        time.tr <- time.tr.original
        time.pretr <- time.pretr.original
      }#End control placebo loop
      
      
      for(i in 1:length(time.reassign)) { #mspe for ctrls
        mspe[1,(i+1)] <- sum(tr.effect[(time.pretr-min(time.pretr)+1),i+1]^2)
      }
      
      
      tr.effect.time <- tr.effect
      mspe.time <- mspe
      
    } #End time placebo study
    
    #data <- as.matrix(data)
    

    if (FALSE %in% placebos){
    
      result <- list(controls=controls, model=ols, counterfactual=results, data=data)
      class(result) <- "pampe"
      return(result) }
    
    if (("controls" %in% placebos)==TRUE&length(placebos)==1) {
      result <-list(controls=controls.estim, model=ols, counterfactual=results,placebo.ctrl=list(mspe=mspe.ctrl, tr.effect=tr.effect.ctrl), data=data)
      class(result) <- "pampe"
      return(result)
    }
    
    if (("time" %in% placebos)==TRUE&length(placebos)==1) {
      result <- list(controls=controls.estim, model=ols, counterfactual=results,placebo.time=list(mspe=mspe.time, tr.effect=tr.effect.time), data=data)
      class(result) <- "pampe"
      return(result)
    }
    
    if (length(placebos)==2) {
      result <- list(controls=controls.estim, model=ols, counterfactual=results,placebo.ctrl=list(mspe=mspe.ctrl, tr.effect=tr.effect.ctrl), placebo.time=list(mspe=mspe.time, tr.effect=tr.effect.time), data=data)
      class(result) <- "pampe"
      return(result)}
    
  }
