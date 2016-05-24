##################################################
##################################################
##												##									
##	GLM / GAM / GBM / ME / RF MODEL FIT LOOP	##
##												##									
##################################################
##################################################

## v1.0 C.RANDIN - UNIBAS - March 2013

###########################################################################################################################################################

# IMPORTING LIBRARY

#library(rms)  
#library(gbm)

#library(gam)
#library(rJava)
#jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
#Sys.setenv(NOAWT=TRUE)
#library(dismo)
#library(randomForest)

# Input dataset format
# Columns 2-3: long and lat geographic coordinates 
# Columns 4-8: Predicting variables
# Columns 9-10: Species presences / absences

ecospat.env <- new.env ()
globalVariables ("ecospat.testData")


######
# Auxiliary function for examples.
######

glm.gam.gbm.me.rf.fit.model  <- function(data,p.pred,p.resp,me.args=c('-J','-P','threshold=false','product=false','quadratic=false','linear=false'),
                                         me.path= outpath,outpath)
{
  
  ###############
  ## SCOPE.GAM ##
  ###############
  
  scope.gam <- function(data.pred,smoother= "s", degree="4")
  {
    scope.list.gam <- as.list(data.pred)
    names(scope.list.gam) <- colnames(data.pred)
    for(i in 1:ncol(data.pred)) 
    {
      if (is.factor(data.pred[[i]])==T)
      {
        junk <- c("1")
        junk <- c(junk, paste(colnames(data.pred[i]), sep = ""))
        junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
        scope.list.gam[[i]] <- junk
      } else if (is.factor(data.pred[[i]])==F)	{
        junk <- c("1")
        junk <- c(junk, paste(smoother,"(", colnames(data.pred[i]),",",degree,")", sep = ""))
        junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
        scope.list.gam[[i]] <- junk	
        
      }	
    }
    return(scope.list.gam)
  }
  
  
  
  sp.names <- colnames(data[p.resp:ncol(data)])
  
  n.pred <- p.resp-p.pred
  rownames(data) <- 1:nrow(data)
  
  # GLM FORMULA
  for (k in (p.pred):(p.resp-1)) 
  { 
    if (k == p.pred && is.factor(data[,k])==FALSE) 
    {
      fmula.pred.glm <- as.vector(paste("pol(",names(data[k]),",2)",sep=""),mode="character")	
    }
    else if (k == p.pred && is.factor(data[,k])==TRUE) 
    {	
      fmula.pred.glm <- as.vector(paste(names(data[k]),sep=""),mode="character")		
    }
    else if (is.factor(data[,k])==FALSE) 
    {
      fmula.pred.glm <- paste(fmula.pred.glm," + pol(",names(data[k]),",2)",sep="")
      
    }
    else if (is.factor(data[,k])==TRUE) 
    {
      
      fmula.pred.glm <- paste(fmula.pred.glm," + ",names(data[k]),sep="")	
    }	
  }		
  
  # GAM FORMULA	
  for (k in (p.pred):(p.resp-1)) 
  { 
    if (k == p.pred && is.factor(data[,k])==FALSE) 
    {
      fmula.pred.gam <- as.vector(paste("s(", names(data[k]),", 4)",sep=""),mode="character")	
    }
    else if (k == p.pred && is.factor(data[,k])==TRUE) 
    {	
      fmula.pred.gam <- as.vector(paste(names(data[k]),sep=""),mode="character")		
    }
    else if (is.factor(data[,k])==FALSE) 
    {
      fmula.pred.gam <- paste(fmula.pred.gam," + s(", names(data[k]),", 4)",sep="")
      
    }
    else if (is.factor(data[,k])==TRUE) 
    {
      
      fmula.pred.gam <- paste(fmula.pred.gam," + ",names(data[k]),sep="")	
    }	
  }		
  
  gam.scope.mod <- scope.gam(data[,p.pred:(p.resp-1)],smoother= "s", degree="4")
  
  # GBM-BOOSTING FORMULA	
  for (k in (p.pred):(p.resp-1)) 
  { 
    if (k == p.pred) 
    {	
      fmula.pred.gbm <- as.vector(paste(names(data[k]),sep=""),mode="character")		
    } else 
    {
      
      fmula.pred.gbm <- paste(fmula.pred.gbm," + ",names(data[k]),sep="")	
    }	
  }		
  
  ##########################################################################################################################################################
  
  ##################################
  ## LOOP FOR SPECIES COMPUTATION ##
  ##################################
  for (i in p.resp:ncol(data)) # LOOP FOR MULTIPLE SPECIES
  {	
    
    #i = p.resp # LOOP FOR MULTIPLE SPECIES
    
    cat("Computations", (i-p.resp+1),"for species ",names(data[i]),"is starting now...", "\n",append = F)
    cat(".............", "\n",append = F)
    
    name.tmp.sp <- names(data[i])
    
    
    ########################
    ## MODELS CALIBRATION ##
    ########################
    
    #df.tmp.input <- na.omit(data[,c(i,p.pred:(p.resp-1))])
    df.tmp.input <<- na.omit(data[,c(i,p.pred:(p.resp-1))])
    
    row.names(df.tmp.input) <- 1:dim(df.tmp.input)[1]
    
    #df.tmp.input <- na.exclude(df.tmp.input)
    df.tmp.input <<- na.exclude(df.tmp.input)
    
    ###########
    # GLM FIT #
    ###########
    cat("GLM", "\n",append = F)
    cat("> calibration", "\n",append = F)
    glm.tmp.step <- step(glm(eval(parse(text = paste(name.tmp.sp, "~",fmula.pred.glm, collapse = ""))),
                             data=df.tmp.input,family=binomial,maxit = 100),trace = F,direction="both")
    
    
    #assign(paste("glm.",name.tmp.sp,sep=""),glm.tmp.step,envir = .GlobalEnv)
    assign(paste("glm.",name.tmp.sp,sep=""),glm.tmp.step,envir = ecospat.env)
    
    
    ###########
    # GAM FIT #
    ###########
    
    cat("GAM", "\n",append = F)
    cat("> calibration", "\n",append = F)	 	
    gam.tmp.step <- step.gam(gam(eval(parse(text = paste(paste(name.tmp.sp), "~ 1", collapse = ""))),
                                 data=df.tmp.input,family=binomial),trace = F,scope=gam.scope.mod,direction="both",control = gam.control(maxit = 50,bf.maxit = 50))
    
    assign(paste("gam.",name.tmp.sp,sep=""),gam.tmp.step,envir = ecospat.env)
    
    
    ###########
    # GBM FIT #
    ###########
    
    cat("GBM", "\n",append = F)
    cat("> calibration", "\n",append = F)
    
    gbm.fmula <- eval(parse(text = paste(paste(name.tmp.sp), "~",fmula.pred.gbm, collapse = "")))
    
    
    
    gbm.tmp <- gbm(gbm.fmula, data=df.tmp.input, var.monotone = rep(0,length =  n.pred), n.trees = 2000, interaction.depth = 3,
                   n.minobsinnode = 10,shrinkage = 0.01, bag.fraction = 0.5, train.fraction = 1, verbose = F, cv.folds = 10)
    
    best.itr.temp <- gbm.perf(gbm.tmp, method="cv", plot.it=F)
    
    gbm.tmp.cal <-gbm.more(gbm.tmp, n.new.trees = best.itr.temp, data=df.tmp.input)
    
    assign(paste("gbm.",name.tmp.sp,sep=""),gbm.tmp.cal,envir = ecospat.env)
    
    
    ##########
    # ME FIT #
    ##########
    
    cat("ME", "\n",append = F)
    cat("> calibration", "\n",append = F)
    
    Sys.setenv(NOAWT=TRUE)
    jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
    #library(dismo)
    
    me.tmp <-maxent(df.tmp.input[,2:ncol(df.tmp.input)],as.vector(df.tmp.input[,1]),path=me.path)
    
    assign(paste("me.",name.tmp.sp,sep=""),me.tmp,envir = ecospat.env)
    
    
    ##########
    # RF FIT #
    ##########
    
    cat("RF", "\n",append = F)
    cat("> calibration", "\n",append = F)
    
    rf.tmp <- randomForest(x=df.tmp.input[,2:ncol(df.tmp.input)], y=as.factor(df.tmp.input[,1]), ntree=1000, importance=TRUE)
    
    assign(paste("rf.",name.tmp.sp,sep=""),rf.tmp,envir = ecospat.env)
    
    #########################################################################################################################################             
    
    cat(".............", "\n",append = F)
    cat(".............", "\n",append = F)	
    
    rm(df.tmp.input)
    # END SPECIES LOOP		
  }
  # END SPECIES LOOP	
  
  # END FUNCTION
}
# END FUNCTION	



##########################
##########################
##########################
##						##
##	ADD-IN FUNCTIONS	##
##						##
##########################
##########################
##########################

########################################################################################################################################################

############
## CV-GLM ##
############

ecospat.cv.glm <- function(glm.obj, K=10, cv.lim = 10, jack.knife = F)
{
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  data.cv <- glm.obj$data
  #data.cv <<- glm.obj$data
  
  
  n <- nrow(data.cv)
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  n.0 <- nrow(data.cv[glm.obj$y==0,])
  n.1 <- nrow(data.cv[glm.obj$y>0,])
  
  
  #################################
  ## Jack-Knife Cross-validation ##
  ################################# 
  
  
  if (jack.knife  == F && n.1 < cv.lim || jack.knife  == F && n.1 < K || jack.knife  == T)
  {
    K <- nrow(data.cv)
    cat("K has been set to ",K," -> Jack-Knife Cross-validation enabled!","\n",append = F)	
    
    for (i in 1:K)
    {
      temp.glm <- update(glm.obj, data = data.cv[-i,])
      if (i == 1)
      {
        vect.id <- i
        vect.predicted <- as.double(predict(temp.glm, data.cv[i,],type="response"))
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, i, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.double(predict(temp.glm, data.cv[i,],type="response")), after=length(vect.predicted))
      }		
      
      
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(glm.obj$y),predictions=df.tmp.res[,2])
    
  } else 
  {	
    id.0 <- as.vector(row.names(data.cv[glm.obj$y==0,]),mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0))
    {
      K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0)
    {
      cat("K.0 has been set to", K.lst.0, "\n",append = F)
    }
    
    id.1 <- as.vector(row.names(data.cv[glm.obj$y>0,]), mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0))
    {
      K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1)
    {
      cat("K.1 has been set to", K.lst.1, "\n",append = F)
    }
    if (K.lst.0 != K.lst.1)
    {
      cat("P/A stratifications have not the same values!", "\n",append = F)
      min.K <- min(K.lst.0,K.lst.1)
      K.lst.0 <- min.K
      K.lst.1 <- min.K
    }
    
    cat("K has been finally set to",K.lst.0, "\n",append = F)
    K.lst <- K.lst.0
    
    # P/A STRATIFICATION
    
    # ABSENCE STRATIFICATION
    f.0 <- ceiling(n.0/K.lst.0)
    s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
    
    # PRESENCE STRATIFICATION
    f.1 <- ceiling(n.1/K.lst.1)
    s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
    
    # RESPONSE PREDICTION
    
    for (i in 1:K.lst)
    {
      j.out <- id.0[(s.0 == i)]
      j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
      j.out <- sort(j.out)
      j.in <- id.0[(s.0 != i)]
      j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
      j.in <- sort(j.in)
      glm.cal <- update(glm.obj, data = data.cv[j.in, , drop = FALSE])
      glm.val <- predict(glm.cal, data.cv[j.out, , drop = FALSE], type = "response")
      
      if (i == 1)
      {
        vect.id <- j.out
        vect.predicted <- as.vector(glm.val)
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, j.out, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
      }	
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(glm.obj$y),predictions=df.tmp.res[,2])
    
  }	
  df.res[df.res[,3]>1,3] <- 1
  df.res[df.res[,3]<0,3] <- 0
  
  return(df.res)
}

#################################################################################################################################
##################################################################################################################################

############
## CV-GAM ##
############



ecospat.cv.gam <- function(gam.obj, K=10, cv.lim = 10, jack.knife = F)
{
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  data.cv <- gam.obj$data
  
  gam.obj <- gam(eval(gam.obj$formula),data=data.cv)
  
  
  n <- nrow(data.cv)
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  n.0 <- nrow(data.cv[gam.obj$y==0,])
  n.1 <- nrow(data.cv[gam.obj$y>0,])
  
  
  
  #################################
  ## Jack-Knife Cross-validation ##
  ################################# 
  
  
  if (jack.knife  == F && n.1 < cv.lim || jack.knife  == F && n.1 < K || jack.knife  == T)
  {
    K <- nrow(data.cv)
    cat("K has been set to ",K," -> Jack-Knife Cross-validation enabled!","\n",append = F)	
    
    for (i in 1:K)
    {
      temp.gam <- update(gam.obj, data = data.cv[-i,])
      
      if (i == 1)
      {
        vect.id <- i
        vect.predicted <- as.double(predict(temp.gam, data.cv[i,],type="response"))
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, i, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.double(predict(temp.gam, data.cv[i,],type="response")), after=length(vect.predicted))
      }		
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(gam.obj$y),predictions=df.tmp.res[,2])
    
  } else 
  {	
    id.0 <- as.vector(row.names(data.cv[gam.obj$y==0,]),mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0))
    {
      K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0)
    {
      cat("K.0 has been set to", K.lst.0, "\n",append = F)
    }
    
    id.1 <- as.vector(row.names(data.cv[gam.obj$y>0,]), mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0))
    {
      K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1)
    {
      cat("K.1 has been set to", K.lst.1, "\n",append = F)
    }
    if (K.lst.0 != K.lst.1)
    {
      cat("P/A stratifications have not the same values!", "\n",append = F)
      min.K <- min(K.lst.0,K.lst.1)
      K.lst.0 <- min.K
      K.lst.1 <- min.K
    }
    
    cat("K has been finally set to",K.lst.0, "\n",append = F)
    K.lst <- K.lst.0
    
    # P/A STRATIFICATION
    
    # ABSENCE STRATIFICATION
    f.0 <- ceiling(n.0/K.lst.0)
    s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
    
    # PRESENCE STRATIFICATION
    f.1 <- ceiling(n.1/K.lst.1)
    s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
    
    # RESPONSE PREDICTION
    
    for (i in 1:K.lst)
    {
      j.out <- id.0[(s.0 == i)]
      j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
      j.out <- sort(j.out)
      j.in <- id.0[(s.0 != i)]
      j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
      j.in <- sort(j.in)
      gam.cal <- update(gam.obj, data = data.cv[j.in, , drop = FALSE])
      gam.val <- predict(gam.cal, data.cv[j.out, , drop = FALSE], type = "response")
      
      if (i == 1)
      {
        vect.id <- j.out
        vect.predicted <- as.vector(gam.val)
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, j.out, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.vector(gam.val), after=length(vect.predicted))
      }	
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(gam.obj$y),predictions=df.tmp.res[,2])
    
    
  }	
  df.res[df.res[,3]>1,3] <- 1
  df.res[df.res[,3]<0,3] <- 0
  
  return(df.res)
}

##################################################################################################################################
##################################################################################################################################


############
## CV-GBM ##
############


ecospat.cv.gbm <- function(gbm.obj,data.cv, K=10, cv.lim = 10, jack.knife = F) 
{
  
  
  best.itr.temp <- gbm.perf(gbm.obj, method="cv", plot.it=F)
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  n <- length(gbm.obj$data$y)
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  n.0 <- length(gbm.obj$data$y[gbm.obj$data$y==0])
  
  n.1 <- length(gbm.obj$data$y[gbm.obj$data$y==1])
  
  
  
  if (jack.knife  == F && n.1 < cv.lim || jack.knife  == F && n.1 < K || jack.knife  == T)
  {
    df.res <- data.frame(id=NA,predicted=NA)
    K <- nrow(data.cv)
    cat("K has been set to ",K," -> Jack-Knife Cross-validation enabled!","\n",append = F)	
    
    for (i in 1:K)
    {	
      gbm.temp <- gbm.more(gbm.obj, n.new.trees =best.itr.temp, data=data.cv[-i,])
      
      
      
      if (i == 1)
      {
        vect.id <- i
        vect.predicted <- as.double(predict.gbm(gbm.temp, data.cv[i,],best.itr.temp,type="response"))
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, i, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.double(predict.gbm(gbm.temp, data.cv[i,],best.itr.temp,type="response")), after=length(vect.predicted))
      }		
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(gbm.obj$data$y),predictions=df.tmp.res[,2])
    
    
  }	else 
  {
    
    id.0 <- as.vector(names(gbm.obj$data$y)[gbm.obj$data$y==0],mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0)) 
    {
      K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0) 
    { 
      cat("K.0 has been set to", K.lst.0,"\n",append=F)
    }
    
    id.1 <- as.vector(names(gbm.obj$data$y)[gbm.obj$data$y==1], mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0)) 
    { 
      K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1) 
    { 
      cat("K.1 has been set to", K.lst.1,"\n",append=F)
    }   
    if (K.lst.0 != K.lst.1)
    {
      cat("P/A stratifications have not the same values","\n",append=F)
      min.K <- min(K.lst.0,K.lst.1)
      K.lst.0 <- min.K
      K.lst.1 <- min.K
    }
    
    cat("K has been finally set to",K.lst.0,"\n",append=F)
    K.lst <- K.lst.0
    
    # P/A STRATIFICATION     
    
    # ABSENCE STRATIFICATION    
    f.0 <- ceiling(n.0/K.lst.0)
    s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
    n.s.0 <- table(s.0)
    ms.0 <- max(s.0)
    
    # PRESENCE STRATIFICATION    
    f.1 <- ceiling(n.1/K.lst.1)
    s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
    n.s.1 <- table(s.1)
    ms.1 <- max(s.1)
    
    # RESPONSE PREDICTION
    df.res <- data.frame(id=NA,predicted=NA)
    for (i in 1:K.lst) 
    {
      j.out <- id.0[(s.0 == i)]
      j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
      j.out <- sort(j.out)
      j.in <- id.0[(s.0 != i)]
      j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
      j.in <- sort(j.in)
      
      
      gbm.cal <-gbm.more(gbm.obj, n.new.trees = best.itr.temp, data=data.cv[j.in,])
      
      
      gbm.val <- predict.gbm(gbm.cal, data.cv[j.out,], best.itr.temp,type = "response")
      if (i == 1)
      {
        vect.id <- j.out
        vect.predicted <- as.vector(gbm.val,mode="numeric")
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, j.out, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.vector(gbm.val), after=length(vect.predicted))
      }	
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(gbm.obj$data$y),predictions=df.tmp.res[,2])
    
    
    
  }
  df.res[df.res[,3]>1,3] <- 1
  df.res[df.res[,3]<0,3] <- 0
  
  return(df.res)
  
}

##################################################################################################################################
##################################################################################################################################

###############
## CV-MaxEnt ##
###############


ecospat.cv.me <- function(data.cv.me,name.sp,names.pred,K=10, cv.lim = 10, jack.knife = F)
{
  jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
  Sys.setenv(NOAWT=TRUE)
  
  rownames(data.cv.me) <- 1:nrow(data.cv.me)
  
  #best.itr.temp <- me.perf(me.obj, method="cv", plot.it=F)
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  n <- nrow(data.cv.me)
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  n.0 <-  nrow(data.cv.me[data.cv.me[,name.sp]==0,])
  n.1 <-  nrow(data.cv.me[data.cv.me[,name.sp]==0,])
  
  
  if (jack.knife  == F && n.1 < cv.lim || jack.knife  == F && n.1 < K || jack.knife  == T)
  {
    df.res <- data.frame(id=NA,predicted=NA)
    K <- nrow(data.cv.me)
    cat("K has been set to ",K," -> Jack-Knife Cross-validation enabled!","\n",append = F)    
    
    
    for (i in 1:K)
    {    
      
      me.cal.jk <-maxent(data.cv.me[-i,names.pred],as.vector(data.cv.me[-i,name.sp]))
      
      if (i == 1)
      {
        vect.id <- i
        vect.predicted <- round(predict(me.cal.jk,data.cv.me[i,names.pred]),3)
      } else if (i > 1)
      {    
        vect.id <- append(vect.id, i, after=length(vect.id))
        vect.predicted <- append(vect.predicted, round(predict(me.cal.jk,data.cv.me[i,names.pred]),3), after=length(vect.predicted))
      }        
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    # To correct
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv.me[,name.sp]),predictions=df.tmp.res[,2])
    
    
  }    else
  {
    
    id.0 <- as.vector(row.names(data.cv.me)[data.cv.me[,name.sp]==0],mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0))
    {
      K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0)
    {
      cat("K.0 has been set to", K.lst.0,"\n",append=F)
    }
    
    id.1 <- as.vector(row.names(data.cv.me)[data.cv.me[,name.sp]==1],mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0))
    {
      K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1)
    {
      cat("K.1 has been set to", K.lst.1,"\n",append=F)
    }   
    if (K.lst.0 != K.lst.1)
    {
      cat("P/A stratifications have not the same values","\n",append=F)
      min.K <- min(K.lst.0,K.lst.1)
      K.lst.0 <- min.K
      K.lst.1 <- min.K
    }
    
    cat("K has been finally set to",K.lst.0,"\n",append=F)
    K.lst <- K.lst.0
    
    # P/A STRATIFICATION     
    
    # ABSENCE STRATIFICATION    
    f.0 <- ceiling(n.0/K.lst.0)
    s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
    n.s.0 <- table(s.0)
    ms.0 <- max(s.0)
    
    # PRESENCE STRATIFICATION    
    f.1 <- ceiling(n.1/K.lst.1)
    s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
    n.s.1 <- table(s.1)
    ms.1 <- max(s.1)
    
    # RESPONSE PREDICTION
    df.res <- data.frame(id=NA,predicted=NA)
    for (i in 1:K.lst)
    {
      j.out <- id.0[(s.0 == i)]
      j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out))
      j.out <- sort(j.out)
      j.in <- id.0[(s.0 != i)]
      j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in))
      j.in <- sort(j.in)
      
      
      me.cal.cv <-maxent(data.cv.me[j.in,names.pred],as.vector(data.cv.me[j.in,name.sp]))
      
      me.val    <- round(predict(me.cal.cv,data.cv.me[j.out,names.pred]),3)
      
      if (i == 1)
      {
        vect.id <- j.out
        vect.predicted <- as.vector(me.val,mode="numeric")
      } else if (i > 1)
      {    
        vect.id <- append(vect.id, j.out, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.vector(me.val), after=length(vect.predicted))
      }    
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv.me[,name.sp]),predictions=df.tmp.res[,2])
    
    
    
  }
  df.res[df.res[,3]>1,3] <- 1
  df.res[df.res[,3]<0,3] <- 0
  
  return(df.res)
  
}

##################################################################################################################################
##################################################################################################################################

###########
## CV-RF ##
###########

ecospat.cv.rf <- function(rf.obj,data.cv, K=10, cv.lim = 10, jack.knife = F)
{
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  n <- nrow(data.cv)
  if ((K > n) || (K <= 1)) {stop("K outside allowable range")}
  
  n.0 <- length(rf.obj$y[rf.obj$y==0])
  
  n.1 <- length(rf.obj$y[rf.obj$y==1])
  
  if (jack.knife  == F && n.1 < cv.lim || jack.knife  == F && n.1 < K || jack.knife  == T)
  {
    df.res <- data.frame(id=NA,predicted=NA)
    K <- nrow(data.cv)
    cat("K has been set to",K,"(leave-one-out CV is enabled!)","\n",append = F)	
    
    
    
    for (i in 1:K)
    {	
      rf.tmp.cal <- randomForest(x=data.cv[-i,2:ncol(data.cv)], y=as.factor(data.cv[-i,1]), ntree=1000, importance=TRUE)
      
      rf.tmp.eval <- predict(rf.tmp.cal, data.cv[i,2:ncol(data.cv)], type="prob")  
      if (i == 1)
      {
        vect.id <- i
        vect.predicted <- as.double(rf.tmp.eval[,2])
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, i, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.double(rf.tmp.eval[,2]), after=length(vect.predicted))
      }		
      
    }
    
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(rf.obj$y),predictions=df.tmp.res[,2])
    
    
  }	else 
  {
    
    id.0 <- as.vector(row.names(data.cv)[data.cv[,1]==0], mode = "numeric")
    K.0 <- K
    K.lst.0 <- round(K)
    kvals.0 <- unique(round(n.0/(1:floor(n.0/2))))
    temp.0 <- abs(kvals.0 - K.lst.0)
    if (!any(temp.0 == 0)) 
    {
      K.lst.0 <- kvals.0[temp.0 == min(temp.0)][1]
    }
    if (K.lst.0 != K.0) 
    { 
      cat("K.0 has been set to", K.lst.0,"\n",append=F)
    }
    
    id.1 <- as.vector(row.names(data.cv)[data.cv[,1]==1], mode = "numeric")
    K.1 <- K
    K.lst.1 <- round(K)
    kvals.1 <- unique(round(n.1/(1:floor(n.1/2))))
    temp.1 <- abs(kvals.1 - K.lst.1)
    if (!any(temp.1 == 0)) 
    { 
      K.lst.1 <- kvals.1[temp.1 == min(temp.1)][1]
    }
    if (K.lst.1 != K.1) 
    { 
      cat("K.1 has been set to", K.lst.1,"\n",append=F)
    }   
    if (K.lst.0 != K.lst.1)
    {
      cat("P/A stratifications have not the same values","\n",append=F)
      min.K <- min(K.lst.0,K.lst.1)
      K.lst.0 <- min.K
      K.lst.1 <- min.K
    }
    
    cat("K has been finally set to",K.lst.0,"\n",append=F)
    K.lst <- K.lst.0
    
    # P/A STRATIFICATION     
    
    # ABSENCE STRATIFICATION    
    f.0 <- ceiling(n.0/K.lst.0)
    s.0 <- sample(rep(1:K.lst.0,f.0), n.0)
    n.s.0 <- table(s.0)
    ms.0 <- max(s.0)
    
    # PRESENCE STRATIFICATION    
    f.1 <- ceiling(n.1/K.lst.1)
    s.1 <- sample(rep(1:K.lst.1,f.1), n.1)
    n.s.1 <- table(s.1)
    ms.1 <- max(s.1)
    
    # RESPONSE PREDICTION
    df.res <- data.frame(id=NA,predicted=NA)
    for (i in 1:K.lst) 
    {
      j.out <- id.0[(s.0 == i)]
      j.out <- append(j.out, id.1[(s.1 == i)], after=length(j.out)) 
      j.out <- sort(j.out)
      j.in <- id.0[(s.0 != i)]
      j.in <- append(j.in, id.1[(s.1 != i)], after=length(j.in)) 
      j.in <- sort(j.in)
      
      set.seed(71)
      rf.tmp.cal <- randomForest(x=data.cv[j.in,2:ncol(data.cv)], y=as.factor(data.cv[j.in,1]), ntree=1000, importance=TRUE)
      
      rf.tmp.val <- predict(rf.tmp.cal, data.cv[j.out,2:ncol(data.cv)], type="prob")
      
      if (i == 1)
      {
        vect.id <- j.out
        vect.predicted <- as.vector(rf.tmp.val[,2],mode="numeric")
      } else if (i > 1)
      {	
        vect.id <- append(vect.id, j.out, after=length(vect.id))
        vect.predicted <- append(vect.predicted, as.vector(rf.tmp.val[,2]), after=length(vect.predicted))
      }	
      
      
    }
    df.tmp.res <-  data.frame(cbind(id=round(as.numeric(vect.id),digits=0),predictions=as.numeric(vect.predicted))[order(vect.id),])
    
    df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(rf.obj$y),predictions=df.tmp.res[,2])
    
    
  }
  df.res[df.res[,3]>1,3] <- 1
  df.res[df.res[,3]<0,3] <- 0
  
  return(df.res)
}

##################################################################################################################################
##################################################################################################################################

################
## Permut-GLM ##
################

ecospat.permut.glm <- function(glm.obj,nperm)  
{
  D2 <- function(glmobj) {
    go <- glmobj
    d2 <- (go$null.deviance - go$deviance)/go$null.deviance
    if (d2 < 0) {
      d2 <- 0
      return(d2)
    }
    return(d2)
  }
  adj.D2 <- function(glmobj) {
    go <- glmobj
    D2 <- (go$null.deviance - go$deviance)/go$null.deviance
    p <- length(go$coefficients)
    n <- length(go$fitted)
    adj.D2 <- 1 - ((n - 1)/(n - p)) * (1 - D2)
    if (adj.D2 < 0) {
      adj.D2 <- 0
      return(adj.D2)
    }
    return(adj.D2)
  }
  mframe <- model.frame(glm.obj)
  #newframe <- model.frame(glm.obj)
  v.coef <- substring(names(glm.obj$coefficients[1]), 2, nchar(names(glm.obj$coefficients[1])) - 
                        1)
  l.coef <- length(glm.obj$coefficients)
  for (l in 2:l.coef) {
    v.coef[l] <- substring(names(glm.obj$coefficients[l]), 
                           regexpr(")", names(glm.obj$coefficients[l]))[1] + 
                             1)
  }
  initial.D2 <- D2(glm.obj)
  initial.adj.D2 <- adj.D2(glm.obj)
  vect.dev <- append(initial.D2, initial.adj.D2, after = length(initial.D2))
  dev.names <- c("D2", "adj.D2")
  l.coeff <- length(glm.obj$coefficients)
  mterms <- attr(mframe, "terms")
  resp <<- attr(mterms, "response")
  resp <- attr(mterms, "response")
  vect.coeff <- c(as.vector(coef(glm.obj), mode = "numeric"), 
                  vect.dev)
  theta <- matrix(0, nperm, l.coeff, dimnames = list(NULL, 
                                                     names(glm.obj$coefficients)))
  dev <- matrix(0, nperm, length(dev.names), dimnames = list(NULL, 
                                                             dev.names))
  newframe <<- mframe
  names.y <- as.vector(glm.obj$formula[2], mode = "character")
  for (i in seq(nperm)) {
    newframe[, resp] <- sample(mframe[, resp])
    mc.glm.obj <- update(glm.obj, newframe[, resp] ~ .)
    theta[i, ] <- as.vector(coef(mc.glm.obj))
    dev[i, 1] <- round(D2(mc.glm.obj), 5)
    dev[i, 2] <- round(adj.D2(mc.glm.obj), 5)
    if (i%%100 == 0) {
      cat("Permutation ", i, "is starting now...", "\n", 
          append = F)
    }
  }
  vect.names <- v.coef
  vect.names <- append(vect.names, dev.names, after = length(vect.names))
  mat.permut <- matrix(0, nperm, length(vect.names), dimnames = list(NULL, 
                                                                     vect.names))
  mat.permut[, 1:l.coeff] <- theta
  mat.permut[, (l.coeff + 1):(l.coeff + 2)] <- dev
  
  
  # CALCULATES P-VALUES  BASED ON RANDOM DISTRIBUTIONS
  mat.pval <- matrix(0, length(vect.names), 1, dimnames = list(vect.names, "pval"))
  cat("...", "\n", append = F)
  cat("Conmputing pvalues by permutations...", "\n", append = F)
  cat("...", "\n", append = F)
  for (j in 1:ncol(mat.permut)) {
    vect.delta.quantile.coeff <- abs(quantile(mat.permut[, 
                                                         j], probs = seq(0, 1, 1e-05)) - vect.coeff[j])
    prob.coeff <- names(vect.delta.quantile.coeff)[vect.delta.quantile.coeff == 
                                                     min(vect.delta.quantile.coeff)]
    temp.val.coeff <- max(as.numeric(substring(prob.coeff, 
                                               1, nchar(prob.coeff) - 1)))
    if (temp.val.coeff == 100 || temp.val.coeff == 0) {
      mat.pval[j, 1] <- 0
    }
    else if (50 <= temp.val.coeff && temp.val.coeff < 100) {
      mat.pval[j, 1] <- 1 - (temp.val.coeff/100)
    }
    else if (0 < temp.val.coeff && temp.val.coeff < 50) {
      mat.pval[j, 1] <- temp.val.coeff/100
    }
  }
  return(mat.pval)
  
}  



##################################################################################################################################
##################################################################################################################################

##############
## Var.part ##
##############
ecospat.varpart<-function (model.1, model.2, model.12) 
{
  vect.estim <- c("D2.model.1", "D2.model.2", "D2.Joined.Contrib.1.2.", 
                  "D2.Unexplained")
  mat.estim <- matrix(0, length(vect.estim), 1, dimnames = list(vect.estim, 
                                                                "Deviance"))
  R2_rsnag <- function(glm.obj) {
    n.s <- length(glm.obj$y)
    r2.rsnag <- (1 - (exp(-glm.obj$null.deviance/2)/exp(-glm.obj$deviance/2))^(2/n.s))/(1 - 
                                                                                          exp(-glm.obj$null.deviance/2)^(2/n.s))
    if (r2.rsnag < 0) 
      r2.rsnag <- 0
    return(r2.rsnag)
  }
  R2.rsnag_12 <- R2_rsnag(model.12)
  R2.rsnag_1 <- R2_rsnag(model.1)
  R2.rsnag_2 <- R2_rsnag(model.2)
  
  # due to 1 alone (among what is explained by 1 and 2)
  R2.rsnag_1al  <- R2.rsnag_12 - R2.rsnag_2
  
  # due to 2 alone (among what is explained by 1 and 2)
  R2.rsnag_2al  <- R2.rsnag_12 - R2.rsnag_1
  
  # common to 1 and 2, note that it could be negative (ie no "common" part, 1 and 2 is larger than 1+2, see Hamilton 87 for linear regression)
  R2.rsnag_12jn  <- R2.rsnag_1 + R2.rsnag_2 - R2.rsnag_12
  
  # unexplained
  R2.rsnag_un  <- 1 - (R2.rsnag_1al + R2.rsnag_2al + R2.rsnag_12jn)
  
  mat.estim[1,1] <- round(R2.rsnag_1al,3)		# Mod1 alone
  mat.estim[2,1] <- round(R2.rsnag_2al,3)			# Mod2 alone
  mat.estim[3,1] <- round(R2.rsnag_12jn,3) 			# Mod1 + Mod2 alone
  mat.estim[4,1] <- round(R2.rsnag_un,3)			# Unexplained
  
  return(mat.estim) 
}

##################################################################################################################################
##################################################################################################################################


########
# CV example function.
########



ecospat.cv.example <- function ()
{
# Position of the first predictor
p.pred <- 4

# Position of the first species (i.e. respone variable)
p.resp <- 9

#Argument for Maxent calibration
me.args=c('-J','-P','threshold=false','product=false','quadratic=false','linear=false')

# Input and Output folders
#inpath.data <- "/Users/crandin/1_Projets/16_BIOASSEMBLE/ECOSPAT/"
#outpath <- "/Users/crandin/1_Projets/16_BIOASSEMBLE/ECOSPAT/"
outpath <- getwd ()

# Path for Maxent output
me.path <- outpath

# IMPORTING DATASET

# Read directly from the package data folder.
#data("ecospat.testData",envir=ecospat.env)
#data <- get("ecospat.testData", envir = ecospat.env)[,1:10]
data <- ecospat.testData[,1:10]

# Run the code with the ECOSPAT example dataset
glm.gam.gbm.me.rf.fit.model(data,p.pred,p.resp,me.args=c('-J',
'-P','threshold=false','product=false','quadratic=false','linear=false'), me.path= outpath,outpath)

#########################################################################################################################################




#########################################################################################################################################

# Examples for running cross-validation function
      
        ########################
		# GLM MODEL EVALUATION #
		########################
		
		# STRATIFIED-CV FOR GLM
		df.out.cv.glm.achatr <- ecospat.cv.glm(get("glm.Achillea_atrata",envir = ecospat.env), K = 10, cv.lim = 10, jack.knife = F)
   
   
		#permut.glm.saopp <- ecospat.permut.glm(glm.SAOPP,100)	
		
			
		########################
		# GAM MODEL EVALUATION #
		########################          
                     
		# STRATIFIED-CV FOR GAM
		df.out.cv.gam.achatr <- ecospat.cv.gam(get("gam.Achillea_atrata",envir = ecospat.env), K = 10, cv.lim = 10, jack.knife = F)

		########################
		# GBM MODEL EVALUATION #
		########################
		
		# STRATIFIED-CV FOR GBM
		df.out.cv.gbm.achatr <- ecospat.cv.gbm(get("gbm.Achillea_atrata",envir = ecospat.env), data, K = 10, cv.lim = 10, jack.knife = F)
		
    	#######################
		# ME MODEL EVALUATION #
		#######################
		
		# STRATIFIED-CV FOR ME
		# df.out.cv.me <- ecospat.cv.me(df.tmp.input,names(df.tmp.input[1]),names(df.tmp.input[2:ncol(df.tmp.input)]),K=10, cv.lim = 10, jack.knife = F)
		
		df.out.cv.me.achatr <- ecospat.cv.me(data, names(data)[9], names(data)[4:8], K = 10, cv.lim = 10, jack.knife = F)
		
	   	#######################
		# RF MODEL EVALUATION #
		#######################
		
		# STRATIFIED-CV FOR RF
		df.out.cv.rf.achatr <- ecospat.cv.rf(get("rf.Achillea_atrata",envir = ecospat.env), data[, c(9, 4:8)], K = 10, cv.lim = 10, jack.knife = F)
}

