ffRanefLMER.fnc <- function(model,
	ran.effects=list(ran.intercepts=as.character(), # or can specify a vector
	slopes=as.character(),# of random effects to consider, e.g.,
	corr=as.character(), # for each slope, if not a factor, specify whether correlations (1) or not (0) should be added. 
	by.vars=as.character()), # c("(0+Length|Subject)","(1+Frequency|Subject)")
    alpha=0.05,
	if.warn.not.add=TRUE,
    log.file=NULL # or other path and file name or FALSE
){
  if(length(alpha)==0){
    stop("please supply a value to the \"alpha\" argument")
  }
 
	if(is.null(log.file)){
		log.file<-file.path(tempdir(),paste("ffRanefLMER_log_",
    		gsub(":","-",gsub(" ","_",date())),".txt",sep=""))
	}
  options(warn=1)

  data<-model@frame

  temp.dir=tempdir()
  tempdir()

  unlink(file.path(temp.dir,"temp.txt"))
  sink(file=NULL,type="message")   
  
  if(log.file!=FALSE)sink(file=log.file,split=TRUE)

  # Get model coefficients and summary
  if(as.vector(model@call[1])=="glmer()"){
    odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
    data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
	wo<-options()$warn
	options(warn=-1)
    temp.lmer <- update(model, . ~ ., family = "gaussian", data = data)
	options(warn=wo)
    coefs=row.names(anova(temp.lmer))
    data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
  } else {
    coefs=row.names(anova(model))
  }
  coefs=unique(unlist(strsplit(coefs,":")))
  coefs=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",coefs)
  coefs=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",coefs)
  coefs=c(names(getME(model,"flist")),coefs)

  if(is.list(ran.effects)){

    if(is.null(ran.effects$ran.intercepts)){
      ran.effects$ran.intercepts=as.character()
    }

    if(is.null(ran.effects$slopes)){
      ran.effects$slopes=as.character()
    }

    if(is.null(ran.effects$corr)){
      ran.effects$corr=as.character()
    }

    if(is.null(ran.effects$by.vars)){
      ran.effects$by.vars=as.character()
    }

    if(length(ran.effects$ran.intercepts)>0){
      intercepts=ran.effects$ran.intercepts
      ### Random intercepts ###
      cat("\t===     random intercepts     ===\n")
      # Determine which ranefs to include as a random effects
      for(intercept in intercepts){
        #if(!intercept%in%coefs){
          #cat("warning:",intercept,"not part of model coefficients\n")
          #cat("\tskipping\n")
        #} else {
          cat("evaluating addition of", paste("by-",intercept,sep=""),"random intercepts to model\n")
          
            wngs=file(file.path(temp.dir,"temp.txt"),open="w+",blocking=TRUE)
            sink(wngs,type="message")
      
          # Fit more complex model
	  model.updated<-NULL
          eval(parse(text=paste("model.updated=update(model,.~.+(1|",intercept,"))",sep="")))
      
        # Should the model with the more complex random-effects structure be kept?
          warn=readLines(wngs)
          sink(file=NULL,type="message")   
          unlink(file.path(temp.dir,"temp.txt"))
      
          if(if.warn.not.add && length(warn)>0){
            if(length(grep("warning",warn,value=TRUE))>0){
              cat(paste("\t",warn,"\n",sep=""))
              cat("\tnot adding",paste("(1 | ",intercept,")",sep=""),"to model\n")
            }
          } else {
              if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
                cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                cat("\tadding",paste("(1 | ",intercept,")",sep=""),"to model\n")
                model=model.updated
              } else {
                cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
                cat("\tnot adding",paste("(1 | ",intercept,")",sep=""),"to model\n")
              }
          }
        #}
      }
    }
  
    ### Random slopes ###
    cat("\t===       random slopes       ===\n")
#    if(length(ran.effects$slopes)==0){
#      if(as.vector(model@call[1])=="glmer()"){
#        odv=data[,as.character(unlist(as.list(model@call))$formula[2])]
#        data[,as.character(unlist(as.list(model@call))$formula[2])]=rnorm(nrow(data),0,1)
#        temp.lmer=update(model,.~.,family="gaussian",data=data)
#        coefs=row.names(anova(temp.lmer))
#        data[,as.character(unlist(as.list(model@call))$formula[2])]=odv
#        slopes=unique(unlist(strsplit(row.names(anova(temp.lmer)),":")))
#        slopes=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",slopes)
#      } else {
#        slopes=unique(unlist(strsplit(row.names(anova(model)),":")))
#        slopes=gsub("^rcs\\((.*), [[:digit:]]*\\)$","\\1",slopes)
#      }
#    } else { 
#      slopes=ran.effects$slopes
#    }

    if(length(ran.effects$slopes)>0){
      slopes<-ran.effects$slopes
      if(length(ran.effects$corr)==0){
		corr<-rep(0,length(ran.effects$slopes))
      }else{
		corr<-ran.effects$corr
      }
      by.vars<-ran.effects$by.vars
  
      for(slope in slopes){
        #if(!slope%in%coefs){
        #  cat("warning:",slope,"not part of model coefficients\n")
        #  cat("\tskipping\n")
        #} else {
          for(variable in by.vars){
	    if(is.factor(model@frame[,slope])){
            	cat("evaluating addition of",paste("(",slope,"|",variable,")",sep=""),"to model\n")
	    }else{
            	cat("evaluating addition of",paste("(",corr[grep(slope,slopes)]," + ",slope,"|",variable,")",sep=""),"to model\n")
	    }
            if(slope!=variable){
                wngs=file(file.path(temp.dir,"temp.txt"),open="w+")
                sink(wngs,type="message")
                #if(!variable%in%coefs){
                #  cat("\twarning:",variable,"not part of model coefficients\n")
                #  cat("\tskipping\n")
                #} else {
                  # Fit more complex model
		  if(is.factor(model@frame[,slope])){
                  	eval(parse(text=paste("model.updated=update(model,.~.+(",slope,"|",variable,"))",sep="")))
		  }else{
                  	eval(parse(text=paste("model.updated=update(model,.~.+(",corr[grep(slope,slopes)],"+",slope,"|",variable,"))",sep="")))
		  }
            
                  # Should the model with the more complex random-effects structure be kept?
                  warn=readLines(wngs)
                  sink(file=NULL,type="message")        
                  unlink(file.path(temp.dir,"temp.txt"))
            
                  if(if.warn.not.add && length(warn)>0){
                    if(length(grep("warning",warn,value=TRUE))>0){
                      cat(paste("\t",warn,"\n",sep=""))
		      if(is.factor(model@frame[,slope])){
                      	cat("\tnot adding",paste("(",slope," | ",variable,")",sep=""),"to model\n")
		      }else{
                      	cat("\tnot adding",paste("(",corr[grep(slope,slopes)],"+",slope," | ",variable,")",sep=""),"to model\n")
		      }
                    }
                  } else {      
                    if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
                      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
		      if(is.factor(model@frame[,slope])){
                      	cat("\tadding",paste("(",slope," | ",variable,")",sep=""),"to model\n")
		      }else{
                      	cat("\tadding",paste("(",corr[grep(slope,slopes)],"+",slope," | ",variable,")",sep=""),"to model\n")
		      }
                      model=model.updated
                    } else {
                      cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
		      if(is.factor(model@frame[,slope])){
                      	cat("\tnot adding",paste("(",slope," | ",variable,")",sep=""),"to model\n")
		      }else{
                      	cat("\tnot adding",paste("(",corr[grep(slope,slopes)],"+",slope," | ",variable,")",sep=""),"to model\n")
		      }
                    }
                 }
               #}
             } # close if(slope!=variable) loop
           } # close for(variable in by.vars) loop
         #}
       } # close for(slope in slopes) loop
     } # close if(length(ran.effects$slopes)>0)
    }else{
    for(ranef in ran.effects){
        ranef<-gsub(" ","",ranef)
        cat("evaluating addition of",ranef,"to model\n")
		model.term<-gsub("\\((.*)\\|.*","\\1",ranef)
		model.term<-gsub(".\\+(.*)","\\1",model.term)
		#model.term<-gsub("(.*)\\|.*","\\1",model.term)
        	#if(!model.term%in%coefs){
          	#	cat("warning: variable",model.term,"not part of model coefficients\n")
          	#	cat("\tskipping\n")
        	#}else{
          		wngs=file(file.path(temp.dir,"temp.txt"),open="w+",blocking=TRUE)
          		sink(wngs,type="message")
          		# Fit more complex model
          		eval(parse(text=paste("model.updated=update(model,.~.+",ranef,")",sep="")))
    
          		# Should the model with the more complex random-effects structure be kept?
          		warn=readLines(wngs)
          		unlink(file.path(temp.dir,"temp.txt"))
          		sink(file=NULL,type="message")        
          		
          		if(if.warn.not.add && length(warn)>0){
            			if(length(grep("Warning",warn,value=TRUE))>0){
              				cat(paste("\t",warn,"\n",sep=""))
              				cat("\tnot adding",ranef,"to model\n")
            			}else{
            				if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
              					cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              					cat("\tadding",ranef,"to model\n")
              					model=model.updated
            				}else{
              					cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              					cat("\tnot adding",ranef,"to model\n")
            				}
				}
          		}else{      
            			if(as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"])<=alpha){
              				cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              				cat("\tadding",ranef,"to model\n")
              				model=model.updated
            			}else{
              				cat("\tlog-likelihood ratio test p-value =",as.vector(anova(model,model.updated)[2,"Pr(>Chisq)"]),"\n")
              				cat("\tnot adding",ranef,"to model\n")
            			}
          		}
      		#}
	}
   }

  options(warn=0)
  sink(file=NULL,type="message")        

  if(log.file!=FALSE){
    cat("log file is",log.file,"\n")
    sink(file=NULL)
  }

   return(model=model)
}

