plot.threshpt <-
function(x, select=NULL, se=T, expdf=4, xlim=NULL, ylim=NULL, xaxt=NULL, yaxt=NULL, 
				col.value=NULL, col.preval=NULL, col.ci=NULL, col.vline=NULL, lwd=NULL, pch=NULL, pch.preval=NULL, 
				main=NULL, xlab=NULL, ylab=NULL, ...)
{
	if(is.null(select)){select<-1}


if(select %in% 1:4)
{	
	myformula=x$input.formula; myfamily=x$family; mydata=x$data
	expvar=x$expvar; resvar=x$resvar; m11=x$glmobj

	if(select %in% 1:2)
	{

      a<-strsplit(as.character(myformula)[3], " ")[[1]]
 	d<-x$varlist
      wc<-x$varclass
      l_fac<-which(wc=="factor")  # l_fac: location of categorical variable

	    if(length(l_fac)==0){
            mean_data<-apply(mydata[,d], 2, mean, na.rm=T)
            mydata_1<-matrix(rep(mean_data, each=dim(mydata)[1]), dim(mydata)[1], length(mean_data))   
            mydata_1<-as.data.frame(mydata_1)
            colnames(mydata_2)<-d
            }

          if(length(l_fac)>0){
            mean_data<-apply(mydata[,d[-l_fac]], 2, mean, na.rm=T)
            mydata_1<-matrix(rep(mean_data, each=dim(mydata)[1]), dim(mydata)[1], length(mean_data))
            mydata_1<-as.data.frame(mydata_1)
            colnames(mydata_1)<-d[-l_fac]
            
            if(length(l_fac)>1){
        		for(i in 1:length(l_fac)){
			mydata_1[,d[l_fac][i]]<-as.factor(matrix(rep(x$glmobj$xlevels[[i]][1], each=dim(mydata)[1]), dim(mydata)[1], 1))
          		} # for
       	  } # if(length(l_fac)>1)
		if(length(l_fac)==1){
		mydata_1[,d[l_fac]]<- as.factor(matrix(rep(x$glmobj$xlevels[[1]][1], each=dim(mydata)[1]), dim(mydata)[1], 1))
              }

          } # if(length(l_fac)>0)

          # For explanatory variable, generate sequence of length N from min(expvar)*0.9 to max(expvar)*1.1  
     	    mydata_1[,expvar]<-seq(min(mydata[,expvar], na.rm=T)*0.9 ,max(mydata[,expvar],na.rm=T)*1.1, length.out=dim(mydata)[1])

		if (select==1) ### nonlinear curve
		{
	
	      a[which(a==expvar)]<-paste("ns(", expvar,",", expdf, ")",collapse="")
      	np_glm_form<-paste(paste(as.character(myformula)[c(2,1)], collapse=""), paste(a, collapse=""), collapse="")

    		m_np<-glm(np_glm_form, data=mydata, family=myfamily)    # np_glm_form: resvar~ ns( expvar , expdf )+others
       	yhat_ns=predict(m_np, se.fit=TRUE, newdata=mydata_1, type="link")
 	
       	#set confidence intervals
    	   	uci_ns=yhat_ns$fit+1.96*yhat_ns$se.fit; lci_ns=yhat_ns$fit-1.96*yhat_ns$se.fit

		if(is.null(main)){main<-paste("Fitted ", x$resvar, " vs ", "spline function of ", x$expvar, sep="")}
		if(is.null(xlab)){xlab<-x$expvar}
		if(is.null(ylab)){
       	 if(x$family=="poisson"){ylab<-paste("log(E(", x$resvar, "))", sep="")}
       	 if(x$family=="gaussian"){ylab<-paste("E(", x$resvar, ")", sep="")}
       	 if(x$family=="binomial"){ylab<-paste("logit(E(", x$resvar, "))", sep="")}
        	}
	
		if(is.null(col.value)){col.value<-"black"}
		if(is.null(col.ci)){col.ci<-"black"}
		if(is.null(lwd)){lwd<-1}

		if(is.null(xlim)){xlim=c(min(mydata_1[,expvar]), max(mydata_1[,expvar]))}
		if(is.null(ylim)){ylim=c(min(lci_ns), max(uci_ns))}
		if(is.null(xaxt)){xaxt<-"s"}
		if(is.null(yaxt)){yaxt<-"s"}


     		plot(mydata_1[,expvar], yhat_ns$fit, xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt,
			xlab=xlab, ylab=ylab ,main=main, col=col.value, lwd=lwd, type='l')
			if(se==TRUE){
     			lines(mydata_1[,expvar],uci_ns, lty=3, col=col.ci)
     			lines(mydata_1[,expvar],lci_ns, lty=3, col=col.ci)
			}
		}

     		if (select==2) # fitted val with mean covs
		{   

     		b<-colnames(mydata)[match(a,colnames(mydata))]
    		para_var<-b[which(is.na(b)==F)]
    		var2<-paste(expvar,"_2", sep="")     
 
           	mydata_1[,var2]<-ifelse(mydata_1[,expvar]-x$best.fit[8]>0, mydata_1[,expvar]-x$best.fit[8], 0)

       	yhat_lm=predict(m11, se.fit=TRUE, newdata=mydata_1, type="link")

       	#set confidence intervals
        	uci_lm=yhat_lm$fit+1.96*yhat_lm$se.fit; lci_lm=yhat_lm$fit-1.96*yhat_lm$se.fit

        	if(myfamily=="poisson"|myfamily=="quasipoisson"){
            predyhat=exp(yhat_lm$fit); predlci=exp(lci_lm); preduci=exp(uci_lm)
       	}
        	if(myfamily=="gaussian"){
            predyhat=yhat_lm$fit; predlci=lci_lm; preduci=uci_lm
         	RR<-NULL
       	}
       	if(myfamily=="binomial"|myfamily=="quasibinomial"){
            predyhat=exp(yhat_lm$fit)/(1+exp(yhat_lm$fit)); 
	      predlci=exp(lci_lm)/(1+exp(lci_lm)); preduci=exp(uci_lm)/(1+exp(uci_lm))
            RR<-NULL
       	}

     		if(is.null(main)){main<-paste("Fitted ", x$resvar, " vs ", x$expvar, " with other mean covariates", sep="")}
      	if(is.null(ylab)){ylab<-paste("Fitted E(", x$resvar, ")", sep="")}
      	if(is.null(xlab)){xlab<-x$expvar}

		if(is.null(col.value)){col.value<-"red"}
		if(is.null(col.ci)){col.ci<-"grey"}
		if(is.null(col.vline)){col.vline<-"green"}

		if(is.null(xlim)){xlim=c(min(mydata_1[,expvar]),max(mydata_1[,expvar]))}
		if(is.null(ylim)){ylim=c(min(predlci), max(preduci))}
		if(is.null(xaxt)){xaxt<-"s"};	if(is.null(yaxt)){yaxt<-"s"}

     		plot(seq(min(mydata_1[,expvar]),max(mydata_1[,expvar]),length=100), seq(min(predlci), max(preduci),length=100)
           		,type="n", ylab=ylab, xlab=xlab, xlim=xlim, ylim=ylim, main=main, xaxt=xaxt, yaxt=yaxt)

        	xx<- mydata_1[,expvar]
        	lines(xx[order(xx)], predyhat[order(xx)], col=col.value, lwd=2, lty=1)
			if(se==TRUE){
        	lines(xx[order(xx)], predlci[order(xx)], col=col.ci, lwd=1, lty=2)
        	lines(xx[order(xx)], preduci[order(xx)], col=col.ci, lwd=1, lty=2)
			}
        	abline(v=x$best.fit[8], col=col.vline)
		}
	} # if 1 or 2

   	if (select==3) # plot of fitted val	
	{
      if(is.null(main)){main<-paste("Observed and fitted ", x$resvar, " vs ", x$expvar, sep="")}
      if(is.null(ylab)){ylab<-x$resvar}
      if(is.null(xlab)){xlab<-x$expvar}

	if(is.null(col.value)){col.value<-"black"}
	if(is.null(col.preval)){col.preval<-"red"}
	if(is.null(col.vline)){col.vline<-"green"}
	if(is.null(pch)){pch<-1}
	if(is.null(pch.preval)){pch.preval<-1}

	if(is.null(xlim)){xlim=c(min(mydata[,x$expvar], na.rm=TRUE), max(mydata[,x$expvar], na.rm=TRUE))}
	if(is.null(ylim)){ylim=c(min(mydata[,x$resvar], na.rm=TRUE), max(mydata[,x$resvar], na.rm=TRUE))}
	if(is.null(xaxt)){xaxt<-"s"};	if(is.null(yaxt)){yaxt<-"s"}

      plot(mydata[,x$expvar], mydata[,x$resvar], xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, col=col.value, pch=pch, xaxt=xaxt, yaxt=yaxt)
      points(mydata[,x$expvar], x$fitted.values, col=col.preval, pch=pch.preval)
      abline(v=x$best.fit[8], col=col.vline)
	}

   	if (select==4) # deviance~threshold point plot
	{
      if(is.null(main)){main<-"Deviance plot"}
      if(is.null(ylab)){ylab<-"deviance"}
      if(is.null(xlab)){xlab<-"threshold"}
	if(is.null(col.value)){col.value<-"black"}
	if(is.null(col.vline)){col.vline<-"green"}
	if(is.null(pch)){pch<-1}

	if(is.null(xaxt)){xaxt<-"s"}
	if(is.null(yaxt)){yaxt<-"s"}

      plot(x$thresholds, x$deviances, main=main, xlab=xlab, ylab=ylab, col=col.value, pch=pch, xaxt=xaxt, yaxt=yaxt)
      abline(v=x$best.fit[8], col=col.vline)
	}

   } ## if select in 1-4
    else{cat("Error: select should be one of 1-4\n")}

}
