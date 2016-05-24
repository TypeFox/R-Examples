threshpt <-
function(formula=formula, family=family, data=data, expvar=expvar, startrng=startrng, endrng=endrng, searchunit=searchunit, ...)
{
   	##library(splines)
	myformula=formula; myfamily=family; mydata=data

    fom<-paste(as.character(myformula)[c(2,1,3)], collapse="")      
    resvar<-as.character(myformula)[2]   #to get the name of response variable
    var2<-paste(expvar,"_2", sep="")     
    myform_glm=as.formula(paste(fom, var2, sep="+")) 

   	a<-strsplit(as.character(myformula)[3], " ")[[1]]
    	b<-colnames(mydata)[match(a,colnames(mydata))]
    	para_var<-b[which(is.na(b)==F)]

	a[substr(a,1,3)=="ns("]<-substr(a[substr(a,1,3)=="ns("], 4, nchar(a[substr(a,1,3)=="ns("])-1)
	a[substr(a,1,7)=="factor("]<-substr(a[substr(a,1,7)=="factor("], 8, nchar(a[substr(a,1,7)=="factor("])-1)
	a[substr(a,1,10)=="as.factor("]<-substr(a[substr(a,1,10)=="as.factor("], 11, nchar(a[substr(a,1,10)=="as.factor("])-1)
    	b2<-colnames(mydata)[match(a, colnames(mydata))]
    	varlist<-b2[which(is.na(b2)==F)]
	

#### Finding threshold

    thresh_best<-array("",c(8))
    mydata<-mydata[which(is.na(mydata[,expvar])==F),]
    if(endrng>max(mydata[,expvar])){endrng=floor(  max(mydata[,expvar])-1)}

    range_pt=c(startrng, endrng)
    diff_pt=endrng-startrng

    if (diff_pt<0) {
    print("error on search range ", range_pt)
    print("Range difference must be greater than 0")
    stop    }

    N=((diff_pt)*(1/searchunit)+1)*8
    M1<-matrix(numeric(N), ncol=8, byrow=TRUE,
    dimnames=list(c(),c("beta1","se1","p1","beta1+2","se1+2","p1+2","deviance","threshold")))
    M1[,8]<-seq(range_pt[1],range_pt[2],searchunit)

    n_for= diff_pt*1/searchunit+1

    	for (i in 1:n_for)
    	{
    	mydata[,var2]<-ifelse(mydata[,expvar]-M1[i,8]>0, mydata[,expvar]-M1[i,8], 0)
    	m11<-glm(myform_glm, data=mydata, family=myfamily)
   	a<-summary(m11)
    	ncoef = dim(vcov(m11))[1]
    	k=which(rownames(a$coefficients)==expvar)
    	c<-rep(0, ncoef); c[k]<-1; c[ncoef]<-1;
    	want = c %*% m11$coef[which(is.na(m11$coef)==F)]
    	se<-sqrt(t(c)%*%vcov(m11)%*%c)
   	Z<- want/se
   	p<-(1-pnorm(abs(Z)))*2
   	m11$aic
    	M1[i,1:7]<-c(c(a$coef[k,1], a$coef[k,2]), a$coef[k,4], c(want,se), p, m11$deviance)
	} # i

    thresh_best=M1[order(M1[,7],decreasing=FALSE)[1],]

    mydata[,var2]<-ifelse(mydata[, expvar]-thresh_best[8]>0,mydata[, expvar]-thresh_best[8],0)

    m11<-glm(myform_glm, data=mydata, family=myfamily)
    a<-summary(m11)
    a$coef


#### Calculate fitted values

      pred_y<-predict(m11, newdata=mydata, type="link")

        if(myfamily=="poisson"|myfamily=="quasipoisson"){
          predy1=exp(pred_y)}
        if(myfamily=="gaussian"){
          predy1=pred_y}
        if(myfamily=="binomial"|myfamily=="quasibinomial"){
          predy1=exp(pred_y)/(1+exp(pred_y))}


#### Return

      wc<-c()
      fn=names(m11$xlevels)
	fn[substr(fn,1,7)=="factor("]<-substr(fn[substr(fn,1,7)=="factor("], 8, nchar(fn[substr(fn,1,7)=="factor("])-1)
	fn[substr(fn,1,10)=="as.factor("]<-substr(fn[substr(fn,1,10)=="as.factor("], 11, nchar(fn[substr(fn,1,10)=="as.factor("])-1)

    	for(i in 1:length(varlist)){
		if(varlist[i] %in% fn)
		{wc[i]<-"factor"       
		} else{
		wc[i]<-class(mydata[, varlist[i]])
		}
	}

      l_fac<-which(wc=="factor")

        aa<-as.data.frame(summary(m11)$coef)
        parm_coef<-rbind(aa[1,-3], aa[which(is.na(match(rownames(aa),para_var[-which(para_var==expvar)]))==F),-3])

        # get coefficients of categorical variables
        cat_dim<-length(l_fac)
        cat_l<-c(); cat_tmp<-NULL

        for(i in 1:cat_dim){
          ll<-length(strsplit(names(m11$xlevels) , "")[[i]])      # get the lengths of categorical variables
          cat_tmp<-rbind(cat_tmp,aa[which(match(substr(rownames(aa), 1, ll), names(m11$xlevels)[i])==1),-3])     #substr: extract substrings in a character vector
          }
          
        parm_coef<-rbind(parm_coef, cat_tmp)
        parm_coef


 	object<-list()

	object$data=mydata	
	object$fitted.values<-predy1
	object$linear.predictors=m11$linear.predictors
	object$residuals=m11$residuals

	object$deviances<-M1[,7]
	object$thresholds<-M1[,8]

	object$best.fit=thresh_best
	object$parm.coef=parm_coef

	object$formula=myform_glm
	object$input.formula=myformula
	object$family=myfamily
	object$aic=m11$aic

	object$expvar=expvar
	object$resvar=resvar
	object$varlist=varlist
	object$varclass=wc

	object$glmobj=m11
	object$call=match.call()

	class(object)<-"threshpt"
	object
}
