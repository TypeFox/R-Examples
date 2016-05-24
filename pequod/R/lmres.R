#(formula=BF_OP2 ~BF_OP1+BF_OP4+BF_OP3+ BF_OP4*BF_OP3,data=esame)

######################################################################
######################################################################

	#########################################################
	#                GENERAL FORMULA                        #
	#########################################################

preditt<-function(formula, data, residual_centering= FALSE, centered="none"){
	
   ## Centering ###################
	
	# dataframe for uncentered continuous predictors
	if("none" %in%centered){ 
		tempdata1<-data
		}
	
	# dataframe for mean-centered continuous predictors	
		else{
			vpos<-0
			options(warn=-1)
			for(i in (1:length(centered))){ 
				vpos[i]<-grep( centered[i], names(data))
				}
	options(warn=0)			
	
	
	centr<-function(k){
		k-mean(k,na.rm=T)
		}
	
	
	tempdata<-data
	tempdata1<-data
	
	if(length(vpos)==1){
		a<-centr((data)[,vpos])		                                     
		tempdata1[,vpos]<-a
		}
		else{
			tempdata<-apply(tempdata[,vpos],2,centr)
			tempdata1[,vpos]<-tempdata
			tempdata1<-as.data.frame(tempdata1)}
			}
	
	
	
	
   ## collinear lm model, model matrix and variable names 
	
	x<-lm( formula,tempdata1,na.action = na.omit)
	xm<-model.matrix(x)
	xmname<-dimnames(xm)[[2]]
	
	## return a list of variable names without ":" 
	jj<-strsplit(xmname,":")
	
	# return number of elements in the list (order)	
	n<-0
	k<-0
	ljj<-length(jj)
	for(i in (1:ljj)){
		n<-n+1
		k[n]<-length(jj[[i]])
		}

    # order frequencies table 
	tab<-table(k) 

	# as data.frame model matrix 
	matrx<-xm[,]
	matrxx<-as.data.frame(matrx)


    # number of terms of order x
	fI<-length(jj[k==1])  #number of first order terms + intercept
	fII<-length(jj[k==2]) #number of second order terms 
	fIII<-length(jj[k==3])#number of third order terms 
	fIV<-length(jj[k==4]) ##number of fourth order terms 

	 # exit and warning if number of fourth order terms >0
	if (fIV>0){
		print("maximum possible interaction <= 3")
		stop()
		}

	# replacing .XX. to :
	nomres<-gsub(":",".XX.",xmname)
	ert<-matrxx
	colnames(ert)<-nomres


	# first order matrix
	matrx1<-ert[,(2:fI)]

	# second order matrix
	
	if(fII>0){
		nres2<-nomres[(fI+1):(fI+fII)]
		matrx2<-ert[,(fI+1):(fI+fII)]
		}

	if(fII==1){
		matrx2<-as.matrix(matrx2,ncol=1)
		}

	if(fII>0){
		colnames(matrx2)<-nres2}


	# third order matrix
	if (fIII>0){
		matrx3<-ert[,(fI+fII+1):(fI+fII+fIII)]
		nres3<-nomres[(fI+fII+1):(fI+fII+fIII)]

	if(fIII==1){
		matrx3<-as.matrix(matrx3,ncol=1)
		}

	colnames(matrx3)<-nres3
		}
	
	
	if(residual_centering== FALSE){
		
		if(fI>0){
		mod2<-matrx1
		
		modfin<-cbind(x$model[1],mod2)	
	ind<-paste(names(mod2)[1:length(names(mod2))],collapse="+")
	
	#  first order regression
	formulax<-paste(names(x$model)[1],"~",ind)
	ortx<-lm(formula=formulax,data=modfin)
    pippox<-(ortx)

	#beta for first order regression (copied from QuantPsych package)
    coe<-pippox$coefficients[2:length(pippox$coefficients)]
    sy<-sapply(pippox$model[1],sd)
    sx<-sapply(pippox$model[2:length(pippox$coefficients)],sd)
    beta.stepI<-coe * sx/sy
	beta.stepI<- c(NA,beta.stepI)
		
		}
		
		if(fII>0){
		
			# new matrix: first order terms + second order residuals
	mod2<-cbind(matrx1,matrx2)
	modfin<-cbind(x$model[1],mod2)	
	# string of indipendent variables (first order + second order residual interactions)
	indxx<-paste(names(mod2)[1:length(names(mod2))],collapse="+")
	
	# second order orthogonalized regression
		formulaxx<-paste(names(x$model)[1],"~",indxx)
		ortxx<-lm(formula=formulaxx,data=modfin)
		pippoxx<-(ortxx)
			
			# beta for second order regression (copied from QuantPsych package)
		coe2<-pippoxx$coefficients[2:length(pippoxx$coefficients)]
		sy2<-sapply(pippoxx$model[1],sd)
		sx2<-sapply(pippoxx$model[2:length(pippoxx$coefficients)],sd)
		beta.stepII<-coe2 * sx2/sy2
		beta.stepII<- c(NA,beta.stepII)
			}
			
			if (fIII>0){
		
	# new matrix: first order terms + second order residuals + third order terms
		mod2<-cbind(matrx1,matrx2,matrx3)
			
			
	modfin<-cbind(x$model[1],mod2)	
		
		# third order orthogonalized regression
		ind3<-paste(names(mod2)[1:length(names(mod2))],collapse="+")
		formula3<-paste(names(x$model)[1],"~",ind3)
		ortxxx<-lm(formula=formula3,data=modfin)
		pippoxxx<-(ortxxx)
	
	# beta for third order regression (copied from QuantPsych package)
		coe3<-pippoxxx$coefficients[2:length(pippoxxx$coefficients)]
		sy3<-sapply(pippoxxx$model[1],sd)
		sx3<-sapply(pippoxxx$model[2:length(pippoxxx$coefficients)],sd)
		beta.stepIII<-coe3 * sx3/sy3
		beta.stepIII<- c(NA,beta.stepIII)
		
		}
		
		}
	
	
	
	
  ## lm with orthogonalized interactions ###################	
	if(residual_centering== TRUE){
		
		
	# string of indipendent variables (first order)
	ind<- paste(colnames(matrx1)[1:(dim(matrx1)[2])], collapse="+")
	
	
	
	# STEP I: FIRST ORDER #
	if(fI>0){
		mod2<-matrx1
		
		modfin<-cbind(x$model[1],mod2)	
		
	#  first order regression
	formulax<-paste(names(x$model)[1],"~",ind)
	ortx<-lm(formula=formulax,data=modfin)
    pippox<-(ortx)

	#beta for first order regression (copied from QuantPsych package)
    coe<-pippox$coefficients[2:length(pippox$coefficients)]
    sy<-sapply(pippox$model[1],sd)
    sx<-sapply(pippox$model[2:length(pippox$coefficients)],sd)
    beta.stepI<-coe * sx/sy
	beta.stepI<- c(NA,beta.stepI)
		
		}
	
	#second order residuals loop
	if(fII>0){
		formula1<-0
		res<-matrix(0,nrow=nrow(matrx2),ncol=fII)
		colnames(res)<-nres2
		n<-0
		for(i in (1:fII)){
			n<-n+1	
			jj<-strsplit(dimnames(matrx2)[[2]][i],".XX.")
			ind <- paste(jj[[1]], collapse="+")
			
			formula1[n]<-paste(dimnames(matrx2)[[2]][i],"~", ind)
			res[,n]<-as.vector(residuals(lm(formula1[n],data=ert)))
			}
			# new matrix: first order terms + second order residuals
	mod2<-cbind(matrx1,res)
	modfin<-cbind(x$model[1],mod2)	
	# string of indipendent variables (first order + second order residual interactions)
	indxx<-paste(names(mod2)[1:length(names(mod2))],collapse="+")
	
	# second order orthogonalized regression
		formulaxx<-paste(names(x$model)[1],"~",indxx)
		ortxx<-lm(formula=formulaxx,data=modfin)
		pippoxx<-(ortxx)
			
			# beta for second order regression (copied from QuantPsych package)
		coe2<-pippoxx$coefficients[2:length(pippoxx$coefficients)]
		sy2<-sapply(pippoxx$model[1],sd)
		sx2<-sapply(pippoxx$model[2:length(pippoxx$coefficients)],sd)
		beta.stepII<-coe2 * sx2/sy2
		beta.stepII<- c(NA,beta.stepII)
			}
			
			if (fIII>0){
		
	# new matrix: first order terms + second order residuals + third order terms
		mod2<-cbind(matrx1,res,matrx3)
			
			#third order residuals loop	
		formula2<-0
		res2<-matrix(0,nrow=nrow(matrx3),ncol=fIII)
		colnames(res2)<-nres3
		n<-0
		for(i in (1:fIII)){
			n<-n+1	
			jj<-strsplit(dimnames(matrx3)[[2]][i],".XX.")
			inx <- paste(jj[[1]], collapse="+")
			inxx<-paste(jj[[1]][1:2],collapse=".XX.")
			inxxx<-paste(jj[[1]][c(1,3)],collapse=".XX.")
			inxxxx<-paste(jj[[1]][2:3],collapse=".XX.")
			indxx<-paste(inx,"+",inxx,"+",inxxx,"+",inxxxx)
			
			formula2[n]<-paste(dimnames(matrx3)[[2]][i],"~", indxx)
			res2[,n]<-as.vector(residuals(lm(formula2[n],data=mod2)))
			}
		# new matrix: first order terms + second order residuals + third order residuals
		mod2<-cbind(matrx1,res,res2)
	modfin<-cbind(x$model[1],mod2)	
		
		# third order orthogonalized regression
		ind3<-paste(names(mod2)[1:length(names(mod2))],collapse="+")
		formula3<-paste(names(x$model)[1],"~",ind3)
		ortxxx<-lm(formula=formula3,data=modfin)
		pippoxxx<-(ortxxx)
	
	# beta for third order regression (copied from QuantPsych package)
		coe3<-pippoxxx$coefficients[2:length(pippoxxx$coefficients)]
		sy3<-sapply(pippoxxx$model[1],sd)
		sx3<-sapply(pippoxxx$model[2:length(pippoxxx$coefficients)],sd)
		beta.stepIII<-coe3 * sx3/sy3
		beta.stepIII<- c(NA,beta.stepIII)
		
		}
		}
		
	
	
		
	
	
	

	
	# max order for the final regression
	mo<-max(k)
	
	
	
 ## output list ################### 

	if(mo==1){
	ort=list(regr.order=mo,formula.Stepfin=formulax,beta.Stepfin= beta.stepI, Stepfin=(pippox))}
	if(mo==2){
	ort=list(regr.order=mo,formula.StepI=formulax,formula.Stepfin=formulaxx,beta.StepI= beta.stepI,  beta.Stepfin= beta.stepII, StepI=(pippox),Stepfin=(pippoxx),F_change=anova(pippox,pippoxx))}
	if(mo==3){
	ort=list(regr.order=mo, formula.StepI=formulax,formula.StepII=formulaxx,formula.Stepfin=formula3,beta.StepI= beta.stepI,beta.StepII= beta.stepII,  beta.Stepfin= beta.stepIII,StepI=(pippox),StepII=(pippoxx),Stepfin=(pippoxxx),F_change=anova(pippox,pippoxx,pippoxxx))}
	if(mo>3){
	ort=list(regr.order=mo,Warning="maximum possible interaction <= 3")}
	
	return(ort)
		}
	
	
	
	
	
	
	#########################################################
	#                OBJECT LMRES                           #
	#########################################################
	
	
	lmres<- function(formula,data,residual_centering,centered, ...) UseMethod("lmres")
	
 ## lmres default ################### 
	lmres.default <- function(formula, data,residual_centering=FALSE,centered="none", ...){
		formula <- as.formula(formula)
		data <- as.data.frame(data)
		centered<-as.character(centered)
		residual_centering<-as.logical(residual_centering)
		est <- preditt(formula, data,residual_centering,centered)
		est$call <- match.call()
		class(est) <- "lmres"
		est
		}
	
 ## print lmres  ################### 	
	print.lmres <- function(x, ...){
		cat("regression order:\n")
		print(x$regr.order)
		cat("\nregr:\n")
		print(x$Stepfin)
		}
	
	
	
	
	
	
	#########################################################
	#                OBJECT SUMMARY LMRES                   #
	#########################################################
	
	summary.lmres <- function(object, type="default",...){
		
		
  ## summary default ##
	
		if(type=="default"){
			type<-"default"
			mm<-object$Stepfin
	
	#formula
			formulafin<-object$formula.Stepfin
	
	
	#collinearity
			coll<- vif(object$Stepfin)
			coll.vif<- coll
			coll.tollerance<- 1/coll.vif
			colli<- cbind(VIF=coll.vif,Tolerance=coll.tollerance)
			colli<- round(colli,4)

	#coefficients
			se <- sqrt(diag(vcov(mm)))
			tval <- coef(mm) / se
			TAB <- cbind(Estimate = coef(mm),StdErr = se,t.value = tval,
			beta=object$beta.Stepfin, p.value = 2*pt(-abs(tval), df=mm$df))
			TAB<-round(TAB,5)

	#residuals
			resi<-summary(mm$residuals)
			resi<-round(resi,4)

	#fstatistic
			su<-summary(mm)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,4)

	#R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square<-su$r.squared
			adj.r.square<-su$adj.r.squared


	#output list
			res<-list(type=type, formulaf=formulafin, R=R, R_squared=R.square, Adjusted_R_squared=adj.r.square,
			fstatistic=fstatistic, residuals=resi,coefficients=TAB,collinearity=colli)
		}
		
			
 ## summary nested	##

	if(type=="nested"){
		regr.order=object$regr.order
		type<-"nested"
	
	
	# SUMMARY STEP I #
	
		if(regr.order==1){
			mm<-object$Stepfin
		
	#formula
			formula1<-as.formula(object$formula.Stepfin)
	
	#coefficients
			se <- sqrt(diag(vcov(mm)))
			tval <- coef(mm) / se
			TAB <- cbind(Estimate = coef(mm),StdErr = se,t.value = tval,
			beta=object$beta.Stepfin, p.value = 2*pt(-abs(tval), df=mm$df))
			TAB<-round(TAB,5)
	
			
	#fstatistic
			su<-summary(mm)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,2)

	#R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square3<-su$r.squared
			adj.r.square<-su$adj.r.squared

			stat<-cbind(R, R.square3, adj.r.square, fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))

			rownames(stat)<-"Model 1:"
	
	
			colnames(stat)<-c("R  ","R^2 "," Adj. R^2", "F  ","df1 ","df2 ","p.value")
			stat<-round(stat,2)
	
	#output list
			res<-list(type=type, stat=stat, regr.order=regr.order, f1=formula1,  coefficients=TAB)
			}	
		
		
		
	# SUMMARY STEP II #
	
		if(regr.order==2){
			mm1<- object$StepI
			mm<-object$Stepfin
	
	#formula step I
			formula1<-as.formula(object$formula.StepI)
	#formula step II
			formulafin<-as.formula(object$formula.Stepfin)
		
	#coefficients step I
			se <- sqrt(diag(vcov(mm1)))
			tval <- coef(mm1) / se
			TAB1 <- cbind(Estimate = coef(mm1),StdErr = se,t.value = tval,
			beta=object$beta.StepI, p.value = 2*pt(-abs(tval), df=mm1$df))
			TAB1<-round(TAB1,5)
	
			MODEL_1<-as.matrix(c(NA,NA,NA,NA,NA))
			MODEL_1<-t(MODEL_1)
			rownames(MODEL_1)="-- Model 1 --"
		
	#coefficients step II
			se <- sqrt(diag(vcov(mm)))
			tval <- coef(mm) / se
			TAB <- cbind(Estimate = coef(mm),StdErr = se,t.value = tval,
			beta=object$beta.Stepfin, p.value = 2*pt(-abs(tval), df=mm$df))
			TAB<-round(TAB,5)
	
			MODEL_2<-as.matrix(c(NA,NA,NA,NA,NA))
			MODEL_2<-t(MODEL_2)
			rownames(MODEL_2)="-- Model 2 --"
	
	# coefficients nested Tab
			TABgg<-rbind(MODEL_1,NA ,TAB1,NA,NA,MODEL_2,NA,TAB)
	
	
	#F_change
			F_change<-object$F_change
	
	
	#fstatistic step I
			su<-summary(mm1)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,4)

	#step I: R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square<-su$r.squared
			adj.r.square<-su$adj.r.squared
			R1<-R.square
			Fs1<-cbind(R, R.square, adj.r.square,R1,fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))

			rownames(Fs1)<-"Model 1"

	#fstatistic step II
			su<-summary(mm)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,2)

	#step II: R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square3<-su$r.squared
			adj.r.square<-su$adj.r.squared
			R3<-R.square3-R.square
			Fs3<-cbind(R, R.square3, adj.r.square,R3,fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))

			rownames(Fs3)<-"Model 2:"
	
			stat<-rbind(Fs1,Fs3)
			colnames(stat)<-c("R  ","R^2 "," Adj. R^2","  Diff.R^2","F  ","df1 ","df2 ","p.value")
			stat<-round(stat,2)
		
		
	#output list
			res<-list(type=type, stat=stat, regr.order=regr.order, f1=formula1,f2=formulafin,F_change=F_change, coefficients=TABgg)
			}	
	
	
	
	
	# SUMMARY STEP III #
	
		if(regr.order==3){
			mm1<- object$StepI
			mm2<-object$StepII
			mm<-object$Stepfin
	
	#formula step I
			formula1<-as.formula(object$formula.StepI)
	#formula step II
			formula2<-as.formula(object$formula.StepII)
	#formula step III
			formulafin<-as.formula(object$formula.Stepfin)
		
		
		
	#coefficients step I
			se <- sqrt(diag(vcov(mm1)))
			tval <- coef(mm1) / se
			TAB1 <- cbind(Estimate = coef(mm1),StdErr = se,t.value = tval,
			beta=object$beta.StepI, p.value = 2*pt(-abs(tval), df=mm1$df))
			TAB1<-round(TAB1,5)
	
			MODEL_1<-as.matrix(c(NA,NA,NA,NA,NA))
			MODEL_1<-t(MODEL_1)
			rownames(MODEL_1)="-- Model 1 --"
	
	#coefficients step II
			se <- sqrt(diag(vcov(mm2)))
			tval <- coef(mm2) / se
			TAB2 <- cbind(Estimate = coef(mm2),StdErr = se,t.value = tval,
			beta=object$beta.StepII, p.value = 2*pt(-abs(tval), df=mm2$df))
			TAB2<-round(TAB2,5)
	
			MODEL_2<-as.matrix(c(NA,NA,NA,NA,NA))
			MODEL_2<-t(MODEL_2)
			rownames(MODEL_2)="-- Model 2 --"
		
	#coefficients step III
			se <- sqrt(diag(vcov(mm)))
			tval <- coef(mm) / se
			TAB <- cbind(Estimate = coef(mm),StdErr = se,t.value = tval,
			beta=object$beta.Stepfin, p.value = 2*pt(-abs(tval), df=mm$df))
			TAB<-round(TAB,5)
	
			MODEL_3<-as.matrix(c(NA,NA,NA,NA,NA))
			MODEL_3<-t(MODEL_3)
			rownames(MODEL_3)="-- Model 3 --"
	
	# coefficients nested Tab
			TABgg<-rbind(MODEL_1,NA ,TAB1,NA,NA,MODEL_2,NA,TAB2,NA,NA,MODEL_3,NA,TAB)
		
	
	#F_change
			F_change<-object$F_change
	
	
	
	#fstatistic step I
			su<-summary(mm1)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,4)
			
			#step I: R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square<-su$r.squared
			adj.r.square<-su$adj.r.squared
			R1<-R.square
			Fs1<-cbind(R, R.square, adj.r.square,R1,fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))
			rownames(Fs1)<-"Model 1"
			
	#fstatistic step II
			su<-summary(mm2)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,4)
			
			#step II: R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square<-su$r.squared
			adj.r.square<-su$adj.r.squared
			R2<-R.square-R1
			Fs2<-cbind(R, R.square, adj.r.square,R2,fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))
			rownames(Fs2)<-"Model 2"
	
	#fstatistic step III
			su<-summary(mm)
			fstatistic<-su$fstatistic
			fstatistic<-round(fstatistic,4)
			



	#step III: R, R square and adjusted R square
			R<-sqrt(su$r.squared)
			R.square3<-su$r.squared
			adj.r.square<-su$adj.r.squared
			R3<-R.square3-R.square
			Fs3<-cbind(R, R.square3, adj.r.square,R3,fstatistic[1],fstatistic[2], fstatistic[3],
			p.value=pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=FALSE))
			rownames(Fs3)<-"Model 3:"
	
			
			stat<-rbind(Fs1,Fs2,Fs3)
			colnames(stat)<-c("R  ","R^2 "," Adj. R^2","  Diff.R^2","F  ","df1 ","df2 ","p.value")
			stat<-round(stat,2)
	
	
	
	#output list
			res<-list(type=type, stat=stat, regr.order=regr.order, f1=formula1,
			f2=formula2,f3=formulafin, F_change=F_change, coefficients=TABgg)
			}	
		}


 		class(res) <- "summary.lmres"
		res
		}


 ## print lsummary.mres  ################### 

	print.summary.lmres <- function(x, ...){
		
	 # default print summary 	
		if(x$type=="default"){
			cat("Formula:\n")
			print(as.formula(x$formulaf))
			cat("\n")
			
			cat("Models\n")	
			Fs<-cbind(x$R, x$R_squared, x$Adjusted_R_squared,value=x$fstatistic[1],
			numdf=x$fstatistic[2], dendf=x$fstatistic[3],
			p.value=pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE))
			colnames(Fs)<-c("R  ","R^2 "," Adj. R^2","F  ","df1 ","df2 ","p.value")
			rownames(Fs)<-"Model"
			printCoefmat(Fs,P.values=TRUE, has.Pvalue=TRUE,justify="centre",digits=3)
			cat("\n")

			cat("Residuals\n")
			print(x$residuals)
			cat("\n")

			cat("Coefficients\n")
			printCoefmat( x$coefficients, P.values=TRUE, has.Pvalue=TRUE,na.print="")
			cat("\n")

			cat("Collinearity\n")
			printCoefmat(x$collinearity,na.print="")
			cat("\n")
			}



  # nested print summary #
		if(x$type=="nested"){
			
	# models for regression order = 1	
			if(x$regr.order==1){
				cat("**Models**\n")	
				cat("\n")
				cat("Model 1: ")
				print(x$f1)
				cat("\n")
				cat("\n")
				}
	
	# models for regression order = 2	
			if(x$regr.order==2){
				cat("**Models**\n")	
				cat("\n")
				cat("Model 1: ")
				print(x$f1)
				cat("\n")
				cat("Model 2: ")
				print(x$f2)
				cat("\n")
				cat("\n")
				}
	# models for regression order = 3
			if(x$regr.order==3){
				cat("**Models**\n")	
				cat("\n")
				cat("Model 1: ")
				print(x$f1)
				cat("\n")
				cat("Model 2: ")
				print(x$f2)
				cat("\n")
				cat("Model 3: ")
				print(x$f3)
				cat("\n")
				cat("\n")	
				}


	# statistics
			cat("**Statistics**\n")
			cat("\n")	
			printCoefmat(x$stat,P.values=TRUE, has.Pvalue=TRUE,justify="centre")
			cat("\n")
			cat("\n")
			
	
	# F_change for regression order > 1
			if(x$regr.order>1){
				cat("**F change**\n")
				cat("\n")	
				Fs<-x$F_change
				printCoefmat(Fs,P.values=TRUE, has.Pvalue=TRUE,justify="centre",digits=3,na.print="")
				cat("\n")
				cat("\n")
				}


	# coefficients
			cat("**Coefficients**\n")
			cat("\n")
			printCoefmat( x$coefficients, P.values=TRUE, has.Pvalue=TRUE,na.print="")
			cat("\n")
	
		}
	}



