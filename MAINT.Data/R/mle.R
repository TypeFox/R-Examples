setMethod("mle",
	signature(Idt = "IData"),
	function(Idt,Model="Normal",Config=1:5,SelCrit=c("AIC","BIC"))
	{
	  Model <- match.arg(Model)
	  SelCrit <- match.arg(SelCrit)
		IdtNmle(Idt,Type="SingDst",Config=Config,SelCrit=SelCrit)
	}
)

setMethod("MANOVA",
	signature(Idt = "IData"),
	function(Idt, grouping, Model="Normal", Config=1:5,
			SelCrit=c("AIC","BIC"), Mxt=c("Hom","Het"), tol=1.0e-4)
	{
		if (!is.factor(grouping)) stop("'grouping' is not a factor\n")
		if ( nrow(Idt) != length(grouping)) stop("The numbers of data and partition observations are different\n")
		Mxt <- match.arg(Mxt)
		Model <- match.arg(Model)
		SelCrit <- match.arg(SelCrit)
		nk <- as.numeric(table(grouping))
		q <- Idt@NIVar 
		p <- 2*q 
		k <- length(nk) 
		if (k==1) stop("The data belongs to one single group. A partition into at least two different groups is required\n")
		H0res <- mle(Idt,Model=Model,Config=Config,SelCrit=SelCrit)
		if (Mxt=="Hom")  H1res <- IdtNmle(Idt,grouping,Type="HomMxt",Config=Config,SelCrit=SelCrit,tol=tol)
		else if (Mxt=="Het") H1res <- IdtHetMxtNmle(Idt,grouping,Config=Config,SelCrit=SelCrit,tol=tol)
		H0Ll <- H0res@logLiks[H1res@BestModel]
		H1Ll <- H1res@logLiks[H1res@BestModel]
		QuiSq <- 2*(H1Ll-H0Ll)
		names(H0Ll) <- names(H1Ll) <- names(QuiSq) <- NULL
		if (Mxt=="Hom") df <- p*(k-1)
		else if (Mxt=="Het") {
			BestConf <- H1res@ModelConfig[H1res@BestModel]
			df <- npar(BestConf,p,q,Ngrps=k,Mxt="Het") - npar(BestConf,p,q,Ngrps=1,Mxt="Hom")
		}
		pvalue <- pchisq(QuiSq,df,lower.tail=FALSE)
		if (Mxt=="Hom") 
			return( new("IdtClMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
							QuiSq=QuiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
		else if (Mxt=="Het") 
			return( new("IdtHetNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
							QuiSq=QuiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
	}
)

setMethod("show",					
	signature(object = "IdtE"),
	function(object)
	{
		cat("Log likelihoods:\n") ; print(object@logLiks)
		if (object@SelCrit=="AIC") { cat("Akaike Information Criteria:\n") ; print(object@AICs) }
		else if (object@SelCrit=="BIC") { cat("Bayesian (Schwartz) Information Criteria:\n") ; print(object@BICs) }
       		cat("Selected Model:\n") ; print(names(object@BestModel)) ; cat("\n")
		cat("Selected Model ML Estimates:\n") ; print(coef(object))
	}
)

setMethod("testMod",					
	signature(Idt = "IdtE"),
	function(Idt,RestMod=Idt@ModelConfig[2]:length(Idt@ModelConfig),FullMod="Next")
	{
		if (is.character(RestMod))  
			RestMod <- sapply(RestMod,function(Rmd) which(Rmd==Idt@ModelNames))
		if (FullMod[1]!="Next" && FullMod[1]!="All" && is.character(FullMod)) 
			FullMod <- sapply(FullMod,function(Fmd) which(Fmd==Idt@ModelNames))
		if (is.element(Idt@ModelConfig[1],RestMod))
			stop("Model",Idt@ModelNames[1],"can not be sepecified as a restricted model\n",
				"since it is the most general model that has been estimated.\n")
		TestRes <- list()
		RestModels <- character()
		FullModels <- character()
		for (RMind in RestMod)  {
			if (FullMod[1]=="Next") FMindices <- NextModel(RMind,ModelType="Normal")
			else if (FullMod[1]=="All") FMindices <- NestedBy(RMind,ModelType="Normal")
			else FMindices <- intersect(FullMod,NestedBy(RMind,ModelType="Normal"))
			for (FMind in FMindices)  {
				H0Ll <- Idt@logLiks[RMind]
				H1Ll <- Idt@logLiks[FMind]
				QuiSq <- 2*(H1Ll-H0Ll)
				q <- Idt@NIVar
				RMnpar <- npar(Idt@ModelConfig[RMind],2*q,q)
				FMnpar <- npar(Idt@ModelConfig[FMind],2*q,q)
				df <-  FMnpar - RMnpar
				pvalue <- pchisq(QuiSq,df,lower.tail=FALSE)
			   	resi <- new("LRTest",H0logLik=H0Ll,H1logLik=H1Ll,QuiSq=QuiSq,df=df,pvalue=pvalue)
			   	TestRes <- c(TestRes,resi)
				RestModels <- c(RestModels,Idt@ModelNames[RMind])
				FullModels <- c(FullModels,Idt@ModelNames[FMind])
			}
		}
		new("ConfTests",TestRes=TestRes,RestModels=RestModels,FullModels=FullModels)
	}
)

setMethod("show",					
	signature(object = "ConfTests"),
	function(object)
	{
		for (i in 1:length(object@TestRes)) if (!is.na(object@TestRes[[i]]@QuiSq)) 
		{
			cat("Testing Model",object@RestModels[i],"against alternative",object@FullModels[i],":\n")
			print(object@TestRes[[i]])
		}
	}
)

setMethod("BestModel",
	signature(Idt = "IdtE"),
	function(Idt,SelCrit=c("IdtCrt","AIC","BIC"))
	{
   		SelCrit <- match.arg(SelCrit)
		if (SelCrit == "IdtCrt") return(Idt@BestModel)
		else if (SelCrit == "AIC") return(which.min(Idt@AICs))
		else if (SelCrit == "BIC") return(which.min(Idt@ABCs))
	}
)

setMethod("coef",
	signature(object = "IdtNDE"),
	function(object,Model=BestModel(object),...)
	{
		list(mu=object@mleNmuE,Sigma=object@Configurations[[Model]]$mleSigE)
	}
)

setMethod("stdEr",
	signature(x = "IdtNDE"),
	function(x,Model=BestModel(x))
	{
		list(mu=x@mleNmuEse,Sigma=x@Configurations[[Model]]$mleSigEse)
	}
)

setMethod("show",					
	signature(object = "LRTest"),
	function(object)
	{
        	cat("Null Model log-likelihood:",object@H0logLik,"\n")
        	cat("Full Model log-likelihood:",object@H1logLik,"\n")
        	cat("Qui-squared statistic:",object@QuiSq,"\n")
        	cat("degrees of freedom:",object@df,"\n")
        	cat("p-value:",object@pvalue,"\n\n")
	}
)

setMethod("H1res",  signature(object = "IdtMANOVA"), function(object) object@H1res)

setMethod("H0res",  signature(object = "IdtMANOVA"), function(object) object@H0res)

setMethod("show",					
	signature(object = "IdtMANOVA"),
	function(object)
	{
		cat("Null Model Log likelihoods:\n") ; print(object@H0res@logLiks)
		cat("Full Model Log likelihoods:\n") ; print(object@H1res@logLiks)
		if (object@H1res@SelCrit=="AIC") { 
			cat("Full Model Akaike Information Criteria:\n")
			print(object@H1res@AICs)
		}
		else if (object@H1res@SelCrit=="BIC") {
			cat("Full Model Bayesian (Schwartz) Information Criteria:\n")
			print(object@H1res@BICs)
		}
       		cat("Selected Model:\n") ; print(names(object@H1res@BestModel)) ; cat("\n")
		callNextMethod()
	}
)

 
IdtNmle <- function(Idt,grouping=NULL,Type=c("SingDst","HomMxt"),Config=1:5,SelCrit=c("AIC","BIC"),tol=1.0e-4)
{
	Type <- match.arg(Type)
	SelCrit <- match.arg(SelCrit)
	q <- Idt@NIVar
	p <- 2*q
	n <- Idt@NObs
	if (Type=="SingDst") k <- 1
	else if (Type=="HomMxt") {
		if (is.null(grouping)) stop("Argument grouping is missing from MANOVA method\n")
		if (!is.factor(grouping)) stop("'grouping' is not a factor")
		nk <- as.numeric(table(grouping))
		if (n != sum(nk)) stop("Number of observations in IData object and grouping factor do not agree with each other\n")
		k <- length(nk) 
		if (k==1) stop("The data belongs to one single group. A partition into at least two different groups is required\n")
	}
  
	X <- cbind(Idt@MidP,Idt@LogR)
	logLiks <- AICs <- BICs <- rep(NA_real_,5)

	Configurations <- vector("list",5)
	names(logLiks)[1:5] <- names(AICs)[1:5] <- names(BICs)[1:5] <- names(Configurations)[1:5]<- paste("NC",1:5,sep="")
	for (model in Config)  
		if (model!=2) Configurations[[model]] <- list(mleSigE=NULL,mleSigEse=NULL,logLik=NULL,AIC=NULL,BIC=NULL)
		else Configurations[[model]] <- list(mleSigE=NULL,mleSigEse=NULL,logLik=NULL,AIC=NULL,BIC=NULL,optres=NULL)
	if (Type=="SingDst") {
		mleNmuE <- colMeans(X)
		mleNmuEse <- sapply(X,sd)/sqrt(n)
		if (is.element(1,Config)) {
			Configurations[[1]]$mleSigE <- var(X) * (n-1)/n
               		XVar <- diag(Configurations[[1]]$mleSigE)
			CnstVar <- XVar[XVar<tol^2]
			if (length(CnstVar)==1)
				stop("Variable ",names(CnstVar)," appears to be constant \n")
			else if (length(CnstVar)>0) 
				stop("Variables ",paste(names(CnstVar),collapse=" ")," appear to be constant\n")
			Configurations[[1]]$mleSigEse  <- sqrt((Configurations[[1]]$mleSigE^2+outer(XVar,XVar))/n)
			dimnames(Configurations[[1]]$mleSigEse) <- dimnames(Configurations[[1]]$mleSigE)		
	  }	
	  else {
		SigNEC1 <- var(X) * (n-1)/n
		XVar <- diag(SigNEC1)
		CnstVar <- XVar[XVar<tol^2]
		if (length(CnstVar)==1)
			stop("Variable ",names(CnstVar)," appears to be constant \n")
		  else if (length(CnstVar)>0)
			stop("Variables ",paste(names(CnstVar),collapse=" ")," appear to be constant\n")
		  SigNEseC1 <- sqrt((SigNEC1^2+outer(XVar,XVar))/n)
		  dimnames(SigNEseC1) <- dimnames(SigNEC1)
	  }
	}
	else if (Type=="HomMxt") {
		mleNmuE <- apply(X,2,function(x) tapply(x,grouping,mean))
		mleNmuEse <- apply(X,2,function(x) tapply(x,grouping,sd))/sqrt(nk)
		rownames(mleNmuE) <- rownames(mleNmuEse) <- levels(grouping)
		if (is.element(1,Config)) {
			Configurations[[1]]$mleSigE <- var(X[grouping==levels(grouping)[1],]) * (nk[1]-1)/n
			for (g in 2:k) Configurations[[1]]$mleSigE <- 
				Configurations[[1]]$mleSigE + var(X[grouping==levels(grouping)[g],]) * (nk[g]-1)/n
			XVar <- diag(Configurations[[1]]$mleSigE)
			CnstVar <- XVar[XVar<tol^2]
			if (length(CnstVar)==1) 
				stop("Variable ",names(CnstVar)," appears to be constant within groups\n")
			else if (length(CnstVar)>0)
				stop("Variables ",paste(names(CnstVar),collapse=" ")," appear to be constant within groups\n")
			Configurations[[1]]$mleSigEse  <- sqrt((Configurations[[1]]$mleSigE^2+outer(XVar,XVar))/n)
			dimnames(Configurations[[1]]$mleSigEse) <- dimnames(Configurations[[1]]$mleSigE)
		}	
		else {
			SigNEC1 <- var(X[grouping==levels(grouping)[1],]) * (nk[1]-1)/n
			for (g in 2:k) SigNEC1 <- SigNEC1 + var(X[grouping==levels(grouping)[g],]) * (nk[g]-1)/n
			XVar <- diag(SigNEC1)
			CnstVar <- XVar[XVar<tol^2]
			if (length(CnstVar)==1)
				stop("Variable ",names(CnstVar)," appears to be constant within groups\n")
			else if (length(CnstVar)>0)
				stop("Variables ",paste(names(CnstVar),collapse=" ")," appear to be constant within groups\n")
			SigNEseC1  <- sqrt((SigNEC1^2+outer(XVar,XVar))/n)
			dimnames(SigNEseC1) <- dimnames(SigNEC1)
		}
	}
	Cnf35 <- Config[Config>2]
	if (length(Cnf35)!=0)  
	for (Conf in Cnf35)  {
		if (is.element(1,Config)) {
			Configurations[[Conf]]$mleSigE <- mleNSigC35(Conf,Configurations[[1]]$mleSigE,q,p,0.)
			Configurations[[Conf]]$mleSigEse <- mleNSigC35(Conf,Configurations[[1]]$mleSigEse,q,p,NA)
		}
		else {
			Configurations[[Conf]]$mleSigE <- mleNSigC35(Conf,SigNEC1,q,p,0.)
			Configurations[[Conf]]$mleSigEse <- mleNSigC35(Conf,SigNEC1,q,p,NA)
		}
	}
	if (is.element(2,Config)) {
		Xsd <- sqrt(XVar)
		lglikdif <- -n*sum(log(Xsd)) 
		if (Type=="SingDst") Xscld <- scale(X,scale=Xsd)
		else if (Type=="HomMxt") {
			Xscld <- scale(X[grouping==levels(grouping)[1],],center=mleNmuE[1,],scale=Xsd)
			for (g in 2:k) Xscld <- 
				rbind(Xscld,scale(X[grouping==levels(grouping)[g],],center=mleNmuE[g,],scale=Xsd))
		}
		C2res <- Cnf2MaxLik(Xscld)  # To do: add possibility of changing the default parameters of Cnf2MaxLik!!
 		Configurations[[2]]$mleSigE <- C2GetCov(C2res$SigmaSr,outer(Xsd,Xsd),q)
		Configurations[[2]]$mleSigEse <- C2GetCovStderr(Configurations[[2]]$mleSigE,X,q,ue=colMeans(X))
		dimnames(Configurations[[2]]$mleSigEse) <- dimnames(Configurations[[2]]$mleSigE)
		logLiks[2] <- Configurations[[2]]$logLik <- C2res$lnLik + lglikdif 
		AICs[2] <- Configurations[[2]]$AIC <- -2*Configurations[[2]]$logLik + 2*npar(2,p,q,Ngrps=k)
		BICs[2] <- Configurations[[2]]$BIC <- -2*Configurations[[2]]$logLik + log(n)*npar(2,p,q,Ngrps=k)
		Configurations[[2]]$optres <- C2res$optres
	}
	if (any(Config!=2)) {
		if (Type=="SingDst") Xdev <- scale(X,scale=FALSE)
		else if (Type=="HomMxt") {
			Xdev <- scale(X[grouping==levels(grouping)[1],],center=mleNmuE[1,],scale=FALSE)
			for (g in 2:k) Xdev <- 
				rbind(Xdev,scale(X[grouping==levels(grouping)[g],],center=mleNmuE[g,],scale=FALSE))
		}
		Conf134 <- Config[Config !=2 & Config!=5]
   		if (length(Conf134)!=0)  for (Conf in Conf134) {
			logLiks[Conf] <- Configurations[[Conf]]$logLik <-  sum( apply(Xdev, 1 , ILogLikNC ,
				SigmaInv=solve(Configurations[[Conf]]$mleSigE) , const=-0.5*(p*log(2*pi) + 
				as.double(determinant(Configurations[[Conf]]$mleSigE,logarithm=TRUE)$modulus))) ) 
			AICs[Conf] <- Configurations[[Conf]]$AIC <- -2*Configurations[[Conf]]$logLik + 2*npar(Conf,p,q,Ngrps=k)
			BICs[Conf] <- Configurations[[Conf]]$BIC <- -2*Configurations[[Conf]]$logLik + log(n)*npar(Conf,p,q,Ngrps=k)
		}
		if (is.element(5,Config)) {
			logLiks[5] <- Configurations[[5]]$logLik <- 
				sum( apply(Xdev , 1 , ILogLikDNC , IVar=1./XVar , const=-0.5*(p*log(2*pi)+sum(log(XVar)))) )
			AICs[5] <- Configurations[[5]]$AIC <- -2*Configurations[[5]]$logLik + 2*npar(5,p,q,Ngrps=k)
			BICs[5] <- Configurations[[5]]$BIC <- -2*Configurations[[5]]$logLik + log(n)*npar(5,p,q,Ngrps=k)
		}
	} 
	if (SelCrit=="AIC") bestmod <- which.min(AICs)
	else if (SelCrit=="BIC") bestmod <- which.min(BICs)

	if (Type=="SingDst")
		new("IdtSngNDE",ModelNames=names(AICs),ModelType=rep("Normal",5),ModelConfig=1:5,
			  NIVar=q,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,Configurations=Configurations,SelCrit=SelCrit,
			  logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod)
	else
		new("IdtMxNDE",ModelNames=names(AICs),ModelType=rep("Normal",5),ModelConfig=1:5,NIVar=q,grouping=grouping,
			  Hmcdt=TRUE,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,Configurations=Configurations,SelCrit=SelCrit,
			  logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod)
}

IdtHetMxtNmle <- function(Idt,grouping,Config,SelCrit,tol=1.0e-4)
{
	n <- Idt@NObs
     	q <- Idt@NIVar
     	p <- 2*q
     	nk <- as.numeric(table(grouping))
     	k <- length(nk) 
	if (k==1) stop("The data belongs to one single group. A partition into at least two different groups is required\n")
	mleNmuE <- matrix(nrow=k,ncol=p)
	mleNmuEse <- matrix(nrow=k,ncol=p)
	rownames(mleNmuE) <- rownames(mleNmuEse) <- levels(grouping)
	colnames(mleNmuE) <- colnames(mleNmuEse) <- names(cbind(Idt@MidP,Idt@LogR))
   	logLiks <- AICs <- BICs <- rep(NA_real_,5)
	Xnams <- colnames(mleNmuE)
        anams <- list(Xnams,Xnams,levels(grouping))
   	Configurations <- vector("list",5)
	names(logLiks)[1:5] <- names(AICs)[1:5] <- names(BICs)[1:5] <- names(Configurations)[1:5]<- paste("NC",1:5,sep="")
	for (model in Config)  { 
		if (model==2) Configurations[[2]] <- list(mleSigE=array(dim=c(p,p,k),dimnames=anams),mleSigEse=array(dim=c(p,p,k),dimnames=anams),
								logLik=0.,AIC=NULL,BIC=NULL,optres=vector("list",k))
		else Configurations[[model]] <- list(mleSigE=array(dim=c(p,p,k),dimnames=anams),mleSigEse=array(dim=c(p,p,k),dimnames=anams),
							logLik=0.,AIC=NULL,BIC=NULL)
	}
        if (is.element(2,Config)) names(Configurations[[2]]$optres) <- levels(grouping)
	for (g in 1:k) {
		Idtg <- Idt[grouping==levels(grouping)[g],]
		Xstdev <- sapply(cbind(Idtg@MidP,Idtg@LogR),sd)
		Cnststdev <- Xstdev[Xstdev<tol]
		if (length(Cnststdev)==1) 
			stop("Variable ",names(Cnststdev)," appears to be constant in group",levels(grouping)[g],"\n")
		else if (length(Cnststdev)>1)
			stop("Variables ",paste(names(Cnststdev),collapse=" ")," appear to be constant in group ",levels(grouping)[g],"\n")
		pres <- IdtNmle(Idtg,grouping,Type="SingDst",Config=Config,SelCrit=SelCrit,tol=tol)
		mleNmuE[g,] <- pres@mleNmuE 
		mleNmuEse[g,] <- pres@mleNmuEse 
	   	for (model in Config)  { 
			Configurations[[model]]$mleSigE[,,g] <- pres@Configurations[[model]]$mleSigE
			Configurations[[model]]$mleSigEse[,,g] <- pres@Configurations[[model]]$mleSigEse
			Configurations[[model]]$logLik <- Configurations[[model]]$logLik + pres@Configurations[[model]]$logLik
	        	if (model==2) Configurations[[2]]$optres[[g]] <- pres@Configurations[[2]]$optres
		}            
	}
	for (model in Config)  {
		logLiks[model] <- Configurations[[model]]$logLik 
		AICs[model] <- Configurations[[model]]$AIC <- -2*Configurations[[model]]$logLik + 2*npar(2,p,q,Ngrps=k,Mxt="Het")
		BICs[model] <- Configurations[[model]]$BIC <- -2*Configurations[[model]]$logLik + log(n)*npar(2,p,q,Ngrps=k,Mxt="Het")
	}
  	if (SelCrit=="AIC") bestmod = which.min(AICs)
  	else if (SelCrit=="BIC") bestmod = which.min(BICs)
 
  	new("IdtMxNDE",ModelNames=names(AICs),ModelType=rep("Normal",5),ModelConfig=1:5,grouping=grouping,
			Hmcdt=FALSE,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,Configurations=Configurations,SelCrit=SelCrit,NIVar=q,
			logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod)
}


mleNSigC35 <- function(Conf,mat,q,p,defval)
{
  Newmat <- matrix(defval,p,p)
  dimnames(Newmat) <- dimnames(mat)
  if (Conf==3) for (i in 1:q) Newmat[c(i,q+i),c(i,q+i)] <- mat[c(i,q+i),c(i,q+i)]  
  if (Conf==4) {
	Newmat[1:q,1:q] <- mat[1:q,1:q]
	Newmat[(q+1):(2*q),(q+1):(2*q)] <- mat[(q+1):(2*q),(q+1):(2*q)]
  }
  if (Conf==5) diag(Newmat) <- diag(mat)
  Newmat  #  return(Newmat)
}

ILogLikNC <- function(Xdev,SigmaInv,const)   const -0.5 * Xdev%*%SigmaInv%*%Xdev
ILogLikNC1 <- function(Xdev,SigmaSrInv,const)   const -0.5 * sum((SigmaSrInv%*%Xdev)^2)
ILogLikDNC <- function(Xdev,IVar,const)  const -0.5 * sum(Xdev^2 * IVar)

npar <- function(Conf,p,q,Ngrps=1,Mxt=c("Hom","Het"))
{
   Mxt <- match.arg(Mxt)
   if (Mxt=="Hom") {
	if (Conf==1)  return(Ngrps*p + p*(p+1)/2)
	if (Conf==2)  return(Ngrps*p + q*(q+1) + p)
	if (Conf==3)  return(Ngrps*p + p + 3*q)
	if (Conf==4)  return(Ngrps*p + q*(q+1))
	if (Conf==5)  return(Ngrps*p + p)
  }
  else if (Mxt=="Het") {
	if (Conf==1)  return(Ngrps*(p + p*(p+1)/2))
	if (Conf==2)  return(Ngrps*(p + q*(q+1) + p))
	if (Conf==3)  return(Ngrps*(p + p + 3*q))
	if (Conf==4)  return(Ngrps*(p + q*(q+1)))
	if (Conf==5)  return(Ngrps*(p + p))
  }
}

NestedBy <- function(Model,ModelType=c("Normal","SKNormal","NrmandSKN"))
{
   ModelType <- match.arg(ModelType)
   if (ModelType=="Normal" || ModelType=="SKNormal")
   {
	if (Model==1)  stop("Configuration 1 is the most general model in this analysis\n")
	if (Model==2)  return(1)
	if (Model==3 || Model==4)  return(1:2)
	if (Model==5)  return(1:4)
  }
  else if (ModelType=="NrmandSKN")
  {
	if (Model==1)  return(6)
	if (Model==2)  return(c(1,6:7))
	if (Model==3)  return(c(1:2,6:8))
	if (Model==4)  return(c(1:2,6:7,9))
	if (Model==5)  return(c(1:4,6:10))
	if (Model==6)  stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n")
	if (Model==7)  return(6)
	if (Model==8 || Model==9)  return(6:7)
	if (Model==10)  return(6:9)
  }
}

NextModel <- function(Model,ModelType=c("Normal","SKNormal","NrmandSKN"))
{
   ModelType <- match.arg(ModelType)
   if (ModelType=="Normal" || ModelType=="SKNormal")
   {
	if (Model==1)  stop("Configuration 1 is the most general model in this analysis\n")
	if (Model==2)  return(1)
	if (Model==3 || Model==4)  return(2)
	if (Model==5)  return(3:4)   # (4:3) ?
  }
  else if (ModelType=="NrmandSKN")
  {
	if (Model==1)  return(6)
	if (Model==2)  return(c(1,7))  # c(7,1) ?
	if (Model==3)  return(c(2,8))
	if (Model==4)  return(c(2,9))
	if (Model==5)  return(c(3:4,10))
	if (Model==6)  stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n")
	if (Model==7)  return(6)
	if (Model==8 || Model==9)  return(7)
	if (Model==10)  return(8:9)
  }
}

