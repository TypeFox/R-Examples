
CPMCGLM<-function(formula,family,link,data,varcod,dicho,nb.dicho,categ,nb.categ,boxcox,nboxcox,N=1000,cutpoint)
{
if(missing(dicho) & missing(nb.dicho) & missing(categ) & missing(nb.categ) & missing(boxcox) & missing(nboxcox) & missing(cutpoint) ){stop("You need to enter at least one transformation in the function")}

#A modifier plus tard
##AD & JR
if(!missing(dicho) & !missing(nb.dicho))stop("Specify only one argument for the dichotomous transformation ('dicho' or 'nb.dicho')")
if(!missing(categ) & !missing(nb.categ))stop("Specify only one argument for the categorical transformation ('categ' or 'nb.categ')")
if(!missing(boxcox) & !missing(nboxcox))stop("Specify only one argument for the boxcox transformation ('boxcox' or 'nboxcox')")

continuous<-"FALSE"
if (missing(boxcox))
{
	if(!missing(nboxcox)){
		if((nboxcox > 5)){stop("The nboxcox argument must not be upper than 5")}
		if(nboxcox<6){
			if (nboxcox==0) {boxcox<-NULL}
			if (nboxcox==1) {boxcox<-1}
			if (nboxcox==2) {boxcox<-c(1,0)}
			if (nboxcox==3) {boxcox<-c(1,0,2)}
			if (nboxcox==4) {boxcox<-c(1,0,2,0.5)}
			if (nboxcox==5) {boxcox<-c(1,0,2,0.5,1.5)}
		}
	}
}


if(missing(dicho) & (!missing(nb.dicho))){
dicho<-c((1:nb.dicho)/(nb.dicho+1))
}


if(missing(categ) & (!missing(nb.categ))){
	categ<-matrix(NA,ncol=nb.categ+1,nrow=nb.categ)
	for(i in 1:nb.categ)
	{
		categ[i,1:(i+1)]<-c(1:(1+i))/(i+2)
	}
}

if(missing(categ) & (!missing(dicho))){
	quantile<-matrix(NA,ncol=1,nrow=length(dicho))
	for( i in 1:length(dicho))
	{
		quantile[i,1]<-dicho[i]
	}
}
if(missing(categ) & (missing(dicho))){quantile<-NULL}

if(!missing(categ) & (missing(dicho))){quantile<-categ}

if(!missing(categ) & (!missing(dicho))){
binc<-matrix(NA,ncol=ncol(categ),nrow=length(dicho))
for( i in 1:length(dicho))
{
binc[i,1]<-dicho[i]
}
quantile<-rbind(binc,categ)
}



	
	call <- match.call()
	if(!(family %in% c("binomial","gaussian","poisson"))) stop("This family of distribution is not taken into account in this function")
	family1<- switch(family,
		"binomial"=binomial,
		"gaussian"=gaussian ,
		"poisson"=poisson
	)

#==========================  14/02/2012 =========================	
	
	if(missing(data)) stop("Missing data argument.")
	if(missing(family)) stop("Missing family argument.")
	if(missing(link)) stop("Missing link argument.")
	if(missing(varcod)) stop("Missing varcod argument.")	
	if(missing(dicho) & missing(boxcox) & missing(categ) & missing(nb.dicho) & missing(nboxcox) & missing(nb.categ) & missing(cutpoint)){
		 stop("No transformation is available")
	}
	
	if(class(formula)!="formula") stop("The argument formula must be a formula.")
	if(class(data)!="data.frame") stop("The argument must be a data frame.")
	if(class(varcod)!="character") stop("The argument varcod must be a character.")
	if(!missing(boxcox)){
		if(is.vector(boxcox)=="FALSE") stop("The argument boxcox must be a vector.")
	}

# Prise en compte de la variable formula
	m <- match.call(expand.dots=FALSE)
	m$family <- m$link <- m$varcod <- m$quantile <- m$dicho <- m$nb.dicho<- m$nb.categ<- m$categ<- m$continuous <- m$boxcox<- m$nboxcox <-m$cutpoint <- m$N <- NULL

	m[[1]] <- as.name("model.frame")
	m <- eval(m) 

# var dependante
	NamesY <- all.names(update(formula,"~1"))[2]
	Y <- as.vector(model.response(model.frame(formula=update(formula,"~1"),data=data)))
	
	
# var expli	
	mat.exp <- if (!is.empty.model(attr(m,"terms")))model.matrix(attr(m,"terms"),m, contrasts)	
	mtt <- model.matrix(attr(m,"terms"),m, contrasts)
	Nom1 <- attributes(mtt)$dimnames[[2]]
	
	type <- paste("factor\\(",varcod,"\\)",sep="")
	if(length(grep(type,Nom1))!= 0){
		stop("The coding variable needs to be continuous")
	}

	if(ncol(mat.exp)==1)stop("the dataset must be contain the varcod variable.")
	mat.exp <- mat.exp[,-1]

# add
	factor.names <- function(x){
		x <- matrix(x,nrow=1)
		Names <- apply(x,MARGIN=2,FUN=function(x){
			if(length(grep("factor",x))!= 0){
				pos1 <- grep("\\(",unlist(strsplit(x,split="")))+1
				pos2 <- grep("\\)",unlist(strsplit(x,split="")))-1
				compris.factor <- substr(x,start=pos1,stop=pos2)
				after.factor <- substr(x,start=(pos2+2),stop=length(unlist(strsplit(x,split=""))))
				paste(compris.factor,after.factor,sep=".")
			}else{
				x
			}
		}
		)
		return(Names)
	}	
	
	if(is.matrix(mat.exp)){
		colnames(mat.exp) <-factor.names(colnames(mat.exp))
	}
	if(is.vector(mat.exp)){
		Z <- NULL
		nb <- 0
		namesZ <- NULL
		var.cod <- as.vector(mat.exp)
		ind.ajust <- 0
	}

# matrice de variable d'ajustement: Z
	if(!is.vector(mat.exp)){
		if(!(varcod %in% colnames(mat.exp)))stop("varcod argument it is not present in the dataset.")
		ind.ajust <- 1
		if(ncol(mat.exp)!=2){
			Z <- mat.exp[,-grep(varcod,colnames(mat.exp))]
			nb <- ncol(Z)
			namesZ <- NULL
			names1 <- unlist(colnames(Z)) 
		
			for(i in 1:nb){
				if(i < nb) namesZ <- paste(namesZ,names1[i],"+")
				else namesZ <- paste(namesZ,colnames(Z)[i])
				
			}
		}else{
			Z <- as.matrix(mat.exp[,-grep(varcod,colnames(mat.exp))])
			nb <- 1
			namesZ <- NULL
			names1 <- as.matrix(colnames(mat.exp))
			posi<-grep(varcod,names1)
			namesZ <-names1[-posi,]	
			colnames(Z)<-namesZ					

		}	
		
		var.cod <- as.vector(mat.exp[,varcod])
	}
# variable a coder

#===============================================================
	n<-nrow(data)
	coup<-NULL
	X3<-Codage(data,quantile,var.cod)[[1]]
	coup2<-Codage(data,quantile,var.cod)[[2]]
	X2<-Boxcox(boxcox,var.cod)[[1]]
	coup1<-Boxcox(boxcox,var.cod)[[2]]
	X4<-codcut(data,cutpoint,var.cod)[[1]]
	coup4<-codcut(data,cutpoint,var.cod)[[2]]
	
	
	if(continuous=="TRUE"){
		X1<-cbind(X3,X4,var.cod,X2)
		coup3<-1
	}else{
		X1<-cbind(X3,X4,X2)
		coup3<-NULL
	}
	
	
	coup<-c(coup2,coup4,coup3,coup1)

	if(ind.ajust==1){
		formula1 <-update(formula,paste("~",namesZ))
	}else{
		formula1 <-update(formula,paste("~",1))
	}

	data <- as.data.frame(cbind(Y,mat.exp))
	colnames(data)[1] <- NamesY

	t.obs<-test.score(formula=formula1,data=data,codage=X1,Z=Z,Y=Y,family1=family1,family=family,link1=link,ind.ajust=ind.ajust)
	
	pval<-p.val<-NULL
	for (i in 1:ncol(X1)){
		pval[i]<-(1-pchisq(t.obs[i],coup[i]))
		p.val[i]<- min(p.val[i-1],pval[i])
	}

	if(sum(is.na(p.val))!=0){stop("At least one of the coding, it isn't adapt for the model specification")}
#Initialisation des param?tres
	pval.exact1<-pval.bonf1<-pval.naive1<-NULL
	control<-rep(1,ncol(X1))


if(all.equal(coup,control)==T){
pval.exact1<-test.liquet(formula=formula1,data=data,codage=X1,Z=Z,Y=Y,family=family,link1=link,family1=family1)}
else {pval.exact1<-"Correction not available for these codings"}
	
	#Reechantillonage
	###Initialisation des parametres
	Y1<-matrix(ncol=N,nrow=n)
	tboot1<-matrix(nrow=N,ncol=ncol(X1))
	alphaboot<-NULL
	alphapermut<-NULL
	tpermut1<-matrix(nrow=N,ncol=ncol(X1))
	#bootstrap
      #Coefficient sous H0
       link0 <- link
	beta<-glm(formula1,family=family1(link0),data)$coeff
	res<-glm(formula1,family=family1(link0),data)$res
#	design <- model.matrix(~data$bmi+data$age)
	design <- cbind(rep(1,n),Z)

	#different pour chaque fonction
	param<- switch(link,
	"logit"= exp(design%*%beta)/(1+exp(design%*%beta)),
	"probit"= pnorm(design%*%beta),
	"log"=exp(design%*%beta),
	"identity"=(design%*%beta)
	)


	tboot<-matrix(ncol=ncol(X1),nrow=N)
	pval.boot<-matrix(ncol=ncol(X1),nrow=N)

	#Simulation des Y1

	progress_bar_tk <- create_progress_bar (title="step 1 :Parametric Bootstrap resampling","tk")
	progress_bar_tk $ init (N)
	
	for (j in 1:N){
		progress_bar_tk $ step ()
		#Pour une regression logistique
		Y1[,j]<- switch(family,
		"gaussian"= rnorm(n,mean=mean(param),sd=sd(res)) ,
		"binomial"= rbinom(n,1,prob=param) ,
		"poisson"= rpois(n,mean(param)),
		)
		#if (family=="binomial"){
		#Y1[,j] <- rbinom(n,1,prob=param)}
		if(ind.ajust==1){
			f1<-paste("~",namesZ)
		}else{
			f1<-"~1"
		}
		
		tboot[j,]<-tmultboot(f1,data,Y1=Y1[,j],Z,X1,family1,family=family,link1=link)
			for (i in 1:ncol(X1)){
			pval.boot[j,i]<-(1-pchisq(tboot[j,i],coup[i]))
			}
	}
	cat("\n")
	progress_bar_tk $ term ()
	#Permutation
	tpermut<-matrix(ncol=ncol(X1),nrow=N)
	pval.permut<-matrix(ncol=ncol(X1),nrow=N)
# barre de progression

		progress_bar_tk <- create_progress_bar (title="step 2 : Permutation resampling","tk")
		progress_bar_tk $ init (N)
		
		for (k in 1:N){
		
			data1<-NULL

			progress_bar_tk $ step ()
			
			ind.permut <- permut(length(Y))
			
			if(ind.ajust==1){
				if(ncol(Z)==1)
				{
					Z11<-as.matrix(Z[ind.permut,])
					colnames(Z11)<-namesZ
					formula2 <-formula(paste("Y","~",namesZ))
					data.permut <- data.frame(Y=Y[ind.permut],varcod=var.cod,Z11)
					tpermut[k,]<-test.score(formula=formula2,data=data.permut,codage=X1,Z=Z11,Y=Y[ind.permut],family1,family=family,link1=link,ind.ajust=ind.ajust)
					}
					else{
					formula2 <-formula(paste("Y","~",namesZ))
					data.permut <- data.frame(Y=Y[ind.permut],varcod=var.cod,Z[ind.permut,])
					tpermut[k,]<-test.score(formula=formula2,data=data.permut,codage=X1,Z=Z[ind.permut,],Y=Y[ind.permut],family1,family=family,link1=link,ind.ajust=ind.ajust)

				}
			}else{
				formula2<-formula(paste("Y","~",1))
				data.permut <- data.frame(Y=Y[ind.permut],varcod=var.cod)

				tpermut[k,]<-test.score(formula=formula2,data=data.permut,codage=X1,Z=Z,Y=Y[ind.permut],family1=family1,family=family,link1=link,ind.ajust=ind.ajust)
			}


			for (i in 1:ncol(X1)){
				pval.permut[k,i]<-(1-pchisq(tpermut[k,i],coup[i]))
			}	
	}
	cat("\n")	
	progress_bar_tk $ term ()
	pval.boot2<-pval.permut2<-matrix(0,ncol=ncol(X1),nrow=N)

		
	for (i in 1:ncol(X1))
	{

		for (j in 1:N){
				pval.boot[is.na(pval.boot)] <- 9999
				pval.permut[is.na(pval.permut)] <- 9999
				if (i>1){
					pval.boot2[j,i]<-min(pval.boot2[j,i-1],min(pval.boot[j,i]))
					pval.permut2[j,i]<-min(pval.permut2[j,i-1],min(pval.permut[j,i]))
				}
				else{
					pval.boot2[j,i]<-min(pval.boot[j,i])
					pval.permut2[j,i]<-min(pval.permut[j,i])
				}
		}

	}

	cat("\n")

	
	for (i in 1:ncol(X1)){
		#Calcul des diff?rentes pvaleur ajust?es
		pval.bonf<-p.val[i]*i
		pval.bonf1<-c(pval.bonf1,pval.bonf)
		pval.naive<-p.val[i]
		pval.naive1<-c(pval.naive1,pval.naive)
	
		#Bootstrap
		#Calcul de la proportion de p-valeur superieur a la p-valeur sur le jeu de donnees initiale
		for (ii in 1:N){
			if(pval.boot2[ii,i]< p.val[i]){tboot1[ii,i]<-1} 		
		}

		alphaboot[i]<-(length(which(tboot1[,i]>0))/N)
		

		#Permutation
		#Calcul de la proportion de p-valeur superieur a la p-valeur sur le jeu de donnees initiale

		for (l in 1:N){
			if(pval.permut2[l,i]< p.val[i]){tpermut1[l,i]<-1}		
		}
		alphapermut[i]<-(length(which(tpermut1[,i]>0))/N)
	}
	cat("\n")
####INFORMATIONS COMPLEMENTAIRES####
#environment(transf)<-environment()

#trans<-transf(continuous=continuous,quantile=quantile,boxcox=boxcox,cutpoint=cutpoint)
	trans<-transf(quantile,continuous,boxcox,cutpoint)
#nombre de variable d'ajustement: A changer plus tard
	adj<-ncol(Z)
#Nombre de transformation par type
	if(continuous=="TRUE"){q<-"TRUE"}else{q<-"FALSE"}
	if(class(quantile)=="NULL"){nbq=0}else{nbq=nrow(quantile)}
	if(missing(cutpoint)){nbc=0}else{nbc=nrow(cutpoint)}
	if(missing(boxcox)){nbb=0}else{nbb=length(boxcox)}
	#if(pval.naive1[1]-alphapermut[1]>0.05){alphapermut="This resampling method it is not adapt for this model"}
#Meilleure transformation
	BC<-trans[which.min(pval.naive1)]
#Meilleure codage

	bestcod<-bestcod(continuous=continuous,quantile=quantile,boxcox=boxcox,pval.naive=pval.naive1,cutpoint=cutpoint)
	bcod<-which(is.na(bestcod)==TRUE)
	if( (is.numeric(bestcod)) & (length(bcod!=0)) ) {bestcod<-bestcod[-bcod]}


	info<-list("link"=link,"nbt"=ncol(X1),"nbb"=nbb,"nbq"=nbq,"nbcont"=q,"adj"=adj,"trans"=trans,"BC"=BC,"bestcod"=bestcod,"nbc"=nbc)
	
	if(substr(BC,1,8)=="Quantile"){
	vq<-quantile(var.cod,bestcod)
	}else{vq<-NULL}

	out <- NULL
	
	out$call <- call
## Complementary information
	out$n <- n
	out$N <- N
	out$family <- family
	out$link <- link
	out$nbt <- ncol(X1)
	out$nbb <- nbb
	out$nbq <- nbq
	out$nbcont <- q
	out$nbc <- nbc
	out$adj <- adj
	out$vq<-vq
	out$trans <- trans
	out$BC <- BC
	out$bestcod <-bestcod 

	##
	out$naive.pvalue <- pval.naive1[ncol(X1)]
	if(is.numeric(pval.exact1)==TRUE)
	{out$exact.pvalue <- pval.exact1[ncol(X1)]}
	else{out$exact.pvalue<-pval.exact1}
	out$bonferroni.adjusted.pvalue <- pval.bonf1[ncol(X1)]
	out$parametric.bootstrap.adjusted.pvalue <- alphaboot[ncol(X1)]
	out$permutation.adjusted.pvalue <- alphapermut[ncol(X1)]
		
	class(out) <- c("CPMCGLM")
	
	out

}
## end

