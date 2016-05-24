# calculates BIC and AIC relative to a saturated loglinear model
# rather than relative to a null model
fitmacro<-function(object) {
	stopifnot(class(object)[1]=="glm",object$family$family=="poisson",object$family$link=="log")
	ncases<-sum(object$y)
	bic<-object$deviance-object$df.residual*log(ncases)
	aic<-object$deviance-object$df.residual*2
	cat("\n","\n")
	cat("deviance:            ",formatC(object$deviance,    width = 12, digits = 3, format = "f"),"\n")
	cat("df:                  ",formatC(object$df.residual, width = 12, digits = 0, format = "f"),"\n")
	cat("bic:                 ",formatC(bic,                width = 12, digits = 3, format = "f"),"\n")
	cat("aic:                 ",formatC(aic,                width = 12, digits = 3, format = "f"),"\n")
	cat("Number of parameters:",formatC(object$rank,        width = 12, digits = 0, format = "f"),"\n")
	cat("Number of cases:     ",formatC(ncases,             width = 12, digits = 0, format = "f"),"\n")
	cat("\n","\n")
}

# utility function to check if the variables are factors
# with the same number of categories
# called by functions for mobility models below
check.square <- function(rowvar,colvar,equal=TRUE) {
	stopifnot(is.factor(rowvar))
	stopifnot(is.factor(colvar))
	if (equal) {
		stopifnot(nlevels(rowvar)==nlevels(colvar))
	}
}

# Quasi-independence
mob.qi <- function(rowvar,colvar,constrained=FALSE,print.labels=FALSE) {
	check.square(rowvar,colvar)
	if (constrained) {
		qi <- ifelse(rowvar==colvar, 1, 0)
		nms<-c("diagonal")
	}
	else {
		qi <- ifelse(rowvar==colvar, rowvar, 0)
		nms<-levels(rowvar)
	}

	qi<-factor(qi)
	qi<-C(qi,contr.treatment,base=1)
	if (print.labels) {
		levels(qi)<-c("offdiag",nms)
	}
	qi
}

# symmetric interaction effects
mob.symint <- function(rowvar,colvar,print.labels=FALSE) {
	check.square(rowvar,colvar)
	# remove factor levels to avoid messy output
	if (!print.labels) {
		attr(rowvar,"levels")<-1:nlevels(rowvar)
		attr(colvar,"levels")<-1:nlevels(colvar)
	}
	mdl<-model.matrix(~rowvar*colvar)
	intrct<-mdl[,attr(mdl,"assign")==3]
	# remove factor names
	colnames(intrct)<-sub("rowvar","",colnames(intrct))
	colnames(intrct)<-sub("colvar","",colnames(intrct))
	w<-ncol(intrct)
	x<-matrix(1:w,sqrt(w),sqrt(w))
	symint<-intrct[,t(x)[lower.tri(x,diag=TRUE)]]+intrct[,x[lower.tri(x,diag=TRUE)]]
	symint
}

# equal main effects, Hope's halfway model
mob.eqmain <- function(rowvar,colvar,print.labels=FALSE) {
	check.square(rowvar,colvar)
	if (!print.labels) {
		attr(rowvar,"levels")<-1:nlevels(rowvar)
		attr(colvar,"levels")<-1:nlevels(colvar)
	}
	rmat<-model.matrix(~rowvar)
	rmat<-rmat[,attr(rmat,"assign")==1]
	cmat<-model.matrix(~colvar)
	cmat<-cmat[,attr(cmat,"assign")==1]
	eqmain<-rmat+cmat
	colnames(eqmain)<-sub("rowvar","",colnames(eqmain))
	eqmain
}
# Crossings parameter
mob.cp <- function(rowvar,colvar) {
	check.square(rowvar,colvar)
	cp<-NULL
	rvar<-as.numeric(rowvar)
	cvar<-as.numeric(colvar)
	for (i in 1:(nlevels(rowvar)-1)) {
		cp<-cbind(cp,as.numeric((rvar <= i & cvar > i) | (rvar > i & cvar <= i)))
	}
	cp
}

# Uniform association
mob.unif <- function(rowvar,colvar) {
	check.square(rowvar,colvar,equal=FALSE)
	as.numeric(rowvar)*as.numeric(colvar)
}

# RC model 1 (unequal row and column effects, page 58)
# Fits a uniform association parameter and row and column effect
# parameters. Row and column effect parameters have the
# restriction that the first and last categories are zero.
mob.rc1 <- function(rowvar,colvar,equal=FALSE,print.labels=FALSE) {
	# the number of row and column categories need not be equal
	# unless an equality restriction is applied
	check.square(rowvar,colvar,equal=equal)
	# use numbers rather than factor levels by default
	if (!print.labels) {
		attr(rowvar,"levels")<-1:nlevels(rowvar)
		attr(colvar,"levels")<-1:nlevels(colvar)
	}

	# row effects, first and last category constrained to 0
	# multiplied by column variable as continuous
	rowvar<-C(rowvar,contr.treatment,base=1)
	rmat<-model.matrix(~rowvar)
	rmat<-rmat[,attr(rmat,"assign")==1]
	rmat<-rmat[,1:(ncol(rmat)-1)]*as.numeric(colvar)

	# column effects, same construction
	colvar<-C(colvar,contr.treatment,base=1)
	cmat<-model.matrix(~colvar)
	cmat<-cmat[,attr(cmat,"assign")==1]
	cmat<-cmat[,1:(ncol(cmat)-1)]*as.numeric(rowvar)

	u<-mob.unif(rowvar,colvar)
	# add a name for the uniform association effect
	dim(u)<-c(length(u),1)
	colnames(u)<-c("U")

	if (equal) {
		# rowmain and colmain are added to impose an equality restriction
		colnames(rmat)<-sub("rowvar","RC",colnames(rmat))
		rc1<-cbind(rmat+cmat,u)
	}
	else {
		colnames(rmat)<-sub("rowvar","R",colnames(rmat))
		colnames(cmat)<-sub("colvar","C",colnames(cmat))
		rc1<-cbind(rmat,u,cmat)
	}
	rc1
}
