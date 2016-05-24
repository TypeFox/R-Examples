
#################### widelly taken from gt{globaltest}
.getXY <- function(Y,X,Z,data){
	call <- match.call()
  # # avoid conflict between "levels" input and "levels" function
  # if (missing(levels)) levels <- NULL
                       
  # data default
  if (missing(data) || is.null(data))
    if(is.data.frame(Y) | (is.matrix(Y)))
	  data <- Y else data <- NULL
  
  if (is.matrix(data))  
    data <- data.frame(data)

  # evaluate Y, which may be one of the colnames of data
  Y <- eval(call$Y, data, parent.frame())

  # settle Z, X and Y if Y is a formula
  if (missing(Z) || is.null(Z)) {
    if (!is.null(dim(Y)[1]))
      Z <- rep(0,dim(Y)[1])  else 	
	  if(is.data.frame(data)) Z <- ~0 else  stop("argument \"Z\" is missing, with no default")  
  }  
  if (missing(X) || is.null(X))  ########livio:  ?
    if (is(Y, "formula")){ 
		#if it is a left+right formula
		if(length(Y)==3) X <- Y[c(1,3)]  else X <- Y
    } else if(is.data.frame(data)) X <- ~1 else  stop("argument \"X\" is missing, with no default")  
	
  if (is(Y, "formula")) {
    name.Y <-  as.character(eval(Y)[[2]]) #serve sta roba?
    Y <- eval(attr(terms(Y, data=data), "variables"), data, environment(Y))[[attr(terms(Y, data=data), "response")]]
  } else {
    name.Y <- deparse(call$Y) #serve sta roba?
  }

  Y =  .makeContrasts(~.,data=data.frame(Y),excludeRefCat=FALSE,excludeIntercept=TRUE)
  # # remove redundant levels from factor Y
  # # and coerce Y to factor in case of levels input
  # if (is.factor(Y) || !is.null(levels))
  #  Y <- factor(Y) 
                                     
  # remove terms from X that are also in Z
  if (is(Z, "formula") && is(X, "formula") && 
        identical(environment(Z), environment(X))) {
	if( !( (length(attr(terms(X, data=data), "term.labels"))==0) & (length(attr(terms(Z, data=data), "term.labels"))==0)  )) {
		dup <- attr(terms(X, data=data), "term.labels") %in% attr(terms(Z, data=data), "term.labels")
		if (all(dup)) stop("all covariates in X also in Z")  
		if (any(dup)) 
			X <- formula(terms(X,data=data)[!dup])
	}
  }
  n <- dim(Y)[1]
  # get Z and X
  Z <- .getNull(Z, data, n)
  offset <- Z$offset   
  Z <- Z$Z
  Z <- Z[,apply(Z,2,function(x) !all(x==0)),drop=FALSE]
  #if(unique(Z[,"(Intercept)"])==0) 
  X <- .getAlternative(X, data, n)
  X <- X[, setdiff(colnames(X),colnames(Z)),drop=FALSE]
  
  # # Adjust input due to levels argument
  # if ((!is.null(levels)) && is.factor(Y)) {
    # if (!all(levels %in% levels(Y)))
      # stop("argument \"levels\" does not match levels(Y)")
    # if (length(levels) > 1) {
      # select <- Y %in% levels
      # Y <- factor(Y[select], levels=levels)
      # X <- X[select,,drop=FALSE]
      # Z$Z <- Z$Z[select,, drop=FALSE]
      # if (!is.null(Z$offset)) 
        # Z$offset <- Z$offset[select]
      # if (length(levels) == 2)
        # model <- "logistic"
    # } else {
      # Y <- factor(Y == levels)
      # levels(Y) <- c("other", levels)
      # model <- "logistic"
    # }
  # }
  
  
  # conservatively impute missing values in X
  all.na <- apply(is.na(X), 2, all)
  some.na <- apply(is.na(X), 2, any) & !all.na
  
  if (ncol(Z) == 0) {
   X[is.na(X)] <- 0
  } else {
	if(any(some.na))
		X[,some.na] <- apply(X[,some.na, drop=FALSE], 2, function(cov) {
			fit <- lm(cov ~ 0 + Z, x = TRUE)
			coefs <- coef(fit)
			coefs[is.na(coefs)] <- 0
			cov[is.na(cov)] <- drop(Z %*% coefs)[is.na(cov)]
			cov
			})
    X[,all.na] <- 0 
  }
  
  return(list(Y=Y,X=X,Z=Z))
}

##########################
# set the contrast for factors
#######################
###TODO : sistemare la funzione .makeContrasts!!!!!! mi pare che non funzioni bene excludeRefCat , non funziona neppure con interazioni & excludeRefCat (non esclude l'interazione di riferimento!)
.makeContrasts <- function(formu, data=data,excludeRefCat=TRUE,excludeIntercept=FALSE){ #excludeRefCat is used only for NOT ordered factors
    # make appropriate contrasts
    mframe <- model.frame(formu, data=data)
	if(length(mframe)>0){
		factors <- names(mframe)[sapply(mframe, is.factor)]
		contrs <- lapply(factors, function(fac) {
		levs <- levels(mframe[[fac]])
		k <- length(levs)
			if (is.ordered(mframe[[fac]])) {
				if (k %% 2 == 1) { 
					contr <- array(0, c(k, k-1), list(levs, paste("[",levs[-k], "<", levs[-1],"]", sep="")))
					contr[outer(1:k,1:(k-1), ">")] <- 1
					contr[,1:((k-1)/2)] <- contr[,1:((k-1)/2)] -1
				} else {
					levsplus <- c(levs[1:(k/2)], "(mid)", levs[k/2+1:(k/2)])
					contr <- array(0, c(k+1, k), list(levsplus, paste("[",levsplus[-(k+1)], "<", levsplus[-1],"]", sep="")))
					contr[outer(1:(k+1),1:k, ">")] <- 1
					contr[,1:(k/2)] <- contr[,1:(k/2)] - 1
					contr <- contr[-(1+k/2),]
					contr[,k/2+c(0,1)] <- contr[,k/2+c(0,1)] / sqrt(2)
				}
				# colnames(contr)=gsub("<",".lower.",colnames(contr))  # inutile
				# colnames(contr)=gsub(">",".greater.",colnames(contr))
			} else {
				contr <- diag(k)
				rownames(contr) <- levs
				colnames(contr) <- paste(".",levs,".",sep="")

				if(excludeRefCat) contr <- contr[,-1]
				
			}
			contr  
		})
		names(contrs) <- factors
	}
    # make the design matrix
    formu <- terms(formu, data=data)
    # if (length(attr(formu, "term.labels")) == 0)
      # stop("empty formu")
    if(excludeIntercept) attr(formu, "intercept") <- 0 else attr(formu, "intercept") <- 1
	# inutile :
#	attributes(attributes(mframe)$terms)$dataClasses[attributes(attributes(mframe)$terms)$dataClasses=="ordered"]="factor"
	# ords=rep(FALSE,dim(data)[2])
	# for(i in 1:length(ords)) ords[i]=is.ordered(data[,i])
	# for(i in which(ords)) data[,i]=factor(data[,i],ordered=FALSE)
	
    formu <- model.matrix(formu, contrasts.arg=contrs, data=data)
#	if(!all(colnames(formu) == "(Intercept)" ) ) { #if only the intercept is present
	#	formu <- formu[,colnames(formu) != "(Intercept)",drop=FALSE]    # ugly, but I've found no other way
 #   }
	formu
  }
############################
# Get the X design matrix
############################
.getAlternative <- function(X, data, n) {

  # coerce X into a matrix
  if (is.data.frame(X) || is.vector(X)) {
    if (all(sapply(X, is.numeric))) {
      X <- as.matrix(X)
    } else {
      stop("argument \"X\" could not be coerced into a matrix")
    }
  }
  
  # if (is(X, "ExpressionSet")) {
    # require("Biobase") || stop("ExpressionSet input but Biobase package not available")
    # X <- t(exprs(X))
  # }
  if (is(X, "formula")) {
    # keep NAs
    old.na.action <- options()$na.action  
    options(na.action="na.pass")

	X=.makeContrasts(X,data=data)
	
    # restore default
    options(na.action = old.na.action)
	}
  #check dimensions and names
  if (nrow(X) != n) {
    stop("the length of \"Y\" (",n, ") does not match the row count of \"X\" (", nrow(X), ")")
  }
  # if (is.null(colnames(X)))
    # stop("colnames missing in X design matrix")
  if(is.null(colnames(X))) colnames(X)=paste("X",1:dim(X)[2],sep="")
  X
}

############################
# Get the Z design matrix
############################
.getNull <- function(Z, data, n) {

  # coerce Z into a matrix and find the offset term
  offset <- NULL
  if (is.data.frame(Z) || is.vector(Z)) {
    if (all(sapply(Z, is.numeric))) {
      Z <- as.matrix(Z)
    } else {
      stop("argument \"Z\" could not be coerced into a matrix")
    }
  }
  if (is(Z, "formula")) {
    if (is.null(data)) {
      tnull <- terms(Z)
      # prevent problems for input ~1 or ~0:
      if (((attr(tnull, "Y") == 0)|| is.null(attr(tnull, "Y")))  && (length(attr(tnull, "term.labels")) == 0)
          && (length(attr(tnull, "offset")) == 0)) {
        if (attr(tnull, "intercept") == 1)
          tnull <- terms(numeric(n) ~ 1)
        else
          tnull <- terms(numeric(n) ~ 0)
      }
      offset <- model.offset(model.frame(tnull))
    } else {
      offset <- model.offset(model.frame(Z, data=data))
      tnull <- terms(Z, data=data)
    }
    data <- model.frame(tnull, data, drop.unused.levels = TRUE)
    Z <- model.matrix(tnull, data)
 
    # # suppress intercept if necessary (can this be done more elegantly?)
    # if (model == "cox") Z <- Z[,names(Z) != "(Intercept)"]
  }

  # check dimensions
  if (nrow(Z) != n) {
    stop("the length of \"Y\" (",n, ") does not match the row count of \"Z\" (", nrow(Z), ")")
  }
  list(Z = Z, offset = offset)
}

##############################################
