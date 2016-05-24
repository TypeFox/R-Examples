#######################
flip.statTest <-
    c("t", "F", "ANOVA",
	"Wilcoxon","Kruskal-Wallis", "kruskal", "rank", "Mann-Whitney",
	"chisq","chisq.separated", "Fisher",
	#"KS", "kolmogorow", "Kolmogorow-Smirnov", "ad",
  "McNemar", "Sign","sum","coeff","NA")

.get.statTest <- function(statTest){ 
	if(is(statTest,"function")) return(statTest) else
	
	statTest <- match.arg(tolower(statTest[1]),tolower(flip.statTest))
	statTest= flip.statTest[which(statTest==tolower(flip.statTest))]
	#synonyms
	if(statTest=="ANOVA") 
		statTest="F" else
	if(statTest=="kruskal") 
		statTest="Kruskal-Wallis" else
	if(statTest=="Mann-Whitney")
		statTest="Wilcoxon" else
# 	if(statTest%in%c("KS", "kolmogorow"))
# 		statTest="Kolmogorow-Smirnov"
# 		
	statTest
}

#########################
flip <- function(Y, X=NULL, Z=NULL, data=NULL, tail = 0, perms = 1000, statTest=NULL, 
                 Strata=NULL, flipReturn, testType=NULL, ...) {

  if(is.null(statTest) ) if(is.null(list(...)$separatedX)   || list(...)$separatedX)   { statTest="t" } else statTest="F"
    statTest <- .get.statTest(statTest)
	
  if(is.null(testType)){
	if(is.null(list(...)$rotationTest) || (!list(...)$rotationTest) ) {testType="permutation"} else { testType="rotation"} 
  } 
  testType=match.arg(testType,c("permutation","rotation","symmetry"))

  if(missing(flipReturn)||is.null(flipReturn)) 
    flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,test=TRUE,permID=TRUE)
                    

  # store the call
  call <- match.call()
  
  if(!is.function(statTest)){
    # get matrices from inputs
    data <- .getXY(Y,X,Z,data,rotationTest=(testType=="rotation"),dummyfy=list(...)$dummyfy,statTest=statTest,Strata=Strata)
    rm(X,Y,Z,Strata)
    
    symmetryTest= is.null(data$X) || (length(unique(data$X))==1)
  
    #check if the problem can be set as one sample problem
    if(!symmetryTest) if(!is.function(statTest))
  	if(statTest%in% c("t","sum","ranks","Wilcoxon","McNemar","Sign"))
  	  if(  !is.null(data$Strata) ){#is.null(data$Z)|| ncol(data$Z)==0)  &
  			keep=setdiff(1:ncol(data$X),.getIntercept(data$X))
  			if( (length(unique(data$X[,keep]))==2) && 
  				(ncol(data$X[,keep,drop=FALSE])==1) )
  					if(all(table(data$X[,keep],unlist(data$Strata))==1)){
  						attrsYassign=attributes(data$Y)$assign
  						attrsYfactors=attributes(data$Y)$factors
              
  						data$X=data$X[,keep,drop=FALSE]
  						levs=unique(data$X)
  						data$Y=t(sapply(unique(as.character(unlist(data$Strata))), function(ids){
                data$Y[(data$Strata==ids)&(data$X==levs[2]),]-
                data$Y[(data$Strata==ids)&(data$X==levs[1]),]}))
              
  						attributes(data$Y)$assign=attrsYassign
  						attributes(data$Y)$factors=attrsYfactors
  						data$X=NULL
  						data$Strata=NULL
  						data$Z=NULL
  						symmetryTest=TRUE
  					}	
  		}
    
    # if symmetry.nptest
    if(symmetryTest){
    		test= .symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,testType=testType,...)
    ##dependence.nptest
    } else 
  	if ( !(any(is.na(data$Y))|| ifelse(is.null(data$X),TRUE,any(is.na(data$X)))) || statTest=="NA"){
      # standard solutions, not missing data
  		test= .dependence.nptest(data, perms=perms,statTest=statTest,  
                               tail = tail,testType=testType,
                               return.permIDs = flipReturn$permID, ...)
  	} else {	stop("Warning: NA values are not allowed unless you use statTest=\"NA\", nothing done.")	}
    #browser()
  } else{
    test= .custom.nptest(Y=Y, X=X, Z=Z, data=data, tail = tail, perms = perms, statTest=statTest, Strata=Strata, flipReturn=flipReturn, testType=testType, ...)
  }
  res <- test$test()
	#build the flip-object
  res$call=call
	res=.getOut(res=res,data=data, call=call, flipReturn=flipReturn,call.env=test)
  return(res)
}

#################

.custom.nptest <- function(Y,X=NULL, Z=NULL, data=NULL, tail = 0, 
          perms = 1000, statTest=NULL, Strata=NULL, flipReturn, 
                           testType=NULL, ...) {
  test<- function() {
    N=nrow(Y)
    if(is.null(N)) N=length(N)
    perms <- make.permSpace(N,perms,return.permIDs=FALSE,Strata=Strata)
    perms$rotFunct=NULL
    digitsK=trunc(log10(perms$B))+1
    
    
    obs=statTest(Y)
    permT=matrix(,perms$B,length(obs))
    colnames(permT)=names(obs)
    permT[1,]=obs
    for(i in 2:(perms$B))
      {permT[i,]=statTest(Y[sample(perms$n),,drop=FALSE] )
              if (i%%10==0) {
                cat(rep("\b", 2*digitsK+3), i, " / ", perms$B, sep="")
                flush.console()
              } 
    }
    flush.console()
    cat("\n")
    if(is.null(colnames(permT))){
      if(ncol(permT)==ncol(Y)) 
        colnames(permT)=.getTNames(Y,,permT=permT,checkUnique=TRUE) else
          colnames(permT)=.getTNames(Y,X,permT=permT,checkUnique=TRUE)
    }
    rownames(permT)=.getTRowNames(permT)		  
    res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="Custom"))
  }
  
  environment(test) <- sys.frame(sys.nframe())
  out <- sys.frame(sys.nframe())
  return(out)
}