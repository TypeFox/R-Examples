#<<BEGIN>>
cornode <- function(...,target, outrank=FALSE, result=FALSE, seed=NULL)
#TITLE Builds a Rank Correlation using the Iman and Connover Method.
#KEYWORDS multivariate
#DESCRIPTION
#This function builds a rank correlation structure between columns of a matrix or between \samp{mcnode} objects
#using the Iman and Connover method (1982).
#INPUTS
#{\ldots}<<A matrix (each of its \samp{n} columns but the first one will be reordered)
#or \samp{n mcnode} objects (each elements but the first one will be reordered).>>
#{target}<<A scalar (only if \samp{n=2}) or a \samp{(n x n)} matrix of correlation.>>
#[INPUTS]
#{outrank}<<Should the order be returned?>>
#{result}<<Should the correlation eventually obtained be printed?>>
#{seed}<<The random seed used for building the correlation. If \samp{NULL} the \samp{seed} is unchanged.>>
#DETAILS
#The arguments should be named.
#
#The function accepts for \samp{data} a matrix or:
#{*}<<some \samp{"V" mcnode} objects separated by a comma;>>
#{*}<<some \samp{"U" mcnode} objects separated by a comma;>>
#{*}<<some \samp{"VU" mcnode} objects separated by a comma. In that case, the structure is built columns by colums (the first column of each \samp{"VU" mcnode}
#will have a correlation structure, the second ones will have a correlation structure, ....).>>
#{*}<<one \samp{"V" mcnode} as a first element and some \samp{"VU" mcnode} objects, separated by a comma.
#In that case, the structure is built between the \samp{"V" mcnode} and each column of the \samp{"VU" mcnode} objects.
#The correlation result (\samp{result = TRUE}) is not provided in that case.>>
#
#The number of variates of the elements should be equal.
#
#\samp{target} should be a scalar (two columns only) or a real symmetric positive-definite square matrix.
#Only the upper triangular part of \samp{target} is used (see \code{\link{chol}}).</>
# The final correlation structure should be
#checked because it is not always possible to build the target correlation structure.</>
#In a Monte-Carlo simulation, note that the order of the values within each \samp{mcnode} will be changed by this function
# (excepted for the first one of the list).
#As a consequence, previous links between variables will be broken.
#The \samp{outrank} option may help to rebuild these links (see the Examples).
#VALUE
#If \samp{rank = FALSE}: the matrix or a list of rearranged \samp{mcnode}s. </>
#If \samp{rank = TRUE}: the order to be used to rearranged the matrix or the \samp{mcnodes} to build the desired correlation structure.
#EXAMPLE
#x1 <- rnorm(1000)
#x2 <- rnorm(1000)
#x3 <- rnorm(1000)
#mat <- cbind(x1,x2,x3)
### Target
#(corr <- matrix(c(1,0.5,0.2,0.5,1,0.2,0.2,0.2,1),ncol=3))
### Before
#cor(mat,method="spearman")
#matc <- cornode(mat,target=corr,result=TRUE)
### The first row is unchanged
#all(matc[,1] == mat[,1])
#
###Using mcnode and outrank
#cook <- mcstoc(rempiricalD, values=c(0,1/5,1/50), prob=c(0.027,0.373,0.600), nsv=1000)
#serving <- mcstoc(rgamma, shape=3.93, rate=0.0806, nsv=1000)
#roundserv <- mcdata(round(serving), nsv=1000)
### Strong relation between roundserv and serving (of course)
#cor(cbind(cook,roundserv,serving),method="spearman")
#
###The classical way to build the correlation structure 
#matcorr <- matrix(c(1,0.5,0.5,1),ncol=2)
#matc <- cornode(cook=cook,roundserv=roundserv,target=matcorr)
### The structure between cook and roundserv is OK but ...
### the structure between roundserv and serving is lost
#cor(cbind(cook=matc$cook,serv=matc$roundserv,serving),method="spearman")
#
###An alternative way to build the correlation structure
#matc <- cornode(cook=cook,roundserv=roundserv,target=matcorr,outrank=TRUE)
### Rebuilding the structure
#roundserv[] <- roundserv[matc$roundserv,,]
#serving[] <- serving[matc$roundserv,,]
### The structure between cook and roundserv is OK and ...
### the structure between roundserv and serving is preserved
#cor(cbind(cook,roundserv,serving),method="spearman")
#CREATED 08-01-08
#REFERENCE
#Connover W., Iman R. (1982). A distribution-free approach to inducing rank correlation among input variables. Technometric, 3, 311-334.
#--------------------------------------------
{
  iman <- function(mat){  #Fonction de base de calcul
	  R <- sapply(1:p,function(i) sample(phi,nvar))
	  Rr <- apply(R,2,rank,ties.method="random")
	  Ret <- Rr %*% P
	  Rret <- apply(Ret,2,rank,ties.method="random")
	  Xr <- apply(mat,2,order)
    rang <- sapply(1:p, function(i) (Xr[,i])[Rret[,i]])
	  rang <- rang[order(rang[,1]),]
	  if(outrank) return(rang)
	  return(sapply(1:p, function(i) mat[rang[,i],i]))
  }

  data <- list(...)
  if(!is.null(seed)) set.seed(seed)
  if(!is.matrix(target)) target <- matrix(c(1,target,target,1),ncol=2)

  if(length(data)==1) { # A matrix
    data <- data[[1]]
    if(!is.matrix(data)) stop("data should be a matrix or a list of mcnodes")
    p <- ncol(data)
    nvar <- nrow(data)

    if(p < 2) stop("the matrix should have at least 2 columns")
    if(!is.matrix(target) || any(dim(target)!=c(p,p))) stop("target should be a matrix of correct dimension")

	  phi <- qnorm((1:nvar)/(nvar+1))
	  P <- chol(target)
    resu <- iman(data)

    colnames(resu) <- colnames(data)
    if(result)
      {cat("Spearman Rank Correlation Post Function\n")
      print(cor(resu,method="spearman"))}
    return(resu)
  }

 # mcnodes
    	list.names <- function(...) {
        	l <- as.list(substitute(list(...)))[-1]
        	nm <- names(l)
        	fixup <- if (is.null(nm)) seq(along = l) else nm == ""
        	dep <- sapply(l[fixup], function(x) if (is.symbol(x)) as.character(x) else "")
        	if (is.null(nm))
            		return(dep)
        	else {
            		nm[fixup] <- dep
            		return(nm)
        	}
    	}

    noms <- list.names(...)
    p <- length(data)
    if(p < 2) stop("the list should have at least 2 mcnodes")
    mcn <- sapply(data,inherits,"mcnode")
    if(!all(mcn))
      stop("the list should be a list of mcnode objects")

    if(!is.matrix(target) || any(dim(target)!=c(p,p)))
      stop("target should be a matrix of dimension l*l where l is the number of mcnodes in the list")

    tmcn <- sapply(data,attr,which="type")
    if(!(all(tmcn=="U") | all(tmcn %in% c("V","VU"))))
      stop("incorrect combination of mcnode: either 'U' or a combinaison of 'V' and 'U'")

    if(tmcn[1]=="U") {
     tmcd <- sapply(data,dim)
      for(i in 1:p) {
        data[[i]][] <- aperm(data[[i]],perm=c(2,1,3))
        dim(data[[i]]) <- tmcd[c(2,1,3),i] }                            # permute without change of structure
      }
    else
      if(any(tmcn=="VU") && any(tmcn=="V")){
        if(result){
          warning("impossible to provide the correlation result")
          result <- FALSE
        }
        if(sum(tmcn=="V")!=1 && tmcn[1]!="V")
          stop("Valid if only one 'V' node in the first position is combined with 'VU' nodes")
      }

    tmcd <- sapply(data,dim)
    nvar <- max(tmcd[1,])
    nunc <- max(tmcd[2,])
    nvariates <- max(tmcd[3,])

    if( !all(tmcd[1,] %in% c(1,nvar)) || !all(tmcd[2,] %in% c(1,nunc)))
      stop("incorrect dimension of mcnodes")
    if( any(tmcd[3,] != nvariates)    )
      stop("incorrect number of variates in mcnodes")

    res   <- data

    datai <- matrix(NA,ncol = p,nrow = nvar)
    phi <- qnorm((1:nvar)/(nvar+1))
    P <- chol(target)

    if(result) {corobs <- array(NA,dim=c(p*p,nunc,nvariates))}

    for(k in 1:nvariates){
      for(i in 1:nunc){
        ir <- ifelse(tmcd[2,1]==1, 1, i)
        datai[,1] <- data[[1]][,ir,k]
        for(j in 2:p) datai[,j] <- data[[j]][,i,k]
          sortim <- iman(datai)
        for(j in 2:p) res[[j]][,i,k] <-  sortim[,j]
        if(result) corobs[,i,k] <- cor(sortim,method="spearman",use="pairwise")
      }
    }

  if(result){
    if(nunc > 1)
      { corobs <- apply(corobs,c(1,3),function(x) c(mean=mean(x),quantile(x,c(0.5,0,1))))
        cat("summary of output Rank Correlation obtained accross the uncertainty dimension for each variates\n")}
      else
        {corobs <- aperm(corobs,perm=c(2,1,3))
        cat("output Rank Correlation per variates\n")}
    for(i in 1:nvariates){
      cat("variates:",i,"\n")
      print(corobs[,,i])}
  }

   if(outrank) res[[1]][] <- 1:nvar

   if(tmcn[1]=="U") {
    for(i in 1:p) {
        res[[i]][] <- aperm(res[[i]],perm=c(2,1,3))
        dim(res[[i]]) <- tmcd[c(2,1,3),i] }                            # repermute without change of structure
    }

   names(res) <- noms

  return(res)}

