#' (Composite) Interval Mapping for QTL detection in multi-parent crosses
#'
#' Interval mapping in multi-parent crosses with options for single-stage mixed
#' model approach; multi-stage approach using predicted means; multi-stage
#' approach including cofactors (CIM)
#'
#' @export
#' @param baseModel Base phenotypic model for analysis
#' @param object Object of class \code{mpcross}
#' @param pheno Phenotypic object
#' @param idname The idname in phenotypic data for which to output predicted means. Should match rownames of the object$finals
#' @param threshold Significance threshold for QTL p-values
#' @param chr Subset of chromosomes for which to compute QTL profile
#' @param step Step size at which to compute the QTL profile. See \code{\link[mpMap]{mpprob}} for further description of default values
#' @param responsename Optional input of response name to look for in object$pheno 
#' @param ncov Number of marker covariates to search for - default is to search for as many as possible using stepAIC (forward/backward selection)
#' @param window Window of cM on each side of markers where we exclude covariates in CIM
#' @param mrkpos Flag for whether to consider both marker positions and step positions or just steps. Is overridden if step=0
#' @param ... Additional arguments
#' @return The original input object with additional component QTLresults containing the following elements:
#' \item{pheno}{Input phenotype data}
#' \item{pvalue}{Each component contains estimated p-values at each position on a given chromosome}
#' \item{wald}{Each component contains Wald statistics at each position on a given chromosome}
#' \item{fndrfx}{Each component contains founder effects estimated at each position on a given chromosome}
#' \item{qtl}{Each component contains the position and effects of a detected QTL}
#' \item{call}{Input arguments to function} 
#' and with attributes describing the number of QTL detected, and the threshold used for detection. Will only return one QTL per chromosome; to find more 
#' QTL see \code{\link[mpMap]{findqtl2}}
#' @details 
#' Depending on the options selected, different models will be fit for QTL
#' detection. If the baseModel input does not include a term matching the 
#' idname input, it will be assumed that a single-stage QTL mapping approach
#' is desired. In this case, no covariates will be added (ncov will be set to
#' 0); all models will be fitted in asreml; and all phenotypic covariates and
#' design factors specified in the baseModel will be fitted along with genetic
#' covariates in mixed model interval mapping. 
#'
#' If the baseModel input does include a term matching the idname, then it will
#' be assumed that a two-stage QTL mapping approach is desired. In this case,
#' the baseModel will be fit using asreml and predicted means will be output
#' to be used as a response in linear model interval mapping. If 
#' \code{ncov>0} additional marker cofactors will be fit; otherwise simple
#' interval mapping will be run. All phenotypic covariates and design factors
#' specified in the baseModel will be fit in the first stage.  
#' 
#' Note that no weights are used in the second stage of analysis which may result in a loss of efficiency compared to a one-stage approach.
#'
#' If no baseModel is input, it will be assumed that predicted means have been
#' included in \code{object} as a phenotypic variable named predmn. In this 
#' case \code{pheno} is not required and asreml does not need to be used. 
#' (composite) Interval mapping will proceed
#' as in the two-stage case depending on the value of \code{ncov}. 
#' @seealso \code{\link[mpMap]{plot.mpqtl}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl")
#' ## Two-stage simple interval mapping 
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")


mpIM <- function(baseModel, object, pheno, idname="id", threshold=1e-3, chr, step=0, responsename="predmn", ncov=1000, window=10, mrkpos=FALSE, ...)
{
  ### Initial setup for all approaches
  lines <- rownames(object$finals)
  n.founders <- nrow(object$founders) 

  if (missing(chr)) chr <- 1:length(object$map) 
  else if (is.character(chr)) chr <- match(chr, names(object$map))

  if (!(inherits(object, "mpprob") && attr(object$prob, "step")==step 
  && attr(object$prob, "mrkpos")==mrkpos))  
  object <- mpprob(object, program="qtl", step=step, chr=chr, mrkpos=mrkpos)

  fmap <- attr(object$prob, "map")

  output <- object

  wald <- list()
  pval <- list()
  fndrfx <- list()
  se <- list()
  degf <- list()

  ### Check for which option should be run
  ### If idname is in the fixed effects of baseModel = lm or CIM
  ### if missing baseModel = lm or CIM and response in object$pheno
  
  if (!missing(baseModel)) {
    require(asreml) 

    if (length(grep(idname, as.character(baseModel$call$fixed)[3]))>0)
	method <- "lm" else method <- "mm"
  } 
  else method <- "lm"

  ### setup for two-stage approach if predicted means have not yet been computed
  if (method=="lm") {
    require(aod)
    if (!missing(baseModel)) {
      ## if baseModel and idname have been input and fixed term 
      ## in baseModel is idname then get out predicted means
      pmmod <- update(baseModel, data=pheno)
    
      ## get out the predicted means
      cf <- pmmod$coefficients$fixed
      pm <- cf[grep(idname, names(cf))]
      names(pm) <- substr(names(pm), nchar(idname)+2, nchar(names(pm)))

      object$pheno[[responsename]] <- pm
      ## otherwise predicted means should already be stored there
    }
   
    ## reset this - this is the only part of the phenotype matrix needed
    pheno <- object$pheno
    if (missing(idname)) idname <- "id"
    pheno[[idname]] <- rownames(object$pheno)

    ## create a vector of responses just for use in finding marker cofactors
    predmn <- pheno[[responsename]][which(!is.na(match(pheno[[idname]], lines)))]

    ### Setting up formulas in the case of covariates
    ### In this case need to just use the markers for selection
    ### set up dataframe with just those values
    if (ncov>0) 
    {
	mrkobj <- mpprob(object, step=0, program="qtl", chr=chr)
	mrkgen <- do.call("cbind", mrkobj$prob)
	mrkgen <- mrkgen[na.omit(match(lines, pheno[[idname]])),]
	## recenter before selecting marker covariates
	for (k in 1:ncol(mrkgen)) mrkgen[,k] <- mrkgen[,k]-mean(mrkgen[,k], na.rm=TRUE)


 	#largest possible formula
    	formula <- as.formula(paste("predmn~", paste(paste("mrkgen[,", seq(1, ncol(mrkgen), n.founders), ":", seq(n.founders, ncol(mrkgen), n.founders), "]", sep=""), collapse="+")))
    	mod <- lm(as.formula(paste("predmn~1", sep="")))
    	modc <- stepAIC(mod, scope=formula, steps=ncov)  

    	covar <- attr(modc$terms, "variables")

    	ncov <- length(covar)-2
    	rmv <- list(length=ncov)
    	# the marker index at which each chromosome begins
	gmap <- attr(mrkobj$prob, "map")
	ngmap <- unlist(lapply(gmap, names))
    	nchrmrk <- unlist(lapply(gmap, length))
    	csnchrmrk <- c(0, cumsum(nchrmrk))

	## format the covariates which will be included in the model
    	if (length(covar)>2)
    	for (k in 3:length(covar))
	{
	    #as.character(covar[[k]][n.founders]) = index into gen3 for this covariate, e.g. 13:16, 21:24
	    #strsplit(as.character(covar[[k]][n.founders]), ":") = split this into two numbers
	    #start = number of the marker that these four covariates represent
	    ###modification, if start is intedend to be the marker index represented by this covariate then it's off by 4, was 7 before. Overflowed
	    start<-(as.numeric(strsplit(as.character(covar[[k]][n.founders]), ":")[[1]][1])+3)/n.founders
	    #which chromosome is this marker on
	    chrn <- which.max(csnchrmrk[csnchrmrk<start])
	    cpos <- grep(ngmap[start], names(fmap[[chrn]]))
  	    # rmvpos is the list of positions for which covariate k should be excluded (with the chromosome stored at the front)
	    rmvpos <- c(chrn, which(fmap[[chrn]] <= fmap[[chrn]][cpos]+window & fmap[[chrn]] >= fmap[[chrn]][cpos] - window))
	    rmv[[k-2]] <- rmvpos
	}
   	    
    	#chromosome numbers of the genetic covariates
    	rmvchr <- unlist(lapply(rmv, function(x) return(x[1])))

    	### set up a list with a matrix for each component representing chr 
    	### one row per genetic covariate, one column per location 
    	### FALSE indicates to EXCLUDE the covariate for that location
    	terms <- fmap
    	for (i in 1:length(terms)) 
    	{
	  terms[[i]] <- matrix(data=TRUE, nrow=ncov, ncol=length(terms[[i]]))
	  ## if cofactors are on this chromosome, set array accordingly
	  if (i %in% rmvchr)
	  for(k in 1:length(rmv))
	  if(rmv[[k]][1] == i) 	terms[[i]][k, rmv[[k]][-1]] <- FALSE
    	}

    	## alter fixed terms in model depending on position
    	formall <- modc$call[[2]]
    	termlab <- attr(modc$terms, "term.labels")
    }
    # matrix has all TRUE values with no covariates
    else 
    {
    	terms <- fmap
    	for (i in 1:length(terms)) 
 	  terms[[i]] <- matrix(data=TRUE, nrow=ncov, ncol=length(terms[[i]]))
    	termlab <- ""
    }
  } # end of "lm" loop

  for (j in chr)
  {
    nam <- names(object$map)[j]
    cat("------Analyzing Chr ",nam,"----------\n")
    gen <- object$prob[[nam]]

    ## All set to NA by default
    wald[[nam]] <- rep(NA, ncol(gen)/n.founders)
    pval[[nam]] <- rep(NA, ncol(gen)/n.founders)
    degf[[nam]] <- rep(NA, length=ncol(gen)/n.founders)
    fndrfx[[nam]] <- matrix(nrow=n.founders, ncol=ncol(gen)/n.founders)
    se[[nam]] <- matrix(nrow=n.founders, ncol=ncol(gen)/n.founders)
    df <- matrix(nrow=nrow(pheno), ncol=ncol(gen))
  
    colnames(gen) <- paste("P", rep(1:(ncol(gen)/n.founders), each=n.founders), "F", LETTERS[1:n.founders], sep="")

    genid <- vector()
    pheid <- vector()
    for (k in 1:length(lines)) {
     matchid <- which(as.character(pheno[[idname]])==as.character(lines[k]))
     genid <- c(genid, rep(k, length(matchid)))
     pheid <- c(pheid, matchid) }
    df[pheid, ] <- as.matrix(gen[genid,])
    df <- cbind(pheno, df)
    df <- as.data.frame(df)

    for (k in 1:ncol(pheno))
    {
     fx <- paste("as.", class(pheno[,k]), sep="")
     df[,k] <- do.call(fx, list(df[,k]))
    }
    names(df) <- c(names(pheno), colnames(gen))

    # recenter probabilities to 0
    for (k in (ncol(pheno)+1):ncol(df))
	df[,k] <- df[,k] - mean(df[,k], na.rm=TRUE)

    for (k in seq(1, ncol(gen), n.founders))
    {
	index <- (k-1)/n.founders+1

	if (method=="mm") {
	  mod <- update(baseModel, 
	  	fixed=eval(as.formula(paste(baseModel$call$fixed[2], 
		baseModel$call$fixed[1], 
		paste(c(as.character(baseModel$call$fixed[3]), 
		names(df)[ncol(pheno)+k:(k+n.founders-1)]), collapse="+"), sep=""))), 
		data=df, Cfixed=TRUE, na.method.X="include")
	
	  summ <- summary(mod, all=TRUE)	
    	  fndrfx[[nam]][,index] <- summ$coef.fixed[n.founders:1,1]
	  se[[nam]][,index] <- summ$coef.fixed[n.founders:1, 2]
    	  man <- list(which(mod$coefficients$fixed[1:n.founders]!=0), "zero")
	  wta <- wald.test.asreml(mod, list(man, "zero"))$zres
    	  wald[[nam]][(k-1)/n.founders+1] <- wta$zwald
    	  degf[[nam]][(k-1)/n.founders+1] <- nrow(wta$ZRows[[1]])
    	  pval[[nam]][(k-1)/n.founders+1] <- wta$zpval
      	}

	if (method=="lm") {
	  if (ncov>0) 
	    # include all necessary covariates
	    form <- as.formula(paste("predmn~", paste(termlab[which(terms[[j]][,index]==1)], collapse="+"), "+", paste(names(df[,2+k:(k+n.founders-1)]), collapse="+"), sep="")) 
	  else form <- as.formula(paste("predmn~", paste(names(df[,2+k:(k+n.founders-1)]), collapse="+"), sep=""))
	  # fit the model
	  mod <- lm(form, data=df)
	  # Collinearity means some coefficients will be NA
	  coe <- coef(mod)[which(!is.na(coef(mod)))]	

	  # if we did not get any non NA values, window was not large enough
	  if(length(grep(paste("P", index, "F", sep=""), names(coe))) != 0)
	  {
	    #test the current location for significance (actually a joint test)
	    wt <- wald.test(vcov(mod), coe, Terms=grep(paste("P", index, "F", sep=""), names(coe)))

	    pval[[nam]][index] <- wt$result$chi2[3] 
	    wald[[nam]][index] <- wt$result$chi2[1]
	    degf[[nam]][index] <- wt$result$chi2[2]

	    a <- summary(mod)$coefficients
	    fndrfx[[nam]][,index] <- c(a[grep(paste("P", index, "F", sep=""), rownames(a)),1], rep(NA, n.founders-length(grep(paste("P", index, "F", sep=""), rownames(a))))) 
	    se[[nam]][,index] <- c(a[grep(paste("P", index, "F", sep=""), rownames(a)),2], rep(NA, n.founders-length(grep(paste("P", index, "F", sep=""), rownames(a))))) 
	  }
	  else
	    cat("Error when testing location ", index, " along chromosome ", j, ". Window may be too small.\n ")
    	} ## end of "lm" loop
     }
  }
  
  minp <- unlist(lapply(pval, function(x) min(x, na.rm=TRUE)))
  maxw <- unlist(lapply(wald, which.max))
  sigchr <- which(minp < threshold)
  map <- attr(object$prob, "map")

  ## select out most significant QTL
  pos <- vector(length=length(fndrfx))
  pos[1:length(pos)] <- NA
  pos[sigchr] <- maxw[sigchr]

  ## positions in cM
  posqtl <- vector(length=length(pos))
  for (i in 1:length(pos))
    posqtl[i] <- (!is.na(pos[i]))*map[[i]][maxw[i]]

  qtl <- list()
  for (i in which(!is.na(pos)))
  {
    qtl[[names(map)[i]]] <- t(as.matrix(c(posqtl[i], fndrfx[[names(map)[i]]][, pos[i]], se[[names(map)[i]]][, pos[i]])))
    attr(qtl[[names(map)[i]]], "index") <- pos[i]
  }
  
  results <- list()
  if (!missing(baseModel))   results$baseModel <- baseModel
  else results$baseModel <- lm(predmn~1, data=pheno)
  results$pheno <- pheno
  results$pvalue <- pval
  results$wald <- wald
  results$se <- se
  results$degf <- degf
  results$fndrfx <- fndrfx
  results$qtl <- qtl
  results$call <- match.call()
  attr(results$qtl, "threshold") <- threshold
  attr(results$qtl, "nqtl") <- sum(!is.na(pos)) 
  attr(results$qtl, "ncov") <- ncov
  attr(results$qtl, "window") <- window
  attr(results, "method") <- method

  output$QTLresults <- results
  class(output) <- c("mpqtl", class(object))

  output
}
