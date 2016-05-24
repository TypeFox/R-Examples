#' Fit a full model including all QTL and effects from base model
#' 
#' Given the output from a scan of chromosomes with significant genetic
#' variation, fits a full mixed model containing all effects in base model and
#' all QTL effects. 
#' @aliases fit fit.mpqtl
#' @usage fit (object, ...)
#' \method{fit}{mpqtl} (object, ...)
#' @export fit
#' @S3method fit mpqtl
#' @method fit mpqtl
#' @param object Object of class \code{mpqtl}
# @param baseModel asreml output from fit of base model
# @param pheno data frame containing phenotypes required to fit base model
# @param effects Flag for whether to include QTL as fixed or random effects
#' @param \dots Additional arguments to be used in \code{asreml}
#' @return An asreml model and summary table of QTL effects, p-values and Wald statistics from fitting the full model.
#' @seealso \code{\link[mpMap]{mpIM}}, \code{\link[mpMap]{summary.mpqtl}}
#' @examples
#' sim.map <- sim.map(len=rep(100, 2), n.mar=11, include.x=FALSE, eq.spacing=TRUE)
#' sim.ped <- sim.mpped(4, 1, 500, 6, 1)
#' sim.dat <- sim.mpcross(map=sim.map, pedigree=sim.ped, qtl=matrix(data=c(1, 10, .4, 0, 0, 0, 1, 70, 0, .35, 0, 0), nrow=2, ncol=6, byrow=TRUE), seed=1)
#' mpp.dat <- mpprob(sim.dat, program="qtl", step=2)
#' mpq.dat <- mpIM(object=mpp.dat, ncov=0, responsename="pheno")
#' fit(mpq.dat)

fit <- function(object, ...)
{
	UseMethod("fit")
}

# note: random effects fitting will eventually be added. Not currently an option
fit.mpqtl <- function(object, ...)
{

  if (!inherits(object, "mpqtl")) stop("Must have object of type mpqtl to
	  fit full model\n")

  output <- list()

  ## given a QTL object want to fit a full model 
  ## containing all the terms in the base model, etc. 
  lines <- rownames(object$finals)
  qtlres <- object$QTLresults$qtl
  nqtl <- attr(qtlres, "nqtl")
  n.founders <- nrow(object$founders)
  map <- object$map
  f3 <- substr(rownames(object$founders), 1, 3)
  fmap <- attr(object$prob, "map")
  baseModel <- object$QTLresults$baseModel
  pheno <- object$QTLresults$pheno
  method <- attr(object$QTLresults, "method")

  ## need to figure out what position they're at
  ## attr(, "index") may be useful for this
  ## is this something updated with findqtl2?

  ## so we need to a) set up a df to include the genetic terms for 
  ## each QTL
  ## b) replace the third term in the fixed model statement
  ## c) extract all the estimates etc. and store them
  index <- unlist(lapply(qtlres, function(x) attr(x, "index")))
  chr <- rep(names(qtlres), unlist(lapply(qtlres, function(x) return(nrow(x)))))
  gen <- list()
  ## check this is extracting the right columns
  for (i in 1:nqtl) 
	gen[[i]] <- object$prob[[chr[i]]][,(index[i]-1)*n.founders+1:n.founders]
  
  gen <- do.call("cbind", gen)
  colnames(gen) <- paste("P", rep(1:(ncol(gen)/n.founders), each=n.founders), "F", LETTERS[1:n.founders], sep="")

  df <- matrix(nrow=nrow(pheno), ncol=ncol(gen))

  genid <- vector()
  pheid <- vector()
  for (k in 1:length(lines)) {
   	  matchid <- which(as.character(pheno$id)==as.character(lines[k]))
     	  genid <- c(genid, rep(k, length(matchid)))
     	  pheid <- c(pheid, matchid) 
  	}
  df[pheid, ] <- as.matrix(gen[genid,])
  df <- cbind(pheno, df)
  df <- as.data.frame(df)

  for (k in 1:ncol(pheno))
  {
    	fx <- paste("as.", class(pheno[,k]), sep="")
     	df[, k] <- do.call(fx, list(df[,k]))
  }
  names(df) <- c(names(pheno), colnames(gen))

  # recenter probabilities to 0
  for (k in (ncol(pheno)+1):ncol(df))
  	df[,k] <- df[,k] - mean(df[,k], na.rm=TRUE)

  wald <- vector(length=nqtl)
  pval <- vector(length=nqtl)
  degf <- vector(length=nqtl)
  fmrkl <- vector(length=nqtl)
  fmrkr <- vector(length=nqtl)
  effects <- "fixed"

  if (effects=="fixed") {

  if (method=="lm") {
    mod <- update(baseModel, 
	formula=eval(as.formula(paste(baseModel$call$formula[2], 
	baseModel$call$formula[1], 
	paste(c(as.character(baseModel$call$formula[3]), 
	names(df)[ncol(pheno)+1:ncol(gen)]), collapse="+"), sep=""))), 
	data=df)

    summ <- summary(mod)$coefficients
    effect <- rep(NA, n.founders)
    se <- rep(NA, n.founders)
    names(effect) <- names(se) <- names(coef(mod))[grep("P", names(coef(mod)))]
    ind <- match(rownames(summ)[grep("P", rownames(summ))], names(effect))
    effect[ind] <- summ[grep("P", rownames(summ)), 1]
    se[ind] <- summ[grep("P", rownames(summ)),2]

    cm <- round(unlist(lapply(qtlres, function(x) return(x[,1]))),2)
    ## these should be done individually for each QTL to test for significance
    for (j in 1:nqtl) {	
	subind <- grep(paste("P", j, "F", sep=""), names(coef(mod)))
  	man <- subind[which(!is.na(coef(mod)[subind]))]
  	wt <- wald.test(vcov(mod), coef(mod)[which(!is.na(coef(mod)))], Terms=man)
	wald[j] <- wt$result$chi2[1]
	degf[j] <- wt$result$chi2[2]
	pval[j] <- wt$result$chi2[3]
	mrkli <- which.max(map[[chr[j]]]*(map[[chr[j]]]<=cm[j]))
	if (length(map[[chr[j]]])>1) {
	if (mrkli==length(map[[chr[j]]])) mrkli <- mrkli-1
	fmrkl[j] <- names(map[[chr[j]]])[mrkli]
	fmrkr[j] <- names(map[[chr[j]]])[mrkli+1]
	} else fmrkl[j] <- fmrkr[j] <- names(map[[chr[j]]])[mrkli]
    }
  } ## end of method=="lm"

  if (method=="mm") {
    require(asreml)
    mod <- update(baseModel, 
	fixed=eval(as.formula(paste(baseModel$call$fixed[2], 
	baseModel$call$fixed[1], 
	paste(c(as.character(baseModel$call$fixed[3]), 
	names(df)[ncol(pheno)+1:ncol(gen)]), collapse="+"), sep=""))), 
	data=df, Cfixed=TRUE, na.method.X="include")

    effect <- summary(mod, all=T)$coef.fixed[ncol(gen):1,1]
    se <- summary(mod, all=T)$coef.fixed[ncol(gen):1, 2]

    cm <- round(unlist(lapply(qtlres, function(x) return(x[,1]))),2)
    ## these should be done individually for each QTL to test for significance
    for (j in 1:nqtl) {	
	subind <- grep(paste("P", j, "F", sep=""), names(mod$coefficients$fixed))
  	man <- subind[which(mod$coefficients$fixed[subind]!=0)]
  	wald[j] <- wald.test.asreml(mod, list(list(man, "zero")))$zres$zwald
    	degf[j] <- nrow(wald.test.asreml(mod, list(list(man, "zero")))$zres$ZRows[[1]])
    	pval[j] <- wald.test.asreml(mod, list(list(man, "zero")))$zres$zpval
	mrkli <- which.max(map[[chr[j]]]*(map[[chr[j]]]<=cm[j]))
	if (length(map[[chr[j]]])>1) {
	if (mrkli==length(map[[chr[j]]])) mrkli <- mrkli-1
	fmrkl[j] <- names(map[[chr[j]]])[mrkli]
	fmrkr[j] <- names(map[[chr[j]]])[mrkli+1]
	} else fmrkl[j] <- fmrkr[j] <- names(map[[chr[j]]])[mrkli]
    }
  } ## end of method=="mm"


   ## these will stay the same - just get the values out separately 
   effect <- t(matrix(round(effect,2), nrow=nqtl, ncol=n.founders, byrow=T))
   se <- t(matrix(round(se,3), nrow=nqtl, ncol=n.founders, byrow=T))
   eff3 <- paste("Effect_",f3,sep="")
   se3 <- paste("SE_",f3,sep="")
   if (n.founders==4)
	table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], "Wald"=round(wald,2), "df"=degf, "pvalue"=signif(pval,3))
   else if (n.founders==8)
	table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], effect[5,], se[5,], effect[6,], se[6,], effect[7,], se[7,], effect[8,], se[8,], "Wald"=round(wald,2), "df"=degf, "pvalue"=signif(pval,3))

   names(table)[seq(5, 4+(2*n.founders), 2)] <- eff3
   names(table)[seq(6, 4+(2*n.founders), 2)] <- se3
  }
  else if (effects=="random") {
    mod <- update(baseModel, 
	random=eval(as.formula(paste(baseModel$call$random[1], 
	paste(c(as.character(baseModel$call$random[3]), 
	names(df)[ncol(pheno)+1:ncol(gen)]), collapse="+"), sep=""))), 
	data=df, Cfixed=TRUE, na.method.X="include")

    wald <- vector(length=nqtl)
    pval <- vector(length=nqtl)
    degf <- vector(length=nqtl)
    fmrkl <- vector(length=nqtl)
    fmrkr <- vector(length=nqtl)

    ## need to think about what summaries make sense with random effects!
#    effect <- summary(mod, all=T)$coef.random[ncol(gen):1,1]
#    se <- summary(mod, all=T)$coef.fixed[ncol(gen):1, 2]

    cm <- round(unlist(lapply(qtlres, function(x) return(x[,1]))),2)
    ## these should be done individually for each QTL to test for significance
    for (j in 1:nqtl) {	
	subind <- grep(paste("P", j, "F", sep=""), names(mod$coefficients$fixed))
  	man <- subind[which(mod$coefficients$fixed[subind]!=0)]
  	wald[j] <- wald.test.asreml(mod, list(list(man, "zero")))$zres$zwald
    	degf[j] <- nrow(wald.test.asreml(mod, list(list(man, "zero")))$zres$ZRows[[1]])
    	pval[j] <- wald.test.asreml(mod, list(list(man, "zero")))$zres$zpval
	mrkli <- which.max(map[[chr[j]]]*(map[[chr[j]]]<=cm[j]))
	if (length(map[[chr[j]]])>1) {
	if (mrkli==length(map[[chr[j]]])) mrkli <- mrkli-1
	fmrkl[j] <- names(map[[chr[j]]])[mrkli]
	fmrkr[j] <- names(map[[chr[j]]])[mrkli+1]
	} else fmrkl[j] <- fmrkr[j] <- names(map[[chr[j]]])[mrkli]
    }

    ## these will stay the same - just get the values out separately 
    effect <- t(matrix(round(effect,2), nrow=nqtl, ncol=n.founders, byrow=T))
    se <- t(matrix(round(se,3), nrow=nqtl, ncol=n.founders, byrow=T))
    eff3 <- paste("Effect_",f3,sep="")
    se3 <- paste("SE_",f3,sep="")
    if (n.founders==4)
	table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], "Wald"=round(wald,2), "df"=degf, "pvalue"=signif(pval,3))
    else if (n.founders==8)
	table <- data.frame("Chr"=chr, "Pos"=cm, "LeftMrk"=fmrkl, "RightMrk"=fmrkr, effect[1,], se[1,], effect[2,], se[2,], effect[3,], se[3,], effect[4,], se[4,], effect[5,], se[5,], effect[6,], se[6,], effect[7,], se[7,], effect[8,], se[8,], "Wald"=round(wald,2), "df"=degf, "pvalue"=signif(pval,3))

    names(table)[seq(5, 4+(2*n.founders), 2)] <- eff3
    names(table)[seq(6, 4+(2*n.founders), 2)] <- se3
  }

  output$table <- table
  output$call <- match.call()
  output$FullModel <- mod
  output
}

