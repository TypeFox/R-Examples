permmodels <- function(model, data, block=NULL, group=NULL, covariate=NULL, nsim=4999, check=FALSE, exact=FALSE, fo=NULL, prt=TRUE) {

  options(scipen=6)
  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  
  if (any(!is.null(block), !is.null(group))) data <- data[do.call(order, data[, c(block, group), drop=FALSE]),]
  rownames(data) <- NULL
  n <- nrow(data)
  tm <- terms(model)
  if (class(model)[1]=="lmerMod") cls <- sapply(all.vars(formula(model, fixed.only=TRUE)),function(x) class(model.frame(model)[[x]])[1]) else cls <- attr(tm,"dataClasses") 
  factn <- names(which(cls=="factor"))
  response <- names(cls)[1]  
  permVar <- c(response, covariate)
  
  if (all(is.null(block), is.null(group))) {
  	permIndex <- replicate(nsim, sample.int(n))
  }else if (all(!is.null(block), is.null(group))) {
    Bsplit <- split(seq_len(n), as.character(data[, block]))    # split 1:n seq into nB blocks
    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) x[sample.int(length(x))]), recursive = FALSE, use.names = FALSE))
  }else if (all(is.null(block), !is.null(group))){   # data=data; group="group"; nsim=5
    permIndex <- replicate(nsim, {Gsplit <- split(seq_len(n), as.character(data[, group]))   # split 1:n seq into separate groups
      nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
      unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)})          # permute the order of the whole groups and repeat process nsim times
  }else{
    ind <- seq_len(n)
    names(ind) <- as.character(data[, group])
    Bsplit <- split(ind, as.character(data[, block]))    # split 1:n seq into blocks

    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) {    # within each block
        Gsplit <- split(x, names(x))  # split each block seq into groups
        nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
        unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)          # permute the order of groups
      }), recursive = FALSE, use.names = FALSE))
  }

  if (all(is.null(block), is.null(group), length(factn) > 0)) {
	colnIndex0 <- as.numeric(interaction(data[, factn], sep=":"))
	colnIndex <- apply(permIndex, 2, function(x) paste(colnIndex0[x], collapse=""))	
  }else{
    colnIndex <- apply(permIndex, 2, function(x) paste(x, collapse=""))
  }
  colnames(permIndex) <- colnIndex
  permIndex <- permIndex[, unique(colnIndex)]
  colnames(permIndex) <- NULL

  if (check) {
    ndata <- cbind(data, data[, permVar, drop=FALSE][permIndex[,1],])
		if (length(permVar)==1) colnames(ndata)[dim(ndata)[2]] <- paste(permVar, "1", sep="")
      	ndata$index1 <- permIndex[,1]
      	ndata <- cbind(ndata, data[, permVar, drop=FALSE][permIndex[,2],])
		if (length(permVar)==1) colnames(ndata)[dim(ndata)[2]] <- paste(permVar, "2", sep="")
      	ndata$index2 <- permIndex[,2]
    rownames(ndata) <- NULL
   	print(ndata)
   	ques <- paste("Is this a right permutation? (Y or N)  ")
    answer <- readline(ques)
   	if (answer=="N") stop("Please create a correct permutation for your data!")
 	} # end of  if (check)

  if (all(is.null(block), is.null(group), length(factn) > 0)) {
  	permdata <- apply(permIndex, 2, function(x) {
      data[, factn] <- data[x, factn]
      return(data)})
  }else{
    permdata <- apply(permIndex, 2, function(x) {
    data[, permVar] <- data[x, permVar]
    return(data)})
  }
    
  permmod <- lapply(permdata, function(x) {
    rfitmodel <- try(update(model, data=x), TRUE)
  	if (class(rfitmodel)[1]==class(model)[1]) {
    	mp <- mymodelparm(rfitmodel)
      aT <- anova(rfitmodel)
      return(list(mp,aT))
      }
  })
  permmod <- permmod[!(sapply(permmod,is.null))]
  nsim <- length(permmod)
  
  if (class(model)[1]!="aov") {
    tTable <- function(x) { # function to construct t-table for the model
    	if (class(model)[1]=="lmerMod"){
    	    mp <- mymodelparm(x)
    		tTable <- cbind(mp$coef, sqrt(base::diag(mp$vcov)), mp$coef/sqrt(base::diag(mp$vcov)))
    		colnames(tTable) <- c("Estimate", "Std. Error", "t value")
    		return(round(tTable, 4))		
    	}else{
    		summ <- summary(x)
    		if (class(x)[1] %in% c("lme", "gls")) {
    			tTable <- summ$tTable
    		}else{
    			tTable <- coef(summ)
    		}
    		return(tTable)
    	}
    }
  
    Tvalue <- colnames(tTable(model))[grep("t.value", colnames(tTable(model)))]  
    
    if (nrow(tTable(model))==1) {
      Tper.p <- (sum(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(tTable(model)[, Tvalue]),6))+!exact)/(nsim+!exact)
     }else{
      Tper.p <- (rowSums(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(tTable(model)[, Tvalue]), 6))+!exact)/(nsim+!exact)
    }
  	if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(round(cbind(tTable(model), "Perm_p_value"=Tper.p),4))
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nsim), "permutations.\n")
  	}
  }
  
  Fvalue <- colnames(anova(model))[grep("F.value",colnames(anova(model)))]
  if (nrow(anova(model))==1) {
    Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(anova(model)[, Fvalue], 6))+!exact)/(nsim+!exact)
  }else{
    Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(anova(model)[, Fvalue], 6))+!exact)/(nsim+!exact)
  }
  ANOVA <- round(cbind(anova(model), "Perm_p_value"=Fper.p),4)
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
  cat("\nANOVA:\n")
  print(ANOVA)
  cat("\nNote: Perm_p_value of F test is obtained using", sQuote(nsim), "permutations.\n\n")
  }
  permlist <- vector("list", nsim)
  for (i in 1:nsim) {
    permlist[[i]] <- permmod[[i]][[1]]
  }
  return(invisible(list("permlist"=permlist, "ANOVA"=ANOVA)))
} 
