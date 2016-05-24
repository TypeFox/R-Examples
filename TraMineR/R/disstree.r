#################
##DisTreeNode
#################
# Donnée accessible
# predictor : liste des prédicteurs
# dissmatrix : matrices des dissimilarités internes
#

#Donnée interne
# split : predicteur choisi(NULL pour noeux terminaux)
# vardis : variabilité interne
# children : noeud enfant (NULL pour noeux terminaux)
# ind: liste des index des individus du noeuds.
# depth: profondeur du noeud
# label: label du noeud (valeurs du predicteur
# R2: R2 du split, (NULL pour noeux terminaux)
#
#source("disstree.dissassoc.R")
.localstuffDissTree <- new.env()
.localstuffDissTree$DTNnodeCounter <- as.integer(1)
DTNInit <- function(ind, vardis, depth, dmat, weights) {
	##cat("Disscenter\n")
	dc <- .Call(TMR_tmrWeightedInertiaContrib, dmat, as.integer(ind),as.double(weights))
	#dc <- .Call("tmrinertiacontrib", dmat, as.integer(ind), PACKAGE="TraMineR")
	medoid <- ind[which.min(dc)]
	info <- list(depth=depth, vardis=vardis, n=sum(weights[ind]), medoid=medoid)
	info[["ind"]] <- ind
	node <- list(id=.localstuffDissTree$DTNnodeCounter, split = NULL, kids = NULL, info = info)
	.localstuffDissTree$DTNnodeCounter <- as.integer(.localstuffDissTree$DTNnodeCounter + 1)
	class(node) <- c("DissTreeNode", class(node))
	return(node)
}

DTNsplit <- function(varindex, index,  prob, info, labels=NULL, breaks=NULL, naGroup=NULL){
	if(is.null(labels)){
		if(is.null(breaks)){
			stop(" [!] Programming error!!!!!")
		}
	}
	retsplit <- list(varindex=varindex, index=index, prob=prob, info=info, breaks=breaks, naGroup=naGroup, labels=labels)
	class(retsplit) <- c("DissTreeSplit", class(retsplit))
	return(retsplit)
}

	

###########################
## Retrieve leaf belonging
###########################

disstreeleaf <- function(tree, label=FALSE) {
	if (inherits(tree, "disstree")) {
		root <- tree$root
	}else{
		root <- tree
		label <- FALSE
	}
	if (!inherits(root, "DissTreeNode")) {
		stop("tree should be a DissTreeNode object")
	}
	DTNnodelabels <- function(dtl){
		split_s <- function(sp){
			formd <- function (x){
				return(format(x, digits =getOption("digits")-2))
			}
			str_split <- character(2)
			vname <- colnames(tree$data)[sp$varindex]
			if (!is.null(sp$breaks)) {
				str_split[1] <- paste("<=", formd(sp$breaks))
				str_split[2] <- paste(">", formd(sp$breaks))
			}
			else {
				str_split[1] <- paste0("[", paste(sp$labels[sp$index==1], collapse=", "),"]")
				str_split[2] <- paste0("[", paste(sp$labels[sp$index==2], collapse=", "),"]")
			}
			if(!is.null(sp$naGroup)){
				str_split[sp$naGroup] <- paste(str_split[sp$naGroup], "with NA")
			}
			return(paste(vname, str_split))
		}
		labelEnv <- new.env()
		labelEnv$label <- list()
		addLabel <- function(n1, n2, val){
			id1 <- as.character(n1$id)
			id2 <- as.character(n2$id)
			labelEnv$label[[id2]] <- c(labelEnv$label[[id1]], val)
		}
		nodeRec <- function(node){
			if(!is.null(node$split)){
				spl <- split_s(node$split)
				addLabel(node, node$kids[[1]], spl[1])
				addLabel(node, node$kids[[2]], spl[2])
				nodeRec(node$kids[[1]])
				nodeRec(node$kids[[2]])
			}
		}
		nodeRec(tree$root)
		l2 <- list()
		for(nn in names(labelEnv$label)){
			l2[[nn]] <- paste0(labelEnv$label[[nn]], collapse=" & ")
		}
		return(factor(factor(dtl, levels=as.numeric(names(l2)), labels=as.character(l2))))
	}
	categorie <- rep(-1, length(root$ind))
	if(label){
		return(DTNnodelabels(DTNdisstreeleaf(root, categorie)))
	}
	return(DTNdisstreeleaf(root, categorie))
}

###########################
## Internal recursion
###########################
DTNdisstreeleaf <- function(node, co) {
	if (is.null(node$kids)) {
		co[node$info$ind] <- node$id
		return(co)
	} else {
		co1 <- DTNdisstreeleaf(node$kids[[1]], co)
		return(DTNdisstreeleaf(node$kids[[2]], co1))
	}
}
DTNaddCovariateSplitschedule <- function(tree) {
	treeSize <- function(node){
		if (is.null(node$kids)) {

			return(1)
		}else{
			return( treeSize(node$kids[[1]])+ treeSize(node$kids[[2]])+1)
		}
	}
	
	trsize <- treeSize(tree$root)
	treeEnv <- environment()
	treeEnv$SCexpl <- numeric(trsize)
	NodeCovariate <- function(node, parent=NULL){
		if(is.null(parent)){ ## root node
			treeEnv$SCexpl[as.character(node$id)] <- node$info$vardis*node$info$n
		}
		else{
			SCtot <- parent$info$vardis*parent$info$n
			SCexpl <- parent$split$info$R2*SCtot
			treeEnv$SCexpl[as.character(node$id)] <- SCexpl
		}
		if (!is.null(node$kids)) {
			NodeCovariate(node$kids[[1]], node)
			NodeCovariate(node$kids[[2]], node)
		}
	}
	NodeCovariate(tree$root)
	treeEnv$depth <- as.integer(factor(1-rank(treeEnv$SCexpl), ordered=TRUE))
	
	NodeDepth <- function(node){
		node$info$splitschedule <- treeEnv$depth[as.character(node$id)==names(treeEnv$SCexpl)]
		if (!is.null(node$kids)) {
			node$kids[[1]] <- NodeDepth(node$kids[[1]])
			node$kids[[2]] <- NodeDepth(node$kids[[2]])
		}
		return(node)
	}
	tree$root <- NodeDepth(tree$root)
	return(tree)
}


###########################
## disstree main function
###########################
disstree <- function(formula, data=NULL, weights=NULL, minSize=0.05, maxdepth=5, R=1000, pval=0.01, object =NULL, weight.permutation="replicate", squared=FALSE, first=NULL) {
	##formula.call <- formula
	tterms <- terms(formula)
	dissmatrix <- eval(formula[[2]], data, parent.frame()) # to force evaluation
	formula[[2]] <- NULL
	## Model matrix from forumla
	predictor <- as.data.frame(model.frame(formula, data, drop.unused.levels = TRUE, na.action=NULL))
	tree <- DTNdisstree(dissmatrix=dissmatrix, predictor=predictor, terms=tterms, 
						weights=weights, minSize=minSize, maxdepth=maxdepth, R=R, 
						pval=pval, object =object, weight.permutation=weight.permutation, 
						squared=squared, first=first)
	return(tree)
	
}
DTNdisstree <- function(dissmatrix, predictor, terms, weights=NULL, minSize=0.05, maxdepth=5, R=1000, pval=0.01, object =NULL, weight.permutation="replicate", squared=FALSE, first=NULL) {
	
	
	if (inherits(dissmatrix, "dist")) {
		dissmatrix <- as.matrix(dissmatrix)
 	}
	if (squared) {
		dissmatrix <- dissmatrix^2
	}
	nobs= nrow(dissmatrix)
	if (nobs!=nrow(predictor)) {
		stop(" [!] dissimilarity matrix and data should be of the same size")
	}
	
	## Allow integer weights for replicates
	if(is.null(weights)) {
		weights <- as.double(rep(1, nobs))
		weight.permutation <- "none"
	}
	if(weight.permutation %in% c("replicate", "rounded-replicate")) {
		rounderror <- sum(abs(round(weights, 0) - weights))
		if(rounderror>0){
			if (weight.permutation=="replicate") {
				stop(" [!] To permute replicate, you should specify integer weights")
			}
			message(" [>] Weigths loss : ", rounderror, " (", (rounderror/sum(weights)), ")")
			weights <- round(weights, 0)
		}
		weights <- as.integer(weights)
	}
	else {
		weights <- as.double(weights)
	}
	
	if(!is.null(first)){
		if(is.character(first)){
			message("[>] Using ", first, " as primary split")
			first <- which(names(predictor)==first)
			if(!(first>0&&first<=ncol(predictor))){
				stop(" [!] Unknow ", first, " variable. It should be a character string appearing in the formula")
			}
		}
		else{
			stop(" [!] Unknow ", first, " variable. It should be a character string appearing in the formula")
		}
	}
	pop <- sum(weights)
	if (minSize<1) {
		minSize <- pop*minSize
	}
	if(R<=1){
		pval <- 1
	}
	else if(pval<(1/sum(R))){
		warning(" [!] Minimum possible p-value using ", R, " permutations is ", 1/sum(R), ". Parameter pval (=", pval, ") changed to ", 1/sum(R))
		pval <- 1/sum(R)
	}
	.localstuffDissTree$DTNnodeCounter <- as.integer(1)
	
	vardis <- dissvar(dissmatrix, weights=weights)
	root <- DTNBuildNode(dmat=dissmatrix, pred=predictor, minSize=minSize, ind=1:nobs,
			vardis=vardis, depth=1, nbperm=R, pval=pval, maxdepth=maxdepth, 
			weights=weights, weight.permutation=weight.permutation, first=first)
	tree <- list()
	tree$fitted <- data.frame(disstreeleaf(root))
	names(tree$fitted) <- "(fitted)"
	tree$info <- list(method="disstree", n=pop, parameters= list(minSize=minSize, maxdepth=maxdepth, R=R, pval=pval), object=object, weight.permutation=weight.permutation)
	if(weight.permutation=="none") {
		tree$info$adjustment <- dissassoc(dissmatrix, tree$fitted[,1], R=R, weights=NULL)
	}
	else {
		tree$info$adjustment <- dissassoc(dissmatrix, tree$fitted[,1], R=R, weights=weights, weight.permutation=weight.permutation)
	}
	tree$data <- predictor
	tree$weights <- weights
	tree$terms <- terms
	##tree <- party(root, data=predictor, fitted =fitted, terms = terms(formula.call),  info = info)
	tree$root <- root
	
	class(tree) <- c("disstree", class(tree))
	tree <- DTNaddCovariateSplitschedule(tree)
	return(tree)
}

###########################
## Building node and finding predictor
###########################
DTNBuildNode <- function(dmat, pred, minSize, ind, vardis,
							depth, nbperm, pval, maxdepth, weights, weight.permutation, first=NULL) {
	node <- DTNInit(ind=ind, vardis=vardis, depth=depth, dmat=dmat, weights=weights)
	SCtot <- vardis*node$info$n
	SCres <- SCtot
	bestSpl <- NULL
	## print(SCtot)
	#varnames <- colnames(pred)
	if (depth>=maxdepth) {
			return(node)
	}
	if(!is.null(first)){
		bestSpl <-  DTNGroupFactorBinary(dissmatrix=dmat, currentSCres=SCres, pred=pred[, first], minSize=minSize, varindex=first, ind=ind, weights=weights)
		SCres <- bestSpl$spl$info$SCres
	}
	else {
		for (p in 1:ncol(pred)) {
			## cat("Checking", names(pred[p]), "...\n")
			spl <- DTNGroupFactorBinary(dissmatrix=dmat, currentSCres=SCres, pred=pred[, p], minSize=minSize, varindex=p, ind=ind, weights=weights)
			## print(str(spl))
			if (!is.null(spl) && (is.null(bestSpl) || spl$spl$info$SCres<bestSpl$spl$info$SCres)) {
				bestSpl <- spl
				SCres <- spl$spl$info$SCres
				## cat(varnames[[p]], " Ok", "\n")
			}
		}
	}
	if (is.null(bestSpl)) {
		return(node)
	}
	if (nbperm>1) {
		spval <- DTNdissassocweighted(dmat=dmat, grp=bestSpl$variable, indiv=ind, weights=weights, R=nbperm, weight.permutation=weight.permutation)
		## print(paste(label, bestSpl$varname, spval))
		if (spval>pval){
			return(node)
		}
		bestSpl$spl$info$pval <- spval
		
	}

	node$split <- bestSpl$spl
	node$split$info$R2 <- 1-(SCres/SCtot)
	node$surrogates <- bestSpl$sur
	left <- DTNBuildNode(dmat=dmat, pred=as.data.frame(pred[bestSpl$variable, ]),
			minSize=minSize, ind=ind[bestSpl$variable], vardis=bestSpl$spl$info$lvar, depth=depth+1,
			nbperm=nbperm, pval=pval, maxdepth=maxdepth, weights=weights, weight.permutation=weight.permutation)
	
	right <- DTNBuildNode(dmat=dmat, pred=as.data.frame(pred[!bestSpl$variable, ]),
			minSize=minSize, ind=ind[!bestSpl$variable], vardis=bestSpl$spl$info$rvar, depth=depth+1,
			nbperm=nbperm, pval=pval, maxdepth=maxdepth, weights=weights, weight.permutation=weight.permutation)
	
	node$kids <- list(left, right)
	return(node)
}

###########################
## Find best binary partition
###########################
DTNGroupFactorBinary <- function(dissmatrix, currentSCres, pred, minSize, varindex, ind, weights) {
	totpop <- sum(weights[ind])
	grp <- factor(pred, ordered=(is.ordered(pred) || is.numeric(pred)))
	lgrp <- levels(grp)
	if (length(lgrp)<=1) {
		#print("Small group return")
		return(NULL)
	}
	nbGrp <- length(lgrp)
	has.na <- FALSE
	llgrp <- lgrp
	## Here we add a group for missing values
	if (sum(is.na(grp))>0) {
		nbGrp <- length(lgrp)+1
		has.na <- TRUE
		llgrp[nbGrp] <- "<Missing>"
	}

	grpCond <- list()
	grpSize <- numeric(length=nbGrp)
	grpSize[] <- 0
	for (i in 1:length(lgrp)) {
		## on crée le groupe en question
		grpCond[[i]] <- (grp==lgrp[i])
		grpCond[[i]][is.na(grpCond[[i]])] <- FALSE
		grpSize[i] <- sum(weights[ind[grpCond[[i]]]])
	}
	## Treating missing values
	if (has.na) {
		grpCond[[nbGrp]] <- is.na(grp)
		grpSize[nbGrp] <- sum(weights[ind[grpCond[[nbGrp]]]])
	}
	inertiaMat <- matrix(0, nrow=nbGrp, ncol=nbGrp)
	for (i in 1:(nbGrp-1)) {
		grpindiv1 <- ind[grpCond[[i]]]
		for (j in (i+1):nbGrp) {
			grpindiv2 <- ind[grpCond[[j]]]
			#cat("Inter Inertia", i, "  ", j, "\n")
			r <- .Call(TMR_tmrWeightedInterInertia, dissmatrix, as.integer(grpindiv1),
						as.integer(grpindiv2), as.double(weights))
				## using only one half of the matrix
			inertiaMat[j, i] <- r
		}
		#cat("Inertia", i,"\n")
		r <- .Call(TMR_tmrWeightedInertiaDist, dissmatrix, as.integer(nrow(dissmatrix)), 
					as.integer(FALSE), as.integer(grpindiv1), as.double(weights), 
					as.integer(FALSE))
		
		inertiaMat[i, i] <- r*grpSize[i]
	}
	## FIXME This step is missing in the loop
	#cat("Inertia", nbGrp,"\n")
	r <- .Call(TMR_tmrWeightedInertiaDist, dissmatrix, as.integer(nrow(dissmatrix)), 
					as.integer(FALSE), as.integer(ind[grpCond[[nbGrp]]]), as.double(weights), 
					as.integer(FALSE))
	inertiaMat[nbGrp, nbGrp] <- r*grpSize[nbGrp]
	## Computing residuals
	#print(inertiaMat)
	SCres <- sum(diag(inertiaMat)/grpSize)
	## An error here might be due to weights=0
	## Check for grpSize==0?
	# if(is.na(SCres)||is.na(currentSCres)){
		# print("SCres")
		# print(SCres)
		# print("inertiaMat")
		# print(inertiaMat)
		# print("grpSize")
		# print(grpSize)
		# print("llgrp")
		# print(llgrp)
		# print("currentSCres")
		# print(currentSCres)
	# }
	## cat("SCres max: ", SCres)
	if (SCres>currentSCres){
		## print("Unsplittable")
		## print(SCres)
		## print(currentSCres)
		
		return(NULL)
	}
	## Fonction to comput inertia based on inertiaMat
	inertiaFunction <- function(inertiaMat, co, pop) {
		## Take care to add one
		## return(((sum(inertiaMat[co, co])+sum(diag(inertiaMat[co, co])))/2)/pop)
		## New way, inertiaMat is triangular -> we can just sum the matrix
		return(sum(inertiaMat[co, co])/pop)
	}
	bestSCres <- currentSCres
	bestRegroup <- NULL
	allgroups <- 1:(nbGrp)
	if (is.ordered(grp)) {
		maxGrp <- nbGrp-1
	} else {
		maxGrp <- ceiling(nbGrp/2)
	}
	for (p in 1:maxGrp) {
		if (is.ordered(grp)) {
			 combi <- list()
			 combi[[1]] <- 1:p
			 if (has.na) {
				 combi[[2]] <- c(1:p, nbGrp)
			 }
		} else {
			 combi <- combn(nbGrp, p, simplify=FALSE)
		}
		for (co in combi) {
			
			 popc <- sum(grpSize[co])
			 popothc <- totpop-popc
			 ## cat("CO:", co, "\n")
			 ## cat("totpop", totpop, popc, popothc, minSize,"\n")
			 if (popc>minSize && popothc>minSize) {
					othc <- allgroups[!(allgroups %in% co)]
					ico <- inertiaFunction(inertiaMat, co, popc)
					iothc <- inertiaFunction(inertiaMat, othc, popothc)
					SCres <- ico+iothc
					## cat("SCres, ", SCres, co, "\n")
					if (SCres<bestSCres) {
						bestSCres <- SCres
						bestRegroup <- list(co=co, othc=othc, ico=ico,
								iothc=iothc, popc=popc, popothc=popothc)
					}
			 }
		}
	}
	if (is.null(bestRegroup)){
		## print("No split found")
		return(NULL)
	}
	prob <- c(bestRegroup$popc, bestRegroup$popothc)/totpop
	
	ret <- list()
	info <- list(
			lpop=bestRegroup$popc,
			rpop=bestRegroup$popothc,
			lvar=bestRegroup$ico/bestRegroup$popc,
			rvar=bestRegroup$iothc/bestRegroup$popothc,
			SCres=bestSCres
	)
	ret$variable <- (grp %in% lgrp[bestRegroup$co])
	if(is.numeric(pred)){
		breaks <- as.numeric(c(llgrp[max(bestRegroup$co[bestRegroup$co<nbGrp])]))
		ret$spl <- DTNsplit(varindex, breaks = breaks, index = as.integer(1:2),  
             prob = prob, info = info)
	}
	else{
		allgroups <- 1:length(lgrp)
		index <- as.integer(ifelse(allgroups %in% bestRegroup$co, 1, 2))
		ret$spl <- DTNsplit(varindex, index = index,  
             prob = prob, info = info, labels=lgrp)
	}
	if (has.na) {
		#index <- ifelse(nbGrp %in% bestRegroup$co, 1, 2)
		#index <- as.integer(c(index, index))
		
		ret$variable[is.na(grp)] <- (nbGrp %in% bestRegroup$co)
		ret$spl$naGroup <- ifelse(nbGrp %in% bestRegroup$co, 1, 2)
	}
	return(ret)
}
