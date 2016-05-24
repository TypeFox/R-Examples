### DList : list of data set matrices; names should be set
### GList : list of gene sets; names should be set // alternatively gmt file also allowed.
### isParallel : if multiple core parallel processing will be used (default: TRUE)
### nCores : how many cores will be used (default: all in unix-like os and 2 in windows)
### useCache : if save GList as cache for next use (default: TRUE)
### filterGenes : whether to use gene filtering (recommended to reduce dimension for fast computation)
MetaQC <- function(DList, GList, isParallel=FALSE, nCores=NULL, useCache=TRUE, filterGenes=TRUE, 
		maxNApctAllowed=.3, cutRatioByMean=.4, cutRatioByVar=.4, minNumGenes=5,
		verbose=FALSE, resp.type=c("Twoclass", "Multiclass", "Survival")) {
	.p <- proto(expr = {
				.verbose <- verbose
				.Names <- names(DList)
				#.DListF0 <- NULL
				.DListF <- NULL 
				.nFeatures <- NULL
				.excluded <- NULL
				.GList <- NULL
				.GListIdx <- NULL
				.cutRatioByMean <- cutRatioByMean
				.cutRatioByVar <- cutRatioByVar
				.minNumGenes <- minNumGenes
				.isParallel <- isParallel
				.useCache <- useCache
				.method.cor <- "pearson"
				.l.norm <- 2
				.PvalOfScores <- NULL #-log pval of scores
				.DistOfStudies <- NULL
				.IScores <- NULL
				.maxNApctAllowed <- maxNApctAllowed
				.PValList <- NULL
				.PValMat <- NULL
				.PValMat0 <- NULL
				.CQCgScores <- NULL
				.PathPValMat <- NULL
				.CQCpScores <- NULL
				.Scores <- NULL
				.summary <- NULL
				.workers <- NULL
				.AQCgScores <- NULL
				.AQCpScores <- NULL
				.resp.type <- "Twoclass"
				
				if(isParallel) {
					if(!is.null(nCores))
						options(cores=nCores)
					if(.Platform$OS.type == "unix") {
						requireAll("doMC")
						registerDoMC()
					} else { #windows
						requireAll("doSNOW")
						if(is.null(getOption('cores')))
							options(cores=2)
						registerDoSNOW(makeCluster(getOption('cores'), type = "SOCK"))
					}
				}
				
				Cleanup <- function(.) {
					warning("This function was deprecated.\n",
							"No more needed to call this function.\n")
				}
				
				.Initialize <- function(., .DList, .GList, .filterGenes, .resp.type) {
#					if(length(.DList)<4)
#						stop("DList should have more than 3 data sets")
					if(any(length(names(.DList))==0))
						stop("DList should have names for all data sets")
					if(any(duplicated(names(.DList))))
						stop("DList should have unique names for all data sets")
					
					.$.resp.type <- .resp.type
					
					.$.DListF <- foreach(d=iter(.DList)) %do% {
						if(!is.list(d)) {
							d <- list(x=d, y=colnames(d))
						} else {
							if(.resp.type %in% c("Twoclass", "Multiclass")) {
								stopifnot(!is.null(d$y) && !is.null(colnames(d$x)))
								colnames(d$x) <- d$y
								if(!is.null(d$geneid))
									rownames(d$x) <- d$geneid
							} else if(.resp.type == "Survival") {
								stopifnot(!is.null(d$y) && !is.null(d$censoring.status))
								if(!is.null(d$geneid))
									rownames(d$x) <- d$geneid
								mode(d$y) <- "numeric"
								mode(d$censoring.status) <- "integer"
							}
						}
						mode(d$x) <- "numeric"
						d
					}
					names(.$.DListF) <- .$.Names
					
					#.$.DListF <- .DList
					if(.filterGenes)
						.$.FilterGenes()
					
					.$.GList <- .GList
				}
				
				.FilterGenes <- function(.) {
					.$.DListF <- foreach(d=iter(.$.DListF)) %dopar% {
						dat <- d$x
						if(any(is.na(rownames(dat))))
							dat <- dat[-which(is.na(rownames(dat))),]
						
						dat <- dat[rowSums(is.na(dat))<=ncol(dat)*.$.maxNApctAllowed,]
						
						.RMrank <- rank(rowMeans(dat, na.rm=T))
						dat <- dat[-order(.RMrank)[1:floor(nrow(dat)*.$.cutRatioByMean)],]
						.RVrank <- rank(rowVars(dat, na.rm=T))
						dat[-order(.RVrank)[1:floor(nrow(dat)*.$.cutRatioByVar)],]
						list(x=dat, y=d$y, censoring.status=d$censoring.status, geneid=d$geneid)
					}
					names(.$.DListF) <- .$.Names
				}
				
				.ConvertToGeneSetIdx <- function(., .DListF=.$.DListF, .GList=.$.GList, .minNumGenes=.$.minNumGenes) {
					.DGList <- foreach(d=iter(.DListF), .packages=c("foreach","iterators")) %dopar% { #each data set wrap each pathway
						.res <- foreach(g=iter(.GList)) %do% {
							.gs <- na.omit(match(g, rownames(d$x)))
							if(length(.gs)<.minNumGenes) return(NA)
							return(as.integer(.gs))
						}
						names(.res) <- names(.GList)
						if(sum(is.na(.res))>0)
							.res <- .res[-which(is.na(.res))]
						return(.res)
					}
					names(.DGList) <- names(.DListF)
					
					stopifnot(length(.DGList)>0)
					
					return(.DGList)
				}
				
				.CalcPval <- function(., .B, .DListF=.$.DListF) {
					
					stopifnot(length(.DListF)==length(.$.GListIdx))
					
					.PvalOfScores <- foreach(i=1:length(.DListF), .combine=c) %do% { 
						.d <- .DListF[[i]]$x
						.isNA <- any(is.na(.d))
						
						.pathMat <- matrix(0, nrow(.d), length(.$.GListIdx[[i]]))
						for(jj in 1:ncol(.pathMat)) {
							.pathMat[.$.GListIdx[[i]][[jj]], jj] <- 1
						}
						rownames(.pathMat) <- 1:nrow(.d) 
						colnames(.pathMat) <- names(.$.GListIdx[[i]])
						
						if(!.isNA) {
							.d <- t(scale(t(.d)))
							.d <- .d %*% t(.d) / (ncol(.d)-1)
						}
						
						.Scores <- foreach(w=1:ncol(.pathMat), .combine=c) %dopar% {
							.g <- as.numeric(rownames(.pathMat)[.pathMat[,w]==1])
							if(!.isNA)
								mean(abs(as.dist(.d[.g,.g]))^.$.l.norm)^(1/.$.l.norm)
							else
								mean(as.dist(abs(cor(t(.d[.g,]), use="pairwise.complete.obs", method=.$.method.cor)^.$.l.norm)))^(1/.$.l.norm)
						}
						
						.pathMat.dupCols <- duplicated(colSums(.pathMat))
						
						.pathList <- foreach(j=iter(.pathMat[,!.pathMat.dupCols],by="col")) %dopar% {
							which(j==1)
						}
						names(.pathList) <- sapply(.pathList,length)
						
						#performance significantly degraded by length of .pathList, which is the number of unique pathway sizes
						.ScoresNullDist <- foreach(b=1:.B, .combine=rbind, .export="printLog") %dopar% {
							
							.g <- sample(nrow(.pathMat))							
							.res <- sapply(.pathList, function(w) {
										w <- .g[w]
										if(!.isNA)
											mean(abs(as.dist(.d[w,w]))^.$.l.norm)^(1/.$.l.norm)
										else
											mean(as.dist(abs(cor(t(.d[w,]), use="pairwise.complete.obs", method=.$.method.cor)^.$.l.norm)))^(1/.$.l.norm)
									})
							
							if(b%%1000==0) printLog(paste("(",.$.Names[i],")",b), .$.verbose)
							
							return(.res)
						}
						
						if(sum(.pathMat.dupCols) > 0) {
							.dupNum <- colSums(.pathMat)[duplicated(colSums(.pathMat))]
							.ScoresNullDist <- cbind(.ScoresNullDist, sapply(.dupNum, function(dn) .ScoresNullDist[,as.character(dn)]))
						} 
						.ScoresNullDist <- rbind(.Scores, .ScoresNullDist)
						
						.ScoresNullDist <- foreach(w=iter(.ScoresNullDist, by="col"), .combine=cbind) %dopar% {
							(.B+2 - rank(w)) / (.B+1)
						}
						colnames(.ScoresNullDist) <- colnames(.pathMat)
						
						.ScoresNullDist <- GetEWPval(.ScoresNullDist)
						
						rank(.ScoresNullDist)[1] / (.B+1)
						
					}
					gc()
					return(.PvalOfScores)
				}
				
				EQC <- function(., nPathCut=NULL, .B=1e4) { 
					printLog("EQC Started", .$.verbose)
					.$.GListIdx <- .$.ConvertToGeneSetIdx()
					
					if(!is.null(nPathCut)) {
						
						.Scores <- foreach(i=1:length(.$.GListIdx)) %do% { 
							foreach(g=iter(.$.GListIdx[[i]]), .combine=c) %dopar% {
								.d <- .$.DListF[[i]]$x
								mean(as.dist(abs(cor(t(.d[g,]), use="pairwise.complete.obs", method=.$.method.cor)^.$.l.norm)))^(1/.$.l.norm)
							}
						}
						
						.allPaths <- unique(unlist(lapply(.$.GListIdx,names)))
						.mat <- matrix(NA,length(.allPaths),length(.$.DListF)); rownames(.mat) <- .allPaths; colnames(.mat) <- names(.$.DListF)
						for(i in 1:ncol(.mat)) {
							.mat[match(names(.$.GListIdx[[i]]),rownames(.mat)),i] <- .Scores[[i]]
						}
						
						.matPval <- foreach(w=iter(.mat, by="col"), .combine=cbind) %dopar% {
							.len <- length(na.omit(w))  
							(.len+1 - rank(w,na.last='keep')) / .len
						}
						.EWPval <- GetEWPval(.matPval)
						
						.mat <- .mat[order(.EWPval),]
						
						nPathCutLimit <- min(apply(.mat,2,function(x)sum(!is.na(x))))
						
						for(i in 1:length(.$.GListIdx)) {
							.selPaths <- rownames(na.omit(.mat[,i,drop=FALSE]))[1:min(nPathCut,nPathCutLimit)]
							.$.GListIdx[[i]] <- .$.GListIdx[[i]][na.omit(match(.selPaths, names(.$.GListIdx[[i]])))]
						}
					}
					.$.PvalOfScores <- .$.CalcPval(.B=.B)
					names(.$.PvalOfScores) <- names(.$.DListF)
					
					printLog("EQC Finished", .$.verbose)
					return(.$.PvalOfScores)
				}
				
				AQCg <- function(., .cutoff=.05, .adjust=TRUE) {
					if(is.null(.$.CQCgScores))
						.$CQCg()
					printLog("AQCg Started", .$.verbose)
					.PValMat <- if(is.null(.$.excluded)) .$.PValMat else .$.PValMat[,-.$.excluded]
					
					.$.AQCgScores <- foreach(i=1:ncol(.PValMat), .combine=c, .export="GetEWPval") %dopar% {
						.dat <- .PValMat[which(!is.na(.PValMat[,i])),]
						.dat <- .dat[rowSums(!is.na(.dat))>=3,]
						#print(paste(i,nrow(.dat)))
						.reduced  <- GetEWPval(.dat[,-i])
						if(.adjust) {
							.reduced <- p.adjust(.reduced, method="BH")
							.obs <- p.adjust(.dat[,i], method="BH")
						} else {
							.obs <- .dat[,i]
						}
						
						.reduced <- ifelse(.reduced<=.cutoff, TRUE, FALSE)
						.obs <- ifelse(.obs<=.cutoff, TRUE, FALSE)						
						.confu <- table(.reduced, .obs)
						if(dim(.confu)[1]<2) .confu <- rbind(.confu,c(0,0))
						if(dim(.confu)[2]<2) .confu <- cbind(.confu,c(0,0))
						
						fisher.test(.confu, alternative="g")$p.value
					}
					names(.$.AQCgScores) <- colnames(.PValMat)
					.$.AQCgScores <- ifelse(.$.AQCgScores < .Machine$double.xmin, .Machine$double.xmin, .$.AQCgScores)
					printLog("AQCg Finished", .$.verbose)
				}
				
				AQCp <- function(., .cutoff=.05, .adjust=TRUE, .GList="c2.all.v3.0.symbols.rda") {
					if(is.null(.$.CQCpScores))
						.$CQCp(.GList=.GList)
					printLog("AQCp Started", .$.verbose)
					.PathPValMat <- if(is.null(.$.excluded)) .$.PathPValMat else .$.PathPValMat[,-.$.excluded]
					
					.$.AQCpScores <- foreach(i=1:ncol(.PathPValMat), .combine=c, .export="GetEWPval") %dopar% {
						.dat <- .PathPValMat[which(!is.na(.PathPValMat[,i])),]
						.dat <- .dat[rowSums(!is.na(.dat))>=3,]
						.reduced  <- GetEWPval(.dat[,-i])
						if(.adjust) {
							.reduced <- p.adjust(.reduced, method="BH")
							.obs <-  p.adjust(.dat[,i], method="BH")
						} else {
							.obs <- .dat[,i]
						}
						
						.reduced <- ifelse(.reduced<=.cutoff, TRUE, FALSE)
						.obs <- ifelse(.obs<=.cutoff, TRUE, FALSE)
						.confu <- table(.reduced, .obs)
						if(dim(.confu)[1]<2) .confu <- rbind(.confu,c(0,0))
						if(dim(.confu)[2]<2) .confu <- cbind(.confu,c(0,0))
						
						fisher.test(.confu, alternative="g")$p.value
					}
					names(.$.AQCpScores) <- colnames(.PathPValMat)
					.$.AQCpScores <- ifelse(.$.AQCpScores < .Machine$double.xmin, .Machine$double.xmin, .$.AQCpScores)
					printLog("AQCp Finished", .$.verbose)
				}
				
				CQCg <- function(.) {
					printLog("CQCg Started", .$.verbose)
					if(is.null(.$.PValMat)) {
						if(is.null(.$.PValList)) {
							.$.PValList <- foreach(dat=iter(.$.DListF), .export="GetPVal") %dopar% {
								GetPVal(dat, .$.resp.type)
							}
							names(.$.PValList) <- .$.Names
						}
						
						.allGNames <- union.rec(lapply(.$.PValList,names))
						
						.$.PValMat <- foreach(pv=iter(.$.PValList), .combine=cbind) %dopar% {
							pv[match(.allGNames, names(pv))]
						}
						colnames(.$.PValMat) <- .$.Names
						rownames(.$.PValMat) <- .allGNames
					}
					
					.PValMat <- if(is.null(.$.excluded)) .$.PValMat else .$.PValMat[,-.$.excluded]
					
					.$.CQCgScores <- foreach(i=1:ncol(.PValMat), .combine=c, .export="GetEWPval") %dopar% {
						.dat <- .PValMat[which(!is.na(.PValMat[,i])),]
						.dat <- .dat[rowSums(!is.na(.dat))>=3,]
						.reduced  <- GetEWPval(.dat[,-i])
						.obs <-  .dat[,i]
						suppressWarnings(cor.test(.reduced, .obs, method="spearman", alternative="g")$p.value)
					}
					names(.$.CQCgScores) <- colnames(.PValMat)
					.$.CQCgScores <- ifelse(.$.CQCgScores < .Machine$double.xmin, .Machine$double.xmin, .$.CQCgScores)
					printLog("CQCg Finished", .$.verbose)
				}
				
				CQCp <- function(., .GList="c2.all.v3.0.symbols.rda") {
					printLog("CQCp Started", .$.verbose)
					if(is.null(.$.PValMat0)) {
						.PValList <- foreach(dat=iter(.$.DListF), .export="GetPVal") %dopar% {
							GetPVal(dat, .$.resp.type)
						}
						names(.PValList) <- .$.Names
						.allGNames <- union.rec(lapply(.PValList,names))
						
						.$.PValMat0 <- foreach(pv=iter(.PValList), .combine=cbind) %dopar% {
							pv[match(.allGNames, names(pv))]
						}
						colnames(.$.PValMat0) <- .$.Names
						rownames(.$.PValMat0) <- .allGNames
					}
					
					if(is.null(.$.PathPValMat)) {
						load(.GList)
						.GListIdx <- .$.ConvertToGeneSetIdx(.GList=GList)
						.PathPValList <- foreach(ii=1:length(.GListIdx), .packages=c("foreach","iterators")) %dopar% {
							.PathPVal <- foreach(jj=iter(.GListIdx[[ii]]), .combine=c) %do% {
								.gnInPath <- rownames(.$.DListF[[ii]]$x)[jj] #gene names in the pathway
								.gMatched <- sort(match(.gnInPath, rownames(.$.PValMat0)))
								.pvInPath <- .$.PValMat0[.gMatched,ii]
								.pvOutPath <- na.omit(.$.PValMat0[-.gMatched,ii])
								suppressWarnings(ks.test(.pvInPath, .pvOutPath, alternative="greater")$p)
							}
							names(.PathPVal) <- names(.GListIdx[[ii]])
							return(.PathPVal)
						}
						names(.PathPValList) <- .$.Names
						
						.allPathNames <- union.rec(lapply(.GListIdx,names))
						
						.$.PathPValMat <- foreach(pv=iter(.PathPValList), .combine=cbind) %dopar% {
							pv[match(.allPathNames, names(pv))]
						}
						colnames(.$.PathPValMat) <- .$.Names
						rownames(.$.PathPValMat) <- .allPathNames
					}
					
					.PathPValMat <- if(is.null(.$.excluded)) .$.PathPValMat else .$.PathPValMat[,-.$.excluded]
					
					.$.CQCpScores <- foreach(i=1:ncol(.PathPValMat), .combine=c, .export="GetEWPval") %dopar% {
						.dat <- .PathPValMat[which(!is.na(.PathPValMat[,i])),]
						.dat <- .dat[rowSums(!is.na(.dat))>=3,]
						.reduced  <- GetEWPval(.dat[,-i])
						.obs <-  .dat[,i]
						suppressWarnings(cor.test(.reduced, .obs, method="spearman")$p.value)
					}
					names(.$.CQCpScores) <- colnames(.PathPValMat)
					.$.CQCpScores <- ifelse(.$.CQCpScores < .Machine$double.xmin, .Machine$double.xmin, .$.CQCpScores)
					printLog("CQCp Finished", .$.verbose)
				} 
				
				IQC <- function(., .excludedN=NULL) { 
					if(is.null(.$.PvalOfScores)) { 
						.$EQC()
					}
					if(!is.null(.excludedN)) {
						.excluded <- .$.GetExcluded(.excludedN)
						if(length(setdiff(.excluded,.$.excluded))>0) { 
							.$.excluded <- c(.$.excluded, setdiff(.excluded,.$.excluded))
						}
					} 
					printLog("IQC Started", .$.verbose)
					.dist <- .$GetDistOfStudies(.$.Names[.$.excluded])
					.distMat <- as.matrix(.dist)
					.IScores0 <- foreach(i=1:attr(.dist,"Size"), .combine=c) %dopar% {
						.own <- .distMat[i,-i]
						.others <- .distMat[-i,-i]; .others <- .others[lower.tri(.others)]
						wilcox.test(.own,.others,alternative="g")$p.value
					}
					names(.IScores0) <- labels(.dist)
					printLog("IQC Finished", .$.verbose)
					
					.$.IScores <- pmax(1-pnorm( qnorm(.IScores0, qnorm(0.95), 1), -qnorm(0.95), 1), 1e-20)
				}
				
				.CalcDistOfStudies <- function(.) {
					stopifnot(!is.null(.$.DListF))
					
					.$.DistOfStudies <- foreach(ii=iter(combinations(length(.$.DListF),2),by="row"), .combine=c, .packages="foreach") %dopar% {
						.gn <- intersect(rownames(.$.DListF[[ii[1]]]$x),rownames(.$.DListF[[ii[2]]]$x))
						
						.DistOfStudies <- foreach(1:100, .combine=c) %do% { #resampling based to fit practical memory limit
							..gn <- sample(.gn,length(.gn)*.1) 
							.CList <- foreach(jj=1:2, .combine=cbind) %do% { 
								as.dist(cor(t(.$.DListF[[ii[jj]]]$x[match(..gn,rownames(.$.DListF[[ii[jj]]]$x)),]), 
												method=.$.method.cor, use="pairwise.complete.obs"))
							}						
							.CList <- na.omit(.CList) #zero variance generate NA
							as.dist((1-cor(.CList, method="spearman"))/2)
						} 
						
						median(.DistOfStudies, na.rm=TRUE) 
					}
					
					attributes(.$.DistOfStudies) <- NULL
					attr(.$.DistOfStudies,"Labels") <- names(.$.DListF)
					attr(.$.DistOfStudies,"Size") <- length(.$.DListF)
					class(.$.DistOfStudies) <- "dist"
					attr(.$.DistOfStudies,"Diag") <- FALSE
					attr(.$.DistOfStudies,"Upper") <- FALSE
				}
				
				RunQC <- function(., nPath=NULL, B=1e4, pvalCut=.05, pvalAdjust=FALSE, fileForCQCp="c2.all.v3.0.symbols.gmt", isCAQC=FALSE) {
					if(!file.exists(fileForCQCp)) {
						res <- Download("MetaQC",fileForCQCp)
						if (inherits(res, "try-error") | res != 0L) {
							file.remove(fileForCQCp)
							stop(gettextf("download of file '%s' failed!\nPlease download gmt files at http://www.broadinstitute.org/gsea/downloads.jsp", fileForCQCp))
						} 
					}
					
					.GList <- paste(sub("(.+)[.][^.]+$", "\\1", basename(fileForCQCp)),".rda",sep="")
					if(!(.$.useCache & file.exists(.GList)))
						GMT2List(fileForCQCp, saveAs=.GList)
					
					.$EQC(nPathCut=nPath, .B=B)
					.$IQC()
					.$AQCg(.cutoff=pvalCut, .adjust=pvalAdjust)
					.$AQCp(.cutoff=pvalCut, .adjust=pvalAdjust, .GList=.GList)
					
					.$.CalcScores(isCAQC)
					
					return(.$.summary)
				} 
				
				.CalcScores <- function(., isCAQC=FALSE) {
					.Scores <- cbind.data.frame(IQC=-log10(.$.IScores), EQC=-log10(.$GetPvalOfScores(.$.Names[.$.excluded])),
							CQCg=-log10(.$.CQCgScores), CQCp=-log10(.$.CQCpScores), AQCg=-log10(.$.AQCgScores), AQCp=-log10(.$.AQCpScores))
					if(isCAQC) {
						.ScoresRankSum <- rank(-.Scores$IQC) + rank(-.Scores$EQC) +  rank(rank(-.Scores$CQCg) + rank(-.Scores$AQCg)) + rank(rank(-.Scores$CQCp) + rank(-.Scores$AQCp))
						.nScores <- 4
					} else {
						.ScoresRankSum <- rank(-.Scores$IQC) + rank(-.Scores$EQC) +  rank(-.Scores$CQCg) + rank(-.Scores$AQCg) + rank(-.Scores$CQCp) + rank(-.Scores$AQCp)
						.nScores <- 6
					}
					
					.$.Scores <- .Scores[order(.ScoresRankSum),]  
					.$.summary <- data.frame(Study=rownames(.$.Scores), round(.$.Scores,2), Rank=round(sort(.ScoresRankSum)/.nScores,2)) 
					
					.tmp <- cbind.data.frame(.$.summary[,1], 
							foreach(d=iter(.$.summary[,-c(1,ncol(.$.summary))],by="row"),.combine=rbind) %do% {ifelse(d < -log10(.05/length(.$.DListF)), paste(d,'*',sep=''), d)},
							.$.summary[,ncol(.$.summary)])
					colnames(.tmp) <- colnames(.$.summary); rownames(.tmp)=1:nrow(.tmp)
					.$.summary <- .tmp
				}
				
				#when need to change gene filter other than default
				SetupGeneFilter <- function(., cutRatioByMean=NULL, cutRatioByVar=NULL, minNumGenes=NULL, maxNApctAllowed=NULL) {
					warning("This function was deprecated.")
#					if(!is.null(cutRatioByMean)) .$.cutRatioByMean <- cutRatioByMean
#					if(!is.null(cutRatioByVar)) .$.cutRatioByVar <- cutRatioByVar
#					if(!is.null(minNumGenes)) .$.minNumGenes <- minNumGenes
#					if(!is.null(maxNApctAllowed)) .$.maxNApctAllowed <- maxNApctAllowed
#					.$.FilterGenes(.$.DListF0)
#					.$.Scores <- NULL
				}
				
				Plot <- function(., .scale.coord.var=4, isCAQC=FALSE) {
					if(is.null(.$.Scores))
						.$RunQC()
					.dat <- apply(.$.Scores, 2, function(s) {
								scale(s)
							})
					
					.dummy <- apply(.$.Scores, 2, function(s) {
								(-log10(.05/length(s)) - mean(s)) / sd(s)
							})
					
					if(isCAQC) {
						.dat <- cbind(.dat[,-match(c('CQCg','AQCg'),colnames(.dat))], CAQCg=rowMeans(.dat[,match(c('CQCg','AQCg'),colnames(.dat))]))
						.dummy <- c(.dummy[-match(c('CQCg','AQCg'),names(.dummy))], CAQCg=mean(.dummy[match(c('CQCg','AQCg'),names(.dummy))]))
						.dat <- cbind(.dat[,-match(c('CQCp','AQCp'),colnames(.dat))], CAQCp=rowMeans(.dat[,match(c('CQCp','AQCp'),colnames(.dat))]))
						.dummy <- c(.dummy[-match(c('CQCp','AQCp'),names(.dummy))], CAQCp=mean(.dummy[match(c('CQCp','AQCp'),names(.dummy))]))
					}
					
					.res <- prcomp(.dat, center=FALSE)
					
					.coord.dummy <- .dummy %*% .res$rotation[,1:2]
					.coord <- .res$x[,1:2]
					.coord <- sweep(.coord, 2, .coord.dummy)
					.coord.var <- sweep(.res$rotation, 2, .res$sdev, "*")[,1:2] #idea from FactoMineR
					
					#force plots be represented to right-upper side
					.sign <- sign(colSums(sign(.coord.var)))
					.sign <- ifelse(.sign>=0,1,-1)
					.coord <- sweep(.coord, 2, .sign, '*') 
					.coord.var <- sweep(.coord.var, 2, .sign, '*')
					
					.pctEig <- (.res$sdev^2/sum(.res$sdev^2)*100)[1:2]
					
					plot(x=.coord[,1],y=.coord[,2],type="n",xlab=bquote(bold(.(sprintf("1st Principal Component (%2.2f%%)",.pctEig[1])))),
							ylab=bquote(bold(.(sprintf("2nd Principal Component (%2.2f%%)",.pctEig[2])))),
							xlim=range(.coord)+c(-1,1)*diff(range(.coord))/4,
							ylim=range(.coord)+c(-1,1)*diff(range(.coord))/4,
							axes=FALSE)
					axis(side=1,lwd=4,tck=-0.02)
					axis(side=2,lwd=4,tck=-0.02)
					box(bty="L",lwd=4)
					
					abline(v=0, lty=2, lwd=2)
					abline(h=0, lty=2, lwd=2)
					
					for (v in 1:nrow(.coord.var)) {
						arrows(0, 0, .coord.var[v, 1]*.scale.coord.var, .coord.var[v, 2]*.scale.coord.var, 
								lwd=3, length = 0.1, angle = 15, code = 2, col=gray(.4)) 
						text(.coord.var[v, 1]*.scale.coord.var, y = .coord.var[v, 2]*.scale.coord.var, 
								labels = bquote(bold(.(rownames(.coord.var)[v]))), pos = 3, cex=1)
					}
					
					points(x=.coord[,1],y=.coord[,2],pch=1,col="black",cex=2.5,lwd=2)
					text(x=.coord[,1],y=.coord[,2],cex=1) 
				}
				
				Print <- function(.) {
					cat("Number of Studies: ", length(.$.DListF), fill=TRUE)
					cat("", fill=TRUE)
					cat("Dimension of Each Study:", fill=TRUE)
					.studies <- sapply(.$.DListF,function(d) dim(d$x)); rownames(.studies) <- c("Genes", "Samples")
					print(.studies)
					cat("", fill=TRUE)
					if(!is.null(.$.summary)) {
						cat("Quality Control Result:", fill=TRUE)
						print(.$.summary)
					}
				}
				
				GetPvalOfScores <- function(., .excludedN=NULL) {
					if(is.null(.$.PvalOfScores))
						.$EQC()
					
					if(is.null(.excludedN) || length(.excludedN)==0)
						.$.PvalOfScores
					else
						.$.PvalOfScores[-.$.GetExcluded(.excludedN)]
				}
				
				GetDistOfStudies <- function(., .excludedN=NULL) {
					if(is.null(.$.DistOfStudies))
						.$.CalcDistOfStudies() 
					if(is.null(.excludedN) || length(.excludedN)==0) {
						.$.DistOfStudies
					} else {
						.excluded <- .$.GetExcluded(.excludedN)
						as.dist(as.matrix(.$.DistOfStudies)[-sort(.excluded), -sort(.excluded)])
					}
				}
			})
	
	if(is.list(GList)) { #GList should be a list of gene sets
		stopifnot(all(length(names(GList))>0) & all(!duplicated(names(GList))))   #must be a pathway name & unique
		stopifnot(all(sapply(GList, is.character)))	#all genes should be a character vector	
	} else {
		if(!file.exists(GList)) {
			res <- Download("MetaQC",GList)
			if (inherits(res, "try-error") | res != 0L) {
				file.remove(GList)	
				stop(gettextf("download of file '%s' failed!\nPlease download gmt files at http://www.broadinstitute.org/gsea/downloads.jsp", GList))
			} 
		}
		
		if(getFileExt(GList)=="gmt") {
			.GList <- paste(getFileName(GList),".rda",sep="")
			if(useCache & file.exists(.GList))
				load(.GList) #loaded as GList
			else
				GList <- GMT2List(GList, saveAs=.GList)
		}
	}
	
	resp.type <- match.arg(resp.type)
	
	.p$.Initialize(.DList=DList, .GList=GList, .filterGenes=filterGenes, .resp.type=resp.type)
	
	return(.p)
}

plot.proto <- function(x, ...) {
	x$Plot(...)
}

print.proto <- function(x, ...) {
	x$Print(...)
}

