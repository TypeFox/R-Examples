#' reconstruct nuclear families
#'
#' Function for reconstruction of nulecar families 
#' or extended pedigrees based on results of 
#' \code{\link{findRelatives}} output (estimatd 
#' meiotic distance matrix). Reconstruction is 
#' based on the fact that parent-offspring pairs 
#' have meiotic distance of '1', and sibs have a 
#' distance '2+2'. If both parents and the offspring 
#' are genotyped, expected distace between offspring 
#' and both parents is '1', and distance between two 
#' parents is >2 (coded as 'NA'). 
#' 
#' @param relationshipGuessMatrix meiotic relationship 
#' matrix, as estimated by \code{\link{findRelatives}}
#' 
#' @param sex Sex, coded with 1 for males and 0 for females
#' 
#' @return A matrix containing reconstructed pedigree(s) 
#' coded in linkage-like format. If "fid" is zero, this 
#' means that a pedigree could not be assigned. 
#' 
#' @author Yurii Aulchenko
#' 
#' @examples 
#' nloci <- 100
#' q <- runif(nloci,min=0.05,max=0.95)
#' # g7---g8
#' #   _|_
#' #  |   |
#' # g9  g10---g11
#' #       __|_
#' #      |   /\
#' #    g12 g13 g14
#' #
#' nids <- 8
#' sex <- c(1,0,0,1,0,0,0,0)
#' names(sex) <- paste("g",c(7:14),sep="")
#' gt <- matrix(ncol=nloci,nrow=nids)
#' rownames(gt) <- paste("g",c(7:14),sep="")
#' gt["g7",] <- rbinom(nloci,2,q)
#' gt["g8",] <- rbinom(nloci,2,q)
#' gt["g11",] <- rbinom(nloci,2,q)
#' gt["g9",] <- generateOffspring(gt["g7",],gt["g8",],q=q)
#' gt["g10",] <- generateOffspring(gt["g7",],gt["g8",],q=q)
#' gt["g12",] <- generateOffspring(gt["g10",],gt["g11",],q=q)
#' gt["g13",] <- generateOffspring(gt["g10",],gt["g11",],q=q)
#' gt["g14",] <- gt["g13",]
#' aa<-findRelatives(gt,q=q,nmei=c(1:2))
#' aa$guess
#' aaPed <- reconstructNPs(aa$guess,sex)
#' aaPed
#' 
#' 
reconstructNPs <- function(relationshipGuessMatrix,sex) 
{
# sanity checks; extracting info in convenient format
	if ( is.null(rownames(relationshipGuessMatrix)) || is.null(colnames(relationshipGuessMatrix)) ) 
		stop("relationshipGuessMatrix should have row and column names")
	if ( any(rownames(relationshipGuessMatrix) != colnames(relationshipGuessMatrix)) )
		stop ("rownames != colnames in relationshipGuessMatrix")
	if (!is(relationshipGuessMatrix,"matrix"))
		stop("relationshipGuessMatrix must be a matrix")
	if (!is(relationshipGuessMatrix[1,1],"character"))
		stop("relationshipGuessMatrix must contain 'character'")
	if (dim(relationshipGuessMatrix)[1] != dim(relationshipGuessMatrix)[2]) 
		stop("ncol != nrows in relationshipGuessMatrix")
	if ( any(relationshipGuessMatrix[upper.tri(relationshipGuessMatrix)] != 
					t(relationshipGuessMatrix)[upper.tri(relationshipGuessMatrix)], na.rm = TRUE) ) {
		warning("lower triangle != upper triangle in relationshipGuessMatrix; using lower.tri")
		
	}
	relationshipGuessMatrix[upper.tri(relationshipGuessMatrix)] <- NA
	save_RelGM <- relationshipGuessMatrix
	save_RelGM[upper.tri(save_RelGM)] <- t(save_RelGM)[upper.tri(save_RelGM)]
	relationshipGuessMatrix <- t(relationshipGuessMatrix)
# also change diag 0 -> NA
	diag(relationshipGuessMatrix) <- NA 
	nids <- dim(relationshipGuessMatrix)[1]
	idnames <- rownames(relationshipGuessMatrix)
	if (length(sex) != nids)
		stop("length(sex) != nids")
	if (!is.null(names(sex))) {
		if (!all( names(sex) %in% idnames)) stop("names(sex) do nor match names in relationshipGuessMatrix")
		sex <- sex[idnames]
	} else {
		warning("no names(sex); assuming the same order as in naming of relationshipGuessMatrix")
		names(sex) <- idnames
	}
# start inferences on Nuclear Pedigrees (NPs)
	outColNames <- c("fid","id","father","mother","sex",
			"MZT","SIB","PO","UNRELATED","PARSEXPRO")
	out <- matrix(nrow=nids,ncol=length(outColNames))
	colnames(out) <- outColNames
	rownames(out) <- idnames
	out[,"id"] <- idnames
	sex[which(sex==0)] <- "female"
	sex[which(sex==1)] <- "male"
	out[,"sex"] <- sex
	out[,"fid"] <- out[,"father"] <- out[,"mother"] <- out[,"MZT"] <- 0 
	out[,"SIB"] <- out[,"UNRELATED"] <- out[,"PARSEXPRO"] <- 0
	out[,"PO"] <- "0"
	#print(relationshipGuessMatrix)
	if ( all(is.na(relationshipGuessMatrix)) ) {
		warning("all non-diag elements of relationshipGuessMatrix are NA; returning 'unrelated'")
		return(out)
	}
# identify specific pairs
	identifySpecificPair <- function(matrixNotion,name,out,notYet="0") {
		for (cid in rownames(out))
			if (any(relationshipGuessMatrix[cid,] == matrixNotion, na.rm = TRUE)) {
				twset <- c(cid,idnames[which(relationshipGuessMatrix[cid,] == matrixNotion)])
				if (all(out[twset,name] == notYet)) {
					out[twset,name] <- paste(twset,collapse="-")
				} else if (all(out[twset,name] != notYet)) {
				} else {
					stop(paste("confused with",name))
				}
			}
		out
	}
# identify twins
	out <- identifySpecificPair("0","MZT",out)
# identify sibs
	out <- identifySpecificPair("2+2","SIB",out)
# identify parent-offspring
	for (cid in idnames) {
		if (any(save_RelGM[cid,] == "1",na.rm=TRUE)) {
			poids <- c(cid,idnames[which(save_RelGM[cid,] == "1")])
			out[cid,"PO"] <- paste(poids,collapse="-")
		} 
	}
# identify unrelated
	curNucFam <- 1
	for (cid in idnames) {
		if (all(is.na(save_RelGM[cid,]) | save_RelGM[cid,] == "0" )) {
			out[cid,"fid"] <- curNucFam
			out[cid,"UNRELATED"] <- "1"
			curNucFam <- curNucFam + 1
		} 
	}
# identify offspring and construct NPs
	marIs.na <- function(matr,MAR) {
		out <- apply(matr,MARGIN=MAR,FUN=function(x) {any(is.na(x))})
		out
	}
	colIs.na <- function(matr) {return(marIs.na(matr,MAR=2))}
	rowIs.na <- function(matr) {return(marIs.na(matr,MAR=1))}
	for (cid in idnames) 
		if (out[cid,"fid"] == "0") {
			ones <- idnames[which(save_RelGM[cid,] == "1")]
			if (length(ones) >= 2) {
				miniRel <- save_RelGM[c(ones),c(ones)]
				diag(miniRel) <- miniRel[upper.tri(miniRel)] <- "0"
				if (sum(is.na(miniRel)) == 1) {
					p1 <- colnames(miniRel)[which(colIs.na(miniRel))]
					p2 <- rownames(miniRel)[which(rowIs.na(miniRel))]
					if (is.na(sex[p1]) || is.na(sex[p2])) {sexprob <- 1}
					else if (sex[p1]==sex[p2]) {sexprob <- 1}
					else {
						if (sex[p1] == "male" && sex[p2] == "female") {
							Fa <- p1; Mo <- p2;
							sexprob <- 0
						} else if (sex[p2] == "male" && sex[p1] == "female") {
							Fa <- p1; Mo <- p2;
							sexprob <- 0
						} else {
							sexprob <- 1
						}
					}
					offspring <- cid
					if (out[cid,"MZT"] != "0") 
						offspring <- c(offspring,out[which(out[,"MZT"] == out[cid,"MZT"]),"id"])
					if (out[cid,"SIB"] != "0") 
						offspring <- c(offspring,out[which(out[,"SIB"] == out[cid,"SIB"]),"id"])
					out[offspring,"father"] <- Fa
					out[offspring,"mother"] <- Mo
					nfam <- c(Fa,Mo,offspring) 
					if (all(out[nfam,"fid"]=="0")) {
						out[nfam,"fid"] <- curNucFam
						curNucFam <- curNucFam + 1
					} else {
						uniqueFamIds <- unique(out[nfam,"fid"])
						uniqueFamIds <- uniqueFamIds[which(uniqueFamIds != "0")]
						minFamId <- min(as.integer(uniqueFamIds))
						for (jj in uniqueFamIds) {out[which(out[,"fid"] == jj),"fid"] <- minFamId;}
						out[nfam,"fid"] <- minFamId;
					}
					if (sexprob) out[nfam,"PARSEXPRO"] <- "1"
				} 
			}
		} 
	
	return(out)
}
