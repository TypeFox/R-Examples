#' Cladistic Data for Dicranograptid Graptolites from Song and Zhang (2014)
#'
#' Character matrix and majority-rule cladogram for 12 dicranograptid
#' (and outgroup) graptoloids, taken from Song and Zhang (2014). Included
#' here for use with functions related to character change.

#' @name SongZhangDicrano
#' @rdname SongZhangDicrano
#' @aliases charMatDicrano cladogramDicrano

#' @details
#' This example dataset is composed of a small cladistic character data for 12 taxa and 24 characters,
#' taken from Song and Zhang (2014). Note that character 22 is a biostratigraphic character, which was
#' not included in all analyses by Song and Zhang.
#'
#' The included cladogram is the majority-rule consensus of a maximum-parsimony analysis on 12 taxa with
#' 24 characters, including a biostratigraphic character. This tree is included (among the four depicted)
#' as it appeared to be the basis for the majority of Song and Zhang's discussion of dicranograptid
#' systematics.
#'
#' Both the matrix and the tree were entered by hand from their flat graphic depiction in Song and Zhang's
#' manuscript.

#' @format 
#' Loading this dataset adds two objects to the R environment.
#' \code{charMatDicrano} is a \code{data.frame} object composed of multiple factors, with \code{NA} values
#' representing missing values (states coded as '?'), read in with \code{readNexus} from package
#' \code{phylobase}. \code{cladogramDicrano} is a cladogram, formatted as a \code{phylo} class object
#' for use with package \code{ape}, without branch-lengths (as this was a consensus tree from a
#' maximum-parsimony analysis).

#' @source 
#' Song, Y., and Y. Zhang. 2014. A preliminary study on the relationship of the
#' early dicranograptids based on cladistic analysis. \emph{GFF} 136(1):243-248.

#' @keywords datasets

#' @docType data

#' @examples
#'
#' \dontrun{
#'
#'  
#'
#' require(ape)
#' require(phylobase)
#' 
#' charMatDicrano<-readNexus(file.choose(),type="data",SYMBOLS = " 0 1 2")
#' 
#' cladogramDicrano<-read.nexus(file.choose())
#' 
#' save(charMatDicrano,cladogramDicrano,file="SongZhangDicrano.rdata")
#' 
#' }
#'
#' data(SongZhangDicrano)
#'
#' # calculate a distance matrix from the character data
#' char<-charMatDicrano
#' charDist<-matrix(,nrow(char),nrow(char))
#' rownames(charDist)<-colnames(charDist)<-rownames(char)
#' for(i in 1:nrow(char)){for(j in 1:nrow(char)){
#' 	charDiff<-logical()
#' 	for(k in 1:ncol(char)){
#' 		selectPair<-char[c(i,j),k]
#' 		if(all(!is.na(selectPair))){
#' 			#drop states that are missing			
#' 			isSame<-identical(selectPair[1],selectPair[2])
#' 			charDiff<-c(charDiff,isSame)
#' 			}
#' 		}
#' 	charDist[i,j]<-1-sum(charDiff)/length(charDiff)
#' 	}}
#' 
#' #####
#' # the following is mostly stolen from the example code for my graptDisparity dataset
#' 
#' #now, is the diagonal zero? (it should be)
#' all(diag(charDist)==0)
#' 
#' #now, is the matrix symmetric? (it should be)
#' isSymmetric(charDist)
#' 
#' #can apply cluster analysis
#' clustRes <- hclust(as.dist(charDist))
#' plot(clustRes)
#' 
#' #can apply PCO (use lingoes correction to account for negative values
#'    #resulting from non-euclidean matrix
#' pco_res <- pcoa(charDist,correction="lingoes")
#' 
#' #relative corrected eigenvalues
#' rel_corr_eig <- pco_res$values$Rel_corr_eig
#' layout(1:2)
#' plot(rel_corr_eig)
#' #cumulative
#' plot(cumsum(rel_corr_eig))
#' 
#' #well let's look at those PCO axes anyway
#' layout(1)
#' pco_axes <- pco_res$vectors
#' plot(pco_axes[,1],pco_axes[,2],pch=16,
#'    xlab=paste("PCO Axis 1, Rel. Corr. Eigenvalue =",round(rel_corr_eig[1],3)),
#'    ylab=paste("PCO Axis 2, Rel. Corr. Eigenvalue =",round(rel_corr_eig[2],3)))
#' 
#' #######
#'
#' # plot majority rule tree from Song and Zhang
#' plot(cladogramDicrano,
#'	main="MajRule_24charX12Taxa_wBiostratChar")
#'
#'
NULL