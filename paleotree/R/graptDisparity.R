#' Morphlogical Character and Range Data for late Ordovician and Early Silurian Graptoloidea
#' 
#' This dataset contains a morphological character matrix (containing both a set of
#' 45 discrete characters, and 4 continuous characters coded as minimum and maximum
#' range values), along with biostratigraphic range data for 183 graptoloid
#' species-level taxa from Bapst et al. (2012, PNAS). Also includes a pre-calculated
#' distance matrix based on the character matrix, using the algorithm applied by Bapst
#' et al (2012). Interval dates for biostratigraphic zones is taken from Sadler et al. 2011.
#'

#' @details 
#' The character matrix contains characters of two differing types with a (very) small but
#' non-neglible amount of missing character data for some taxa. This required the
#' use of an unconventional ad hoc distance metric for the published analysis, resulting
#' in a (very slightly) non-Euclidean distance matrix. This breaks some assumptions of
#' some statistical analyses or requires special corrections, such as with PCO.
#'
#' Note that taxonomic data were collected only for species present within an interval
#' defined by the base of the Uncinatus biozone (~448.57 Ma) to the end of the cyphus
#' biozone (~439.37 Ma). Many taxa have first and last appearance dates listed in \code{graptRanges}
#' which are outside of this interval (see examples).

#' @name graptDisparity

#' @rdname graptDisparity

#' @aliases graptDisparity graptCharMatrix graptRanges graptDistMat

#' @docType data

#' @format This dataset is composed of three objects: 
#'
#' \describe{
#' \item{graptCharMatrix}{A matrix composed of mixed character data
#' and a group code for 183 graptoloid taxa, with rows named with species
#' names and columns named with character names.}
#'
#' \item{graptRanges}{A list containing two matrices: the first matrix
#' describes the first and last interval times for 20 graptolite biozones
#' and the second matrix contains the first and last appearances of 183
#' graptolite species in those same biozones. (In other words, graptRanges
#' has the 'timeList' format called by some paleotree functions).}
#'
#' \item{graptDistMat}{A 183x183 matrix of pair-wise distances (dissimilarities)
#' for the 183 graptolite species, using the algorithm for discrete characters
#' and min-max range values described in Bapst et al.}}
#'

#' @source 
#' Source for stratigraphic ranges and character data: 
#'
#' Bapst, D. W., P. C. Bullock, M. J. Melchin, H. D. Sheets, and C. E. Mitchell.
#' 2012. Graptoloid diversity and disparity became decoupled during the Ordovician
#' mass extinction. \emph{Proceedings of the National Academy of Sciences}
#' 109(9):3428-3433.
#' 
#' Source for interval dates for graptolite zones: 
#'
#' Sadler, P. M., R. A. Cooper,
#' and M. Melchin. 2009. High-resolution, early Paleozoic (Ordovician-Silurian)
#' time scales. \emph{Geological Society of America Bulletin} 121(5-6):887-906.

#' @seealso 
#' For more example graptolite datasets, see \code{\link{retiolitinae}}
#'
#' This data was added mainly to provide an example dataset
#' for \code{\link{nearestNeighborDist}}

#' @keywords datasets

#' @examples
#' 
#' #load data
#' data(graptDisparity)
#' 
#' #separate out two components of character matrix
#' 
#' #45 discrete characters
#' discChar <- graptCharMatrix[,1:45]
#' 
#' #min ranges for 4 continuous characters
#' cMinChar <- graptCharMatrix[,c(46,48,50,52)]
#' #max ranges for 4 continuous characters
#' cMaxChar <- graptCharMatrix[,c(47,49,51,53)]
#' 
#' #group (clade/paraclade) coding
#' groupID <- graptCharMatrix[,54]
#' 
#' #number of species
#' nspec <- nrow(graptCharMatrix)
#' 
#' #some plotting information from Bapst et al.'s plotting scripts
#' grpLabel <- c("Normalo.","Monogr.","Climaco.",
#' 		"Dicrano.","Lasiogr.","Diplogr.","Retiol.")
#' grpColor <- c("red","purple",colors()[257],colors()[614],
#' 		colors()[124],"blue",colors()[556])
#' 
#' ##########
#'
#' #plot diversity curve of taxa
#' taxicDivDisc(graptRanges)
#'
#' #but the actual study interval for the data is much smaller
#' abline(v=448.57,lwd=3) #start of study interval
#' abline(v=439.37,lwd=3) #end of study interval
#'
#' #plot diversity curve just for study interval
#' taxicDivDisc(graptRanges, timelims=c(448.57,439.37))
#'
#' ############
#'
#' #distance matrix is given as graptDistMat
#'    #to calculate yourself, see code below in DoNotRun section
#' 
#' #now, is the diagonal zero? (it should be)
#' all(diag(graptDistMat)==0)
#' 
#' #now, is the matrix symmetric? (it should be)
#' isSymmetric(graptDistMat)
#' 
#' #can apply cluster analysis
#' clustRes <- hclust(as.dist(graptDistMat))
#' plot(clustRes,labels=FALSE)
#' 
#' #use ape to plot with colors at the tips
#' dev.new(width=15) 	# for a prettier plot
#' plot.phylo(as.phylo(clustRes),show.tip.label=FALSE,
#' 		no.margin=TRUE,direction="upwards")
#' tiplabels(pch=16,col=grpColor[groupID+1])
#' legend("bottomright",legend=grpLabel,col=grpColor,pch=16)
#' dev.set(2)
#' 
#' #can apply PCO (use lingoes correction to account for negative values
#'    #resulting from non-euclidean matrix
#' pco_res <- pcoa(graptDistMat,correction="lingoes")
#' 
#' #relative corrected eigenvalues
#' rel_corr_eig <- pco_res$values$Rel_corr_eig
#' layout(1:2)
#' plot(rel_corr_eig)
#' #cumulative
#' plot(cumsum(rel_corr_eig))
#' 
#' #first few axes account for very little variance!!
#' 
#' 
#' 
#' #well let's look at those PCO axes anyway
#' layout(1)
#' pco_axes <- pco_res$vectors
#' plot(pco_axes[,1],pco_axes[,2],pch=16,col=grpColor[groupID+1],
#'    xlab=paste("PCO Axis 1, Rel. Corr. Eigenvalue =",round(rel_corr_eig[1],3)),
#'    ylab=paste("PCO Axis 2, Rel. Corr. Eigenvalue =",round(rel_corr_eig[2],3)))
#' legend("bottomright",legend=grpLabel,col=grpColor,pch=16,ncol=2,cex=0.8)
#' 
#' 
#' ##########m##############
#'
#' \dontrun{ 
#'
#' #calculate a distance matrix (very slow!)
#' #Bapst et al. calculated as # char diffs / total # of chars
#'    #but both calculated for only non-missing characters for both taxa
#' #non-identical discrete states = difference for discrete traits
#' #non-overlapping ranges for continuous characters = difference for cont traits
#' 
#' distMat<-matrix(,nspec,nspec)
#' rownames(distMat)<-colnames(distMat)<-rownames(graptCharMatrix)
#' for(i in 1:nspec){ for(j in 1:nspec){    #calculate for each pair of species
#'    #discrete characters
#'    di <- discChar[i,]     #discrete character vector for species i
#'    dj <- discChar[j,]     #discrete character vector for species j
#'    #now calculate pair-wise differences for non-missing characters
#'    discDiff <- (di!=dj)[!is.na(di)&!is.na(dj)] #logical vector
#'    #
#'    #continuous characters: need another for() loop
#'    contDiff <- numeric()
#'    for(ct in 1:4){
#'        #if they do not overlap, a min must be greater than a max value
#'        contDiff[ct] <- cMinChar[i,ct]>cMaxChar[j,ct] | cMinChar[j,ct]>cMaxChar[i,ct]
#'        }
#'    #remove NAs
#'    contDiff <- contDiff[!is.na(contDiff)]
#'    #combine
#'    totalDiff <- c(discDiff,contDiff)
#'    #divide total difference 
#'    distMat[i,j] <- sum(totalDiff)/length(totalDiff)
#'    }}
#' 
#' #but is it identical to the distance matrix already provided?
#' identical(distMat,graptDistMat)
#' #ehh, numerical rounding issues...
#' 
#' #A somewhat speeder alternative to calculate a distance matrix
#' distMat<-matrix(,nspec,nspec)
#' rownames(distMat)<-colnames(distMat)<-rownames(graptCharMatrix)
#' for(i in 1:(nspec-1)){ for(j in (i+1):nspec){    #calculate for each pair of species
#'    #now calculate pair-wise differences for non-missing characters
#'    discDiff <- (discChar[i,]!=discChar[j,])[
#'       !is.na(discChar[i,])&!is.na(discChar[j,])] #logical vector
#'    #continuous characters: if they do not overlap, a min must be greater than a max value
#'       contDiff<-sapply(1:4,function(ct)
#'             cMinChar[i,ct]>cMaxChar[j,ct] | cMinChar[j,ct]>cMaxChar[i,ct])
#'    #remove NAs, combine, divide total difference 
#'    distMat[i,j] <- distMat[j,i] <- sum(c(discDiff,contDiff[!is.na(contDiff)]))/length(
#'         c(discDiff,contDiff[!is.na(contDiff)]))
#'    }}
#' diag(distMat)<-0
#' 
#' #but is it identical to the distance matrix already provided?
#' identical(distMat,graptDistMat)
#' #ehh, MORE numerical rounding issues...
#' 
#' }
#' 
NULL