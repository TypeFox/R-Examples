#' Binary classification with Rotation Forest (Rodriguez en Kuncheva, 2006)
#'
#' \code{rotationForest} implements an ensemble method where each base classifier (tree) is fit on the principal components of the variables of random partitions of the feature set.
#' 
#' @param x A data frame of predictors (numeric, or integer). Categorical variables need to be transformed to indicator (dummy) variables.
#' @param y A factor containing the response vector. Only \{0,1\} is allowed.
#' @param K The number of variable subsets. The default is the value \code{K} that results in three features per subset.
#' @param L The number of base classfiers (trees using the \code{rpart} package). The default is 10.
#' @param ... Arguments to \code{rpart.control} of the \code{rpart} package. Firt run \code{library(rpart)}.
#' @examples
#' data(iris)
#' y <- as.factor(ifelse(iris$Species[1:100]=="setosa",0,1))
#' x <- iris[1:100,-5]
#' rF <- rotationForest(x,y)
#' predict(object=rF,newdata=x)
#' @references Rodriguez, J.J., Kuncheva, L.I., 2006. Rotation forest: A new classifier ensemble method. IEEE Trans. Pattern Anal. Mach. Intell. 28, 1619-1630. doi:10.1109/TPAMI.2006.211
#' @seealso \code{\link{predict.rotationForest}}
#' @return An object of class \code{rotationForest}, which is a list with the following elements:
#'   \item{models}{A list of trees.}
#'   \item{loadings}{A list of loadings.}
#'   \item{columnnames}{Column names of x.}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @keywords classification
rotationForest <- function(x,y,K=round(ncol(x)/3,0),L=10,...){

    x <- data.frame(sapply(x,as.numeric))
		#reset K as to the remainder of the number of features and K should is 0
    #%% moduls
		while (ncol(x) %% K != 0) {
									K <- K - 1
									}

		#calculate number of features M=n/K per subset
		M <- round(ncol(x)/K)


	predicted <- list()
	fit <- numeric()
	Ri <- list()
	Ria <- list()
	fit <- list()
  predicted <- matrix(NA,nrow=nrow(x),ncol=L)
			
	subsets <- list()

	SelectedClass<- list()
	IndependentsClassSubset<- list()
	IndependentsClassSubsetBoot<- list()
	pcdata <- list()
	loadings <- list()



#For i=1..L=number of classifiers
for (i in 1:L) {

	#############prepare the rotation matrix RAi
		#Split F(the feature set) randomly into K subsets: Fi,j (for j= 1...K)
      
			#randomize order of variables
			Independents <- x[,sample(1:ncol(x),ncol(x))]
			
			n <- 0

			subsets[[i]] <- list()

			SelectedClass[[i]]<- list()
			IndependentsClassSubset[[i]]<- list()
			IndependentsClassSubsetBoot[[i]]<- list()
			pcdata[[i]] <- list()
			loadings[[i]] <- list()

		#For j=1...K (this is n in subsets[[n]]
			for (j in seq(1,K)) {
				
				
				

				#let Xi,j be the dataset X for the feature set Fi,j
					
					n <- n + M

					subsets[[i]][[j]] <- data.frame(Independents[,(n-(M-1)):n],y)

									
				#Eliminate from Xi,j a random subset of classes

          SelectedClass[[i]][[j]]<- as.integer(sample(levels(as.factor(y)),1))
					IndependentsClassSubset[[i]][[j]] <- subsets[[i]][[j]][subsets[[i]][[j]]$y == SelectedClass[[i]][[j]],] 


				#Select a bootstrap sample from Xi,j of size 75% of the number of objects in Xi,j. Denote the new set by XijBoot
					IndependentsClassSubsetBoot[[i]][[j]] <- IndependentsClassSubset[[i]][[j]][sample(1:dim(IndependentsClassSubset[[i]][[j]])[1], round(0.75*nrow(IndependentsClassSubset[[i]][[j]])), replace = TRUE),]
	
				#Apply PCA on XijBoot to obtain the coefficients in a matrix Ci,j (here called loadings)
					pcdata[[i]][[j]] <- princomp(IndependentsClassSubsetBoot[[i]][[j]][,!colnames(IndependentsClassSubsetBoot[[i]][[j]]) %in% "y"])
          loadings[[i]][[j]] <- pcdata[[i]][[j]]$loadings[,]

					#rename columns with coefficients so they can be merged properly
          #colnames(loadings[[i]][[j]]) <- gsub("V", "a",dimnames(loadings[[i]][[j]])[[1]])
	        colnames(loadings[[i]][[j]]) <- dimnames(loadings[[i]][[j]])[[1]]
					#add unique ID and name it rowID (for merging)
					loadings[[i]][[j]] <- data.frame(dimnames(loadings[[i]][[j]])[[1]],loadings[[i]][[j]])
					colnames(loadings[[i]][[j]])[1] <- "rowID"
	
			}



		#Arrange the Ci,j, for j=1...K in a Rotation matrix Ri
			#pregrow data frame C (Coefficients = loadings)
			#because the number of columns and rows can differ we should do a loop across all loading subsets and sum the nrow's and ncol's
			#Ri <- data.frame(matrix(nrow=(n*dim(loadings[[j]])[[1]]), ncol=(2*dim(loadings[[j]])[[2]]-(j-1))) )
			


    Ri[[i]] <- Reduce(function(x, y) merge(x, y, by='rowID',all=TRUE),loadings[[i]])

		#Construct Ria by rearranging the columns of Ri so as to match the order of features in F (Feature set)
			#replace NA's with 0
			Ri[[i]][is.na(Ri[[i]])] <- 0


			#at this point we have to sort everything (rows and columns) for subsequent multiplication
			
			#sort rows and columns
      Ria[[i]] <- Ri[[i]][order(match(Ri[[i]]$rowID,colnames(x))),order(match(colnames(Ri[[i]]),colnames(x)))]

			rownames(Ria[[i]]) <- Ria[[i]]$rowID
      Ria[[i]]$rowID <- NULL


      final <- data.frame(as.matrix(x) %*% as.matrix(Ria[[i]]),y)




	#############build classifier Di using (X Ria, y) as the training set

		fit[[i]] <- rpart(y~ ., method="class", data=final,...)


}
res <- list(models=fit, loadings=Ria, columnnames=colnames(x))
class(res) <- "rotationForest"
res
}