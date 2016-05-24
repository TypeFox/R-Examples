## RcppShark -- An R interface to the Shark machine learning library
##
## Copyright (C) 2015  Aydin Demircioglu
##
## This file is part of the RcppShark library for GNU R.
## It is made available under the terms of the GNU General Public
## License, version 3, or at your option, any later version,
## incorporated herein by reference.
##
## This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public
## License along with this program; if not, write to the Free
## Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
## MA 02111-1307, USA


#' Budgeted SGD Train.
#'
#' This will train a support vector machine by an SGD (Pegasos-like).
#' It will limit the number of support vector machines by applying the given
#' budget maintenance strategy.
#' See more on http://image.diku.dk/shark/sphinx_pages/build/html/rest_sources/tutorials/algorithms/kernelBudgetedSGD.html
#' 
#' @param x		matrix with input data
#' @param y		vector with labels
#' @param verbose		verbose output?
#' @param budget		size of budget 
#' @param strategy		strategy to use to maintain the budget size. 
#' 		choises are 'Merge', 'RemoveSmallest', 'RemoveRandom, 'Project'
#' @param C		regularization constant
#' @param gamma		kernel bandwidth for RBF kernel
#' @param epochs		number of iterations through data set 
#'
#' @note		Currently works only for binary classification. Uses only RBF kernel.
#'
#' @export
SharkBSGDTrain <- function(x, y = NULL, 
                 verbose = FALSE, budget = 500,
                 strategy = "Merge",
                 C = 1, gamma=1,  epochs = 1)
{
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertVector(y)
	checkmate::checkNumber(C, lower = 0)
	checkmate::checkNumber(gamma, lower = 0)
	checkmate::assertCount(budget)
	checkmate::assertCount(epochs)
	checkmate::assertString (strategy)
	checkmate::assertFlag (verbose)

	# train
	model <- .Call("RcppShark_BSGDWrapperTrain",
				X=x, y=y,
				C=C,
				budget = budget,
				gamma=gamma,
				epochs = epochs,
				budgetMaintenanceStrategy = "Merge",
				verbose = verbose,
				PACKAGE="RcppShark")
	
	# add for prediction
	model$gamma = gamma
	class(model) <- c("RcppShark.BSGD")
	return (model)
}



	
#' Budgeted SGD Predict.
#'
#' This will do prediction with a support vector machine trained by BSGD.
#' See more on http://image.diku.dk/shark/sphinx_pages/build/html/rest_sources/tutorials/algorithms/kernelBudgetedSGD.html
#'
#' @param x		matrix with input data
#' @param model		a model trained with BSGD Train
#' @param verbose		verbose output?
#'
#' @note		Currently works only for binary classification. Uses only RBF kernel.
#'
#' @examples
#'		x = as.matrix(iris[,1:4])
#'		y = as.vector(as.numeric(iris[,5]))
#'		y = replace(y, y == 2, 0)
#'		y = replace(y, y == 3, 0)
#'		model = SharkBSGDTrain (x, y, C = 0.0001, 
#'		      budget = 5, gamma = 1, epochs = 1, strategy = "Merge")
#'		results = SharkBSGDPredict (x, model)
#'		cat ("BSGD training error is ", sum(abs(y - results$predictions))/length(y), "\n")
#' @export
SharkBSGDPredict <- function(x, model, verbose = FALSE)
{
	# checkmate checks
	checkmate::assertClass (model, "RcppShark.BSGD")
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertFlag (verbose)	
	
	# predict
	res <- .Call("RcppShark_BSGDWrapperPredict",
				X=x, alpha = model$alpha,
				SV = model$SV,
				offset = model$offset,
				gamma = model$gamma,
				verbose = verbose,
				PACKAGE="RcppShark")
				
	return (res)
}

