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



#' Training a simple deep network.
#' 
#' This will train a simple 'deep' neural network with two hidden layers.
#' It will use an autoencoder for pretraining. 
#' For more information refer to the Shark tutorial at http://image.diku.dk/shark/sphinx_pages/build/html/rest_sources/tutorials/algorithms/deep_denoising_autoencoder_network.html
#'
#' @param x		matrix with input data
#' @param y		vector with labels
#' @param nHidden1		number of nodes of first hidden layer (part of network model)
#' @param nHidden2		number of nodes of second hidden layer (part of network model)
#' @param unsupRegularisation 	regularization factor of supervised training
#' @param noiseStrength		noise strength for unsupervised training
#' @param unsupIterations			iteration number for unsupervised training
#' @param regularisation		regularisation factor for supervised training
#' @param iterations		iteration number for supervised training
#' @param verbose		print extra information?
#'
#' @examples
#'		x = as.matrix(iris[,1:4])
#'		y = as.vector(as.numeric(iris[,5]))
#'		y = replace(y, y == 2, 0)
#'		y = replace(y, y == 3, 0)
#'		model = DeepNetworkTrain (x, y, nHidden1 = 32, nHidden2 = 32)
#'		results = DeepNetworkPredict (x, model)
#'		networkPrediction = apply (results$prediction, 1, which.max) - 1
#'		errors = sum(abs(y - networkPrediction))/length(y)
#'		cat("Network produced ", errors, "errors.\n")
#' @export
DeepNetworkTrain <- function(x, y, nHidden1 = 8L, nHidden2 = 8L, unsupRegularisation = 0.001, noiseStrength = 0.3, unsupIterations = 100L, regularisation = 0.0001, iterations = 200L, verbose = FALSE) 
{
	checkmate::checkMatrix (x)
	checkmate::checkVector (y)
	checkmate::checkCount (nHidden1, positive = TRUE)
	checkmate::checkCount (nHidden2, positive = TRUE)
	checkmate::checkNumber (unsupRegularisation, lower = 0)
	checkmate::checkNumber (noiseStrength, lower = 0)
	checkmate::checkCount (unsupIterations, positive = TRUE)
	checkmate::checkNumber (regularisation, lower = 0)
	checkmate::checkCount (iterations, positive = TRUE)
	checkmate::checkFlag (verbose)
	
	model = .Call('RcppShark_DeepNetworkWrapperTrain', PACKAGE = 'RcppShark', 
		x, y, nHidden1, nHidden2, unsupRegularisation, noiseStrength, unsupIterations, regularisation, iterations, verbose)
	
	class (model) <- c("RcppShark.DeepNetwork")
	return (model)
}



#' Predictions from a simple deep network.
#'
#' This will do prediction using a trained simple 'deep' neural network with two hidden layers.
#' For more information refer to the Shark tutorial at http://image.diku.dk/shark/sphinx_pages/build/html/rest_sources/tutorials/algorithms/deep_denoising_autoencoder_network.html
#'
#' @param x		matrix with input data
#' @param model		a model trained with the deep network trainer.
#' @param verbose		verbose output?
#'
#' @export
DeepNetworkPredict <- function (x, model, verbose = FALSE) {
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertClass (model, "RcppShark.DeepNetwork")
	checkmate::assertFlag (verbose)	

	# extract information
	weights = model$weights
	nHidden1 = model$nHidden1
	nHidden2 = model$nHidden2
	inputSize = model$inputSize 
	outputSize = model$outputSize

	res = .Call('RcppShark_DeepNetworkWrapperPredict', PACKAGE = 'RcppShark', x, weights, nHidden1, nHidden2, inputSize, outputSize, verbose)
	return (res)
}

