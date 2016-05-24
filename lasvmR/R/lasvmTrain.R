#
# lasvmRrain.R
#
# Copyright (C) 2015  Aydin Demircioglu, aydin.demircioglu /at/ ini.rub.de
#
# This file is part of the lasvmR library for GNU R.
# It is made available under the terms of the GNU General Public
# License, version 2, or at your option, any later version,
# incorporated herein by reference.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
# MA 02111-1307, USA
#


#' lasvmTrain
#' 
#' Use lasvm to train a given problem.
#'
#'  @param	x		data matrix 
#'  @param y		labels
#'  @param  gamma		RBF kernel parameter
#'  @param  cost		regularization parameter
#'  @param  degree 	degree for poly kernel
#'  @param  coef0 	coefficient for poly kernel
#'  @param  optimizer 	type of optimizer
#'  @param  kernel 	kernel type
#'  @param  selection 	selection strategy
#'  @param  termination 	criterion for stopping
#'  @param  sample		time for stopping/number of iterations tec
#'  @param  cachesize 	size of kernel cache
#'  @param  bias 	use  bias?
#'  @param  epochs 	number of epochs
#'  @param epsilon 	stopping criterion parameter
#'  @param	verbose		verbose output?
#'
#'  @return	a list consisting of
#'	alpha		alpha for SVs as vector
#'	SV		support vectors as matrix 
#'
#' @examples
#' model = lasvmR::lasvmTrain (x = as.matrix(iris[seq(1,150,2),1:4]),
#' 	y = (as.numeric(iris[seq(1,150,2),5]) %% 2)*2-1,
#' 	gamma = 1, 
#' 	cost = 1, 
#' 	kernel = 2)
#' ytrue = (as.numeric(iris[seq(2,150,2),5]) %% 2)*2-1
#' result = lasvmPredict (x = as.matrix(iris[seq(2,150,2),1:4]), model)
#' ypred = result$predictions
#' error = sum(abs(ypred - ytrue))/length(ytrue)
#' cat ("Error rate =", error*100)
#' @export
lasvmTrain = function (x, y, 
	gamma = 1,
	cost = 1,
	 degree = 3,
	 coef0 = 0,
	 optimizer = 1,
	 kernel = 2,
	 selection = 0,
	 termination = 0,
	 sample = 100000000, # from source code
	 cachesize = 256,
	 bias = 1,
	 epochs = 1,
	epsilon = 0.001,
	verbose = FALSE)
{
	# check arguments TODO: addl imits
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertVector(y)
	checkmate::checkNumber(cost, lower = 0)
	checkmate::checkNumber(gamma, lower = 0)
	checkmate::checkNumber(degree, lower = 0)
	checkmate::checkNumber(coef0)
	checkmate::checkNumber(sample)
	checkmate::assertCount(optimizer)
	checkmate::assertInt(kernel, lower = 0, upper = 3)
	checkmate::assertInt(selection, lower = 0, upper = 2)
	checkmate::assertCount(termination)
	checkmate::assertCount(cachesize)
	checkmate::assertCount(bias)
	checkmate::assertCount(epochs)
	checkmate::checkNumber(epsilon)
	checkmate::assertFlag (verbose)

	model = lasvmTrainWrapper (x, y, 
		gamma = gamma, 
		cost = cost,
		degree = 	 degree,
		coef0 = 	 coef0,
		optimizer = 	 optimizer,
		kernel = 	 kernel,
		selection = 	 selection,
		termination = 	 termination,
		sample = sample,
		cachesize = 	 cachesize,
		bias = 	 bias,
		epochs = 	 epochs,
		epsilon = 	epsilon,
		verbose = verbose)

	# move other needed data into model
	model$gamma = gamma
	model$degree = degree
	model$coef0 = coef0
	model$cost = cost
	model$kernel = kernel
		
	class(model) <- c("lasvmR.model")
	return (model)
}

