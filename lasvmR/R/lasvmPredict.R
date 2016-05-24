#
# lasvmPredict.R
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


#' lasvmPredict
#' 
#' Use lasvm to train a given problem.
#'
#'  @param	x		data matrix 
#'  @param model		trained model
#'  @param	verbose		verbose output?
#'
#'  @return	a list consisting of
#'	predictions		the predicted labels
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
lasvmPredict = function (x, model, verbose = FALSE)
{
	# check arguments
	checkmate::assertMatrix(x, min.rows = 1)
	checkmate::assertClass (model, "lasvmR.model")
	checkmate::assertFlag (verbose)
	
	results = lasvmPredictWrapper (x, model$SV, model$alpha, 
		gamma = model$gamma,
		kdegree = model$degree,
		kcoef0 = model$coef0,
		bias = 	 model$bias,
		kerneltype = model$kernel,
		verbose = verbose)

	return (results);
}

