#' PNN
#' 
#' Probabilistic neural network.
#' 
#' The package \pkg{PNN} implements the algorithm proposed by Specht (1990). It is written in the statistical langage R. It solves a common problem in automatic learning. Knowing a set of observations described by a vector of quantitative variables, we classify them in a given number of groups. Then, the algorithm is trained with this datasets and should guess afterwards the group of any new observation.  This neural network has the main advantage to begin generalization instantaneously even with a small set of known observations.
#' 
#' The package \pkg{PNN} exports four functions. These funtions are documented with examples and provided with unit tests:
#' \itemize{
#'  \item \code{\link{learn}}: Create a new Probabilist neural network with a new training set or update an existing one with new known observations.
#'  \item \code{\link{smooth}}: Set the smoothing parameter. If the value is not known, a genetic algorithm search the best value.
#'  \item \code{\link{perf}}: Compute the performance of the Probabilist neural network.
#'  \item \code{\link{guess}}: Guess the category of a new observation.
#' }
#' 
#' To help the use of \pkg{PNN}, the package contains a dataset \code{\link{norms}}. You could find more documentation at the package website: \url{http://flow.chasset.net/pnn/}.
#' 
#' The Probabilist neural network ist the main object used by the four functions. It is a \code{list} with several description fields:
#' \itemize{
#'  \item \code{model}: A name of the model ("Probabilistic neural network" by default).
#'  \item \code{set}: The raw training set.
#'  \item \code{category.column}: See above.
#'  \item \code{categories}: The categories found in the \code{category.column} field.
#'  \item \code{k}: The number of variables.
#'  \item \code{n}: The number of observations.
#'  \item \code{sigma}: The smoothing parameter.
#'  \item \code{observed}: A \code{list} of observed categories.
#'  \item \code{guessed}: A \code{list} of guessed categories.
#'  \item \code{success}: The number of times that the neural network chooses the right category.
#'  \item \code{fails}: The number of times that the neural network fails to guess the right category.
#'  \item \code{success_rate}: The rate of sucess over all observations in training set.
#'  \item \code{bic}: It is an adapted version of the Bayesian Information Criterion helping to compare different version of Probabilist neural networks.
#' }
#' 
#' @author Pierre-Olivier Chasset
#' @docType package
#' @seealso \code{\link{learn}}, \code{\link{smooth}}, \code{\link{perf}}, \code{\link{guess}}, \code{\link{norms}}
#' @keywords Neural network, Probability
#' @references Specht D.F. (1990). Probabilistic neural networks. Neural networks, 3(1):109-118.
#' @examples
#' library(pnn)
#' data(norms)
#' 
#' # The long way
#' pnn <- learn(norms)
#' pnn <- smooth(pnn, sigma=0.9)
#' pnn$sigma
#' \dontrun{pnn <- perf(pnn) # Optional}
#' \dontrun{pnn$success_rate # Optional}
#' guess(pnn, c(1,1))
#' guess(pnn, c(2,1))
#' guess(pnn, c(1.5,1))
#' 
#' # The short way
#' guess(smooth(learn(norms), sigma=0.8), c(1,1))
#' guess(smooth(learn(norms), sigma=0.8), c(2,1))
#' guess(smooth(learn(norms), sigma=0.8), c(1.5,1))
#' 
#' # Demonstrations
#' \dontrun{demo("norms-trainingset", "pnn")}
#' \dontrun{demo("small-trainingset", "pnn")}
#' @aliases pnn
#' @name pnn-package
NULL
