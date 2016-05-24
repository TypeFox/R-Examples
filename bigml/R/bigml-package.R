#' R bindings for BigML API
#'
#' \tabular{ll}{
#' Package: \tab bigml\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1-1\cr
#' Date: \tab 20012-04-30\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' A set of methods that enable straightforward usage of the BigML API.
#' The methods use R idioms and native datatypes where appropriate, while
#' also providing access to more conventional API usage.
#'
#'
#' @name bigml-package
#' @aliases bigml
#' @docType package
#' @template author
#' @keywords package
#' @examples
#' \dontrun{
#' 	# set default credentials
#' 	setCredentials('username', 'key')
#' 	model = quickModel(iris, 'Species')
#' 	quickPrediction(model, c(Petal.Width=0.2, Petal.Length=1.4))
#'
#' 	# use specific credentials
#' 	quickPrediction(model, c(Petal.Width=0.2, Petal.Length=1.4),
#'     username='someuser', api_key='somekey')
#'
#'	# list most recent sources
#'  listSources()
#'
#' 	# specify limit and offset
#'  listModels(limit=15,offset=300)
#'
#'	# specify filter criteria
#'  listDatasets(size__gt=1048576)
#' }
NULL