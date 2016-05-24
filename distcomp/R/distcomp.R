#' Distributed Computing with R
#'
#' \code{distcomp} is a collection of methods to fit models to data that may be
#' distributed at various sites. The package arose as a way of addressing the
#' issues regarding data aggregation; by allowing sites to have control over
#' local data and transmitting only summaries, some privacy controls can be
#' maintained. Even when participants have no objections in principle to data
#' aggregation, it may still be useful to keep data local and expose just the
#' computations. For further details, please see the arxiv paper cited below.
#'
#' The initial implementation consists of a stratified Cox model fit with
#' distributed survival data and a Singular Value Decomposition
#' of a distributed matrix. General Linear Models will soon be added.
#' Although some sanity checks and balances are present, many more are needed
#' to make this truly robust. We also hope that other methods will be added by users.
#'
#' We make the following assumptions in the implementation:
#' (a) the aggregate data is logically a stacking of data at each site, i.e.,
#' the full data is row-partitioned into sites where the rows are observations;
#' (b) Each site has the package \code{distcomp} installed and a workspace setup
#' for (writeable) use by the \code{opencpu} server
#' (see \code{\link{distcompSetup}}); and (c) each site is exposing \code{distcomp}
#' via an \code{opencpu} server.
#'
#' The main computation happens via a master process, a script of R code,
#' that makes calls to \code{distcomp} functions at worker sites via \code{opencpu}.
#' The use of \code{opencpu} allows developers to prototype their distributed implementations
#' on a local machine using the \code{opencpu} package that runs such a server locally
#' using \code{localhost} ports.
#'
#' Note that \code{distcomp} computations are not intended for speed/efficiency;
#' indeed, they are orders of magnitude slower. However, the models that are fit are
#' not meant to be recomputed often. These and other details are discussed in the
#' paper mentioned above.
#'
#' The current implementation, particularly the Stratified Cox Model, makes direct use of
#' code from \code{\link[survival]{coxph}}. That is, the underlying Cox model code is
#' derived from that in the R \code{survival} survival package.
#'
#' For an understanding of how this package is meant to be used, please see the documented
#' examples and the reference.
#' @seealso The examples in \code{system.file("doc", "examples.html", package="distcomp")}
#' @seealso The source for the examples: \code{system.file("doc_src", "examples.Rmd", package="distcomp")}.
#' @docType package
#' @references Software for Distributed Computation on Medical Databases:
#' A Demonstration Project \url{http://arxiv.org/abs/1412.6890}
#' @references Appendix E of Modeling Survival Data: Extending the Cox Model by
#' Terry M. Therneau and Patricia Grambsch. Springer Verlag, 2000.
#' @name distcomp
NULL

#' Make an appropriate opencpu URL for a specified function and url prefix for the
#' distcomp package
#'
#' @description .makeOpencpuURL returns an appropriate URL to call a function in the distcomp
#' package given the name of the function and a url prefix.
#'
#' @param fn is the name of the function in the distcomp package
#' @param urlPrefix is the URL of the opencpu server with the distcomp package installed
#' @return the formatted url as a string
#'
#' @rdname distcomp-internal
#' @import utils
#' @importFrom stringr str_trim
#' @importFrom httr POST
#' @importFrom httr headers
#' @importFrom httr add_headers
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite toJSON
#'
#' @examples
#' distcomp:::.makeOpencpuURL("foo", "http://localhost:9999/ocpu")
#'
#' @keywords internal
.makeOpencpuURL <- function(fn, urlPrefix, package="distcomp") {
  paste(urlPrefix, "library", package, "R", fn, "json", sep="/")
}

#' Check that a definition object meets minimal requirements
#'
#' @description .defnOK returns TRUE or FALSE depending on whether the definition object
#' meets minimimal requirements.
#'
#' @param defn is the definition object passed
#' @return TRUE or FALSE depending on the result
#'
#' @rdname distcomp-internal
#'
#' @examples
#' distcomp:::.defnOK(data.frame()) ## FALSE
#' distcomp:::.defnOK(data.frame(id = "ABC", stringsAsFactors=FALSE)) ## TRUE
#'
#' @keywords internal
.defnOK <- function(defn) {
    ! is.null(defn$id) && nchar(defn$id) > 0
}

#' Setup a workspace and configuration for a distributed computation
#' @description The function \code{discompsetup} sets up a distributed computation
#' and configures some global parameters such as definition file names,
#' data file names, instance object file names, and ssl configuration parameters. The
#' function creates some of necessary subdirectories if not already present and throws
#' an error if the workspace areas are not writeable
#' @seealso \code{getConfig}
#'
#' @param workspacePath a folder specifying the workspace path. This has to be
#' writable by the opencpu process. On a cloud opencpu server on Ubuntu, for example,
#' this requires a one-time modification of apparmor profiles to enable write
#' permissions to this path
#' @param defnPath the path where definition files will reside, organized by
#' computation identifiers
#' @param instancePath the path where instance objects will reside
#' @param defnFileName the name for the compdef definition files
#' @param dataFileName the name for the data files
#' @param instanceFileName the name for the instance files
#' @param ssl_verifyhost integer value, usually \code{1L}, but for testing with
#' snake-oil certs, one might set this to \code{0L}
#' @param ssl_verifypeer integer value, usually \code{1L}, but for testing with
#' snake-oil certs, one might set this to \code{0L}
#' @return TRUE if all is well
#'
#' @importFrom httr config
#'
#' @examples
#' \dontrun{
#' distcompSetup(workspacePath="./workspace")
#' }
#' @export
distcompSetup <- function(workspacePath = "",
                          defnPath = paste(workspacePath, "defn", sep=.Platform$file.sep),
                          instancePath = paste(workspacePath, "instances", sep=.Platform$file.sep),
                          defnFileName = "defn.rds",
                          dataFileName = "data.rds",
                          instanceFileName = "instance.rds",
                          ssl_verifyhost = 1L,
                          ssl_verifypeer = 1L) {
  ## TODO: In the next version, this should be stuffed in an R6 class
  testFileName <- paste(sample(letters, 15), collapse="")
  if (!file.exists(defnPath)) {
    defnOk <- dir.create(defnPath)
    if (!defnOk) {
      stop("distcompSetup: workspace permissions issue: not writable!")
    }
  } else {
    ## defnPath exists; check if it is writable
    testFile <- paste(defnPath, testFileName, sep=.Platform$file.sep)
    createOk <- file.create(testFile)
    if (!createOk) {
      stop("distcompSetup: workspace permissions issue: not writable!")
    }
    file.remove(testFile)
  }

  if (!file.exists(instancePath)) {
    defnOk <- dir.create(instancePath)
    if (!defnOk) {
      stop("distcompSetup: workspace permissions issue: not writable!")
    }
  } else {
    ## instancePath exists; check if it is writable
    testFile <- paste(instancePath, testFileName, sep=.Platform$file.sep)
    createOk <- file.create(testFile)
    if (!createOk) {
      stop("distcompSetup: workspace permissions issue: not writable!")
    }
    file.remove(testFile)
  }

  distcompEnv <- getOption("distcompEnv")
  distcompEnv[["config"]] <- list(workspacePath = workspacePath,
                                  defnPath = defnPath,
                                  instancePath = instancePath,
                                  defnFileName = defnFileName,
                                  dataFileName = dataFileName,
                                  instanceFileName = instanceFileName,
                                  sslConfig = config(ssl_verifyhost = ssl_verifyhost,
                                    ssl_verifypeer = ssl_verifypeer))
  distcompEnv[["computationInfo"]] <- list()
  TRUE
}


#' Return the workspace and configuration setup values
#' @description The function \code{getConfig} returns the values of the
#' configuration parameters set up by \code{distcompSetup}
#' @seealso \code{distcompSetup}
#' @param ... any further arguments
#' @return a list consisting of
#' \item{workspacePath}{a folder specifying the workspace path. This has to be
#' writable by the opencpu process. On a cloud opencpu server on Ubuntu, for example,
#' this requires a one-time modification of apparmor profiles to enable write
#' permissions to this path}
#' \item{defnPath}{the path where definition files will reside, organized by
#' computation identifiers}
#' \item{instancePath}{the path where instance objects will reside}
#' \item{defnFileName}{the name for the compdef definition files}
#' \item{dataFileName}{the name for the data files}
#' \item{instanceFileName}{the name for the instance files}
#' \item{ssl_verifyhost}{integer value, usually \code{1L}, but for testing with
#' snake-oil certs, one might set this to \code{0L}}
#' \item{ssl_verifypeer}{integer value, usually \code{1L}, but for testing with
#' snake-oil certs, one might set this to \code{0L}}
#'
#' @examples
#' \dontrun{
#' getConfig()
#' }
#' @export
getConfig <- function(...) {
  ## TODO: In the next version, this should be stuffed in an R6 class
  getOption("distcompEnv")[["config"]]
}


#' Make a worker object given a definition and data
#' @description The function \code{makeWorker} returns an object of the
#' appropriate type based on a computation definition and sets the data for
#' the object. The types of objects that can be created depend upon the
#' available computations
#' @seealso \code{\link{availableComputations}}
#' @param defn the computation definition
#' @param data the data for the computation
#' @return a worker object of the appropriate class based on the definition
#'
#' @export
makeWorker <- function (defn, data) {
  compType <- defn$compType
  available <- availableComputations()
  k <- match(compType, names(available))
  if (is.na(k)) {
    stop(sprintf("No such computation: %s", compType))
  } else {
    available[[k]]$makeWorker(defn, data)
  }
}

#' Make a master object given a definition
#' @description The function \code{makeMaster} returns a master object
#' corresponding to the definition. The types of master objects that can
#' be created depend upon the available computations
#' @seealso \code{\link{availableComputations}}
#' @param defn the computation definition
#' @return a master object of the appropriate class based on the definition
#'
#' @export
makeMaster <- function(defn) {
  compType <- defn$compType
  available <- availableComputations()
  k <- match(compType, names(available))
  if (is.na(k)) {
    stop(sprintf("No such computation: %s", compType))
  } else {
    available[[k]]$makeMaster(defn)
  }
}

#' Return the currently available (implemented) computations
#' @description The function \code{availableComputations} returns a list
#' of available computations with various components. The names of this list
#' (with no spaces) are unique canonical tags that are used throughout the
#' package to unambiguously refer to the type of computation; web applications
#' particularly rely on this list to instantiate objects. As more computations
#' are implemented, this list is augmented.
#' @seealso \code{\link{getComputationInfo}}
#' @return a list with the components corresponding to a computation
#' \item{desc}{a textual description (25 chars at most)}
#' \item{definitionApp}{the name of a function that will fire up a shiny webapp
#' for defining the particular computation}
#' \item{workerApp}{the name of a function that will fire up a shiny webapp
#' for setting up a worker site for the particular computation}
#' \item{masterApp}{the name of a function that will fire up a shiny webapp
#' for setting up a master for the particular computation}
#' \item{makeDefinition}{the name of a function that will return a data frame
#' with appropriate fields needed to define the particular computation assuming
#' that they are populated in a global variable. This function is used by web
#' applications to construct a definition object based on inputs specified
#' by the users. Since the full information is often gathered incrementally by
#' several web applications, the inputs are set in a global variable and
#' therefore retrieved here using the function \code{getComputationInfo}
#' designed for the purpose}
#' \item{makeMaster}{a function that will construct a master object for the
#' computation given the definition and a logical flag indicating
#' if debugging is desired}
#' \item{makeWorker}{a function that will construct
#' a worker object for that computation given the definition and data}
#'
#' @examples
#' availableComputations()
#' @export
availableComputations <- function() {
  list(
    StratifiedCoxModel = list(
      desc = "Stratified Cox Model",
      definitionApp = "defineNewCoxModel",
      setupWorkerApp = "setupCoxWorker",
      setupMasterApp = "setupCoxMaster",
      makeDefinition = function() {
        data.frame(id = getComputationInfo("id"),
                   compType = getComputationInfo("compType"),
                   projectName = getComputationInfo("projectName"),
                   projectDesc = getComputationInfo("projectDesc"),
                   formula = getComputationInfo("formula"),
                   stringsAsFactors=FALSE)
      },
      makeMaster = function(defn, debug = FALSE) CoxMaster$new(defnId = defn$id, formula = defn$formula, debug=debug),
      makeWorker = function(defn, data) CoxWorker$new(data = data, formula = defn$formula)
    ),
    RankKSVD = list(
      desc = "Rank K SVD",
      definitionApp = "defineNewSVDModel",
      setupWorkerApp = "setupSVDWorker",
      setupMasterApp = "setupSVDMaster",
      makeDefinition = function() {
        data.frame(id = getComputationInfo("id"),
                   compType = getComputationInfo("compType"),
                   projectName = getComputationInfo("projectName"),
                   projectDesc = getComputationInfo("projectDesc"),
                   rank = getComputationInfo("rank"),
                   ncol = getComputationInfo("ncol"),
                   stringsAsFactors=FALSE)
      },
      makeMaster = function(defn, debug = FALSE) SVDMaster$new(defnId = defn$id, k = defn$rank, debug = debug),
      makeWorker = function(defn, data) SVDWorker$new(x = data)
    )
  )
}

#' Return currently implemented data sources
#' @description The function \code{availableDataSources} returns the
#' currently implemented data sources such as CSV files, Redcap etc.
#'
#' @return a list of named arguments, each of which is another list, with
#' required fields named \code{desc}, a textual description and
#' \code{requiredPackages}
#'
#' @examples
#' availableDataSources()
#' @export
availableDataSources <- function() {
  list(
    CSVFile = list(desc = "CSV File", requiredPackages = list())
    , Redcap = list(desc = "Redcap API", requiredPackages = list("redcapAPI"))
    ##, Postgres = list(desc = "Postgres", requiredPackages = list("RPostgreSQL"))
  )
}


#' Make a computation definition given the computation type
#' @description The function \code{makeDefinition} returns a computational
#' definition based on current inputs (from the global store) given a
#' canonical computation type tag. This is a utility function for web
#' applications to use as input is being gathered
#'
#' @seealso \code{\link{availableComputations}}
#' @param compType the canonical computation type tag
#' @return a data frame corresponding to the computation type
#' @examples
#' \dontrun{
#' makeDefinition(names(availableComputations())[1])
#' }
#' @export
makeDefinition <- function(compType) {
  available <- availableComputations()
  k <- match(compType, names(available))
  if (is.na(k)) {
    stop(paste("makeDefintion: No such computation:", compType))
  } else {
    available[[k]]$makeDefinition()
  }
}

#' Given the id of a serialized object, invoke a method on the object with arguments
#' @description The function \code{executeMethod} is really the heart of distcomp.
#' It executes an arbitrary method on an object that has been serialized to the
#' distcomp workspace with any specified arguments. The result, which is dependent
#' on the computation that is executed, is returned. If the object needs to save
#' state between iterations on it, it is automatically serialized back for the ensuing
#' iterations
#'
#' @param objectId the (instance) identifier of the object on which to invoke a method
#' @param method the name of the method to invoke
#' @param ... further arguments as appropriate for the method
#' @return a result that depends on the computation being executed
#' @export
executeMethod <- function(objectId, method, ...) {
  config <- getConfig()
  filePath <- paste(config$instancePath, objectId, config$instanceFileName, sep=.Platform$file.sep)
  object <- readRDS(file=filePath)
  call <- substitute(object$METHOD(...), list(METHOD = as.name(method)))
  result <- eval(call)
  if (object$getStateful()) {
    saveRDS(object, file=filePath)
  }
  result
}

#' Deserialize the result of a http response
#' @description .deSerialize will convert the JSON result of a http response as needed,
#' else the raw content is returned.
#'
#' @rdname distcomp-internal
#' @import utils
#' @importFrom stringr str_trim
#' @importFrom httr POST
#' @importFrom httr headers
#' @importFrom httr add_headers
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite toJSON
#' @param q the result of a httr response
#' @return the converted result, if JSON, or the raw content
#'
#' @keywords internal
## Deserialize URL result, only handle JSON at the moment
.deSerialize <- function(q) {
  cType <- headers(q)['content-type']
  if (cType == "application/json")
    fromJSON(rawToChar(q$content))
  else
    q$content
}


#' Given the definition identifier of an object, instantiate and store object in workspace
#' @description The function \code{createInstanceObject} uses a definition identified by
#' defnId to create the appropriate object instance. The instantiated object is assigned
#' the instanceId and saved under the dataFileName if the latter is specified.
#' This instantiated object may change state between iterations when a computation executes
#' @seealso \code{\link{availableComputations}}
#' @param defnId the identifier of an already defined computation
#' @param instanceId an indentifier to use for the created instance
#' @param dataFileName a file name to use for saving the data. Typically \code{NULL}, this
#' is only needed when one is using a single opencpu server to behave like multiple
#' sites in which case the data file name serves to distinguish the site-specific data files.
#' When it is \code{NULL}, the data file name is taken from the configuration settings
#' @import utils
#' @return TRUE if everything goes well
#' @export
createInstanceObject <- function (defnId, instanceId, dataFileName=NULL) {
  config <- getConfig()
  defn <- readRDS(paste(config$defnPath, defnId, config$defnFileName, sep=.Platform$file.sep))
  compType <- defn$compType
  available <- availableComputations()
  if (!(compType %in% names(available))) {
    stop(paste("createInstanceObject: No such computation:", compType))
  }
  data <- readRDS(paste(config$defnPath, defnId,
                        if (is.null(dataFileName)) config$dataFileName else dataFileName,
                        sep=.Platform$file.sep))

  object <- makeWorker(defn, data)

  ## Check if the instance folder exists
  thisInstancePath <- paste(config$instancePath, instanceId, sep=.Platform$file.sep)
  dir.create(thisInstancePath)
  ## Save it under the instance id to find it.
  saveRDS(object, file=paste(thisInstancePath, config$instanceFileName, sep=.Platform$file.sep))
  TRUE
}

#' Destroy an instance object given its identifier
#' @description The function \code{destroyInstanceObject} deletes an object associated
#' with the instanceId. This is typically done after a computation completes and results
#' have been obtained.
#' @param instanceId the id of the object to destroy
#' @seealso \code{\link{createInstanceObject}}
#' @import utils
#' @return TRUE if everything goes well
#' @export
destroyInstanceObject <- function (instanceId) {
  config <- getConfig()
  file.remove(paste(config$instancePath, instanceId, config$instanceFileName, sep=.Platform$file.sep))
  file.remove(paste(config$instancePath, instanceId, sep=.Platform$file.sep))
  TRUE
}

#' Save a computation instance, given the computation definition, associated data and
#' possibly a data file name to use
#' @description The function \code{saveNewComputation} uses the computation definition to save
#' a new computation instance. This is typically done for every site that wants to participate
#' in a computation with its own local data. The function examines the computation definition
#' and uses the identifier therein to uniquely refer to the computation instance at the site.
#' This function is invoked (maybe remotely) on the opencpu server by
#' \code{\link{uploadNewComputation}} when a worker site is being set up
#' @seealso \code{\link{uploadNewComputation}}
#' @param defn the identifier of an already defined computation
#' @param data the (local) data to use
#' @param dataFileName a file name to use for saving the data. Typically \code{NULL}, this
#' is only needed when one is using a single opencpu server to behave like multiple
#' sites in which case the data file name serves to distinguish the site-specific data files.
#' When it is \code{NULL}, the data file name is taken from the configuration settings
#' @return TRUE if everything goes well
#' @export
saveNewComputation <- function(defn, data, dataFileName=NULL) {
  config <- getConfig()
  defnId <- defn$id
  thisDefnPath <- paste(config$defnPath, defnId, sep=.Platform$file.sep)
  ## Should tryCatch this
  if (!file.exists(thisDefnPath)) {
    dir.create(thisDefnPath)
  }
  saveRDS(object=defn, file=paste(thisDefnPath, config$defnFileName, sep=.Platform$file.sep))
  if (is.null(dataFileName)) {
    saveRDS(object=data, file=paste(thisDefnPath, config$dataFileName, sep=.Platform$file.sep))
  } else {
    saveRDS(object=data, file=paste(thisDefnPath, dataFileName, sep=.Platform$file.sep))
  }
  TRUE
}

#' Upload a new computation and data to an opencpu server
#' @description The function \code{uploadNewComputation} is really a remote version
#' of \code{\link{saveNewComputation}}, invoking that function on an opencpu server.
#' This is typically done for every site that wants to participate in a computation
#' with its own local data. Note that a site is always a list of at least a unique
#' name element (distinguishing the site from others) and a url element.
#' @seealso \code{\link{saveNewComputation}}
#' @param site a list of two items, a unique \code{name} and a \code{url}
#' @param defn the identifier of an already defined computation
#' @param data the (local) data to use
#' @importFrom stringr str_trim
#' @importFrom httr POST
#' @importFrom httr headers
#' @importFrom httr add_headers
#' @importFrom jsonlite fromJSON
#' @importFrom jsonlite toJSON
#'
#' @return TRUE if everything goes well
#' @export
uploadNewComputation <- function(site, defn, data) {
  if (! .defnOK(defn)) {
      stop("uploadNewComputation: Improper definition")
  }
  localhost <- (grepl("^http://localhost", site$url) || grepl("^http://127.0.0.1", site$url))
  payload <- if (localhost) {
    list(defn = defn, data = data, dataFileName = paste0(site$name, ".rds"))
  } else {
    list(defn = defn, data = data)
  }
  q <- POST(.makeOpencpuURL(urlPrefix = site$url, fn = "saveNewComputation"),
            body = toJSON(payload),
            add_headers("Content-Type" = "application/json"),
            config=getConfig()$sslConfig
            )
  .deSerialize(q)
}


#' Generate an identifier for an object
#' @description A hash is generated based on the contents of the object
#'
#' @seealso \code{\link[digest]{digest}}
#' @param object the object for which a hash is desired
#' @param algo the algorithm to use, default is "xxhash64" from
#' \code{\link[digest]{digest}}
#' @importFrom digest digest
#' @return the hash as a string
#'
#' @export
generateId <- function(object, algo="xxhash64") digest::digest(object, algo=algo)

##
## We have many apps that need to communicate values between each other.
## So a global area is needed to share variables.
## The function saveModelInfo will save model information in a named list.
## The function getModelInfo will return the value of an item in that list.
##


#' Set a name to a value in a global variable
#' @description In distcomp, several web applications need to communicate
#' between themselves. Since only one application is expected to be
#' active at any time, they do so via a global store, essentially a hash table.
#' This function sets a name to a value
#'
#' @seealso \code{\link{getComputationInfo}}
#' @param name the name for the object
#' @param value the value for the object
#' @return invisibly returns the all the name value pairs
#'
#' @export
setComputationInfo <- function(name, value) {
  ## TODO: In the next version, this should be stuffed in an R6 class
  distcompEnv <- getOption("distcompEnv")
  currentValue <- distcompEnv[["computationInfo"]]
  ## This is guaranteed to be non-null because of distcompSetup
  ## Set the value
  currentValue[[name]] <- value
  ## Then save it back
  distcompEnv[["computationInfo"]] <- currentValue
  invisible(currentValue)
}

#' Get the value of a variable from the global store
#' @description In distcomp, several web applications need to communicate
#' between themselves. Since only one application is expected to be
#' active at any time, they do so via a global store, essentially a hash table.
#' This function retrieves the value of a name
#'
#' @seealso \code{\link{setComputationInfo}}
#' @param name the name for the object
#' @return the value for the variable, \code{NULL} if not set
#'
#' @export
getComputationInfo <- function(name) {
  ## TODO: In the next version, this should be stuffed in an R6 class
  distcompEnv <- getOption("distcompEnv")
  ## The statement below is guaranteed to work because of distcompSetup
  currentValue <- distcompEnv[["computationInfo"]]
  ## Note return value could be null
  currentValue[[name]]
}

#' Clear the contents of the global store
#' @description In distcomp, several web applications need to communicate
#' between themselves. Since only one application is expected to be
#' active at any time, they do so via a global store, essentially a hash table.
#' This function clears the store, except for the working directory.
#'
#' @seealso \code{\link{setComputationInfo}} \code{\link{getComputationInfo}}
#' @return an empty list
#'
#' @export
resetComputationInfo <- function() {
  ## TODO: In the next version, this should be stuffed in an R6 class
  workingDir <- getComputationInfo("workingDir")
  distcompEnv <- getOption("distcompEnv")
  if (is.null(workingDir)) {
    invisible(distcompEnv[["computationInfo"]] <- list())
  } else {
    invisible(distcompEnv[["computationInfo"]] <- list(workingDir = workingDir))
  }
}

#' Run a specified distcomp web application
#' @description Web applications can define computation, setup worker sites or masters.
#' This function invokes the appropriate web application depending on the task
#' @seealso \code{\link{defineNewComputation}}, \code{\link{setupWorker}},
#' \code{\link{setupMaster}}
#' @import shiny
#' @param appType one of three values: "definition", "setupWorker", "setupMaster"
#' @return the results of running the web application
#'
#' @export
runDistcompApp <- function(appType = c("definition", "setupWorker", "setupMaster")) {
  resetComputationInfo()
  appType <- match.arg(appType)
  app <- paste0(appType, "App")
  appPath <- system.file("webApps", app, package="distcomp")
  compType <- shiny::runApp(appPath, launch.browser=TRUE)
  available <- availableComputations()
  availableNames <- names(available)
  index <- match(compType, availableNames)
  subApp <- available[[index]][[app]]
  appPath <- system.file("webApps", app, subApp, package="distcomp")
  ##print(appPath)
  shiny::runApp(appPath, launch.browser=TRUE)
}

#' Define a new computation
#' @description This function just calls \code{\link{runDistcompApp}} with the
#' parameter "definition"
#' @seealso \code{\link{runDistcompApp}}
#' @return the results of running the web application
#'
#' @export
defineNewComputation <- function() {
  setComputationInfo("workingDir", getwd())
  runDistcompApp(appType = "definition")
}

#' Setup a worker site
#' @description This function just calls \code{\link{runDistcompApp}} with the
#' parameter "setupWorker"
#' @seealso \code{\link{runDistcompApp}}
#' @return the results of running the web application
#'
#' @export
setupWorker <- function() {
  setComputationInfo("workingDir", getwd())
  runDistcompApp(appType = "setupWorker")
}

#' Setup a computation master
#' @description This function just calls \code{\link{runDistcompApp}} with the
#' parameter "setupMaster"
#' @seealso \code{\link{runDistcompApp}}
#' @return the results of running the web application
#'
#' @export
setupMaster <- function() {
  resetComputationInfo()
  setComputationInfo("workingDir", getwd())
  appPath <- system.file("webApps", "setupMasterApp", package="distcomp")
  shiny::runApp(appPath, launch.browser=TRUE)
}

#' Write the code necessary to run a master process
#' @description Once a computation is defined, worker sites are set up, the master process
#' code is written by this function. The current implementation does not allow one to mix
#' localhost URLs with non-localhost URLs
#' @seealso \code{\link{setupMaster}}
#' @param defn the computation definition
#' @param sites a named list of site URLs participating in the computation
#' @param outputFileName the name of the output file to which code will be written
#' @return the value \code{TRUE} if all goes well
#' @export
writeCode <- function(defn, sites, outputFileName) {
  compType <- defn$compType
  siteNames <- names(sites)
  wd <- getComputationInfo("workingDir")
  f <- file(paste(wd, outputFileName, sep=.Platform$file.sep), open="w")
  writeLines("library(distcomp)", con=f)
  dump(c("defn", "sites"), f)
  writeLines("master <- makeMaster(defn)", f)
  writeLines("for (site in sites) {", f)
  writeLines("   master$addSite(name = site$name, url = site$url)", f)
  writeLines("}", f)
  writeLines("result <- master$run()", f)
  writeLines("print(master$summary())", f)
  close(f)
  TRUE
}
