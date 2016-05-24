#' Class providing object with methods for communication with lightning-viz server
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom RCurl postForm
#' @importFrom RJSONIO fromJSON toJSON
#' @importFrom httr POST
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods for communication with lightning-viz server.
#' @format \code{\link{R6Class}} object.
#' @examples
#' Lightning$new("http://localhost:3000/")
#' Lightning$new("http://your-lightning.herokuapp.com/")
#' @field serveraddress Stores address of your lightning server.
#' @field sessionid Stores id of your current session on the server.
#' @field url Stores url of the last visualization created by this object.
#' @field autoopen Checks if the server is automatically opening the visualizations.
#' @field notebook Checks if the server is in the jupyter notebook mode.
#' #' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/lightning-viz/lightining-r/}
#'   \item{\code{new(serveraddress)}}{This method is used to create object of this class with \code{serveraddress} as address of the server object is connecting to.}
#'
#'   \item{\code{sethost(serveraddress)}}{This method changes server that you are contacting with to \code{serveraddress}.}
#'   \item{\code{createsession(sessionname = "")}}{This method creates new session on the server with optionally given name in \code{sessionname}.}
#'   \item{\code{usesession(sessionid)}}{This method changes currently used session on the server to the one with id given in \code{sessionid} parameter.}
#'   \item{\code{openviz(vizid = NA)}}{This method by default opens most recently created by this object visualization. If \code{vizid} parameter is given, it opens a visualization with given id instead.}
#'   \item{\code{enableautoopening()}}{This method enables auto opening of every visualisation that you create since that moment. Disabled by default.}
#'   \item{\code{disableautoopening()}}{This method disables auto opening of every visualisation that you create since that moment. Disabled by default.}
#'   \item{\code{line(series, index = NA, color = NA, label = NA, size = NA, xaxis = NA, yaxis = NA, logScaleX = "false", logScaleY = "false")}}{This method creates a line visualization for vector/matrix with each row representing a line, given in \code{series}.}
#'   \item{\code{scatter(x, y, color = NA, label = NA, size = NA, alpha = NA, xaxis = NA, yaxis = NA)}}{This method creates a scatterplot for points with coordinates given in vectors \code{x, y}.}
#'   \item{\code{linestacked(series, color = NA, label = NA, size = NA)}}{This method creates a plot of multiple lines given in matrix \code{series}, with an ability to hide and show every one of them.}
#'   \item{\code{force(matrix, color = NA, label = NA, size = NA)}}{This method creates a force plot for matrix given in \code{matrix}.}
#'   \item{\code{graph(x, y, matrix, color = NA, label = NA, size = NA)}}{This method creates a graph of points with coordinates given in \code{x, y} vectors, with connection given in \code{matrix} connectivity matrix.}
#'   \item{\code{map(regions, weights, colormap)}}{This method creates a world (or USA) map, marking regions given as a vector of abbreviations (3-char for countries, 2-char for states) in \code{regions} with weights given in \code{weights} vector and with \code{colormap} color (string from colorbrewer).}
#'   \item{\code{graphbundled(x, y, matrix, color = NA, label = NA, size = NA)}}{This method creates a bundled graph of points with coordinates given in \code{x, y} vectors, with connection given in \code{matrix} connectivity matrix. Lines on this graph are stacked a bit more than in the \code{graph} function.}
#'   \item{\code{matrix(matrix, colormap)}}{This method creates a visualization of matrix given in \code{matrix} parameter, with its contents used as weights for the colormap given in \code{colormap} (string from colorbrewer).}
#'   \item{\code{adjacency(matrix, label = NA)}}{This method creates a visualization for adjacency matrix given in \code{matrix} parameter.}
#'   \item{\code{scatterline(x, y, t, color = NA, label = NA, size = NA)}}{This method creates a scatterplot for coordinates in vectors \code{x, y} and assignes a line plot to every point on that plot. Each line is given as a row in \code{t} matrix.}
#'   \item{\code{scatter3(x, y, z, color = NA, label = NA, size = NA, alpha = NA)}}{This method creates a 3D scatterplot for coordinates given in vectors \code{x, y, z}.}
#'   \item{\code{image(imgpath)}}{This method uploads image from file \code{imgpath} to the server and creates a visualisation of it.}
#'   \item{\code{gallery(imgpathvector)}}{This method uploads images from vector of file paths \code{imgpathvector} to the server and creates a gallery of these images.}}

Lightning <- R6Class("Lightning",
   public = list(
      serveraddress = NA,
      sessionid = NA,
      notebook = NA,
      url = NA,
      autoopen = FALSE,
      initialize = function(serveraddress, notebook = F) {
         if(!missing(serveraddress)){
            self$serveraddress <- serveraddress
         }

         self$notebook = notebook
      },
      line = function(series, index = NA, color = NA, label = NA, size = NA, xaxis = NA, yaxis = NA, logScaleX = "false", logScaleY = "false") {
         listbuilder <- list(type = "line", options = list(logScaleX = logScaleX, logScaleY = logScaleY))
         features <- list(series = series)

         if (!is.na(index)) {
            features$index <- index
         }
         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }
         if (!is.na(xaxis)) {
            features$xaxis <- xaxis
         }
         if (!is.na(yaxis)) {
            features$yaxis <- yaxis
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }



         jsonbody <- toJSON(listbuilder, .na = "")
         oldbody = "placeholder"

         while (!(jsonbody == oldbody)) {
            oldbody <- jsonbody
            jsonbody <- gsub(", ,", ",", jsonbody)
         }
         jsonbody <- gsub(",  ]", " ]", jsonbody)
         jsonbody <- gsub('[ ,', '[', jsonbody, fixed = TRUE)
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")

         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      scatter = function(x, y, color = NA, label = NA, size = NA, alpha = NA, xaxis = NA, yaxis = NA){
         listbuilder <- list(type = "scatter", options = list())
         points <- matrix(ncol = 2, nrow = length(x))
         points[,1] <- x
         points[,2] <- y
         features <- list(points = points)

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }
         if (!is.na(alpha)) {
            features$alpha <- alpha
         }
         if (!is.na(xaxis)) {
            features$xaxis <- xaxis
         }
         if (!is.na(yaxis)) {
            features$yaxis <- yaxis
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         if (self$autoopen) {
            browseURL(url)
         }
         self$url <- url
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      linestacked = function(series, color = NA, label = NA, size = NA){
         listbuilder <- list(type = "line-stacked", options = list())
         features <- list(series = series)

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "")

         oldbody = "placeholder"

         while (!(jsonbody == oldbody)) {
            oldbody <- jsonbody
            jsonbody <- gsub(", ,", ",", jsonbody)
         }
         jsonbody <- gsub(",  ]", " ]", jsonbody)
         jsonbody <- gsub('[ ,', '[', jsonbody, fixed = TRUE)

         jsonbody <- gsub('"options"', "{}", jsonbody)
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      force = function(matrix, color = NA, label = NA, size = NA){
         ##matrix- connectivity matrix n by n, where value is the weight of the edge

         listbuilder <- list(type = "force", options = list())
         nodes <- vector(mode = "numeric", length = nrow(matrix))
         for (i in 0:(length(nodes)-1)) {
            nodes[i] <- i
         }
         features <- list(nodes = nodes)

         ##convert matrix to links
         links <- matrix(ncol = 3, byrow = TRUE)
         for (i in 1:ncol(matrix)) {
            for (j in 1:nrow(matrix)) {
               if (!(matrix[i,j] == 0)) {
                  links <- rbind(links, c((i-1), (j-1), (matrix[i,j])))
               }
            }
         }
         links <- links[-1,]
         features$links <- links

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }


         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      graph = function(x, y, matrix, color = NA, label = NA, size = NA) {
         listbuilder <- list(type = "graph", options = list())
         points <- matrix(ncol = 2, nrow = length(x))
         points[,1] <- x
         points[,2] <- y
         features <- list(nodes = points)

         links <- matrix(ncol = 3, byrow = TRUE)
         for (i in 1:ncol(matrix)) {
            for (j in 1:nrow(matrix)) {
               if (!(matrix[i,j] == 0)) {
                  links <- rbind(links, c((i-1), (j-1), (matrix[i,j])))
               }
            }
         }
         links <- links[-1,]
         features$links <- links

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }
         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      map = function(regions, weights, colormap) {
         listbuilder <- list(type = "map", options = list())
         features = list(regions = regions, values = weights, colormap = colormap)
         listbuilder$data <- features

         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      graphbundled = function(x, y, matrix, color = NA, label = NA, size = NA){
         listbuilder <- list(type = "graph-bundled", options = list())
         points <- matrix(ncol = 2, nrow = length(x))
         points[,1] <- x
         points[,2] <- y
         features <- list(nodes = points)

         links <- matrix(ncol = 3, byrow = TRUE)
         for (i in 1:ncol(matrix)) {
            for (j in 1:nrow(matrix)) {
               if (!(matrix[i,j] == 0)) {
                  links <- rbind(links, c((i-1), (j-1), (matrix[i,j])))
               }
            }
         }
         links <- links[-1,]
         features$links <- links

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }

         listbuilder$data <- features

         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      matrix = function(matrix, colormap){
         listbuilder <- list(type = "matrix", options = list())
         features = list(matrix = matrix, colormap = colormap)
         listbuilder$data <- features

         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         return(list(url = url, id = response$id))
      },
      adjacency = function (matrix, label = NA) {
         listbuilder <- list(type = "adjacency", options = list())
         nodes <- c(0:(nrow(matrix) - 1))
         links <- matrix(ncol = 3, byrow = TRUE)
         for (i in 1:ncol(matrix)) {
            for (j in 1:nrow(matrix)) {
               if (!(matrix[i,j] == 0)) {
                  links <- rbind(links, c((i-1), (j-1), (matrix[i,j])))
               }
            }
         }
         links <- links[-1,]
         features <- list(links = links, nodes = nodes)

         if (!is.na(label)) {
            features$label <- label
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      scatterline = function(x, y, t, color = NA, label = NA, size = NA){
         listbuilder <- list(type = "scatter-line", options=list())
         points <- matrix(ncol = 2, nrow = length(x))
         points[,1] <- x
         points[,2] <- y
         features <- list(points = points, series = t)

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }

         listbuilder$data <- features
         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "")
         oldbody = "placeholder"

         while (!(jsonbody == oldbody)) {
            oldbody <- jsonbody
            jsonbody <- gsub(", ,", ",", jsonbody)
         }
         jsonbody <- gsub(",  ]", " ]", jsonbody)
         jsonbody <- gsub('[ ,', '[', jsonbody, fixed = TRUE)

         jsonbody <- gsub('"options"', "{}", jsonbody)
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      scatter3 = function(x, y, z, color = NA, label = NA, size = NA, alpha = NA) {
         listbuilder <- list(type = "scatter-3", options=list())
         points <- matrix(ncol = 3, nrow = length(x))
         points[,1] <- x
         points[,2] <- y
         points[,3] <- z
         features <- list(points = points)

         if (!is.na(color)) {
            features$color <- color
         }
         if (!is.na(label)) {
            features$label <- label
         }
         if (!is.na(size)) {
            features$size <- size
         }
         if (!is.na(alpha)) {
            features$alpha <- alpha
         }

         listbuilder$data <- features

         if(self$notebook) {
            listbuilder$options$width = 600
         }

         jsonbody <- toJSON(listbuilder, .na = "{}")
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      image = function(imgpath) {
         body = list(file = upload_file(imgpath), type = "image", options = list())
         if(self$notebook) {
            body$options$width = 600
         }

         rawresponse <- POST(url = paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""), encode = 'multipart', body = body)
         jsonstring <- rawToChar(rawresponse$content)
         response <- fromJSON(jsonstring)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      gallery = function(imgpathvector) {
         firstpath <- imgpathvector[1]
         otherpaths <- imgpathvector[-1]

         body = list(file = upload_file(firstpath), type = "gallery", options = list())
         if(self$notebook) {
            body$options$width = 600
         }

         rawresponse <- POST(url = paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", sep=""), encode = 'multipart', body = body)
         jsonstring <- rawToChar(rawresponse$content)
         response <- fromJSON(jsonstring)
         url <- paste(self$serveraddress, "visualizations/", response$id, "/", sep="")
         self$url <- url
         for (var in otherpaths) {
            POST(url = paste(self$serveraddress, "sessions/", self$sessionid, "/visualizations/", response$id, "/data/images/", sep=""), encode = 'multipart', body = list(file = upload_file(var), type = "image"))
         }
         if (self$autoopen) {
            browseURL(url)
         }
         if(self$notebook) {
            return(display_html(getURLContent(paste(url, "embed/", sep=""))))
         }
         return(list(url = url, id = response$id))
      },
      createsession = function(sessionname = ""){
         jsonbody <- toJSON(list(name=sessionname))
         httpheader<- c(Accept = "text/plain", "Content-Type" = "application/json")
         response = postForm(paste(self$serveraddress, "sessions/", sep=""),
                             .opts = list(httpheader = httpheader, postfields=jsonbody))
         response = fromJSON(response)
         self$sessionid <- response["id"]
         return(response)
      },
      usesession = function(sessionid){
         self$sessionid <- sessionid
      },
      sethost = function(serveraddress){
         self$url <- NA
         self$sessionid <- NA
         self$serveraddress <- serveraddress
      },
      enableautoopening = function(){
         self$autoopen <- TRUE
      },
      disableautoopening = function(){
         self$autoopen <- FALSE
      },
      openviz = function(vizid = NA){
         if (is.na(vizid)) {
            if (is.na(self$url)) {
               print("No vizualisation to show")
            } else{
               browseURL(self$url)
            }
         } else{
            url <- paste(self$serveraddress, "visualizations/", vizid, "/", sep="")
            browseURL(url)
         }
      }
   )
)
