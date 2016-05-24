library(methods)

#' @import igraph
#' @importFrom ggplot2 scale_size_continuous scale_size_identity geom_point aes ggplot geom_tile scale_fill_brewer scale_alpha_identity scale_colour_brewer geom_abline
NULL

#' Use first two spatial coordinates as position in layout
#' 
#' @note You can pass any igraph layout algorithm to plot
#' 
#' @param g GNG object
#' 
#' @export
gng.plot.layout.v2d <- function(g){
  cbind(V(g)$v0, V(g)$v1)
}

gng.plot.color.label <- 'label'

gng.plot.color.fast.cluster <- 'fast.cluster'

gng.plot.color.cluster <- 'cluster'

gng.plot.color.none <- 'none'

gng.plot.layout.igraph.fruchterman.fast <- layout.fruchterman.reingold

gng.plot.layout.igraph.auto <- layout.auto

gng.plot.2d <- "2d"

gng.plot.2d.errors <- "2d.errors"

gng.type.default <- function(){
	c(2)
}

gng.type.optimized <- function(minimum=0, maximum=10){
  c(0, minimum, maximum)
}

gng.type.utility<- function(k=1.3){
  c(1, k)
}

.gng.dataset.bagging.prob <- 3
.gng.dataset.bagging <- 2
.gng.dataset.sequential <-1

.GNG <- NULL

#' Plot GNG
#'
#' @title plot GNG object
#' @description Plot resulting graph using igraph plotting
#' @rdname plot.gng
#' @export
#' @method plot Rcpp_GNGServer
#'
#' @param x GNG object
#' @param mode \code{"2d"} (igraph plot)
#' \code{"2d.errors"} (igraph plot with mean error log plot)
#' 
#' @param layout igraph layout to be used when plotting. Defaults to \code{layout.fruchterman.reingold}. 
#' Other good choice is using \code{gng.plot.layout.v2d}, which returns two first spatial coordinates.
#' 
#' @param vertex.color How to color vertexes. Possible values: \code{"fast.cluster"} (vertex color is set to fastgreedy.community clustering),
#' \code{"label"} (rounds to integer label if present), \code{list of integers} (colors vertices according to provided list), \code{"none"} (every node is white),
#' 
#' @param vertex.size Size of plotted vertices
#' @param ... other arguments not used by this method.
#' 
#' @note If you want to "power-use" plotting and plot for instance a subgraph, you might be interested in
#' exporting igraph with convertToIGraph function 
#' 
#' @examples
#' \dontrun{
#' gng <- GNG(scaled.wine)
#' # Plots igraph using first 2 coordinates and colors according to clusters
#' plot(gng, mode=gng.plot.2d.errors, layout=gng.plot.layout.v2d, vertex.color=gng.plot.color.cluster)
#' 
#' # For more possibilities see gng.plot.* constants
#' }
plot.Rcpp_GNGServer <- NULL

#' Save model to binary format
#' @title gngSave
#' @description Writes model to a disk space efficient binary format. 
#' @export
#'  
#' @param object GNG object
#' @param filename File where binary will be saved
gngSave <- NULL


#' Load model from binary format
#'
#' @title gngLoad
#' @description Writes model to a disk space efficient binary format. 
#' @export
#' 
#' @param filename Binary file location
gngLoad <- NULL

#' Get centroids
#'
#' @title calculateCentroids
#' @description Using passed community.detection finds communities and for each community pick node with biggest betweenness score
#' @export
#' 
#' @param object GNG object
#' @param community.detection.algorithm Used algorithm from igraph package, by default spinglass.community
#' 
#' @examples
#' \dontrun{
#' gng <- GNG(gng.preset.sphere(100))
#' print(node(gng, calculateCentroids(gng)[1])$pos)
#' }
calculateCentroids <- NULL

#' Find closest node
#'
#' @title findClosests
#' @description Finds closest node from given list to vector. Often used together with calculateCentroids
#' @export

#' @param object GNG object
#' @param node.ids List of indexes of nodes in gng. 
#' @param x Can be either \code{vector} or \code{data.frame.}
#' 
#' @examples
#' \dontrun{
#' gng <- GNG(gng.preset.sphere(100))
#' # Find closest centroid to c(1,1,1)
#' found.centroids <- calculateCentroids(gng)
#' findClosests(gng, found.centroids, c(1,1,1))
#' }
findClosests <- NULL

#' Check if GNG is running
#'
#' @title isRunning
#' @description Returns TRUE if GNG object is training
#' @export
#' @param object GNG object
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' # FALSE, because did not pass train.online to constructor
#' print(isRunning(gng))
#' 
isRunning <- function(object) {
  return(object$.isRunning())
}

#' Find closest component
#' @name predictComponent
#' @title predictComponent
#' @description Finds connected component closest to given vector(s). On the first
#' execution of function strongly connected components are calculated using igraph::cluster function.
#' @export
#' @rdname predictComponent-methods
#' @docType methods
#'
#' @param object GNG object
#' @param x Can be either \code{vector} or \code{data.frame}.
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' # Find closest component to c(1,1,1)
#' predictComponent(gng,  c(1,1,1))
#' 
#' @aliases predictComponent
predictComponent <- NULL

#' Get GNG node
#' @name node
#' @title node
#' @description Retrieves node from resulting graph
#' @rdname node-methods
#' @export
#' 
#' @param x GNG object
#' @param gng_id Id of the node to retrieve. This is the id returned by functions like predict, or centroids
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' print(node(gng, 10)$pos)
#' 
#' @aliases node
#' 
node <- function(x, gng_id) UseMethod("node")

#' Predict 
#' @name predict.gng
#' @title predict
#' @description Retrieves prediction from trained GNG model
#' @rdname predict.gng
#' @export
#' 
#' @param object Trained model
#' @param x Vector or matrix of examples
#' @param ... other arguments not used by this method
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' predict(gng, c(1,2,2))
predict.Rcpp_GNGServer <- function(object, x, ...){
  if( is.vector(x)){
    object$.predict(x)
  }else{
    if ( !is(x, "data.frame") && !is(x, "matrix") && !is(x,"numeric")  ) {
      stop(gmum.error(GMUM_WRONG_PARAMS, "Wrong target class, please provide data.frame, matrix or numeric vector"))
    }
    
    if (!is(x, "matrix")) {
      x <- data.matrix(x)
    }
    
    y <- rep(NA, nrow(x))
    
    for(i in 1:nrow(x)){
      y[i] <- object$.predict(x[i,])
    }
    
    y
  }
}

#' @export 
node.Rcpp_GNGServer <- NULL

#' @name run
#' @title run
#' @rdname run-methods
#' @description Run algorithm (in parallel)
#' @export
#' 
#' @param object GNG object
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' run(gng)
#' print(isRunning(gng))
#' 
run <- function(object) UseMethod("run")

#' @export
run.Rcpp_GNGServer <- NULL

#' @title pause
#' @description Pause algorithm
#' @export
#'
#' @param object GNG object
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' pause(gng)
#' print(isRunning(gng))
pause <- function(object) UseMethod("pause")

#' @export
pause.Rcpp_GNGServer <- NULL

#' @title terminate
#' @name terminate
#' @description Terminate algorithm
#' @export
#' @rdname terminate-methods
#' @docType methods
#' 
#' @param object GNG object
#'
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' terminate(gng)
#' 
#' @aliases terminate
#'
terminate <- function(object) UseMethod("terminate")

#' @export
terminate.Rcpp_GNGServer <- NULL

#' @title meanError
#' @description Gets mean error of the graph (note: blocks the execution, O(n))
#' @param object GNG object
#' 
#' @export
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' meanError(gng)
meanError <- NULL

#' @title errorStatistics
#' @description Gets vector with errors for every second of execution
#' @export
#' 
#' @param object GNG object
#' 
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' errorStatistics(gng)
errorStatistics <- NULL

#' @title Constructor of Optimized GrowingNeuralGas object. 
#' @rdname optimized-gng
#' 
#' @export 
#' 
#' @description Construct simplified and optimized GNG object. Can be used to train offline, or online. Data dimensionality shouldn't be too big, if
#' it is consider using dimensionality reduction techniques.
#'
#' @param beta Decrease the error variables of all node 
#' nodes by this fraction (forgetting rate). Default 0.99
#' 
#' @param alpha Decrease the error variables of the nodes neighboring to 
#' the newly inserted node by this fraction. Default 0.5
#' 
#' @param lambda New vertex is added every lambda iterations. Default 200
#' 
#' @param max.nodes Maximum number of nodes 
#' (after reaching this size it will continue running, but new noes won't be added)
#' 
#' @param eps.n Strength of adaptation of neighbour node. Default \code{0.0006}
#' 
#' @param eps.w Strength of adaptation of winning node. Default \code{0.05}
#' 
#' @param max.iter If training offline will stop if exceedes max.iter iterations. Default \code{200}
#'
#' @param train.online If used will run in online fashion. Default \code{FALSE}
#'
#' @param min.improvement Used for offline (default) training. 
#' Controls stopping criterion, decrease if training stops too early. Default \code{1e-3}
#'
#' @param dim Used for training online, specifies dataset example dimensionality
#'
#' @param value.range All example features should be in this range, required for optimized version of the algorithm. Default \code{(0,1)} 
#' 
#' @param x Passed data (matrix of data.frame) for offline training
#' 
#' @param labels Every example can be associated with labels that are added to nodes later. By default empty
#' 
#' @param max.edge.age Maximum edge age. Decrease to increase speed of change of graph topology. Default \code{200}
#' 
#' @param verbosity How verbose should the process be, as integer from \eqn{[0,6]}, default: \code{0}
#' 
#' @param seed Seed for internal randomization
#' 
#' @examples
#' \dontrun{
#' # Train online optimizedGNG. All values in this dataset are in the range (-4.3, 4.3)
#' X <- gng.preset.sphere(100)
#' gng <- OptimizedGNG(train.online = TRUE, value.range=c(min(X), max(X)), dim=3, max.nodes=20)
#' insertExamples(gng, X)
#' run(gng)
#' Sys.sleep(10)
#' pause(gng)
#' }
OptimizedGNG <- NULL

#' @name clustering
#' @title clustering
#' 
#' @description Gets vector with node indexes assigned to examples in the dataset
#' 
#' @method clustering Rcpp_GNGServer 
#' @export
#' 
#' @rdname clustering-methods
#' 
#' @docType methods
#'
#' @examples
#' gng <- GNG(gng.preset.sphere(100))
#' clustering(gng)
clustering.Rcpp_GNGServer <- NULL

#' @title Constructor of GrowingNeuralGas object. 
#' 
#' @rdname gng
#' 
#' @export 
#' 
#' @description Construct GNG object. Can be used to train offline, or online.
#' 
#' @param beta Decrease the error variables of all node 
#' nodes by this fraction (forgetting rate). Default 0.99
#' 
#' @param alpha Decrease the error variables of the nodes neighboring to 
#' the newly inserted node by this fraction. Default 0.5
#' 
#' @param lambda Every lambda iteration is added new vertex. Default 200
#' 
#' @param max.nodes Maximum number of nodes 
#' (after reaching this size it will continue running, but won't add new nodes)
#' 
#' @param eps.n How strongly adapt neighbour node. Default \code{0.0006}
#' 
#' @param eps.w How strongly adapt winning node. Default \code{0.05}
#' 
#' @param max.iter Uf training offline will stop if exceedes max.iter iterations. Default \code{200}
#'
#' @param train.online default FALSE. If used will run in online fashion
#'
#' @param min.improvement Used for offline (default) training. 
#' Controls stopping criterion, decrease if training stops too early. Default \code{1e-3}
#'
#' @param dim Used for training online, specifies training example size
#'
#' @param k Utility constant, by default turned off. Good value is 1.3. Constant controlling speed of erasing obsolete nodes, 
#' see \url{http://sund.de/netze/applets/gng/full/tex/DemoGNG/node20.html}
#' 
#' @param x Passed data (matrix of data.frame) for offline training
#' 
#' @param labels Every example can be associated with labels that are added to nodes later. By default empty
#' 
#' @param max.edge.age Maximum edge age. Decrease to increase speed of change of graph topology. Default \code{200}
#' 
#' @param verbosity How verbose should the process be, as integer from \eqn{[0,6]}, default: \code{0}
#' 
#' @param seed Seed for internal randomization
#' 
#' @examples
#' \dontrun{
#' X <- gng.preset.sphere(100)
#' y <- round(runif(100))
#' # Train in an offline manner
#' gng <- GNG(X, labels=y, max.nodes=20)
#' # Plot
#' plot(gng)
#'
#' # Train in an online manner with utility (erasing obsolete nodes)
#' gng <- GNG(max.nodes=20, train.online=TRUE, k=1.3, dim=3)
#' insertExamples(gng, X, labels=y)
#' run(gng)
#' Sys.sleep(10)
#' terminate(gng)
#' # Plot
#' plot(gng)
#' }
GNG <- NULL

#' @title convertToIGraph
#' @description Converts GNG to igraph object, where every vertex contains attributes gng.index, error, data.label and 3 first spatial coordinates (as attributes v0, v1, v2).
#' Additionally utility attribute is present if utility GNG is used.
#' 
#' @param object GNG object
#' @param calculate.dist If true will calculate all \code{n^2} distances in the graph
#' 
#' @export
convertToIGraph <- NULL

#' @title numberNodes
#' @description Get current number of nodes in the graph
#' 
#' @param object GNG object
#' 
#' @export
numberNodes <- function(object){
  object$getNumberNodes()
}

#' @name insertExamples
#' @title insertExamples
#' @description Insert examples with optional labels.
#' 
#' @export
#' 
#' @param object GNG object
#' @param examples \code{matrix} or \code{data.frame} with rows as examples. Note: if training online make sure
#' number of columns matches dim parameter passed to GNG constructor.
#' @param labels \code{vector} of labels, that will be associated with nodes in the graph. GNG will assign to each
#' node a mean of labels of closest examples.
#'
#' @examples
#' X <- gng.preset.sphere(100)
#' gng <- GNG(X, train.online=TRUE)
#' # Add more examples
#' X = gng.preset.sphere(100)
#' insertExamples(gng, X)
#' 
#' @note It copies your examples twice in RAM. You might want to use object$.insertExamples.
insertExamples <- NULL

.GNG <- function(x=NULL, labels=c(),
                  beta=0.99, 
                  alpha=0.5, 
                  max.nodes=100, 
                  eps.n=0.0006, 
                  eps.w= 0.05, 
                  max.edge.age = 200, 
                  type = gng.type.default(),
                  max.iter=200,
                  train.online=FALSE,
                  min.improvement=1e-3,
                  lambda=200,
                  dim=-1,
                  verbosity=0,
                  seed=-1
){
  
  
  config <- new(GNGConfiguration)
  
  config$seed = seed
  
  if(is.data.frame(x)){
    warning("Converting data.frame to matrix. Please make sure you pass only numerics to GNG.")
    x <- as.matrix(x)
  }
  
  # Fill in configuration
  if(train.online){
    if(is.null(x)){
      if (dim == -1) {
        stop(gmum.error(GMUM_WRONG_PARAMS, "To train online, please pass desired dimensionality in dim parameter"))
      }
      config$dim = dim
    }else{
      config$dim = ncol(x)
    }
    config$max_iter = -1
  }else{
    config$dim = ncol(x)
    config$max_iter = max.iter
  }
  
  if(type[1] == gng.type.optimized()[1]){
    config$.uniformgrid_optimization = TRUE
    config$.lazyheap_optimization = TRUE  
    config$.set_bounding_box(type[2], type[3])
    
    if(!train.online){
      if(!max(x) <= type[3] && !min(x) >= type[2]){
        stop(gmum.error(GMUM_WRONG_PARAMS, "Passed incorrect parameters. The dataset is not in the defined range"))
      }
    }
    
  }else{
    config$.uniformgrid_optimization = FALSE
    config$.lazyheap_optimization = FALSE
  }
  
  if(type[1] == gng.type.utility()[1]){
    config$.experimental_utility_k = type[2]
    config$.experimental_utility_option = 1
  }
  else{
    config$.experimental_utility_option = 0
  }
  
  
  config$.dataset_type=.gng.dataset.bagging
  config$beta = beta
  config$max_edge_age = max.edge.age
  config$alpha = alpha  
  config$max_nodes = max.nodes
  config$eps_n = eps.n
  config$eps_w = eps.w
  
  config$lambda = lambda
  config$verbosity = verbosity
  
  if(!config$.check_correctness()){
    stop(gmum.error(GMUM_WRONG_PARAMS, "Passed incorrect parameters."))
  }
  
  # Construct server
  server = new(GNGServer, config)
  
 
  if(train.online){
    if(!is.null(x)){
      insertExamples(server, x, labels)
      run(server)
    }
  }
  if(!train.online){
    
    print("Training offline")
    if(is.null(x)){
      stop(gmum.error(GMUM_ERROR, "Passed null data and requested training offline"))
    }else{
      insertExamples(server, x, labels)
      run(server)
      
      max_iter = max.iter
      min_relative_dif = min.improvement
      iter = 0
      previous_iter = -1
      best_so_far = 1e10
      initial_patience = 3
      error_index = -1 # always bigger than 0
      patience = initial_patience
      
      tryCatch({
        
        while(server$getCurrentIteration() == 0 || server$.isRunning()){}
        
        # max_iter is checked in GNG
        while(iter == 0 || server$.isRunning()){
          Sys.sleep(0.1)
          iter = server$getCurrentIteration()
          
          if(previous_iter != iter && iter %% (max_iter/100) == 0){    
            print(paste("Iteration", iter))
          }
          
          
          
          if(length(server$getErrorStatistics()) > 5){
            errors = server$getErrorStatistics()
            
            best_previously = min(errors[(length(errors)-5):length(errors)])
            
            # this is same as (best_so_far-best_previously)/best_so_far < min_relative_di
            # we get minimum of window 5 and look at the history
            if( (error_index - server$.getGNGErrorIndex()) > 4 && 
                  (best_so_far - best_previously) < best_so_far*min_relative_dif){
              patience = patience - 1
              if(patience <= 0){
                print(sprintf("Best error during training: %f", best_so_far))
                print(sprintf("Best error in 5 previous iterations %f", best_previously))
                print(errors[(length(errors)-5):length(errors)])
                print("Patience (which you can control) elapsed, bailing out")
                break
              }
            }else{
              patience = initial_patience
            }
            
            
            error_index = server$.getGNGErrorIndex()
            best_so_far = min(best_previously, best_so_far)
          }
          
        }
        
        print(paste("Iteration", iter))
        previous_iter = iter
        
        
        
        if(server$.isRunning()){
          terminate(server)
        }
        
        server$.updateClustering()
        
      }, interrupt=
        function(interrupt){
          if(server$.isRunning()){
            terminate(server)
          }
          
        })
      
    }
  }else{
  }
  
  
  
  server
}


GNG <- function(x=NULL, labels=c(),
                 beta=0.99, 
                 alpha=0.5, 
                 max.nodes=1000, 
                 eps.n=0.0006, 
                 eps.w=0.05, 
                 max.edge.age=200,
                 train.online=FALSE,
                 max.iter=200,
                 dim=-1,
                 min.improvement=1e-3,
                 lambda=200,
                 verbosity=0,
                 seed=-1,
                 k=NULL
){
  gng <- NULL
  call <- match.call(expand.dots = TRUE)
  if(is.null(k)){
    gng <- .GNG(x=x, labels=labels, beta=beta, alpha=alpha, max.nodes=max.nodes, 
                eps.n=eps.n, eps.w=eps.w, max.edge.age=max.edge.age, type=gng.type.default(), train.online=train.online, max.iter=max.iter, dim=dim, min.improvement=min.improvement, lambda=lambda, verbosity=verbosity, seed=seed)
  }else{
    gng <- .GNG(x=x, labels=labels, beta=beta, alpha=alpha, max.nodes=max.nodes, 
                eps.n=eps.n, eps.w=eps.w, max.edge.age=max.edge.age, type=gng.type.utility(k=k), train.online=train.online, max.iter=max.iter, dim=dim, min.improvement=min.improvement, lambda=lambda, verbosity=verbosity, seed=seed)		
  }
  assign("call", call, gng)
  gng
}

OptimizedGNG <- function(x=NULL, labels=c(),
                          beta=0.99, 
                          alpha=0.5, 
                          max.nodes=1000, 
                          eps.n=0.0006, 
                          eps.w= 0.05, 
                          max.edge.age = 200,
                          train.online=FALSE,
                          max.iter=200,
                          dim=0,
                          min.improvement=1e-3,
                          lambda=200,
                          verbosity=0,
                          seed=-1,
                          value.range=c(0,1)
){
  if(value.range[1] >= value.range[2]){
    stop(gmum.error(GMUM_ERROR, "Incorrect range"))
    return()
  }
  call <- match.call(expand.dots = TRUE)
  gng <- .GNG(x=x, labels=labels, beta=beta, alpha=alpha, max.nodes=max.nodes, 
              eps.n=eps.n, eps.w=eps.w, max.edge.age=max.edge.age, type=gng.type.optimized(minimum=value.range[1], maximum=value.range[2]), train.online=train.online, max.iter=max.iter, dim=dim, min.improvement=min.improvement, lambda=lambda, verbosity=verbosity, seed=seed)
  assign("call", call, gng)
  gng
}    

predictComponent <- function(object, x){
  tryCatch(if(is.null(object$components.membership)){
    assign("components.membership", clusters(convertToIGraph(object))$membership, object)
  }, error=function(...) 
    assign("components.membership", clusters(convertToIGraph(object))$membership, object))
  
  object$components.membership[predict(object, x)]
}

plot.Rcpp_GNGServer <- function(x, vertex.color=gng.plot.color.cluster, 
                      layout=layout.fruchterman.reingold, mode=gng.plot.2d, 
                      vertex.size=3, ...){
  if(vertex.size <= 0){
    stop("Please pass positivie vertex.size")
  }
  
  if(!(is.list(vertex.color) || is.vector(vertex.color) || vertex.color %in% c(gng.plot.color.cluster, 
                           gng.plot.color.fast.cluster, gng.plot.color.label, gng.plot.color.none))){
    stop("Please pass correct vertex.color")
  }
  
  
  if(x$getNumberNodes() > 4000){
    warning("Trying to plot very large graph (>4000 nodes). It might take a long time especially if using layout function.")
  }
  
  if(x$getNumberNodes() == 0){
    warning("Empty graph")
    return()
  }
  
  if(mode == gng.plot.2d){
    .gng.plot2d(x, vertex.color, layout, vertex.size=vertex.size)
  }
  else if(mode == gng.plot.2d.errors){
    .gng.plot2d.errors(x, vertex.color, layout, vertex.size=vertex.size)
  }else{
    stop("Unrecognized mode")
  }
}


node.Rcpp_GNGServer  <- function(x, gng_id){
  x$getNode(gng_id)
}

run.Rcpp_GNGServer  <- function(object){
  # Invalidate components
  assign("components.membership", NULL, object)
  object$.run()
}

pause.Rcpp_GNGServer  <- function(object){
  object$.pause()
  n = 0.0
  sleep = 0.1
  while(object$.isRunning()){
    Sys.sleep(sleep)  
    n = n + 1
    if(n > 2/sleep){
      print("Warning: GNG has not paused! Check status with gng$.isRunning(). Something is wrong.")
      return()
    }
  }
}

terminate.Rcpp_GNGServer <- function(object){
  object$.terminate()
}

meanError <- function(object){
  object$getMeanError()
}  

errorStatistics <- function(object){
  object$getErrorStatistics()
}  

clustering.Rcpp_GNGServer <- function(c){
  c$getClustering()
}  

gngSave <- function(object, filename){
  warning("Saving does not preserve training history")
  object$.save(filename)
}

gngLoad <- function(filename){
  warning("Saving does not preserve training history")
  fromFileGNG(filename)
}

calculateCentroids  <- function(object, community.detection.algorithm=spinglass.community){
  ig <- convertToIGraph(object)
  
  cl = clusters(ig)
  components = lapply(levels(as.factor(cl$membership)), function(x) induced.subgraph(ig, cl$membership==as.numeric(x)))
  
  centroids <- c()
  for(cc in components){
    communities <- community.detection.algorithm(cc)
    for(i in 1:length(communities)){
      #Get subcommunity
      community_graph <- induced.subgraph(cc, which(membership(communities)==i))
      #Get index of centroid (which is ordered by betwenness)
      centroid_index = which(order(betweenness(community_graph))==1)
      # Append
      centroids<- c(centroids, V(community_graph)$gng.index[centroid_index])
    }
  }
  centroids
}


convertToIGraph <- function(object, calculate.dist=TRUE){
  was_running = object$.isRunning()
  if(was_running){
    pause(object)
  }
  
  if(object$getNumberNodes() == 0){
    return(graph.empty(n=0, directed=FALSE))
  }
  
  #Prepare index map. Rarely there is a difference in indexing
  #due to a hole in memory representation of GNG graph (i.e.
  #indexing in gng can be non-continuous)
  
  # Warning: This is a hack. If there is a bug look for it here
  indexesGNGToIGraph <- 1:(2*object$.getLastNodeIndex()) 
  indexesIGraphToGNG <- 1:object$getNumberNodes()
  
  if(object$.getLastNodeIndex() != object$getNumberNodes()){
    igraph_index = 1
    for(i in (1:object$.getLastNodeIndex())){
      node <- node(object, i)
      if(length(node) != 0){
        indexesGNGToIGraph[i] = igraph_index
        indexesIGraphToGNG[igraph_index] = i
        igraph_index = igraph_index + 1
      }
    }
  }
  
  adjlist<-list()
  for(i in 1:object$.getLastNodeIndex()){
    node <- node(object, i)
    if(length(node) != 0){
      
      igraph_index = indexesGNGToIGraph[i]
      #print(paste(igraph_index, node$neighbours))
      neighbours = node$neighbours[node$neighbours > i]
      adjlist[[igraph_index]] <- sapply(neighbours, function(x){ indexesGNGToIGraph[x] })
    } else{
      print("Empty node")
    }
  }
  
  
  g <- graph.adjlist(adjlist, mode = "all", duplicate=FALSE)
  for(i in 1:object$.getLastNodeIndex()){
    node <- node(object, i)
    if(length(node) != 0){
      igraph_index = indexesGNGToIGraph[i]
      #TODO: it is more efficient to assign whole vectors
      #TODO: refactor in whole code v0 v1 v2 to pos_1 pos_2 pos_3
      V(g)[igraph_index]$v0 <- node$pos[1]
      V(g)[igraph_index]$v1 <- node$pos[2]
      V(g)[igraph_index]$v2 <- node$pos[3]
      V(g)[igraph_index]$label <- node$index
      V(g)[igraph_index]$data.label <- node$label
      V(g)[igraph_index]$error <- node$error
      V(g)[igraph_index]$gng.index <- node$index
      if(!is.null(node$utility)){
        V(g)[igraph_index]$utility = node$utility
      }
    } 
  }
  
  if(calculate.dist){
    # Add distance information
    dists <- apply(get.edges(g, E(g)), 1, function(x){ 
      object$.nodeDistance(indexesIGraphToGNG[x[1]], indexesIGraphToGNG[x[2]])
    })
    E(g)$dists = dists
  }
  
  if(was_running){
    run(object)
  }
  
  g
}

findClosests <- function(object, node.ids, x){
  .findClosests <- function(object, node.ids, x){
    # Returns all dists from given pos to given nodes
    get_all_dists <- function(pos, nodes, gng){
      sapply(nodes, function(node_index) sqrt(sum((pos-node(gng, node_index)$pos)^2)))
    }
    
    which.min(get_all_dists(x, node.ids, object))
  }
  if( is.vector(x)){
    .findClosests(object, node.ids, x)
  }else{
    if ( !is(x, "data.frame") && !is(x, "matrix") && !is(x,"numeric")  ) {
      stop(gmum.error(GMUM_WRONG_PARAMS, "Wrong target class, please provide data.frame, matrix or numeric vector"))
    }
    
    if (!is(x, "matrix")) {
      x <- data.matrix(x)
    }
    
    y <- rep(NA, nrow(x))
    
    for(i in 1:nrow(x)){
      y[i] <- .findClosests(object, node.ids, x[i,])
    }
    
    y
  }
}

insertExamples <- function(object, examples, labels=c()){    
  if(length(labels) == 0){
    object$.insertExamples(examples)
  }else if(typeof(labels) == "character"){
    if(typeof(labels) == "list"){
      if(is.null(examples$labels)){
        stop(gmum.error(GMUM_WRONG_PARAMS, "Empty labels column"))
      }else{
        label.column <- examples$labels
        examples$labels <- NULL
        object$.insertLabeledExamples(examples, label.column)
      }
    }else{
      stop(gmum.error(GMUM_WRONG_PARAMS, "Please pass data frame"))
    }
  }else{
    object$.insertLabeledExamples(examples, labels)
  }    
}

loadModule('gng_module', TRUE)

#' Class Rcpp_GNGServer.
#'
#' Class \code{Rcpp_GNGServer} defines a GNGServer class. 
#'
#' @name Rcpp_GNGServer-class
#' @exportClass Rcpp_GNGServer
setClass(Class = "Rcpp_GNGServer")

# Lazy loading to allow for discovery of all files
evalqOnLoad( {
  .wine <<- NULL
})

