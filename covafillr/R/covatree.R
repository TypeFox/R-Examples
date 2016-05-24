#' A Reference Class for Search Tree Approximated Local Polynomial Regression with covatree.
#' 
#' @field ptr External pointer to the covatree C++ object
#'
#' @examples
#' getRefClass('covatree')
#' fn <- function(x) x ^ 4 - x ^ 2
#' x <- runif(2000,-3,3)
#' y <- fn(x) + rnorm(2000,0,0.1)
#' ct <- covatree(coord = x,obs = y,h = 0.5,p = 3L, minLeft = 50)
#' ct$getDim()
#' x0 <- seq(-1,1,0.1)
#' y0 <- ct$predict(x0)
#' par(mfrow=c(2,1))
#' plot(x0,y0[,1], main = "Function")
#' lines(x0,fn(x0))
#' plot(x0, y0[,2], main = "First derivative")
#' lines(x0, 4 * x0 ^ 3 - 2 * x0)
#' 
#' @export covatree
#' @importFrom methods setRefClass new
#' @exportClass covatree
covatree <- setRefClass("covatree",
                        fields = list(ptr = "externalptr"),
                        methods = list(

                            initialize = function(.Object,
                                                  coord,
                                                  obs,
                                                  h = 1.0, 
                                                  p = 2L,
                                                  minLeft = length(obs)/10,
                                                  ...){
                                "Method to initialize the covafill. coord is a matrix of coordinates, obs is a vector of corresponding observations, h is a vector of bandwiths, p is the polynomial degree, and minLeft is the minimum number of observations that will create a sub tree."
                                ## Check input
                                if(missing(coord) || missing(obs))
                                    stop("coord and obs must be specified")
                                if(!is.numvec(coord) & !is.nummat(coord))
                                    stop("coord must be a numeric vector or matrix.")

                                if(!is.numvec(obs))
                                    stop("Obs must be a numeric vector")
                                
                                if(!is.samelength(coord,obs))
                                    stop("obs and coord must have matching dimensions.")
                                if(is.numvec(coord))
                                    coord <- matrix(coord,ncol=1)
                                if(length(h) > dim(coord)[2])
                                    warning("Additional bandwiths are ignored.")
                                if(dim(coord)[2] > 3)
                                    stop("coord must have between 1 and 3 columns.")
                                h <- rep(h,len = dim(coord)[2])

                                if(!is.numsca(p) & !is.numint(p))
                                    stop("p must be an integer.")

                                if(!is.numsca(minLeft))
                                    stop("minLeft must be a scalar")
                                
                                if(!is.numint(p))
                                    warning("p is coerced to an integer.")

                                p <- as.integer(p)

                                if(p < 1)
                                    stop("p must be 1 or greater")

                                ## Create pointer
                                ptr0 <- .Call("MakeTree",coord,obs,h,p,minLeft,
                                              PACKAGE="covafillr")

                                initFields(ptr = ptr0)

                            },
                            getDim = function(){
                                "Get the dimension of the coordinates."
                                return(.Call("getTreeDim",.self$ptr,
                                              PACKAGE="covafillr"))
                            },
                            predict = function(coord){
                                "Predict function value and first order derivatives with search tree approximated local polynomial regression at coord."
                                d <- .self$getDim()
                                if(!is.numvec(coord) & !is.nummat(coord))
                                    stop("coord must be a numeric vector or matrix.")
                                if(is.numvec(coord))
                                    coord <- matrix(coord,ncol=1)
                                if(dim(coord)[2] != d)
                                    stop(paste("coord must have",d,"columns."))

                                val <- .Call("predictTree",.self$ptr,coord,
                                             PACKAGE="covafillr")

                                if(is.null(colnames(coord))){
                                    cnam <- 1:d
                                }else{
                                    cnam <- colnames(coord)
                                }

                                if(is.null(rownames(coord))){
                                    rnam <- 1:dim(coord)[1]
                                }else{
                                    rnam <- rownames(coord)
                                }
                                
                                colnames(val) <- c('fn',
                                                   paste('gr',cnam,sep='_'))
                                rownames(val) <- rnam

                                
                                return(val)
                            }                            
                        )
                        )
covatree$lock("ptr")
