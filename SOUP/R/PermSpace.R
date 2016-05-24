#-----------------------------------#
#    'PermSpace' CLASS DEFINITION    #
#-----------------------------------#
setClassUnion("arrayOrList", c("array", "list"))

#' Permutation Space
#' 
#' Contains the permutation space of the test statistic, useful for 
#' reproducibility of the analyses
#' 
#' @title Class \code{PermSpace}
#' @aliases PermSpace
#' @aliases PermSpace-class
#' @aliases show,PermSpace,PermSpace-method
#' @aliases print,PermSpace,PermSpace-method
#' @aliases plot,PermSpace,PermSpace-method
#' @rdname PermSpace-class
#' @keywords classes
#' @section Objects from the Class:
#'     Objects can be created by calls of the form 
#'     \code{new("PermSpace", ...)}. 
#'     It contains information of permutation spaces 
#'     used in the analysis, the combined test statistics and \emph{p}-values, 
#'     \code{IDs} (row indexes) and the \code{seed} for the \code{RNG}, 
#'     the \code{rawStats} (non-combined test statistics) and 
#'     \code{comb.funct} (the nonparametric combining function). 
#'     But objects of the class are principally supposed to be created and 
#'     used internally for storing results of \code{\link{SOUP}}.
#' @section Slots:
#'     \describe{
#'     \item{\code{seed}:}{
#'         \code{integer} seed for the Random Number Generator}
#'     \item{\code{T.H0Low}:}{
#'         \code{matrix} containing the permutation space of \emph{combined} 
#'         test statistics with null hypothesis 
#'         \eqn{H_0: x_i \ge x_h}{H_0: x_i <= x_h}, \eqn{i < h}, 
#'         \eqn{i,h = 1,\ldots,G}}
#'     \item{\code{T.H0Gre}:}{
#'         \code{matrix} containing the permutation space of \emph{combined} 
#'         test statistics with null hypothesis 
#'         \eqn{H_0: x_i \le x_h}{H_0: x_i <= x_h}, \eqn{i < h}, 
#'         \eqn{i,h = 1,\ldots,G}}
#'     \item{\code{P.H0Low}:}{
#'         \code{matrix} containing the permutation space of \emph{combined} 
#'         \emph{p}-values with null hypothesis 
#'         \eqn{H_0: x_i \ge x_h}{H_0: x_i <= x_h}, \eqn{i < h}, 
#'         \eqn{i,h = 1,\ldots,G}}
#'     \item{\code{P.H0Gre}:}{
#'         \code{matrix} containing the permutation space of \emph{combined} 
#'         \emph{p}-values with null hypothesis 
#'         \eqn{H_0: x_i \le x_h}{H_0: x_i <= x_h}, \eqn{i < h}, 
#'         \eqn{i,h = 1,\ldots,G}}
#'     \item{\code{IDs}:}{
#'         \code{matrix} permutation space of row indexes}
#'     \item{\code{rawStats}:}{
#'         \eqn{3}-way \code{array} containing the permutation space of 
#'         \emph{non-combined} test statistics}
#'     \item{\code{comb.funct}:}{
#'         nonparametric combining function used for \code{\link{NPC}} of 
#'         \code{rawStats}}
#'     }
#' @section Methods:
#'     \describe{
#'     \item{initialize}{
#'         constructor used when calling \code{new(PermSpace, ...)}}
#'     \item{show}{
#'         \code{signature(object = "PermSpace")}: shows only the main 
#'         information (on screen) for the object}
#'     \item{print}{
#'         \code{signature(x = "PermSpace")}: It prints the whole object 
#'         on screen (mostly useful for external saving)}
#'     \item{\code{signature(x = "PermSpace")}}{
#'         Plots a bivariate representation of the permutation space, when 
#'         there are more than \eqn{2} (original) variables then a Principal 
#'         Component Analysis is performed and the first \eqn{2} variables in 
#'         the transformed space are shown}
#'     }
#' @examples 
#'     showClass("PermSpace")
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @exportClass PermSpace
#' @import tensor
#' @import methods

setClass("PermSpace",
    representation = representation(
        seed       = "integer",
        T.H0Low    = "arrayOrList",
        T.H0Gre    = "arrayOrList",
        P.H0Low    = "arrayOrList",
        P.H0Gre    = "arrayOrList",
        IDs        = "array",
        rawStats   = "array",
        comb.funct = "character"
    )
)#-END-

#---------------------------#
#    'PermSpace' METHODS        #
#---------------------------#
##- inizializator (constructor)
#' Methods for function \code{initialize} in package \pkg{SOUP}: constructor
#' 
#' @name initialize
#' @aliases initialize,PermSpace-method
#' @docType methods
#' @rdname initialize-methods
#' @section Methods:
#'     \describe{
#'     \item{\code{signature(object = "PermSpace")}}{
#'          constructor of the object (\code{initialize})}
#'     }
#' @exportMethod initialize
setMethod(
    f = 'initialize',
    signature = "PermSpace",
    definition = function(.Object, seed, T.H0Low, T.H0Gre, P.H0Low, P.H0Gre,          
            IDs, rawStats, comb.funct)
    {
        if(missing(seed))
        {
            .Object@seed <- integer()
        } else
        {
            .Object@seed <- as.integer(seed)
        }# END:if-seed
        if(missing(T.H0Low))
        {
            .Object@T.H0Low <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@T.H0Low <- T.H0Low
        }# END:if-T.H0Low
        if(missing(T.H0Gre)) 
        {
            .Object@T.H0Gre <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@T.H0Gre <- T.H0Gre
        }# END:if-T.H0Low
        if(missing(P.H0Low))
        {
            .Object@P.H0Low <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@P.H0Low   <- P.H0Low
        }# END:if-P.H0Low
        if(missing(P.H0Gre))
        {
            .Object@P.H0Gre <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@P.H0Gre <- P.H0Gre
        }# END:if-P.H0Gre
        if(missing(IDs))
        {
            .Object@IDs <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@IDs <- IDs
        }# END:if-IDs
        if(missing(rawStats))
        {
            .Object@rawStats <- array(NA, c(0, 0, 0))
        } else
        {
            .Object@rawStats <- rawStats
        }# END:if-rawStats
        if(missing(comb.funct))
        {
            .Object@comb.funct <- character()
        } else
        {
            .Object@comb.funct <- comb.funct
        }# END:if-comb.funct
        return(.Object)
    }# END:definition
)# END:initialize

##- show
#' \describe{
#' \item{\code{signature(object = "PermSpace")}}{
#'      Shows only the main information (on screen) for the object}
#' }
#' 
#' @name show
#' @aliases show,PermSpace-method
#' @docType methods
#' @rdname show-methods
#' @exportMethod show

setMethod(
    f = 'show',
    signature = "PermSpace",
    definition = function(object) {
        cat("*** \"PermSpace\" object of package \"SOUP\" ***\n")
        
        ##  statistics H0:'<='
        cat("* Permutation Space of the combined test-statistics with H0:'<=' (limited to 21x10 matrices) *\n")
        if(is.list(object@T.H0Low))
        {
            nm <- names(object@T.H0Low)
            for(i in seq_along(object@T.H0Low)) {
                dimStats <- c(min(21, NROW(object@T.H0Low[[i]])), min(10, NCOL(object@T.H0Low[[i]])))
                print(nm[i])
                print(object@T.H0Low[[i]][seq_len(dimStats[1]), seq_len(dimStats[2])])
                if((NROW(object@T.H0Low[[i]]) > 21) || (NCOL(object@T.H0Low[[i]]) > 10)) {
                    cat("...\t...\t...\n")
                } else { }# END:if-high-dims
            }# END:for-i
        } else
        {
            dimStats <- c(min(21, NROW(object@T.H0Low)), min(10, NCOL(object@T.H0Low)))
            print(object@T.H0Low[seq_len(dimStats[1]), seq_len(dimStats[2])])
            if((NROW(object@T.H0Low) > 21) || (NCOL(object@T.H0Low) > 10)) {
                cat("...\t...\t...\n")
            } else {}# END:if-high-dims
        }# END:if-list
        ## END:T.H0Low
        
        ##  statistics H0:'>='
        cat("* Permutation Space of the combined test-statistics with H0:'>=' (limited to 21x10 matrices) *\n")
        if(is.list(object@T.H0Gre))
        {
            for(i in seq_along(object@T.H0Gre)) {
                dimStats <- c(min(21, NROW(object@T.H0Gre[[i]])), min(10, NCOL(object@T.H0Gre[[i]])))
                print(nm[i])
                print(object@T.H0Gre[[i]][seq_len(dimStats[1]), seq_len(dimStats[2])])
                if((NROW(object@T.H0Gre[[i]]) > 21) || (NCOL(object@T.H0Gre[[i]]) > 10)) {
                    cat("...\t...\t...\n")
                } else { }# END:if-high-dims
            }# END:for-i
        } else
        {
            dimStats <- c(min(21, NROW(object@T.H0Gre)), min(10, NCOL(object@T.H0Gre)))
            print(object@T.H0Gre[seq_len(dimStats[1]), seq_len(dimStats[2])])
            if((NROW(object@T.H0Gre) > 21) || (NCOL(object@T.H0Gre) > 10)) {
                cat("...\t...\t...\n")
            } else { }# END:if-high-dims
        }# END:if-list
        ## END:T.H0Gre
        
        ## p.values H0:'<='
        cat("* Permutation Space of the p.values with H0:'<=' (limited to 21x10 matrices) *\n")
        if(is.list(object@P.H0Low))
        {
            for(i in seq_along(object@P.H0Low)) {
                dimP.values <- c(min(21, NROW(object@P.H0Low[[i]])), min(10, NCOL(object@P.H0Low[[i]])))
                print(nm[i])
                print(object@P.H0Low[[i]][seq_len(dimP.values[1]), seq_len(dimP.values[2])])
                if((NROW(object@P.H0Low[[i]]) > 21) || (NCOL(object@P.H0Low[[i]]) > 10)){
                    cat("...\t...\t...\n")
                } else { }# END:if-high-dims
            }# END:for-i
        } else
        {
            dimP.values <- c(min(21, NROW(object@P.H0Low)), min(10, NCOL(object@P.H0Low)))
            print(object@P.H0Low[seq_len(dimP.values[1]), seq_len(dimP.values[2])])
            if((NROW(object@P.H0Low) > 21) || (NCOL(object@P.H0Low) > 10)){
                cat("...\t...\t...\n")
            } else { }# END:if-high-dims
        }# END:if-list
        ## END:P.H0Low
        
        ## p.values H0:'>='
        cat("* Permutation Space of the p.values with H0:'>=' (limited to 21x10 matrices) *\n")
        if(is.list(object@P.H0Gre))
        {
            for(i in seq_along(object@P.H0Gre)) {
                dimP.values <- c(min(21, NROW(object@P.H0Gre[[i]])), min(10, NCOL(object@P.H0Gre[[i]])))
                print(nm[i])
                print(object@P.H0Gre[[i]][seq_len(dimP.values[1]), seq_len(dimP.values[2])])
                if((NROW(object@P.H0Gre[[i]]) > 21) || (NCOL(object@P.H0Gre[[i]]) > 10)){
                    cat("...\t...\t...\n")
                } else { }# END:if-high-dims
            }# END:for-i
        } else
        {
            dimP.values <- c(min(21, NROW(object@P.H0Gre)), min(10, NCOL(object@P.H0Gre)))
            print(object@P.H0Gre[seq_len(dimP.values[1]), seq_len(dimP.values[2])])
            if((NROW(object@P.H0Gre) > 21) || (NCOL(object@P.H0Gre) > 10)){
                cat("...\t...\t...\n")
            } else { }# END:if-high-dims
        }# END:if-list
        ## END:P.H0Gre
        
        ## comb.funct
        cat("* Combining function: ", object@comb.funct, ", seed for the R.N.G.: ", object@seed, " *\n", sep = "")
        cat("*** End \"PermSpace\" ***\n")
    }
)# END:show


##- print
#' \describe{
#' \item{\code{signature(x = "PermSpace")}}{
#'      It prints the whole object on screen (mostly useful for 
#'      external saving)}
#' }
#' 
#' @name print
#' @docType methods
#' @aliases print,PermSpace-method
#' @rdname print-methods
#' @exportMethod print

setGeneric("print")
setMethod(
    f = 'print', 
    signature = "PermSpace",
    definition = function(x, ...) {
        cat("*** \"PermSpace\" object of package \"SOUP\" ***\n")
        ##  statistics H0:'<='        
        cat("* Permutation Space of the combined test-statistics with H0:'<=' *\n")
        print(x@T.H0Low)
        ##  statistics H0:'>='
        cat("* Permutation Space of the combined test-statistics with H0:'>=' *\n")
        print(x@T.H0Gre)
        ## p.values H0:'<='
        cat("* Permutation Space of the p.values with H0:'<=' *\n")
        print(x@P.H0Low)
        ## p.values H0:'>='
        cat("* Permutation Space of the p.values with H0:'>=' *\n")
        print(x@P.H0Gre)
        ## comb.funct
        cat("* Combining function:", x@comb.funct, ", seed for the R.N.G.: ", x@seed, " *\n", sep = "")
        cat("*** End \"PermSpace\" ***\n")
    }
)# END:print


##- plot (using PCA)
#' Methods for function \code{plot} in package \pkg{SOUP}
#' 
#' @name plot
#' @docType methods
#' @aliases plot,PermSpace-method
#' @rdname PermSpace-class
#' @exportMethod plot
setGeneric("plot")
setMethod(
    f = "plot",
    signature = c("PermSpace"),
    definition = function(x, main = NULL, 
            xlab = NULL, ylab = NULL , ...)
    {
        plotStats <- function(stats)
        {
            switch(as.character(NCOL(stats)),
                # '1' = {
                    # hist(stats, freq = FALSE)
                    # lines(density(stats), col = 4)
                    # rug(stats)
                    # title("Univariate Permutation Space")
                # },
                
                '2' = {
                    plot(stats[, 1], stats[, 2], lwd = 1, pty = "o", xlab = colnames(stats)[1],
                        ylab = colnames(stats)[2])
                    points(stats[1, 1], stats[1, 2], col =  "red", lwd = 3)
                    title("Bivariate Permutation Space") 
                },
                
                {
                    pc <- prcomp(stats[, apply(stats, 2, var) > 0], scale. = TRUE, center = FALSE)
                    pc$rotation[1,] <- pc$rotation[1,] * sign(pc$x[1, 1])
                    pc$rotation[2,] <- pc$rotation[2,] * sign(pc$x[1, 2]) 
                    pc$rotation[3,] <- pc$rotation[3,] * sign(pc$x[1, 3]) 
                    pc$x[, 1] <- pc$x[, 1] * sign(pc$x[1, 1])
                    pc$x[, 2] <- pc$x[, 2] * sign(pc$x[1, 2]) 
                    lam <- pc$sdev[1:2]
                    pc$x[, 1:2] <- pc$x[, 1:2] / lam
                    pc$rotation[, 1:2] <- pc$rotation[, 1:2] * lam
                    plot(pc$x[, 1], pc$x[, 2], lwd = 1, pty = "o", xlim = range(pc$x[, 1]) * 1.2,
                        ylim = range(pc$x[, 2]) * 1.2,
                        xlab = paste("PC1 (", round(pc$sdev[1]^2 / sum(pc$sdev^2) * 100,2), " %)", sep = ""),
                        ylab = paste("PC2 (", round(pc$sdev[2]^2 / sum(pc$sdev^2) * 100,2), " %)", sep = ""),
                        col = "gray", pch = 21, bg = "gray"
                    )
                    points(pc$x[1, 1], pc$x[1, 2], col = "red", lwd = 3, pch = 21, bg = "red")
                    text(pc$x[1, 1] * 1.1, pc$x[1, 2] * 1.1, col = "red", "Obs")
                    arrows(0, 0, 2 * pc$rotation[, 1], 2 * pc$rotation[, 2], lwd = 1, col = 1)
                    text(2.1 * pc$rotation[, 1], 2.1 * pc$rotation[, 2], rownames(pc$rotation), cex = 1.5, col = "black")
                    title("PCA of Permutation Space") 
                }
            )
        }# END:plotStats
        
        if(is.list(x@T.H0Low))
        {
            for(i in seq_along(x@T.H0Low))
            {
#                windows()
                par(mfrow = c(2, 1))
                plotStats(stats = x@T.H0Low[[i]])
                title(main = "\n \n statistics with H0:'<='", cex.main = .8, font.main = 3)
                plotStats(stats = x@T.H0Gre[[i]])
                title(main = "\n \n statistics with H0:'>='", cex.main = .8, font.main = 3)
            }# END:for
        } else
        {
#            windows()
            par(mfrow = c(2, 1))
            plotStats(stats = x@T.H0Low)
            title(main = "\n \n statistics with H0:'<='", cex = .8, font.main = 3)
            plotStats(stats = x@T.H0Gre)
            title(main = "\n \n statistics with H0:'>='", cex = .8, font.main = 3)
        }# END:if-list
    }# END:definition
)

##- documenting it
# promptClass("PermSpace", file = "man/PermSpace.Rd")
# promptMethods("plot", file = "man/PermSpace-plot.Rd")

