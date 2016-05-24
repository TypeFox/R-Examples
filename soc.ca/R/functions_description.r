#### Functions for describtion and printing

#' Print soc.ca objects
#' 
#' Prints commonly used measures used in the analysis of multiple correspondence
#' analysis
#' @param x is a soc.ca class object
#' @param ... further arguments are ignored
#' @return Active dimensions is the number of dimensions remaining after the
#'   reduction of the dimensionality of the analysis.
#' @return Active modalities is the number of modalities that are not set as
#'   passive.
#' @return Share of passive mass is the percentage of the total mass that is
#'   represented by the passive modalities.
#' @return The values represented in the scree plot are the adjusted inertias,
#'   see \link{variance}
#' @return The active variables are represented with their number of active
#'   modalities and their share of the total variance/inertia.
#' @export print.soc.mca
#' @seealso \link{soc.mca}, \link{contribution}
#'   
#' @examples
#' example(soc.ca)
#' print(result)
#' @export

print.soc.mca  <- function(x, ...){
  
  # Help functions
  scree       <- function(x, dim = 6){
    set.dim   <- dim
    dim       <- ifelse((nrow(x$adj.inertia)<dim) == TRUE, nrow(x$adj.inertia), dim)
    adj       <- round(x$adj.inertia[1:dim, 4], digits = 1)
    stars     <- round(round(adj)/2)
    starscree <- vector("list", set.dim)
    for (i in 1:length(stars)){
      starscree[[i]] <- noquote(rep("*", stars[i]))
    }
    return(starscree)
    # x is a soc.ca class object
    # Dim is the number of dimensions included in the plot
  }
    
  Nmodal       <- x$n.mod
  Nsup         <- sum(x$freq.sup != 0)
  Nid          <- x$n.ind
  Share.of.var <- round((x$modal[, 3] - 1)/ (length(x$names.passive) + Nmodal - nrow(x$modal)), 2) * 100
  Vnames       <- paste(rownames(x$modal), " [", x$modal[, 3], " - ", format(Share.of.var), "%]", sep = "")
  Vnames       <- Vnames[order(Share.of.var, decreasing = TRUE)]
  Submass 	   <- 1 - round(sum(x$mass.mod), digits = 2) 
  act.dim 	   <- nrow(x$adj.inertia)
  dim80 		   <- which.min(x$adj.inertia[, 5] < 80)
  scree.dim	   <- 7
  N.pas.mod    <- length(x$names.passive)
  stars 		   <- scree(x, scree.dim)
  adj.dim      <- 1:scree.dim
  
  dim.a        <- ifelse((scree.dim < nrow(x$adj.inertia)), scree.dim, nrow(x$adj.inertia))
  adj          <- vector(mode = "numeric", length = scree.dim)
  adj[1:dim.a] <- x$adj.inertia[1:dim.a, 4]
  adj	         <- paste(formatC(adj,format = "f",digits = 1), sep = "", collide = "%")
  
  ## Output
  # Soc.csa title
  if (inherits(x, "soc.csa") == TRUE) cat(format("Class Specific Multiple Correspondence Analysis:", 	width = 90, justify = "centre"),"\n", "\n")
  # Soc.mca title
  if (inherits(x, "soc.csa") == FALSE) cat(format("Specific Multiple Correspondence Analysis:",   width = 90, justify = "centre"),"\n", "\n")
  
  cat(format("Statistics",  width = 50, justify = "centre"), format("Scree plot", width = 40, justify = "centre"),"\n",
      
      format("	Active dimensions: ", 			width = 40,), format(act.dim, 	width = 10, justify = "right"),
      format("|  1.", width = 10, justify = "centre" ), format(adj[1], width = 10, justify = "centre"), format(paste(stars[[1]]), width = 1), "\n",

      format("	Dimensions explaining 80% of inertia: ",width = 40,), format(dim80, 	width = 10, justify = "right"), 
      format("|  2.", width = 10, justify = "centre" ), format(adj[2], width = 10, justify = "centre"), format(paste(stars[[2]]), width = 1), "\n",
      
      format("	Active modalities: ", 			width = 40,), format(Nmodal, 	width = 10, justify = "right"), 
      format("|  3.", width = 10, justify = "centre" ), format(adj[3], width = 10, justify = "centre"), format(paste(stars[[3]]), width = 1), "\n",
      
      format("	Supplementary modalities: ",		width = 40,), format(Nsup, 	width = 10, justify = "right"), 
      format("|  4.", width = 10, justify = "centre" ), format(adj[4], width = 10, justify = "centre"), format(paste(stars[[4]]), width = 1), "\n",
      
      format("	Individuals: ",		 		width = 40,), format(Nid, 	width = 10, justify = "right"), 
      format("|  5.", width = 10, justify = "centre" ), format(adj[5], width = 10, justify = "centre"), format(paste(stars[[5]]), width = 1), "\n",
      
      format("	Share of passive mass:",	 		width = 40,), format(Submass, 	width = 10, justify = "right"), 
      format("|  6.", width = 10, justify = "centre" ), format(adj[6], width = 10, justify = "centre"), format(paste(stars[[6]]), width = 1), "\n",
      
      format("\tNumber of passive modalities:",	 		width = 40, ), format(N.pas.mod, 	width = 10, justify = "right"), 
      format("|  7.", width = 10, justify = "centre" ), format(adj[7], width = 10, justify = "centre"), format(paste(stars[[7]]), width = 1), "\n",
      
      "\n",
      format(paste("The", length(Vnames),"active variables: [No. modalities - share of variance]"), 			width = 100, justify = "centre" ),
      "\n",
      "\n",
      sep = "")
  cat(format(Vnames, width = 25, justify = "right"), fill = 100)
   
}

#' Contribution balance
#' 
#' Calculates the balance of the contribution of each dimension. This measure
#' indicates whether too much of a dimensions contribution is placed on either
#' the + or - side of the dimension.
#' @param object is a soc.ca class object
#' @param act.dim is the number of active dimensions to be measured
#' @return A matrix with the share of contribution on each side of 0 and their
#'   balance (+/-)
#' @seealso \link{soc.mca}, \link{contribution}
#' @examples
#' example(soc.ca)
#' balance(result)
#' balance(result, act.dim = 3)
#' @export

balance   <- function(object, act.dim = object$nd){
  coord   <- object$coord.mod[, 1:act.dim]
  contrib <- object$ctr.mod[, 1:act.dim]
  pm      <- matrix(, nrow = act.dim, ncol = 3)
  for (i in 1:act.dim){
    temp  <- cbind(coord[, i], contrib[, i])
    temp  <- temp[order(temp[, 1]), ]
    plus  <- temp[which(temp[, 1] >= 0), ]
    minus <- temp[which(temp[, 1] <= 0), ]
    pcontrib <- sum(plus[, 2])
    mcontrib <- sum(minus[, 2])
    pm[i, 1] <- pcontrib
    pm[i, 2] <- mcontrib
    pm[i, 3] <- pcontrib/mcontrib
  }
  colnames(pm) <- c("+ Contrib.", "- Contrib.", "Balance (+/-)")
  return(pm)
}

#' Summaries of contribution values
#' 
#' Different forms of contribution summaries for \link{soc.ca} objects. Results
#' are presented according to the specified \code{mode}
#' @param object a \link{soc.ca} object
#' @param dim the included dimensions
#' @param all If TRUE returns all modalities instead of just those that
#'   contribute above average
#' @param indices If TRUE; returns a vector with the row indices of the
#'   modalities or individuals
#' @param mode indicates which form of output. Possible values: \code{"sort"},
#'   \code{"mod"}, \code{"ind"}, \code{"variable"}. If the mode is
#'   \code{"variable"}, \code{dim} can be a sequence of dimensions: \code{1:5}
#' @return Each mode prints different results:
#' @return   \item{"mod"}{Ranks all modalities according to their contribution}
#'   \item{"sort"}{Ranks all modalities according to their contribution and then sorts them according to their coordinates}
#'   \item{"ind"}{Ranks all individuals according to their contribution}
#'   \item{"variable"}{Sorts all modalities according to their variable and sums the contributions per variable}
#' @return The values reported:
#' \item{Ctr}{Contribution values in percentage. Contribution values for individuals are reported in permille}
#' \item{Coord}{Principal coordinates}
#' \item{Cor}{The correlation with the dimension}
#' @seealso \link{map.ctr}
#' @examples
#' 
#' example(soc.ca)
#' contribution(result)
#' contribution(result, 2)
#' contribution(result, dim = 3, all = TRUE)
#' contribution(result, indices = TRUE)
#' contribution(result, 1:2, mode = "variable")
#' @export

contribution <- function(object, dim = 1, all = FALSE, indices = FALSE, mode = "sort"){
  
  if (indices == TRUE & mode == "mod"){
    ctr     <- object$ctr.mod[,dim]
    av.ctr  <- as.vector(apply(as.matrix(ctr), 2, function(x) which(x >= mean(x, na.rm = TRUE))))
    if (is.list(av.ctr) == TRUE) av.ctr  <- unlist(av.ctr[dim], use.names = FALSE)
    av.ctr  <- av.ctr[duplicated(av.ctr) == FALSE]    
    return(av.ctr)
  }
    
  # Modalities
  if (identical(mode, "mod")){
    if (length(dim) > 1 ) stop("This mode does not support more than 1 dimension")
    ctr     <- round(100 * object$ctr.mod[, dim], digits = 1)
    cor     <- round(100 * object$cor.mod[, dim], digits = 1)
    coord   <- round(object$coord.mod[, dim], digits = 2)
    names   <- object$names.mod
    if (identical(all, FALSE) == TRUE){
      av.ctr <- contribution(object, dim = dim, indices = TRUE, mode = mode)    
      header <- paste("The modalities contributing above average to dimension: ", dim, ".", sep = "")
    }
    if (identical(all, TRUE) == TRUE){
      av.ctr <- 1:length(ctr)
      header <- paste("The contribution of all modalities to dimension: ", dim, ".", sep = "")
    }
    
    out           <- data.frame(ctr[av.ctr], cor[av.ctr], coord[av.ctr])
    rownames(out) <- names[av.ctr]
    colnames(out) <- c("   Ctr.", "   Cor." , "   Coord")
    out           <- out[order(-out[, 1]), ]
    maxwidth      <- max(nchar(names)) + sum(nchar(colnames(out)))
    cat("\n", format(header, width = maxwidth, justify = "centre"), "\n", "\n")
    print(out)
  }
  # Individuals  
  if (identical(mode, "ind")){
    individuals(object, dim, indices = indices, all = FALSE)
  }
  # Side sorted modalities
  if (identical(mode, "sort")){
    if(length(dim) > 1 ) stop("Sort mode does not support more than 1 dimension")
    tab.dim(object, dim)
  }
  # Variables
  if (identical(mode, "variable")){
    tab.variable(object, dim)
  }
}

# ' The most contributing individuals
# ' 
# ' Returns the individuals with above average contribution to the selected dimension
# ' @param object is a soc.ca object
# ' @param dim is the included dimensions
# ' @param all: If TRUE returns all individuals instead of just those that contribute above average
# ' @param ind.indices: If TRUE returns a vector with the row indices of the individuals
# ' @return Ctr is the contribution in 1000
# ' @return Cor is the correlation with the dimension
# ' @return Coord is the principal coordinate
# ' @seealso \link{tab.dim}, \link{soc.mca}, \link{contribution}, \link{p.id}
# ' @export

individuals <- function(object, dim = 1, all = FALSE, indices = FALSE){
  
  if (identical(indices, TRUE) == TRUE){
    ctr        <- object$ctr.ind[,dim]
    av.ctr     <- as.vector(apply(as.matrix(ctr), 2, function(x) which(x >= mean(x, na.rm = TRUE))))
    if(is.list(av.ctr) == TRUE) av.ctr  <- unlist(av.ctr[dim], use.names = FALSE)
    av.ctr     <- unique(av.ctr)
    return(av.ctr)
  }else{
    ctr        <- object$ctr.ind[, dim]
    ctr.round  <- round(100 * object$ctr.ind[, dim], 2)
    coord      <- round(object$coord.ind[, dim], 2)
    names      <- object$names.ind
    if (identical(all, FALSE) == TRUE){
      av.ctr   <- individuals(object, dim = dim, indices = TRUE, all = FALSE)    
      header   <- paste("The individuals contributing above average to dimension: ",
                        paste(dim, collapse = ", "), ".", sep = "")
    }
    if (identical(all, TRUE) == TRUE){
      av.ctr   <- 1:length(ctr)
      header   <- paste("The contribution of all individuals to dimension: ",
                        paste(dim, collapse = ", "), ".", sep = "")
    }
    
    out        <- data.frame(ctr.round, coord)[av.ctr, ]
    rownames(out) <- names[av.ctr]
    colnames(out) <- c(paste("   Ctr.", dim), paste("   Coord.", dim))
    out        <- out[order(-out[, 1]), ]
    
    maxwidth   <- max(nchar(names)) + sum(nchar(colnames(out)))
    cat("\n", format(header, width = maxwidth, justify = "centre"), "\n", "\n")
    print(out)
  }

}

# ' The most contributing modalities according to direction on dimension
# ' 
# ' Gives the most contributing modalities sorted according to direction on dimension
# ' @param x is a soc.ca object
# ' @param dim is the dimension
# ' @param label.plus is the label of the dimensions plus side
# ' @param label.minus is the label of the dimensions minus side
# ' @param all defines whether all modalities are to be printed
# ' @seealso \link{contribution}, \link{soc.mca}, \link{p.ctr}
# ' @examples
# ' example(soc.ca)
# ' tab.dim(result, 2)
# ' tab.dim(result, 2, label.plus = "Technical capital", label.minus = "Organizational capital")
# ' @export


tab.dim <- function(x, dim = 1, label.plus = NULL, label.minus = NULL, all = FALSE){
  
  if (identical(label.plus, NULL) == TRUE){
    label.plus    <- paste("Dimension ", dim ,". (+)", sep = "")
  }
  
  if (identical(label.minus, NULL) == TRUE){
    label.minus   <- paste("Dimension ", dim ,". (-)", sep = "")
  }
  
  ctr             <- round(100 * x$ctr.mod[, dim], digits = 1)
  coord           <- round(x$coord.mod[, dim], digits = 2)
  names           <- x$names.mod
  
  if (identical(all, FALSE) == TRUE){
    av.ctr        <- contribution(x, dim = dim, indices = TRUE, mode = "mod")    
  }
  if (identical(all, TRUE) == TRUE){
    av.ctr        <- seq(x$n.mod)
  }
  
  out             <- data.frame(ctr[av.ctr], coord[av.ctr])
  names           <- names[av.ctr]
  maxwidth        <- max(nchar(names))
  
  for (i in seq(names)){
    width         <- maxwidth-nchar(names[i])
    fill          <- paste(rep(" ", width), sep = "", collapse = "")
    names[i]      <- paste(names[i], fill, sep = "", collapse = "")
  }
  rownames(out)   <- names
  ctr.lab         <- paste("Ctr")
  coord.lab       <- paste("Coord")
  colnames(out)   <- c(ctr.lab, coord.lab)
  out             <- out[order(-out[, 1]), ]
  out.label       <- c(ctr.lab, coord.lab)
  outminus        <- out[which(out[, 2] <= 0), ]
  outplus         <- out[which(out[, 2] >= 0), ]
  
  cat("\n", format(label.plus, width = maxwidth, justify = "centre"), "\n")
  print(format(outplus, justify = "centre", width = 8))
  cat("\n", format(label.minus, width = maxwidth, justify = "centre"), "\n")
  print(format(outminus, justify = "centre", width = 8))
  
  }

# ' Contribution per variabel
# ' 
# ' tab.variable returns the contribution values of all modalities ordered by variable
# ' 
# ' @param object is a soc.ca object
# ' @param dim is the included dimensions. The default is 1:3
# ' @param If sup = TRUE the coordinates of the supplementary variables are given instead
# ' @return If assigned using <- tab.variable returns a list of matrixes with the contribution values
# ' @return The returned list is a tab.variable class object and can be exported with the \link{export} function included in the soc.mca package.  
# ' @seealso \link{export}, \link{contribution}
# ' @examples
# ' example(soc.ca)
# ' tab.variable(result)
# ' tab.variable(result, dim = c(1, 3))
# ' tab.variable(result, sup = TRUE)
# ' @export

tab.variable    <- function(object, dim = 1:3, sup = FALSE){
    variable    <- as.factor(object$variable)
    ctr.mod     <- as.matrix(object$ctr.mod[, dim])
    lev.var     <- levels(variable)
    names.mod   <- object$names.mod
    freq.mod    <- object$freq.mod
    
    var.list    <- list()
    for (i in seq(length(lev.var))){
        var.ctr           <- round(ctr.mod[variable == lev.var[i], ] * 100, digits = 1)
        var.ctr           <- cbind(var.ctr, freq.mod[variable == lev.var[i]])
        var.ctr           <- rbind(var.ctr, colSums(var.ctr))
        rownames(var.ctr) <- c(names.mod[variable == lev.var[i]], "Total")
        colnames(var.ctr) <- c(paste(" Dim.", dim, sep = ""), "  Freq")
        
        var.list[[paste(lev.var[i])]] <- var.ctr
    }
    
    ### Supplementary modalities
    
    if (identical(sup, TRUE)){
      
      variable    <- as.factor(object$variable.sup)
      coord.sup   <- object$coord.sup[, dim]
      lev.var     <- levels(variable)
      names.mod   <- object$names.sup
      freq.mod    <- object$freq.sup
      
      var.list    <- list()
      for (i in seq(length(lev.var))){
        var.ctr           <- round(coord.sup[variable == lev.var[i], ], digits = 2)
        var.ctr           <- cbind(var.ctr, freq.mod[variable == lev.var[i]])
        rownames(var.ctr) <- c(names.mod[variable == lev.var[i]])
        colnames(var.ctr) <- c(paste(" Dim.", dim, sep = ""), "  Freq")
        
        var.list[[paste(lev.var[i])]] <- var.ctr
      }
          
      
    }
    
    # The printing
    
    av.ctr      <- round(100/object$n.mod, digits = 1)
    maxwidth    <- max(nchar(names.mod))
    l           <- ncol(var.ctr)
        
    if (identical(sup, FALSE)) cat("The contribution of the active variables")
    if (identical(sup, TRUE))  cat("The coordinates of the supplementary variables")
    
    # Beautiful printing!
    for (i in seq(length(lev.var))){
        var.ctr <- var.list[[i]]
        cat("\n", "\n", format(lev.var[i], width = maxwidth), colnames(var.ctr))
        
        for (q in seq(nrow(var.ctr))){
            cat("\n", format(rownames(var.ctr)[q], width = maxwidth),
                format(var.ctr[q, -l], width = 6), format(var.ctr[q, length(dim) + 1], width = 6, drop0trailing = TRUE))
        }
        
    }
    
    
    if (identical(sup, FALSE)) cat("\n", "Average contribution per modality: ", av.ctr, sep = "")
    cat("\n", "Total number of individuals: ", object$n.ind, sep = "")
    
    class(var.list) <- "tab.variable"
    invisible(var.list)
}


#' Variance table
#' 
#' variance returns a table of variance for the selected dimensions.
#' 
#' @param object is a soc.ca object
#' @param dim is the included dimensions, if set to NULL, then only the
#'   dimensions explaining approx. more than 0.90 of the adjusted variance are included
#' @return If assigned variance returns a matrix version of the table
#'   of variance.
#' @seealso \link{soc.mca}, \link{print.soc.mca}
#' @examples
#' example(soc.ca)
#' variance(result)
#' variance(result, dim = 1:4)
#' @export

variance    <- function(object, dim = NULL){
    
    variance <- object$adj.inertia
    if (identical(dim, NULL) == TRUE){
        dim  <- variance[, 5] <= 91
    }
    variance <- t(variance[dim, ])
    line.dim <- paste(1:ncol(variance), ".", sep = "")
    cat("\n", "Dim        ", format(line.dim, width = 6), sep = "")
    cat("\n", "Eigen    ", format(round(variance[2, ], 2), width = 6), sep = "")
    cat("\n", "Var     ", format(round(variance[3, ], 2), width = 6), sep = "")
    cat("\n", "Adj.Var ", format(round(variance[4, ], 1), width = 6), sep = "")
    cat("\n", "Cum %   ", format(round(variance[5, ], 1), width = 6), sep = "")
    
    invisible(variance)
}



#' Average coordinates
#' 
#' Find the average coordinates for each category in a variable on two dimensions.
#' 
#' @param object is soc.ca result object
#' @param x is a variable of the same length and order as the active variables used to construct the soc.ca object
#' @param dim is the two dimensions used
#' @return a matrix with the mean points and frequencies of the given variable
#' @export
#' @examples
#' example(soc.ca)
#' average.coord(result, sup$Income)

average.coord <- function(object, x, dim = c(1, 2)){
  
  coord             <- object$coord.ind[, dim]
  point.coord       <- aggregate(coord, list(x), mean)
  point.coord[, 2]  <- point.coord[, 2] / sqrt(object$eigen[dim[1]])
  point.coord[, 3]  <- point.coord[, 3] / sqrt(object$eigen[dim[2]])
  point.coord       <- cbind(point.coord, table(x)) 
  
  colnames(point.coord) <- c("label", "X", "Y", "group", "Freq")
  point.coord
}

#' CSA measures
#' 
#' Several measures for the evaluation of the relations between the dimensions of the CSA and the dimensions the of original MCA
#' 
#' @param csa.object is a "soc.csa" class object created by the \link{soc.csa} function
#' @param correlations if TRUE correlations calculated by the \link{cor} function is returned
#' @param cosines if TRUE cosine similarities are returned
#' @param cosine.angles if TRUE angles are calculated in the basis of the cosine values
#' @param dim is the dimensions included
#' @param format if TRUE results are formatted, rounded and printed for screen reading, if FALSE the raw numbers are returned
#' @param ... furhter arguments are send to the \link{cor} function
#' @return A list of measures in either formatted or raw form.
#' @export
#' @examples
#' example(soc.csa)
#' csa.measures(res.csa)
#' csa.measures(res.csa, correlations = FALSE, cosine.angles = FALSE, dim = 1:10, format = FALSE)

csa.measures <- function(csa.object, correlations = TRUE, cosines = TRUE, cosine.angles = TRUE,
                         dim = 1:5, format = TRUE, ...){
  
  csca.coord      <- csa.object$coord.ind
  object          <- csa.object$original.result
  class.indicator <- csa.object$original.class.indicator
  ca.coord        <- object$coord.ind[class.indicator, ]
  
  csca.coord      <- csca.coord[, dim]
  ca.coord        <- ca.coord[, dim]
  
  ##################################
  # Correlations
  cor.mat             <- cor(csca.coord, ca.coord, ...)
  rownames(cor.mat)   <- paste("CSA:", dim)
  colnames(cor.mat)   <- paste("MCA:", dim)
  
  ####################################
  # Cosines similarity
  cosine.similarity      <- function(x, y) x %*% y / sqrt(x %*% x * y %*% y)
  cosine.mat             <- matrix(ncol = ncol(ca.coord), nrow = ncol(ca.coord))
  rownames(cosine.mat)   <- paste("CSA:", dim)
  colnames(cosine.mat)   <- paste("MCA:", dim)
    
    for (i in 1:ncol(ca.coord)){
      cosine.mat[i, ]       <- apply(ca.coord, 2, cosine.similarity, csca.coord[, i])
    }

  #####################################
  # Angles
  cosine.to.angle       <- function(x) acos(x)/pi * 180 
  angles                <- cosine.to.angle(cosine.mat)
  
  ####################################
  # Out.list
  out.list              <- list() 
  if (cosines == TRUE)        out.list$cosines      <- cosine.mat
  if (cosine.angles == TRUE)  out.list$angles       <- angles
  if (correlations == TRUE)   out.list$cor          <- cor.mat
  
  if (identical(format, FALSE)) return(out.list)
  
  ###################################
  # Formatted output
  if (identical(format, TRUE)){
    
    cat("\n", format("Measures for Class Specific Multiple Correspondence Analysis:",
                     width = 90, justify = "centre"), "\n", "\n")
    
  #############
  # Cosines
  if (cosines == TRUE){
   
    cat("\n", format("Cosine similarity:",   width = 10, justify = "right"), "\n", "\n")
    cosine.fat           <- round(cosine.mat, 2) 
    rownames(cosine.fat) <- format(rownames(cosine.fat), width = 10, justify = "centre")
    colnames(cosine.fat) <- format(colnames(cosine.fat), width = 8, justify = "right")
    cosine.fat           <- format(cosine.fat, width = 8, justify = "left")
    print(noquote(cosine.fat))
    cat("\n", "\n")
    }
    
  #############
  # Angles
  if (cosine.angles == TRUE){
    cat("\n", format("Cosine angles:",   width = 10, justify = "right"), "\n", "\n")
    
    angles.fat           <- round(angles, 1) 
    rownames(angles.fat) <- format(rownames(angles.fat), width = 10, justify = "centre")
    colnames(angles.fat) <- format(colnames(angles.fat), width = 8, justify = "right")
    angles.fat           <- format(angles.fat, width = 8, justify = "left")
    print(noquote(angles.fat))
    cat("\n", "\n") 
    }

    ##############
    # Correlations
    if(correlations == TRUE){
      cat("\n", format("Correlations:",   width = 10, justify = "right"), "\n", "\n")
      cor.fat           <- round(cor.mat, 2) 
      rownames(cor.fat) <- format(rownames(cor.mat), width = 10, justify = "centre")
      colnames(cor.fat) <- format(colnames(cor.mat), width = 8, justify = "right")
      cor.fat           <- format(cor.fat, width = 8, justify = "left")
      print(noquote(cor.fat))
      cat("\n", "\n")
    }
  }
}