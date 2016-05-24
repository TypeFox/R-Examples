#######################################################################
##                       S2. Class                                   ##
#######################################################################

multi.ax.S2 <- function(thecall) {
  ## function to plot ADEgS when an s.* function is called and 'xax/yax' arguments are vectors of length > 1
  
  listGraph <- list()
  thenewcall <- thecall
 
  ## update some arguments
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$plot <- FALSE

  if(thenewcall[[1]] == "s.value") {
    if(is.null(thenewcall$psub.position))
      thenewcall$psub.position <- "topleft"
  }

  ## evaluate some arguments in the correct frame
  xax <- eval(thecall$xax, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
  yax <- eval(thecall$yax, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)) 
  ## create ADEg plots
  for(i in yax) {
    for(j in xax) {
      thenewcall$xax <- j
      thenewcall$yax <- i
      thenewcall$psub.text <- paste("xax=", j, ", yax=", i, collapse = "", sep = "")
      listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
    }
  }

  ## create the multiplot ADEgS
  names(listGraph) <- paste("x", apply(expand.grid(xax, yax), 1, paste, collapse = "y"), sep = "")
  posmatrix <- layout2position(c(length(yax), length(xax)), ng = length(listGraph), square = FALSE)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = length(listGraph), nrow = length(listGraph)), Call = as.call(thecall))
  return(object)
}


##
##
##

multi.facets.S2 <- function(thecall, adepar, samelimits = TRUE) {
  ## function to plot ADEgS when the 'facets' argument is used
  listGraph <- list()
  
  oldparamadeg <- adegpar()
  on.exit(adegpar(oldparamadeg))
  adegtot <- adegpar(adepar)
  
  ## update some arguments in the newcall
  thenewcall <- thecall
  thenewcall$plot <- FALSE
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$facets <- NULL
  
  ## evaluate some arguments in the correct frame
  if(thecall[[1]] != "s.match")
    dfxy <- eval(thecall$dfxy, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
  else
    dfxy <- do.call("rbind", list(thecall$dfxy1, thecall$dfxy2), envir = sys.frame(sys.nframe() + thenewcall$pos + 2))

  facets <- factor(eval(thecall$facets, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)))
  
  ## same limits for all sub-graphics
  if((isTRUE(samelimits) | is.null(samelimits)) & (thecall[[1]] != "s.corcircle")) {
    xax <- thenewcall$xax
    yax <- thenewcall$yax
    if(is.null(thenewcall$Sp))
      lim.global <- setlimits2D(minX = min(dfxy[, xax]), maxX = max(dfxy[, xax]), minY = min(dfxy[, yax]), maxY = max(dfxy[, yax]),
        origin = adegtot$porigin$origin, aspect.ratio = adegtot$paxes$aspectratio, includeOr = adegtot$porigin$include)
    else { ## Sp: ex map, alors par defaut on prend la bbox
      limsSp <- bbox(eval(thenewcall$Sp))
      lim.global <- setlimits2D(minX = limsSp[1, 1], maxX = limsSp[1, 2], minY = limsSp[2, 1], maxY = limsSp[2, 2], origin = rep(adegtot$porigin$origin, le = 2), aspect.ratio = adegtot$paxes$aspectratio, includeOr = adegtot$porigin$include) 
    }
    if(is.null(thecall$xlim))
      thenewcall$xlim <- lim.global$xlim
    if(is.null(thecall$ylim))
      thenewcall$ylim <- lim.global$ylim
  }
    
  ## creation of the plots (ADEg objects)
  for(i in 1:nlevels(facets)) {
    thenewcall$psub.text <- levels(facets)[i]

    ## specific arguments for the different functions
    if(thecall[[1]] == "s.match") {
      thenewcall$dfxy1 <- call("[[", call("split", call("as.data.frame", thecall$dfxy1), thecall$facets), i)
      thenewcall$dfxy2 <- call("[[", call("split", call("as.data.frame", thecall$dfxy2), thecall$facets), i)
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)
    } else {
      thenewcall$dfxy <- call("[[", call("split", call("as.data.frame", thecall$dfxy), thecall$facets), i)
    }
    
    if(thecall[[1]] == "s.class") {
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
      thenewcall$wt <- call("[[", call("split", thecall$wt, thecall$facets), i)
    }
    
    if(thecall[[1]] == "s.distri")
      thenewcall$dfdistri <- call("[[", call("split", thecall$dfdistri, thecall$facets), i)
    
    if(thecall[[1]] == "s.image")
      thenewcall$z <- call("[[", call("split", thecall$z, thecall$facets), i)
    
    if(thecall[[1]] == "s.label" || thecall[[1]] == "s.corcircle"|| thecall[[1]] == "s.arrow")
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)
    
    if(thecall[[1]] == "s.logo")
      thenewcall$logos <- call("[[", call("split", thecall$logos, thecall$facets), i)
    
    if(thecall[[1]] == "s.traject") {
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
      if(!is.null(thecall$order))
        thenewcall$order <- call("[[", call("split", thecall$order, thecall$facets), i)
    }
    
    if(thecall[[1]] == "s.value") {
      thenewcall$z <- call("[[", call("split", thecall$z, thecall$facets), i)
      if(is.null(thenewcall$breaks)) {
        ## same breaks for all groups
        z <- eval(thecall$z, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
        breaks <- pretty(z, thenewcall$n)
        thenewcall$breaks <- breakstest(breaks, z, n = length(breaks))
      }
      if(is.null(thenewcall$psub.position))
        thenewcall$psub.position <- "topleft"
    }
   
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }

  ## creation of the multi-plot (ADEgS object)
  names(listGraph) <- levels(facets)
  posmatrix <- layout2position(.n2mfrow(nlevels(facets)), ng = nlevels(facets), square = FALSE)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = nlevels(facets), nrow = nlevels(facets)), Call = as.call(thecall))
  ## change pos et frame a posteriori ??
  return(object)
}

##
##
##

multi.variables.S2 <- function(thecall, arg.vary) {
  ## function to plot ADEgS when an s.* function is called and an argument is multivariable (e.g., z in s.value, fac in s.class, etc)
  ## the name of the varying argument is in name.vary
  
  listGraph <- list()
  thenewcall <- thecall
  ## update some arguments
  thenewcall$pos <- eval(thecall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  name.vary <- thenewcall[[arg.vary]]
  dfvary <- eval(name.vary, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)) 
  
  ## create ADEg plots
  for(j in 1:ncol(dfvary)) {
    thenewcall[[arg.vary]] <- call("[", name.vary, substitute(1:nrow(name.vary)), j)
    thenewcall$psub.text <- colnames(dfvary)[j]
    if(thenewcall[[1]] == "s.class" || thenewcall[[1]] == "s.traject") {
      thenewcall$labels <- call("levels", call("as.factor", thenewcall[[arg.vary]]))
    }

    if(thenewcall[[1]] == "s.value") {
      if(is.null(thenewcall$psub.position))
        thenewcall$psub.position <- "topleft"
    }
    
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }
  
  ## create the multiplot ADEgS
  names(listGraph) <- colnames(dfvary)
  posmatrix <- layout2position(.n2mfrow(ncol(dfvary)), ng = ncol(dfvary), square = FALSE)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = ncol(dfvary), nrow = ncol(dfvary)), Call = as.call(thecall))
  return(object)
}


#######################################################################
##                       C1. Class                                   ##
#######################################################################

multi.score.C1 <- function(thecall) {
  ## function to plot ADEgS when an s1d.* function is called and score is a data.frame with multiple columns
  listGraph <- list()
  thenewcall <- thecall
 
  ## update some arguments
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  if(thenewcall[[1]] != "s1d.interval") {
    score <- eval(thecall$score, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    name.score <- thecall$score
  } else {
    score <- eval(thecall$score1, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    name.score <- thecall$score1
  }
  nc <- ncol(score)
 
  ## create ADEg plots
  for(i in 1:nc) {
    thenewcall$psub.text <- colnames(score)[i]
    ## specific arguments for the different functions
    if(thenewcall[[1]] != "s1d.interval") {
      thenewcall$score <- call("[", thecall$score, substitute(1:nrow(name.score)), i)
   } else {
      thenewcall$score1 <- call("[", thecall$score1, substitute(1:nrow(name.score)), i)
      thenewcall$score2 <- call("[", thecall$score2, substitute(1:nrow(name.score)), i)
    }

    if(thenewcall[[1]] == "s1d.barchart") {
      if(is.null(thenewcall$labels))
        thenewcall$labels <- call("rownames", thecall$score)
    }
 
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }

  ## create the multiplot ADEgS
  names(listGraph) <- colnames(score)
  posmatrix <- layout2position(.n2mfrow(nc), ng = nc)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = length(listGraph), nrow = length(listGraph)), Call = as.call(thecall))
  return(object)
}


##
##
##

multi.facets.C1 <- function(thecall, adepar, samelimits = TRUE) {
  ## function to plot ADEgS when the 'facets' argument is used
  listGraph <- list()
  
  oldparamadeg <- adegpar()
  on.exit(adegpar(oldparamadeg))
  adegtot <- adegpar(adepar)
  
  ## update some arguments in the newcall
  thenewcall <- thecall
  thenewcall$plot <- FALSE
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$facets <- NULL
  
  ## evaluate some arguments in the correct frame
  if(thenewcall[[1]] != "s1d.interval") {
    score <- eval(thecall$score, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
  } else {
    score1 <- eval(thecall$score1, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    score2 <- eval(thecall$score2, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    score <- c(score1, score2)
  }
  
  facets <- factor(eval(thecall$facets, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)))
  
  ## same limits for all graphics
  if(isTRUE(samelimits) | is.null(samelimits)) {
    lim.axe1 <- setlimits1D(min(score), max(score), origin = adegtot$porigin$origin[1], includeOr = adegtot$porigin$include)
    if(adegtot$p1d$horizontal & is.null(thecall$xlim)) {
      thenewcall$xlim <- lim.axe1
    }
    if(!adegtot$p1d$horizontal & is.null(thecall$ylim)) {
      thenewcall$ylim <- lim.axe1
    }
  }
  
  ## creation of the plots (ADEg objects)
  for(i in 1:nlevels(facets)) {
    thenewcall$psub.text <- levels(facets)[i]
    if(thecall[[1]] == "s1d.interval") {
      thenewcall$score1 <- call("[[", call("split", thecall$score1, thecall$facets), i)
      thenewcall$score2 <- call("[[", call("split", thecall$score2, thecall$facets), i)
      thenewcall$at <- call("[[", call("split", thecall$at, thecall$facets), i)
      } else {
      thenewcall$score <- call("[[", call("split", thecall$score, thecall$facets), i)
    }
    
    
    if(thecall[[1]] == "s1d.barchart" & !is.null(thecall$labels))
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)

     if(thecall[[1]] == "s1d.dotplot" | thecall[[1]] == "s1d.curve")
      thenewcall$at <- call("[[", call("split", thecall$at, thecall$facets), i)
     
    if(thecall[[1]] == "s1d.density")
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
    
    if(thecall[[1]] == "s1d.gauss") {
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
      thenewcall$wt <- call("[[", call("split", thecall$wt, thecall$facets), i)
    }
    
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }

  ## creation of the multi-plot (ADEgS object)
  names(listGraph) <- levels(facets)
  posmatrix <- layout2position(.n2mfrow(nlevels(facets)), ng = nlevels(facets))
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = nlevels(facets), nrow = nlevels(facets)), Call = as.call(thecall))
  ## change pos et frame a posteriori ??
  return(object)
}

##
##
##

multi.variables.C1 <- function(thecall, arg.vary) {
  ## function to plot ADEgS when an s1d.* function is called and an argument is multivariable (e.g., fac in s1d.density)
  ## the name of the varying argument is in name.vary
  
  listGraph <- list()
  thenewcall <- thecall
 
  ## update some arguments
  thenewcall$pos <- eval(thecall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  name.vary <- thenewcall[[arg.vary]]
  dfvary <- eval(name.vary, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)) 
  
  ## create ADEg plots
  for(j in 1:ncol(dfvary)) {
    thenewcall[[arg.vary]] <- call("[", name.vary, substitute(1:nrow(name.vary)), j)
    thenewcall$psub.text <- colnames(dfvary)[j]
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }
  
  ## create the multiplot ADEgS
  names(listGraph) <- colnames(dfvary)
  posmatrix <- layout2position(.n2mfrow(ncol(dfvary)), ng = ncol(dfvary))
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = ncol(dfvary), nrow = ncol(dfvary)), Call = as.call(thecall))
  return(object)
}



#######################################################################
##                       S1. Class                                   ##
#######################################################################

multi.score.S1 <- function(thecall) {
  ## function to plot ADEgS when an s1d.* function is called and score is a data.frame with multiple columns
  
  listGraph <- list()
  thenewcall <- thecall
 
  ## update some arguments
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  if(thenewcall[[1]] != "s1d.match") {
    score <- eval(thecall$score, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    name.score <- thecall$score
  } else {
    score <- eval(thecall$score1, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    name.score <- thecall$score1
  }
  
  ## create ADEg plots
  nc <- ncol(score)
  for(i in 1:nc) {
    ## specific arguments for the different functions
    if(thenewcall[[1]] != "s1d.match") {
      thenewcall$score <- call("[", thecall$score, substitute(1:nrow(name.score)), i)
    } else {
      thenewcall$score1 <- call("[", thecall$score1, substitute(1:nrow(name.score)), i)
      thenewcall$score2 <- call("[", thecall$score2, substitute(1:nrow(name.score)), i)
    }
    
    thenewcall$psub.text <- colnames(score)[i]
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }
 
  ## create the multiplot ADEgS
  names(listGraph) <- colnames(score)
  posmatrix <- layout2position(.n2mfrow(nc), ng = nc)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = length(listGraph), nrow = length(listGraph)), Call = as.call(thecall))
  return(object)
}


##
##
##

multi.facets.S1 <- function(thecall, adepar, samelimits = TRUE) {
  ## function to plot ADEgS when the 'facets' argument is used
  listGraph <- list()
  
  oldparamadeg <- adegpar()
  on.exit(adegpar(oldparamadeg))
  adegtot <- adegpar(adepar)
  
  ## update some arguments in the newcall
  thenewcall <- thecall
  thenewcall$plot <- FALSE
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$facets <- NULL
  
  ## evaluate some arguments in the correct frame
  if(thenewcall[[1]] != "s1d.match") {
    score <- eval(thecall$score, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
  } else {
    score1 <- eval(thecall$score1, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    score2 <- eval(thecall$score2, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
    score <- c(score1, score2)
  }
  facets <- factor(eval(thecall$facets, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)))
  
  ## same limits for all graphics
  if(isTRUE(samelimits) | is.null(samelimits)) {
    lim.global <- setlimits1D(min(score), max(score), origin = adegtot$porigin$origin[1], includeOr = adegtot$porigin$include)
    if(adegtot$p1d$horizontal & is.null(thecall$xlim))
      thenewcall$xlim <- lim.global
    if(!adegtot$p1d$horizontal & is.null(thecall$ylim))
      thenewcall$ylim <- lim.global
  }
  
  ## creation of the plots (ADEg objects)
  for(i in 1:nlevels(facets)) {
    thenewcall$psub.text <- levels(facets)[i]

    if(thecall[[1]] == "s1d.match") {
      thenewcall$score1 <- call("[[", call("split", thecall$score1, thecall$facets), i)
      thenewcall$score2 <- call("[[", call("split", thecall$score2, thecall$facets), i)
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)
    } else {
      thenewcall$score <- call("[[", call("split", thecall$score, thecall$facets), i)
    }
    
    if(thecall[[1]] == "s1d.label" & !is.null(thecall$labels))
       thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)

    if(thecall[[1]] == "s1d.class") {
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
      thenewcall$wt <- call("[[", call("split", thecall$wt, thecall$facets), i)
    }

    if(thecall[[1]] == "s1d.distri")
      thenewcall$dfdistri <- call("[[", call("split", thecall$dfdistri, thecall$facets), i)

    if(thecall[[1]] == "s1d.boxplot")
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
         
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }

  ## creation of the multi-plot (ADEgS object)
  names(listGraph) <- levels(facets)
  posmatrix <- layout2position(.n2mfrow(nlevels(facets)), ng = nlevels(facets))
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = nlevels(facets), nrow = nlevels(facets)), Call = as.call(thecall))
  ## change pos et frame a posteriori ??
  return(object)
}

##
##
##

multi.variables.S1 <- function(thecall, arg.vary) {
  ## function to plot ADEgS when an s1d.* function is called and an argument is multivariable (e.g., z in fac in s1d.class)
  ## the name of the varying argument is in name.vary
  
  listGraph <- list()
  thenewcall <- thecall
  
  ## update some arguments
  thenewcall$pos <- eval(thecall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  name.vary <- thenewcall[[arg.vary]]
  dfvary <- eval(name.vary, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)) 
  
  ## create ADEg plots
  for(j in 1:ncol(dfvary)) {
    thenewcall[[arg.vary]] <- call("[", name.vary, substitute(1:nrow(name.vary)), j)
    thenewcall$psub.text <- colnames(dfvary)[j]
    
    if(thenewcall[[1]] == "s1d.class")
      thenewcall$labels <- call("levels", call("as.factor", thenewcall[[arg.vary]]))
    
    if(thenewcall[[1]] == "s1d.boxplot" || thenewcall[[1]] == "s1d.distri")
      thenewcall$at <- call("seq", 1, call("nlevels", call("as.factor", thenewcall[[arg.vary]])))
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }
  
  ## create the multiplot ADEgS
  names(listGraph) <- colnames(dfvary)
  posmatrix <- layout2position(.n2mfrow(ncol(dfvary)), ng = ncol(dfvary), square = FALSE)
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = ncol(dfvary), nrow = ncol(dfvary)), Call = as.call(thecall))
  return(object)
}


#######################################################################
##                       Tr. Class                                   ##
#######################################################################


multi.facets.Tr <- function(thecall, samelimits = TRUE) {
  ## function to plot ADEgS when the 'facets' argument is used
  listGraph <- list()
  
  ## update some arguments in the newcall
  thenewcall <- thecall
  thenewcall$plot <- FALSE
  thenewcall$pos <- eval(thenewcall$pos) - 3
  thenewcall$facets <- NULL
  
  ## evaluate some arguments in the correct frame
	if(thecall[[1]] != "triangle.match")
  	dfxyz <- eval(thecall$dfxyz, envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
	else
  	dfxyz <- do.call("rbind", list(thecall$dfxyz1, thecall$dfxyz2), envir = sys.frame(sys.nframe() + thenewcall$pos + 2))
	
  facets <- factor(eval(thecall$facets, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)))
 
	## same limits for all graphics
  if(isTRUE(samelimits) | is.null(samelimits)) {
    #lim.global <- .trranges(df = dfxyz, scale = thecall$scale, min3 = NULL, max3 = NULL)
    lim.global <- .trranges(df = dfxyz, adjust = TRUE, min3 = NULL, max3 = NULL) 
    if(is.null(thecall$min3d))
      thenewcall$min3d <- lim.global$mini
    if(is.null(thecall$max3d))
      thenewcall$max3d <- lim.global$maxi
  }
  
  ## creation of the plots (ADEg objects)
  for(i in 1:nlevels(facets)) {
    thenewcall$psub.text <- levels(facets)[i]

    ## specific arguments for the different functions
    if(thecall[[1]] == "triangle.match") {
      thenewcall$dfxyz1 <- call("[[", call("split", call("as.data.frame", thecall$dfxyz1), thecall$facets), i)
      thenewcall$dfxyz2 <- call("[[", call("split", call("as.data.frame", thecall$dfxyz2), thecall$facets), i)
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)
    } else {
      thenewcall$dfxyz <- call("[[", call("split", call("as.data.frame", thecall$dfxyz), thecall$facets), i)
    }
    
    if(thecall[[1]] == "triangle.class") {
      thenewcall$fac <- call("[[", call("split", thecall$fac, thecall$facets), i)
      thenewcall$wt <- call("[[", call("split", thecall$wt, thecall$facets), i)
    }
    
    if(thecall[[1]] == "triangle.label")
      thenewcall$labels <- call("[[", call("split", thecall$labels, thecall$facets), i)
 
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }

  ## creation of the multi-plot (ADEgS object)
  names(listGraph) <- levels(facets)
  posmatrix <- layout2position(.n2mfrow(nlevels(facets)), ng = nlevels(facets))
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = nlevels(facets), nrow = nlevels(facets)), Call = as.call(thecall))
  ## change pos et frame a posteriori ??
  return(object)
}

##
##
##

multi.variables.Tr <- function(thecall, arg.vary) {
  ## function to plot ADEgS when an triangle.* function is called and an argument is multivariable (e.g., fac in triangle.class, etc)
  ## the name of the varying argument is in name.vary
  
  listGraph <- list()
  thenewcall <- thecall
  ## update some arguments
  thenewcall$pos <- eval(thecall$pos) - 3
  thenewcall$plot <- FALSE

  ## evaluate some arguments in the correct frame
  name.vary <- thenewcall[[arg.vary]]
  dfvary <- eval(name.vary, envir = sys.frame(sys.nframe() + thenewcall$pos + 2)) 
  
  ## create ADEg plots
  for(j in 1:ncol(dfvary)) {
    thenewcall[[arg.vary]] <- call("[", name.vary, substitute(1:nrow(name.vary)), j)
    thenewcall$psub.text <- colnames(dfvary)[j]
    if(thenewcall[[1]] == "triangle.class") {
      thenewcall$labels <- call("levels", call("as.factor", thenewcall[[arg.vary]]))
    }
    listGraph <- c(listGraph, do.call(as.character(thenewcall[[1]]), thenewcall[-1]))
  }
  
  ## create the multiplot ADEgS
  names(listGraph) <- colnames(dfvary)
  posmatrix <- layout2position(.n2mfrow(ncol(dfvary)), ng = ncol(dfvary))
  object <- new(Class = "ADEgS", ADEglist = listGraph, positions = posmatrix, add = matrix(0, ncol = ncol(dfvary), nrow = ncol(dfvary)), Call = as.call(thecall))
  return(object)
}
