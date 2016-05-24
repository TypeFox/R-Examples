##############################################
##               general class              ##
##############################################

setClass(
  Class = "ADEgS",
  slots = c(
    ADEglist = "list",
    positions = "matrix",
    add = "matrix",  ## n*n, if xij = 1, j superposed to i
    Call = "call"),
  ## slots checking
  validity = function(object) {
    ng <- length(object@ADEglist)
    add <- object@add
    if(ncol(object@positions) != 4)
      stop("wrong positions matrix, only 4 columns expected (x0, y0, x1, y1)")
    if(nrow(object@positions) != ng)
      stop("not enough positions: rows number of the positions matrix should be equal to the number of graphics in the ADEglist")
    ## checking add:
    if((NROW(add) != NCOL(add)) | (NCOL(add) != ng))
      stop("add matrix dimensions are not equal to the number of graphics in ADEglist")
    if(any(add != 0 & add != 1))
      stop("add matrix can only contain 0/1 values")
    for(i in 1:ng) {
      j <- 1:i
      if(any(add[i, j] != 0))
        stop("upper diagonal matrix expected for add, only 0 are allowed for xij, when j > = i")
    }
    return(TRUE)
  })

##############################################          
##                 initialize               ##
##############################################       

setMethod(
  f = "initialize",
  signature = "ADEgS",
  function(.Object, ADEglist, positions, add, Call) {
    ## add linking
    superpose <- list()
    ng <- length(ADEglist)
    for(i in 1:ng) {
      superpose <- c(superpose, list(which(add[, i] == 1)))  ## where i is superposed to 1
      if(length((superpose[[i]]))) {
        for(j in superpose[[i]]) {
          add[superpose[[j]], i] <- 1
          superpose[[i]] <- c(superpose[[i]], superpose[[j]])
        }}}
    .Object@add <- add
    
    ## check names of the list ADEglist
    if(is.null(names(ADEglist)))
      names(ADEglist) <- paste("g", lapply(1:length(ADEglist), function(i) i), sep="")
    else
      names(ADEglist) <- make.names(names(ADEglist),unique = TRUE)
    
    ## assignation
    .Object@ADEglist <- ADEglist
    .Object@positions <- positions
    .Object@add <- add
    .Object@Call <- Call
    ## checking validations
    validObject(.Object)
    return(.Object)
  })


setClassUnion(name = "ADEgORADEgSORtrellis", members = c("ADEg", "ADEgS", "trellis"))


##############################################          
##      Get elements/information            ##
##############################################

setMethod(
  f = "getcall",
  signature = "ADEgS",
  definition = function(object) {
    return(object@Call)
  })


setMethod(
  f = "getgraphics",
  signature = "ADEgS",
  definition = function(object) {
    return(object@ADEglist)
  })

setMethod(
  f = "getpositions",
  signature = "ADEgS",
  definition = function(object) {	
    return(object@positions)
  })


setMethod(
  f = "length",
  signature = "ADEgS",
  definition = function(x) {
    return(length(x@ADEglist))
  })


setMethod(
  f = "names",
  signature = "ADEgS",
  definition = function(x) {
    return(names(x@ADEglist))
  })


setMethod(
  f = "names<-",
  signature = c("ADEgS", "character"),
  definition = function(x, value) {
    nameobj <- deparse(substitute(x))
    names(x@ADEglist) <- value
    x
  })

##############################################          
##             Extract graphics             ##
##############################################

## [: if drop =TRUE can return a ADEg or a length-1 ADEgS
## else return a ADEgS no manipulation is made on the positions

setMethod(
  f = "[",
  signature = c("ADEgS", "numeric", "missing", "logical"),
  definition = function(x, i, j, drop = TRUE) {
    if(drop && length(i) == 1)
      return(x@ADEglist[[i]])   ## return(adeG)
    else
      return(new(Class = "ADEgS", ADEglist = x@ADEglist[i], positions = x@positions[i, , drop = drop], add = x@add[i, i, drop = drop], Call = match.call()))
  })


setMethod(
  f = "[",
  signature = c("ADEgS", "numeric", "missing", "missing"),
  definition = function(x, i, j, drop) {
    object <- x[i, drop = FALSE]
    object@Call <- match.call()
    return(object)
  })


setMethod(
  f = "$",
  signature = "ADEgS",
  definition = function(x, name) {
    invisible(x@ADEglist[[name]])
  })

setMethod(
  f = "[[",
  signature = c("ADEgS", "numeric", "missing"),
  definition = function(x, i, j, ...) {
    invisible(x@ADEglist[[i]])
  })


setMethod(
  f = "[[",
  signature = c("ADEgS", "character", "missing"),
  definition = function(x, i, j, ...) {
    invisible(x@ADEglist[[i]])
  })

setMethod(
  f = "[[<-",
  signature = c("ADEgS", "numeric", "missing", "ADEg"),
  definition = function(x, i, j, ..., value) {
    x@ADEglist[[i]] <- value
    invisible(x)
  })

setMethod(
  f = "[[<-",
  signature = c("ADEgS", "numeric", "missing", "ADEgS"),
  definition = function(x, i, j, ..., value) {
    x@ADEglist[[i]] <- value
    invisible(x)
  })

##############################################          
##               superposition              ##
##############################################

setMethod(
  f = "superpose",
  signature = c("ADEgS", "ADEgORtrellis", "numeric", "logical"),
  definition = function(g1, g2, which, plot) {
    ## new ADEgS
    ngraph <- length(g1)
    if(which > ngraph)
        stop("Values in 'which' should be lower than the length of g1")
    
    if(!inherits(g1[[which]], "ADEg")) 
      stop("superposition is only available between two ADEg")
    addi <- cbind(rbind(g1@add, rep(0, ngraph)), rep(0, ngraph + 1))
    addi[which, ngraph + 1] <- 1  ## new graph superpose to which
    ADEglist <- g1@ADEglist
    ADEglist[[ngraph + 1]] <- g2
    ADEgS <- new(Class = "ADEgS", ADEglist = ADEglist, positions = rbind(g1@positions, g1@positions[which,]), add = addi, Call = match.call())
    if(plot) 
      print(ADEgS)
    invisible(ADEgS)
  })


setMethod(
  f = "superpose",
  signature = c("ADEgS", "ADEgORtrellis", "numeric", "ANY"),
  definition = function(g1, g2, which, plot) {
    objectnew <- superpose(g1, g2, which = which, plot = FALSE)
    objectnew@Call <- match.call()
    if(plot)
      print(objectnew)
    invisible(objectnew)            
  })


setMethod(
  f = "superpose",
  signature = c("ADEgS", "ADEgORtrellis", "missing", "ANY"),
  definition = function(g1, g2, which, plot) {
    if(!inherits(g1[[length(g1)]], "ADEg"))
      stop("superposition is only available between two ADEg")
    objectnew <- superpose(g1, g2, which = length(g1), plot = FALSE)
    objectnew@Call <- match.call()
    if(plot)
      print(objectnew)
    invisible(objectnew)
  })


setMethod(
  f = "superpose",
  signature = c("ADEgS", "ADEgS", "missing", "ANY"),
  definition = function(g1, g2, which, plot) {
    ## superpose two ADEgS which have the same number of graphics and the same positions 
    if(length(g1) != length(g2))
      stop("The two ADEgS objects should contain the same number of graphics")
    if(!isTRUE(all.equal(g1@positions, g2@positions, check.attributes = FALSE)))
      stop("The two ADEgS objects should have the same 'positions' slot")
    
    f1 <- function(x, y) {
      if(inherits(x, "ADEg")) {
        addi <- matrix(0, 2, 2)
        addi[1,2] <- 1
        thecall <- call("superpose", x@Call, y@Call)
        obj <- new(Class = "ADEgS", ADEglist = list(x, y), positions = matrix(rep(c(0, 1), each = 4), 2, 4), add = addi, Call = thecall)
      } else if(inherits(x, "ADEgS")) {
        addi <- x@add
        ng <- ncol(addi)
        posi <- x@positions
        
        ## check that positions in posi are all equal and one 1 in each column of addi (i.e. graphs are still superposed)
        checkadd <- all(colSums(addi[, -1, drop = FALSE]) > 0) 
        checkpos <- isTRUE(all.equal(matrix(posi[1, ], nrow = nrow(posi), ncol = ncol(posi), byrow = TRUE), posi))
        if(!checkpos | !checkadd)
          stop("ADEgS object should contain only superposition")
        
        ## superpose
        addi <- rbind(addi, rep(0, ng))
        addi <- cbind(addi, rep(0, ng + 1))
        addi[ng, ng + 1] <- 1
        
        posi <- rbind(posi, posi[1,])
        thecall <- call("superpose", x@Call, y@Call)
        obj <- new(Class = "ADEgS", ADEglist = c(x@ADEglist, list(y)), positions = posi, add = addi, Call = thecall)
      }
      invisible(obj)
    }
    
    res <- lapply(1:length(g1), FUN = function(i) {f1(g1[[i]], g2[[i]])})
    obj <- new(Class = "ADEgS", ADEglist = res, positions = g1@positions, add = g1@add, Call = match.call())
    
    if(plot)
      print(obj)
    invisible(obj)
  })

  
setMethod(
  f = "+",
  signature = c("ADEgS", "ADEg"),
  definition = function(e1, e2) {
    newobj <- superpose(e1, e2, plot = TRUE)
    newobj@Call <- match.call()
    return(newobj)
  })


setMethod(
  f = "+",
  signature = c("ADEg", "ADEgS"),
  definition = function(e1, e2) {
    newobj <- superpose(e2, e1, plot = TRUE)
    warning("the second graph is below the first one ; the reverse situation is not yet implemented", call. = FALSE)
    newobj@Call <- match.call()
    return(newobj)
  })


setMethod(
  f = "cbindADEg", 
  signature = c("ADEgORADEgSORtrellis", "ADEgORADEgSORtrellis"),
  definition = function(g1, g2, ..., plot = FALSE) {
  	if(try(is.list(...), silent = TRUE) == TRUE)
      glist <- as.list(c(g1, g2, ...))
   
	  else
	    glist <- list(g1, g2, ...)
    
    nbg <- length(glist)
    obj <- ADEgS(adeglist = glist, layout = c(1, nbg), add = matrix(0, ncol = nbg, nrow = nbg), plot = FALSE)
    obj@Call <- match.call()

    if(plot)
      print(obj)
    invisible(obj)
  })


setMethod(
  f = "rbindADEg", 
  signature = c("ADEgORADEgSORtrellis", "ADEgORADEgSORtrellis"),
  definition = function(g1, g2, ..., plot = FALSE) {
  	if(try(is.list(...), silent = TRUE) == TRUE)
      glist <- as.list(c(g1, g2, ...))
	  else
	    glist <- list(g1, g2, ...)
    
    nbg <- length(glist)
    obj <- ADEgS(adeglist = glist, layout = c(nbg, 1), add = matrix(0, ncol = nbg, nrow = nbg), plot = FALSE)
    obj@Call <- match.call()
    
    if(plot)
      print(obj)
    invisible(obj)
  })


##############################################        
##                   insertion              ##
##############################################

setMethod(
  f = "insert", 
  signature = c("ADEgS", "missing"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot, which, dispatch) {
    positions <- .getposition(posi, w = ratio, h = ratio) + inset
    currentgraphic <- get("currentadeg", envir = .ADEgEnv)
    if(!(length(currentgraphic)))
      stop("no existing graphics")
    else
      newADEgS <- insert(graphics = graphics, oldgraphics = currentgraphic, posi = posi, ratio = ratio, inset = inset, plot = plot, which = which, dispatch = dispatch)
    
    if(plot)
      print(newADEgS[length(newADEgS)], newpage = FALSE)
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })


setMethod(
  f = "insert", 
  signature = c("ADEgS", "ADEg"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot) {
    positions <- .getposition(posi, w = ratio, h = ratio) + inset
    thecall <- call("insert", graphics@Call, oldgraphics@Call)
    newADEgS <- new(Class = "ADEgS", ADEglist = list(oldgraphics, graphics), positions = rbind(c(0, 0, 1, 1), positions), add = matrix(0, ncol = 2, nrow = 2), thecall)
    if(plot)      
      print(newADEgS)
    
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })


setMethod(
  f = "insert", 
  signature = c("ADEgORtrellis", "ADEgS"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot, which) {
    thecall <- call("insert", graphics@Call, oldgraphics@Call)
    if(missing(which)) {
      positions <- .getposition(posi, w = ratio, h = ratio) + inset
      newADEgS <- new(Class = "ADEgS", ADEglist = c(oldgraphics@ADEglist, list(graphics)), positions = rbind(oldgraphics@positions, positions), add = rbind(cbind(oldgraphics@add, rep(0, length.out = nrow(oldgraphics@add))), rep(0, length.out = ncol(oldgraphics@add) + 1)), Call = thecall)
    } else {
      l <- sapply(1:length(oldgraphics), FUN = function(i) {if(i %in% which) {insert(graphics, oldgraphics@ADEglist[[i]], posi = posi, ratio = ratio, inset = inset, plot = FALSE)} else oldgraphics@ADEglist[[i]]})
      newADEgS <- new(Class = "ADEgS", ADEglist = l, positions = oldgraphics@positions, add = oldgraphics@add, Call = thecall)
    }
    
    if(plot)
      print(newADEgS)
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })


setMethod(
  f = "insert", 
  signature = c("ADEgS", "ADEgS"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot, which, dispatch) {
    thecall <- call("insert", graphics@Call, oldgraphics@Call)
    if(!dispatch){
      if(missing(which)) {
        positions <- .getposition(posi, w = ratio, h = ratio) + inset
        newADEgS <- new(Class = "ADEgS", ADEglist = c(oldgraphics@ADEglist, list(graphics)), positions = rbind(oldgraphics@positions, positions), add = rbind(cbind(oldgraphics@add, rep(0, length.out = nrow(oldgraphics@add))), rep(0, length.out = ncol(oldgraphics@add) + 1)), Call = thecall)
      } else {
        l <- sapply(1:length(oldgraphics), FUN = function(i) {if(i %in% which) {insert(graphics, oldgraphics@ADEglist[[i]], posi = posi, ratio = ratio, inset = inset, plot = FALSE)} else oldgraphics@ADEglist[[i]]})
        newADEgS <- new(Class = "ADEgS", ADEglist = l, positions = oldgraphics@positions, add = oldgraphics@add, Call = thecall)
      }
    } else {
      if(length(graphics) != length(oldgraphics))
        stop("dispatch option is not allowed with ADEgS object of different length")
      else {
        l <- sapply(1:length(oldgraphics), FUN = function(i) {insert(graphics@ADEglist[[i]], oldgraphics@ADEglist[[i]], posi = posi, ratio = ratio, inset = inset, plot = FALSE)})
        newADEgS <- new(Class = "ADEgS", ADEglist = l, positions = oldgraphics@positions, add = oldgraphics@add, Call = thecall)
      }
    }
    
    if(plot)
      print(newADEgS)
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })

##############################################          
##                  Update                  ##
##############################################

## update the modified parameters
setMethod(
  f = "update",
  signature = "ADEgS",
  definition = function(object, ..., plot = TRUE) {
    nameobj <- deparse(substitute(object, env = parent.frame()))
    ## object is in parent.frame() because 'update' method pattern is different with 'update' generic method pattern
    ## see https://stat.ethz.ch/pipermail/r-help/2008-January/152296.html
    
    slots <- list()
    slots$names <- names(object)
    slots$positions <- object@positions
    
    ## extract specific slots used in function call
    pattern <- c("names", "positions")
    lpattern <- as.list(rep("", length(pattern)))
    names(lpattern) <- pattern
    
    ## sort parameters
    sep <- separation(..., pattern = lpattern)
    slots <- modifyList(slots, sep[[1]], keep.null = TRUE)
    sep[[2]] <- sortparamADEgS(sep[[2]], graphsnames = slots$names)

    ADEglist <- sapply(1:length(object@ADEglist), FUN = function(x) {if(inherits(object@ADEglist[[x]], "ADEg") | inherits(object@ADEglist[[x]], "ADEgS")) update(object@ADEglist[[x]], plot = FALSE, sep[[2]][[x]]) else do.call("update", c(list(object = object@ADEglist[[x]]), sep[[2]][[x]]))})
    object <- new("ADEgS", ADEglist = ADEglist, positions = slots$positions, add = object@add, Call = match.call())
    names(object) <- slots$names
    
    if(plot)
      print(object)
    assign(nameobj, object, envir = parent.frame(n = 2))
    ## see also https://stat.ethz.ch/pipermail/r-help/2008-January/152296.html
    
    assign("currentadeg", object, envir = .ADEgEnv)
  })

##############################################          
##                   Display                ##
##############################################

setMethod(
  f = "show",
  signature = "ADEgS",
  definition = function(object) {
    print(object)
  })


setMethod(
  f = "plot",
  signature = c("ADEgS", "ANY"),
  definition = function(x, y) {
    print(x)
  })


setMethod(
  f = "print",
  signature = "ADEgS",
  definition = function(x, closeViewport = TRUE, square = NULL) {
    oldtextcex <- trellis.par.get("fontsize")$text
    oldpointcex <- trellis.par.get("fontsize")$points
    oldmarginH <- trellis.par.get("layout.heights")
    oldmarginW <- trellis.par.get("layout.widths")
    
    trellis.par.set(layout.heights = list(top.padding = .2 + oldmarginH$top.padding, bottom.padding = .2 + oldmarginH$bottom.padding), layout.widths = list(left.padding = .2 + oldmarginW$left.padding, right.padding = .2 + oldmarginW$right.padding))
    on.exit(trellis.par.set(list("fontsize" = list("text" = oldtextcex, "points" = oldpointcex), "layout.widths" = list("left.padding" = oldmarginW$left.padding, "right.padding" = oldmarginW$right.padding), "layout.heights" = list("top.padding" = oldmarginH$top.padding, "bottom.padding" = oldmarginH$bottom.padding))))
    
    gettextsize <- function(widG, heigG) {
      ## Adjust text size to viewport size
      if(widG < 1 / 2 || heigG < 1 / 2)
        return(0.66 / 1.25)
      if(widG == 1 / 2 && heigG == 1 / 2)
        return(0.83 / 1.25)
      if(widG == 1 && heigG == 1)
        return(1)
      else return(1 / 1.25)              
    }
    
    getxscale <- function(object) {
      ## Obtain limits for x
      if(inherits(object, "ADEg"))
        object <- gettrellis(object)
      if(class(object) == "trellis") {
        res <- object$x.limits
      } else {
        res <- c(0, 1)
      }
      return(res)
    }
    
    getyscale <- function(object) {
      ## Obtain limits for y
      if(inherits(object, "ADEg"))
        object <- gettrellis(object)
      if(class(object) == "trellis") {
        res <- object$y.limits
      } else {
        res <- c(0, 1)
      }
      return(res)
    }
    
    printADEGs <- function(adegobject, closeViewport, square) {
      if(closeViewport)
        grid.newpage()
      
      positions <- adegobject@positions
      listG <- adegobject@ADEglist
      
      ## create the list of viewport and push it
      unit.vpL <- "npc"
      if(isTRUE(square))
        unit.vpL <- "snpc"
      
      vpL <- do.call("vpList", lapply(1:length(listG), function(i) do.call("viewport", args = list(x = positions[i, 1], y = positions[i, 2], width = positions[i, 3] - positions[i, 1], height = positions[i, 4] - positions[i, 2], just = c(0, 0), name = names(listG)[i], xscale = getxscale(listG[[i]]), yscale = getyscale(listG[[i]]), default.units = unit.vpL))))
      pushViewport(vpL)
      
      upViewport(0)
      width.root <- convertWidth(unit(1, unit.vpL), "inches", valueOnly = TRUE)
      height.root <- convertHeight(unit(1, unit.vpL), "inches", valueOnly = TRUE)
      
      for(i in 1:length(listG)) {
        object <- listG[[i]]
        seekViewport(names(listG)[i])
        
        if(inherits(object, "ADEg") | class(object) == "trellis") {
          if(inherits(object, "ADEg"))
            trobject <- gettrellis(object)
          else
            trobject <- object
          
          square.i <- ifelse(is.null(square), !trobject$aspect.fill, square)
          unit.vpi <- "npc"
          if(isTRUE(square.i))
            unit.vpi <- "snpc"
          
          vp <- viewport(x = 0, y = 0, width = 1, height = 1, just = c(0, 0), name = "current", xscale = getxscale(listG[[i]]), yscale = getyscale(listG[[i]]), default.units = unit.vpi)
          pushViewport(vp)

          width.current <- convertWidth(unit(1, unit.vpi), "inches", valueOnly = TRUE)
          height.current <- convertHeight(unit(1, unit.vpi), "inches", valueOnly = TRUE)
          ratio.width <- width.current / width.root
          ratio.height <- height.current / height.root
          
          cst <- gettextsize(ratio.width, ratio.height)
          sup <- adegobject@add[, i]
          trellis.par.set(list("fontsize" = list("text" = oldtextcex * cst, "points" = oldpointcex * cst)))
          if(any(sup == 1))
            printSuperpose(g1 = object, refg = listG[[which(adegobject@add[, i] == 1)[1]]])
          else
            print(object, newpage = FALSE)
          
          popViewport()
          
        } else if(inherits(object, "ADEgS")) {
          names(object) <- paste(names(listG)[i], names(object), sep = ".")
          printADEGs(object, closeViewport = FALSE, square = square)
        } else {
          stop(paste("Not implemented for class:", class(object), sep = " "))
        }
        popViewport()
      }
    }
    
    printADEGs(x, closeViewport = closeViewport, square = square)
    assign("currentadeg", x, envir = .ADEgEnv)
  })


##############################################          
##                   Creation               ##
##############################################

ADEgS <- function(adeglist, positions, layout, add = NULL, plot = TRUE) {
  
  m <- matrix(0, length(adeglist), length(adeglist))
  
  if(missing(layout) & (is.null(add) | identical(add, m)) & missing(positions))
    layout <- .n2mfrow(length(adeglist))
  
  if(missing(positions) & !missing(layout)) {
    if(is.list(layout)) ## in layout: width and heights informations, layout is a list
      positions <- do.call("layout2position", layout)
    else
      positions <- layout2position(layout, ng = length(adeglist))
  }
  
  if(missing(positions)) 
    positions <- matrix(rep(c(0, 0, 1, 1), length.out = length(adeglist) * 4), byrow = TRUE, ncol = 4)
  
  if(is.null(add))
    add <- m
  
  ADEgObject <- new(Class = "ADEgS", ADEglist = adeglist, positions = positions, add = add, Call = match.call())
  if(plot)
    print(ADEgObject)
  invisible(ADEgObject)
}
