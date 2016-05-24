##############################################
##               general class              ##
##############################################

setClass(
  Class = "ADEg",
  contains = "VIRTUAL",
  slots = c(
    trellis.par = "list",
    adeg.par = "list",
    lattice.call = "list",
    g.args = "list",
    stats = "list",
    s.misc = "list",
    Call = "call")
	## validity
	## validity = function(object){return(TRUE)}
)

##############################################          
##                 initialize               ##
##############################################       

setMethod(
  f = "initialize",
  signature = "ADEg",
  definition = function(.Object, trellis.par = list(), adeg.par = list(), lattice.call = list(), g.args = list(), stats = list(), s.misc = list(), Call = call("emptycall"), ...) {
    .Object@trellis.par <- trellis.par
    .Object@adeg.par <- adeg.par
    .Object@lattice.call <- lattice.call
    .Object@g.args <- g.args
    .Object@stats <- stats
    .Object@s.misc <- s.misc
    .Object@Call <- Call
    return(.Object)
  })


setOldClass("trellis")
setClassUnion(name = "ADEgORtrellis", members = c("ADEg", "trellis"))


setMethod(
  f = "panelbase",
  signature = "ADEg",
  definition = function(object, x, y) {
    
    whichpanel <- packet.number()
    sub <- lapply(object@adeg.par$psub, FUN = function(x) {rep(x, length.out = max(whichpanel, 2))}) ## repeat at least twice for position coordinates
    if(is.numeric(object@adeg.par$psub$position)) 
      num <- 1
    else 
      num <- 0
    
    if(sub$cex[whichpanel] && (sub$text[whichpanel] != "")) {
      if(!num)
        posit <- .setposition(sub$position[whichpanel])     
      else
        posit <- .setposition(sub$position[(2 * whichpanel - 1):(2 * whichpanel)])
      text <- textGrob(label = sub$text[whichpanel], x = posit$posi[1], y = posit$posi[2], gp = gpar(cex = sub$cex[whichpanel], col = sub$col[whichpanel]), just = posit$just, name = paste("subtitle_", whichpanel, sep = ""))            
      grid.rect(x = posit$posi[1], y = posit$posi[2], width = grobWidth(text), height = grobHeight(text), just = posit$just, gp = gpar(fill =  ifelse(class(object) == "S2.corcircle" | inherits(object, "ADEg.Tr"), "transparent", object@adeg.par$pbackground$col), alpha = 1, col = "transparent"))
      grid.draw(text)
    }
  })

##############################################          
##      Get elements/information            ##
##############################################

setMethod(
  f = "getparameters",
  signature = "ADEg",
  definition = function(object, number = 0) {	
    if(number == 0)
      return(list(trellis.par = object@trellis.par, adeg.par = object@adeg.par, g.args = object@g.args))
    if(number == 1)
      return(object@trellis.par)
    if(number == 2)
      return(object@adeg.par)
    stop("wrong number for getparameters")
  })


setMethod(
  f = "getlatticecall",
  signature = "ADEg",
  definition = function(object) {
    return(object@lattice.call)
  })


setMethod(
  f = "getcall",
  signature = "ADEg",
  definition = function(object) {
    return(object@Call)
  })


setMethod(
  f = "getstats",
  signature = "ADEg",
  definition = function(object) {
    return(object@stats)
  })


##############################################          
##               superposition              ##
##############################################

## g1 superpose on refg
## settings about margin limits ect... taken from refg
## modified object only displayed (not save), original limits ect are kept 
setMethod(
  f = "printSuperpose",
  signature = c("ADEgORtrellis", "ADEgORtrellis"),
  definition = function(g1, refg) {
    ## to respect axis, limits, etc., we work directly on the trellis objects.
    if(inherits(refg, "ADEg")) {
      trelref <- gettrellis(refg)
      if(inherits(g1, "ADEg")) {
      	g1@adeg.par$pgrid$draw <- FALSE
      	g1@g.args$xlim <- refg@g.args$xlim
      	g1@g.args$ylim <- refg@g.args$ylim
      	g1@adeg.par$paxes$draw <- refg@adeg.par$paxes$draw
        g1@adeg.par$pbackground$col <- "transparent"	## useful for S2.corcircle
      	g1@adeg.par$porigin$draw <- FALSE
      	g1@s.misc$scales <- refg@s.misc$scales
        
      	if(inherits(g1, "ADEg.Tr") & inherits(refg, "ADEg.Tr")) {
  	      g1@g.args$min3d <- refg@g.args$min3d
        	g1@g.args$max3d <- refg@g.args$max3d
        	g1@adeg.par$pgrid$text$cex <- 0 ## no text corner for g1
        	g1@lattice.call$arguments$par.settings$axis.text$cex <- 0
      	}
      	setlatticecall(g1)
      	trel1 <- gettrellis(g1)
      } else {
      	trel1 <- g1
      }
    } else {  ## refg is a trellis
      trelref <- refg
      if(inherits(g1, "ADEg")) {
        g1@adeg.par$pgrid$draw <- FALSE
      	g1@g.args$xlim <- refg$x.limits
      	g1@g.args$ylim <- refg$y.limits
      	g1@adeg.par$paxes$draw <- refg$x.scales$draw * refg$y.scales$draw
      	g1@adeg.par$porigin$draw <- FALSE
        g1@g.args$xlab <- ""
        g1@g.args$ylab <- ""
        
      	setlatticecall(g1)
      	trel1 <- gettrellis(g1)
      } else {
        trel1 <- g1
      }
    }
    
    trel1$par.settings$panel.background$col <- "transparent"
    trel1$par.settings$axis.text$alpha <- 0
    trel1$par.settings$axis.line$col <- "transparent"
    
    names <- c("x.scales", "y.scales", "xlab", "ylab", "main", "sub", "x.between", "y.between", "as.table", "x.limits", "y.limits", "aspect.ratio")            
    for(i in names)
      trel1[[i]] <- trelref[[i]]

    print(trel1, newpage = FALSE)
  })

  
setMethod(
  f = "superpose",
  signature = c("ADEgORtrellis", "ADEgORtrellis", "ANY", "ANY"),
  definition = function(g1, g2, which, plot) {
    addi <- matrix(0, 2, 2)
    addi[1, 2] <- 1
    obj <- new(Class = "ADEgS", ADEglist = list(g1, g2), positions = matrix(rep(c(0, 1), each = 4), 2, 4), add = addi, Call = match.call())
    if(plot)
      print(obj)
    invisible(obj)
  })

  
##############################################          
##                   insertion              ##
##############################################

setMethod(
  f = "insert", 
  signature = c("ADEgORtrellis", "missing"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot, which) {
    positions <- .getposition(posi, w = ratio, h = ratio) + inset
    currentgraphic <- get("currentadeg", envir = .ADEgEnv)
    if(!(length(currentgraphic)))
      stop("no existing graphics")
    else
      newADEgS <- insert(graphics = graphics, oldgraphics = currentgraphic, posi = posi, ratio = ratio, inset = inset, plot = plot, which = which)
    
    if(plot)
      print(newADEgS[length(newADEgS)], closeViewport = FALSE)
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })


setMethod(
  f = "insert", 
  signature = c("ADEgORtrellis", "ADEg"),
  definition = function(graphics, oldgraphics, posi, ratio, inset, plot) {
    positions <- .getposition(posi, w = ratio, h = ratio) + inset
    thecall <- call("insert", graphics@Call, oldgraphics@Call)
    newADEgS <- new(Class = "ADEgS", ADEglist = list(oldgraphics, graphics), positions = rbind(c(0, 0, 1, 1), positions), add = matrix(0, ncol = 2, nrow = 2), Call = thecall)
    if(plot)     
      print(newADEgS)
    assign("currentadeg", newADEgS, envir = .ADEgEnv)
    invisible(newADEgS)
  })

  
##############################################          
##                   Add                    ##
##############################################

setMethod(
  f = "+",
  signature = c("ADEg", "ADEg"),
  definition = function(e1, e2) {
    newobj <- superpose(e1, e2)
    newobj@Call <- match.call()
    return(newobj)
  })


setMethod(
  f = "add.ADEg",
  signature = c("ADEg"),
  definition = function(object) {
    previous <- get("currentadeg", envir = .ADEgEnv)
    if(!(length(previous))) 
      stop("no graph to add to")
    
    objects <- superpose(previous, object)
    if(inherits(previous, "ADEg"))              
      printSuperpose(object, previous, position = c(0, 0, 1, 1))
    
    else if(inherits(previous, "ADEgS"))
      printSuperpose(object, previous[[length(previous)]], position = previous@positions[length(previous), ])
#    lattice:::lattice.setStatus(print.more = FALSE)
    assign("currentadeg", objects, envir = .ADEgEnv)
    invisible(objects)
  })

##############################################          
##                  Update                  ##
##############################################

## update the modified parameters
setMethod(
  f = "update",
  signature = "ADEg",
  definition = function(object, ..., plot = TRUE) {
    nameobj <- deparse(substitute(object, env = parent.frame(n = 1)))
    ## object is in parent.frame() because 'update' method pattern is different with 'update' generic method pattern
    ## see https://stat.ethz.ch/pipermail/r-help/2008-January/152296.html
    
    ## extract specific slots used in function call
    pattern <- names(object@g.args)
    lpattern <- as.list(rep("", length(pattern)))
    names(lpattern) <- pattern
    
    ## sort parameters
    sep <- separation(..., pattern = lpattern)
    selection <- sortparamADEg(sep[[2]])
    selection$g.args <- c(selection$g.args, sep[[1]])
    
    if(length(selection$rest))
      warning(c("Unused parameters: ", paste(unique(names(selection$rest)), " ", sep = "")), call. = FALSE)
    
    object@adeg.par <- modifyList(object@adeg.par, selection$adepar, keep.null = TRUE)
    object@trellis.par <- modifyList(object@trellis.par, selection$trellis, keep.null = TRUE)
    object@g.args <- modifyList(object@g.args, selection$g.args, keep.null = TRUE)
    
    prepare(object)
    setlatticecall(object)
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
  signature = "ADEg",
  definition = function(object) {
    print(object)
  })


setMethod(
  f = "plot",
  signature = c("ADEg", "ANY"),
  definition = function(x, y, adjust = FALSE) {
    print(x, adjust = adjust)
  })


setMethod(
  f = "print",
  signature = c("ADEg"),
  definition = function(x, adjust = FALSE, newpage = TRUE) {
		## if adjust, graphic limits are readjust according to the device size.
		## for now it is only available if only an ADEg is drawn (not in ADEgS)
    if(adjust) {
      aspp <- dev.size()  ## device size (in inches)
      ratid <- aspp[1] / aspp[2]
      oxlim <- x@lattice.call$arguments$xlim  ## old xlim
      oylim <- x@lattice.call$arguments$ylim
      ratig <- diff(oxlim) / diff(oylim)
      ## if not mandatory ...?
      if((ratid / ratig) > 1) { 
        ## width device bigger (relative) than width graphic        
        centerx <- oxlim[1] + diff(oxlim) / 2
        nxlim <- rep(centerx, 2) + c(-1, 1) * ((ratid * diff(oylim)) / 2)
        nylim <- oylim
      }
      else if((ratid / ratig) < 1) {  ## then relative device height bigger than relative graphic height
        centery <- oylim[1] + diff(oylim) / 2
        nylim <- rep(centery, 2) + c(-1, 1) * (1 / ratid * diff(oxlim)) / 2
        nxlim <- oxlim
      }
      x@s.misc$backgrid <- .getgrid(xlim = nxlim, ylim = nylim, x@adeg.par$pgrid$nint, rep(x@adeg.par$porigin$origin, le = 2), asp = x@adeg.par$paxes$aspectratio)
      setlatticecall(x) ## passing backgrid
      ## changing limits
      x@lattice.call$arguments$xlim <- nxlim
      x@lattice.call$arguments$ylim <- nylim
    }
    object <- x
    
    if(!length(object@lattice.call))
      stop("no graphics instruction")
    else { 
      tmp_object <- gettrellis(x)
      print(tmp_object, newpage = newpage)
      assign("currentadeg", x, envir = .ADEgEnv)
    }
  })

