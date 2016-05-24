

#######################
# "animValue" stuff
#######################

# An animValue is a vector PLUS a timeid PLUS an id
animValue <- function(x, timeid=NULL, id=NULL) {
    if (!is.atomic(x))
        stop("'x' must be a atomic")
    if (!is.null(timeid))
        timeid <- rep(timeid, length.out=length(x))
    if (!is.null(id))
        id <- rep(id, length.out=length(x))
    tu <- list(values=x, timeid=timeid, id=id)
    class(tu) <- "animValue"
    tu
}

is.animValue <- function(x) inherits(x, "animValue")

as.animValue <- function(x, ...) {
    UseMethod("as.animValue")
}

as.animValue.animValue <- function(x, ...) x

as.animValue.numeric <- function(x, ...) {
    animValue(x)
}

as.animValue.character <- function(x, ...) {
    animValue(x)
}

# 'multVal' controls whether columns of the matrix are used as
# 'timeid' or 'id'
as.animValue.matrix <- function(x, multVal=FALSE, ...) {
    if (multVal) {
        animValue(x, timeid=rep(1:ncol(x), each=nrow(x)))
    } else {
        animValue(x, id=rep(1:ncol(x), each=nrow(x)))
    }
}

as.animValue.list<- function(x, multVal=FALSE, ...) {
    if (!all(sapply(x, is.atomic)))
        stop("All components of list must be atomic")
    if (multVal) {
        animValue(unlist(x),
                  timeid=rep(1:length(x), sapply(x, length)))
    } else {
        animValue(unlist(x),
                  id=rep(1:length(x), sapply(x, length)))
    }
}

listFromAnimValue <- function(x) {
    if (is.null(x$id)) {
        if (is.null(x$timeid)) {
            n <- length(x$values)
            animValueList <- as.list(x$values)
        } else {
            times <- unique(x$timeid)
            n <- length(times)
            animValueList <- split(x$values, x$timeid)
        }
        names(animValueList) <- paste("t", 1:n, sep="")
    } else {
        shapes <- unique(x$id)
        ns <- length(shapes)
        animValueList <- vector("list", ns)
        for (i in 1:ns) {
            animValueList[[i]] <-
                listFromAnimValue(animValue(x$values[x$id == i],
                                            x$timeid[x$id == i]))
        }
        names(animValueList) <- paste("id", 1:ns, sep="")
    }
    animValueList
}

print.animValue <- function(x, ...) {
    # Generate list from animValue and then print the list
    print(listFromAnimValue(x))
}

#######################
# "animUnit" stuff

# An animUnit is a unit PLUS a timeid PLUS an id
# The timeid breaks the values in the unit into different time
# periods, and the id breaks the values into different shapes

animUnit <- function(x, timeid=NULL, id=NULL) {
    if (!is.unit(x))
        stop("'x' must be a unit object")
    if (!is.null(timeid))
        timeid <- rep(timeid, length.out=length(x))
    if (!is.null(id))
        id <- rep(id, length.out=length(x))
    tu <- list(values=x, timeid=timeid, id=id)
    class(tu) <- "animUnit"
    tu
}

is.animUnit <- function(x) inherits(x, "animUnit")

as.animUnit <- function(x, ...) {
    UseMethod("as.animUnit")
}

as.animUnit.animUnit <- function(x, ...) x

as.animUnit.numeric <- function(x, unit=NULL, ...) {
    if (is.null(unit))
        stop("Require 'unit' to convert numeric vector")
    animUnit(unit(x, unit))
}

as.animUnit.unit <- function(x, ...) {
    animUnit(x)
}

# 'multVal' controls whether columns of the matrix are used as
# 'timeid' or 'id'
as.animUnit.matrix <- function(x, unit=NULL, multVal=FALSE, ...) {
    if (is.null(unit))
        stop("Require 'unit' to convert matrix")
    if (multVal) {
        animUnit(unit(x, unit), timeid=rep(1:ncol(x), each=nrow(x)))
    } else {
        animUnit(unit(x, unit), id=rep(1:ncol(x), each=nrow(x)))
    }
}

as.animUnit.list<- function(x, multVal=FALSE, ...) {
    if (!all(sapply(x, is.unit)))
        stop("All components of list must be units")
    if (multVal) {
        animUnit(do.call("unit.c", x),
                 timeid=rep(1:length(x), sapply(x, length)))
    } else {
        animUnit(do.call("unit.c", x),
                 id=rep(1:length(x), sapply(x, length)))
    }
}

listFromAnimUnit <- function(x) {
    if (is.null(x$id)) {
        if (is.null(x$timeid)) {
            n <- length(x$values)
            animUnitList <- vector("list", n)
            for (i in 1:n)
                animUnitList[[i]] <- x$values[i]
        } else {
            times <- unique(x$timeid)
            n <- length(times)
            animUnitList <- vector("list", n)
            for (i in 1:n)
                animUnitList[[i]] <- x$values[x$timeid == i]
        }
        names(animUnitList) <- paste("t", 1:n, sep="")
    } else {
        shapes <- unique(x$id)
        ns <- length(shapes)
        animUnitList <- vector("list", ns)
        for (i in 1:ns) {
            animUnitList[[i]] <-
                listFromAnimUnit(animUnit(x$values[x$id == i],
                                          x$timeid[x$id == i]))
        }
        names(animUnitList) <- paste("id", 1:ns, sep="")
    }
    animUnitList
}

print.animUnit <- function(x, ...) {
    # Generate list from animUnit and then print the list
    print(listFromAnimUnit(x))
}

# duration says how many SECONDS the animation lasts for
# id indicates the identity of multiple animated values
#   (i.e., allows a vector of animated values)
#   If "auto" then it depends on the number and size
#     of the elements being animated.  If there is
#     only one element, it is NULL.
# rep says how many times to repeat the animation
#   (TRUE means indefinitely;  FALSE means once)
# revert says whether to revert to the start value of the
#   animation upon completion

autoid <- function(id) {
    if (!is.numeric(id))
        if (id == "auto")
            TRUE
        else
            stop("Invalid id")
    else
        FALSE
}

animationSet <- function(...,
                         duration=1, rep=FALSE, revert=FALSE,
                         begin=0, interp="linear") {
    animations <- list(...)
    if (is.null(animations[[1]]))
        stop("need argument to animate")
    list(animations=animations,
         begin=begin, interp=interp, duration=duration, rep=rep, revert=revert)
}

animateGrob <- function(grob, ...,
                        duration=1, 
                        rep=FALSE, revert=FALSE,
                        begin=0, interpolate="linear", group=FALSE) {
    if (!interpolate %in% c("linear", "discrete"))
        stop("Invalid interpolation method")
    as <- animationSet(...,
                       duration=duration, rep=rep, revert=revert,
                       begin=begin, interp=interpolate)
    cl <- class(grob)
    if (group) {
        grob$groupAnimationSets <- c(grob$groupAnimationSets, list(as))
    } else {
        grob$animationSets <- c(grob$animationSets, list(as))
    }
    class(grob) <- unique(c("animated.grob", cl))
    grob
}
  
grid.animate <- function(path, ..., group=FALSE, redraw = FALSE,
                         strict=FALSE, grep=FALSE, global=FALSE) {
    grobApply(path,
              function(path) {
                  grid.set(path, animateGrob(grid.get(path), ...,
                                             group=group),
                           redraw = redraw)
              }, strict = strict, grep = grep, global = global)
    invisible()
}

applyAnimation <- function(x, ...) {
  UseMethod("applyAnimation")
}

# Convert to animValue then take value(s) for shape "i"
# This function is designed for animValues where timeid is NULL
# so that each time period has only ONE value
# RETURN a VECTOR
ithValue <- function(animValues, i) {
    av <- as.animValue(animValues)
    if (!is.null(av$timeid))
        stop("Expecting only one value per time point")
    if (is.null(av$id))
        av$values
    else
        av$values[av$id == i]
}

# Convert to animUnit then take unit(s) for shape "i"
# This function is designed for animUnits where timeid is NULL
# so that each time period has only ONE value
# RETURN a UNIT
ithUnit <- function(animValues, origValue, i) {
    au <- as.animUnit(animValues,
                      # Only take the first "unit" value
                      unit=attr(origValue, "unit")[1])
    if (!is.null(au$timeid))
        stop("Expecting only one value per time point")
    if (is.null(au$id))
        au$values
    else
        au$values[au$id == i]
}

# Convert to animValue then take value(s) for shape "i"
# This function is designed for animValues where there is a timeid 
# so that each time period has MULTIPLE values
# RETURN an ANIMVALUE
ithAnimValue <- function(animValues, i) {
    av <- as.animValue(animValues, multVal=TRUE)
    if (is.null(av$timeid))
        stop("Expecting multiple values per time point")
    if (is.null(av$id))
        av
    else
        animValue(av$values[av$id == i],
                  av$timeid[av$id == i])
}

# Convert to animUnit then take unit(s) for shape "i"
# This function is designed for animUnit where there is a timeid 
# so that each time period has MULTIPLE values
# RETURN an ANIMUNIT
ithAnimUnit <- function(animValues, origValue, i) {
    au <- as.animUnit(animValues,
                      # Only take the first "unit" value
                      unit=attr(origValue, "unit")[1],
                      multVal=TRUE)
    if (is.null(au$timeid))
        stop("Expecting multiple values per time point")
    if (is.null(au$id))
        au
    else
        animUnit(au$values[au$id == i],
                 au$timeid[au$id == i])
}

######################
# applyAnimation methods
#
# There is one of these for each primitive, but they all have similar
# structure:

# if (group)
#     animate the <g> element

# else

#     some sets of values (e.g., x/y) are animated together
#     so bail out if this combination has already been animated

#     recycle animation values to full length

#     for each shape ...

#         select anim values
#         animate sets of values
#         animate anything else

######################

######################
# FIXME:
# When animating some cobination x/y/width/height/size AT THE SAME TIME
# the code below only makes sense if the number of time periods
# is the same for all of x/y/width/height/size (that are being animated)
######################

applyAnimation.rect <- function(x, animSet, animation, group, dev) {

    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
  # We may be dealing with multiple rects that need animating
  n <- max(length(x$x), length(x$y), length(x$width), length(x$height))

  # Rep the original x/y/width/height out to be the same length
  x$x <- rep(x$x, length.out=n)
  x$y <- rep(x$y, length.out=n)
  x$width <- rep(x$width, length.out=n)
  x$height <- rep(x$height, length.out=n)
  
  # Repeating animation parameters so that each element can have
  # distinct values
  begin <- rep(animSet$begin, length.out = n)
  interp <- rep(animSet$interp, length.out = n)
  dur <- rep(animSet$duration, length.out = n)
  rep <- rep(animSet$rep, length.out = n)
  rev <- rep(animSet$revert, length.out = n)

  angle <- current.rotation()
  
  for (i in 1:n) {
    subName <- subGrobName(x$name, i)

    # If x AND y change, need to transform together
    # If width/height changes, have to animate x/y as well
    # because SVG <rect> does not have justification
    if ("x" %in% names(animSet$animations))
        xi <- ithUnit(animSet$animations$x, x$x, i)
    else
        xi <- x$x[i]
    if ("y" %in% names(animSet$animations))
        yi <- ithUnit(animSet$animations$y, x$y, i)
    else
        yi <- x$y[i]
    if ("width" %in% names(animSet$animations))
        widthi <- ithUnit(animSet$animations$width, x$width, i)
    else
        widthi <- x$width[i]
    if ("height" %in% names(animSet$animations))
        heighti <- ithUnit(animSet$animations$height, x$height, i)
    else
        heighti <- x$height[i]
    lb <- leftbottom(xi, yi, widthi, heighti, x$just, x$hjust, x$vjust, dev)
    
    switch(animation,
           x={
               svgAnimateXYWH("x", cx(lb$x, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i],
                              subName, dev@dev)
               if (angle != 0) {
                   if (!("y" %in% names(animSet$animations))) {
                       svgAnimateXYWH("y", cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      subName, dev@dev)
                   }
                   svgAnimateRotation(angle, cx(lb$x, dev), cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      id=subName, svgdev=dev@dev)
               }
           },
           y={
               svgAnimateXYWH("y", cy(lb$y, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i],
                              subName, dev@dev)
               if (angle != 0) {
                   if (!("x" %in% names(animSet$animations))) {
                       svgAnimateXYWH("x", cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      subName, dev@dev)
                   }
                   svgAnimateRotation(angle, cx(lb$x, dev), cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      id=subName, svgdev=dev@dev)
               }
           },
           width={
               # If x is also animated, this has already been handled above
               if (!("x" %in% names(animSet$animations))) {
                   svgAnimateXYWH("x", cx(lb$x, dev),
                                  begin[i], interp[i], dur[i], rep[i], rev[i],
                                  subName, dev@dev)
               } 
               if (angle != 0) {
                   if (!("y" %in% names(animSet$animations))) {
                       svgAnimateXYWH("y", cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      subName, dev@dev)
                   }
                   svgAnimateRotation(angle, cx(lb$x, dev), cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      id=subName, svgdev=dev@dev)
               }
               dim <- dimToInches(ithUnit(animSet$animations$width,
                                          x$width, i),
                                  x$height[i], dev)
               svgAnimateXYWH("width", cw(dim$w, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i],
                              subName, dev@dev)
           },
           height={
               if (!("y" %in% names(animSet$animations))) {
                   svgAnimateXYWH("y", cy(lb$y, dev),
                                  begin[i], interp[i], dur[i], rep[i], rev[i],
                                  subName, dev@dev)
               } 
               if (angle != 0) {
                   if (!("x" %in% names(animSet$animations))) {
                       svgAnimateXYWH("x", cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      subName, dev@dev)
                   }
                   svgAnimateRotation(angle, cx(lb$x, dev), cy(lb$y, dev),
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      id=subName, svgdev=dev@dev)
               }
               dim <- dimToInches(x$width[i],
                                  ithUnit(animSet$animations$height,
                                          x$height, i),
                                  dev)
               svgAnimateXYWH("height", ch(dim$h, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i],
                              subName, dev@dev)
           },
           # Any other attribute
           {
               svgAnimate(animation,
                          paste(ithValue(animSet$animations[[animation]], i),
                                collapse=";"),
                          begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
           })
  }
    }
}

applyAnimation.circle <- function(x, animSet, animation, group, dev) {

    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
  # We may be dealing with multiple circles that need animating
  n <- max(length(x$x), length(x$y), length(x$r))

  # Rep the original x/y/width/height out to be the same length
  x$x <- rep(x$x, length.out=n)
  x$y <- rep(x$y, length.out=n)
  x$r <- rep(x$r, length.out=n)
  
  # Repeating animation parameters so that each element can have
  # distinct values
  begin <- rep(animSet$begin, length.out = n)
  interp <- rep(animSet$interp, length.out = n)
  dur <- rep(animSet$duration, length.out = n)
  rep <- rep(animSet$rep, length.out = n)
  rev <- rep(animSet$revert, length.out = n)

  # Because grobs can produce multiple elements, if animation is to
  # occur on a grob it is assumed to occur on all elements, but
  # elements may simply have their properties assigned to the same
  # value multiple times.
  #
  # Also note that when casting to a matrix, units lose their "unit"
  # attribute, we have to set this to the same unit as the grob
  # attribute that is being animated, for this reason, attributes should
  # be in the same unit prior to calling grid.animate()
  for (i in 1:n) {
    subName <- subGrobName(x$name, i)

    if ("x" %in% names(animSet$animations))
        xi <- ithUnit(animSet$animations$x, x$x, i)
    else
        xi <- x$x[i]
    if ("y" %in% names(animSet$animations))
        yi <- ithUnit(animSet$animations$y, x$y, i)
    else
        yi <- x$y[i]

    switch(animation,
           x={
               loc <- locToInches(xi, yi, dev)
               svgAnimateXYWH("cx", cx(loc$x, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
           },
           y={
               loc <- locToInches(xi, yi, dev)
               svgAnimateXYWH("cy", cy(loc$y, dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
           },
           r={
               svgAnimateXYWH("r", cd(ithUnit(animSet$animations$r, x$r, i), dev),
                              begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
           },
           # Any other attribute
           {
               svgAnimate(animation,
                          paste(ithValue(animSet$animations[[animation]], i),
                                collapse=";"),
                          begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
           })
  }
     }
}

applyAnimation.points <- function(x, animSet, animation, group, dev) {

    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
  # We may be dealing with multiple points that need animating
  n <- max(length(x$x), length(x$y), length(x$size))

  # Rep the original x/y/width/height out to be the same length
  x$x <- rep(x$x, length.out=n)
  x$y <- rep(x$y, length.out=n)
  x$pch <- rep(x$pch, length.out = n)
  x$size <- rep(x$size, length.out = n)

  # Need to grab the lwd so that we can keep line thickness the same
  # as we change the size of a point
  sw <- if (! is.null(x$gp$lwd)) x$gp$lwd
        else get.gpar()$lwd
  sw <- rep(as.numeric(devLwdToSVG(sw, dev)), length.out = n)

  # Repeating animation parameters so that each element can have
  # distinct values
  begin <- rep(animSet$begin, length.out = n)
  interp <- rep(animSet$interp, length.out = n)
  dur <- rep(animSet$duration, length.out = n)
  rep <- rep(animSet$rep, length.out = n)
  rev <- rep(animSet$revert, length.out = n)

  # Because grobs can produce multiple elements, if animation is to
  # occur on a grob it is assumed to occur on all elements, but
  # elements may simply have their properties assigned to the same
  # value multiple times.
  for (i in 1:n) {
      subName <- subGrobName(x$name, i)
      
      if ("x" %in% names(animSet$animations))
          xi <- ithUnit(animSet$animations$x, x$x, i)
      else
          xi <- x$x[i]
      if ("y" %in% names(animSet$animations))
          yi <- ithUnit(animSet$animations$y, x$y, i)
      else
          yi <- x$y[i]
      if ("size" %in% names(animSet$animations))
          pointsize <- ithUnit(animSet$animations$size, x$size, i)
      else
          pointsize <- x$size[i]

      # Enforce gp$cex or gp$fontsize
      pointsize <- adjustSymbolSize(pointsize, x$gp)
      
      angle <- current.rotation()
      
      switch(animation,
             x={
                 loc <- locToInches(xi, yi, dev)
                 svgAnimateXYWH("x", cx(loc$x, dev),
                                begin[i], interp[i], dur[i],
                                rep[i], rev[i],
                                id=subName, svgdev=dev@dev)
                 if (angle != 0) {
                     if (!"y" %in% names(animSet$animations)) {
                         svgAnimateXYWH("y", cy(loc$y, dev),
                                        begin[i], interp[i], dur[i],
                                        rep[i], rev[i], 
                                        id=subName, svgdev=dev@dev)
                     }
                     svgAnimateRotation(angle, cx(loc$x, dev), cy(loc$y, dev),
                                        begin[i], interp[i], dur[i],
                                        rep[i], rev[i],
                                        id=subName, svgdev=dev@dev)
                     if (!"size" %in% names(animSet$animations)) {
                         svgAnimateTranslation(-cd(pointsize, dev)/2,
                                               -cd(pointsize, dev)/2,
                                               begin[i], interp[i], dur[i],
                                               rep[i], rev[i],
                                               additive="sum",
                                               id=subName, svgdev=dev@dev)
                     }
                 }
             },
             y={
                 loc <- locToInches(xi, yi, dev)
                 svgAnimateXYWH("y", cy(loc$y, dev),
                                begin[i], interp[i], dur[i],
                                rep[i], rev[i], subName, dev@dev)
                 if (angle != 0) {
                     if (!"x" %in% names(animSet$animations)) {
                         svgAnimateXYWH("x", cx(loc$x, dev),
                                        begin[i], interp[i], dur[i],
                                        rep[i], rev[i], 
                                        id=subName, svgdev=dev@dev)
                         svgAnimateRotation(angle,
                                            cx(loc$x, dev), cy(loc$y, dev),
                                            begin[i], interp[i], dur[i],
                                            rep[i], rev[i],
                                            id=subName, svgdev=dev@dev)
                         if (!"size" %in% names(animSet$animations)) {
                             svgAnimateTranslation(-cd(pointsize, dev)/2,
                                                   -cd(pointsize, dev)/2,
                                                   begin[i], interp[i], dur[i],
                                                   rep[i], rev[i],
                                                   additive="sum",
                                                   id=subName, svgdev=dev@dev)
                         }
                     }
                 }
             },
             size={
                 pchi <- x$pch[i]
                 docharanim <- (is.character(pchi) && pchi != ".") ||
                               (is.numeric(pchi) && (pchi >= 32 && pchi != 46))
                 donumanim <- is.numeric(pchi) && pchi <= 25

                 # If we don't have a good pch, don't bother
                 if (! any(c(donumanim, docharanim)))
                     return()

                 dimsize <- cd(pointsize, dev)
                 svgAnimateXYWH("width", dimsize,
                                begin[i], interp[i], dur[i],
                                rep[i], rev[i],
                                id=subName, svgdev=dev@dev)
                 svgAnimateXYWH("height", cd(pointsize, dev),
                                begin[i], interp[i], dur[i],
                                rep[i], rev[i],
                                id=subName, svgdev=dev@dev)
                 # Centering the point
                 trdimsize <- -dimsize / 2
                 additive <- "replace"
                 if (angle != 0) {
                     loc <- locToInches(xi, yi, dev)
                     svgAnimateRotation(angle,
                                        cx(loc$x, dev), cy(loc$y, dev),
                                        begin[i], interp[i], dur[i],
                                        rep[i], rev[i],
                                        id=subName, svgdev=dev@dev)
                     additive <- "sum"
                 }
                 svgAnimateTranslation(trdimsize, trdimsize,
                                       begin[i], interp[i], dur[i],
                                       rep[i], rev[i],
                                       additive,
                                       id=subName, svgdev=dev@dev)
                 # Ensuring that stroke-width stays the same.
                 # Only do this with low numeric pchs because 
                 if (donumanim & length(dimsize) > 1) {
                    swi <- sw[i]
                    scalef <- dimsize / 10
                    swi <- swi / scalef
                    svgAnimatePointSW(swi,
                                      begin[i], interp[i], dur[i],
                                      rep[i], rev[i],
                                      id=subName, svgdev=dev@dev)
                 }
             },
             # Any other attribute
             {
                 svgAnimate(animation,
                            paste(ithValue(animSet$animations[[animation]], i),
                                  collapse=";"),
                            begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
             })
  }
    }
}

applyAnimation.text <- function(x, animSet, animation, group, dev) {
    
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
    # We may be dealing with multiple points that need animating
    n <- max(length(x$x), length(x$y), length(x$label))

    # Rep the original x/y/width/height out to be the same length
    x$x <- rep(x$x, length.out=n)
    x$y <- rep(x$y, length.out=n)
    
    # Repeating animation parameters so that each element can have
    # distinct values
    begin <- rep(animSet$begin, length.out = n)
    interp <- rep(animSet$interp, length.out = n)
    dur <- rep(animSet$duration, length.out = n)
    rep <- rep(animSet$rep, length.out = n)
    rev <- rep(animSet$revert, length.out = n)

    angle <- current.rotation()
    
    for (i in 1:n) {
        subName <- subGrobName(x$name, i)

        if ("x" %in% names(animSet$animations))
            xi <- ithUnit(animSet$animations$x, x$x, i)
        else
            xi <- x$x[i]
        if ("y" %in% names(animSet$animations))
            yi <- ithUnit(animSet$animations$y, x$y, i)
        else
            yi <- x$y[i]

        switch(animation,
               x={
                   loc <- locToInches(xi, yi, dev)
                   additive <- "replace"
                   if (angle != 0) {
                       svgAnimateRotation(angle, cx(loc$x, dev), cy(loc$y, dev),
                                          begin[i], interp[i], dur[i],
                                          rep[i], rev[i],
                                          id=subName, svgdev=dev@dev)
                       additive <- "sum"
                   }
                   svgAnimateTranslation(cx(loc$x, dev), cy(loc$y, dev),
                                         begin[i], interp[i], dur[i],
                                         rep[i], rev[i],
                                         additive, subName, dev@dev)
               },
               y={
                   if (!("x" %in% names(animSet$animations))) {
                       loc <- locToInches(xi, yi, dev)
                       additive <- "replace"
                       if (angle != 0) {
                           svgAnimateRotation(angle,
                                              cx(loc$x, dev), cy(loc$y, dev),
                                              begin[i], interp[i], dur[i],
                                              rep[i], rev[i],
                                              id=subName, svgdev=dev@dev)
                           additive <- "sum"
                       }
                       svgAnimateTranslation(cx(loc$x, dev), cy(loc$y, dev),
                                             begin[i], interp[i], dur[i],
                                             rep[i], rev[i],
                                             additive, subName, dev@dev)
                   }
               },
               # Any other attribute
               {
                   svgAnimate(animation,
                              paste(ithValue(animSet$animations[[animation]], i),
                                    collapse=";"),
                              begin[i], interp[i], dur[i], rep[i], rev[i],
                              # Apply these to child <text> element rather
                              # than parent <g>
                              paste(subName, "text",
                                    sep=getSVGoption("id.sep")),
                              dev@dev)
               })
    }
    }
}

doNotAnimate <- function(animSet, animation) {
    # Avoid doing BOTH x and y if BOTH animated
    animNames <- names(animSet$animations)
    if ((all(c("x", "y") %in% animNames) &&
         animation %in% c("x", "y") &&
         match(animation, animNames) ==
         max(match(c("x", "y"), animNames))) ||
        (sum(c("x0", "y0", "x1", "y1") %in% animNames) > 1 &&
         animation %in% c("x0", "y0", "x1", "y1") &&
         match(animation, animNames) >
         min(match(c("x0", "y0", "x1", "y1"), animNames), na.rm=TRUE)))
        TRUE
    else
        FALSE
}

applyAnimation.lines <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
    if (doNotAnimate(animSet, animation))
        return()
    
    # NOTE:  only ever drawing ONE line
    begin <- animSet$begin
    interp <- animSet$interp
    dur <- animSet$duration
    rep <- animSet$rep
    rev <- animSet$revert
    
    subName <- subGrobName(x$name, 1)
    
    if ("x" %in% names(animSet$animations)) {
        au <- ithAnimUnit(animSet$animations$x, x$x, 1)
        xx <- au$values
        timeid <- au$timeid
    } else {
        xx <- x$x
    }
    if ("y" %in% names(animSet$animations)) {
        au <- ithAnimUnit(animSet$animations$y, x$y, 1)
        yy <- au$values
        timeid <- au$timeid
    } else {
        yy <- x$y
    }
    
    if (any(c("x", "y") %in% names(animSet$animations))) {
        loc <- locToInches(xx, yy, dev)
        svgAnimatePoints(cx(loc$x, dev), cy(loc$y, dev), timeid,
                         begin, interp, dur, rep, rev, subName, dev@dev)
    }
    # Any other attribute
    if (!(animation %in% c("x", "y"))) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   begin, interp, dur, rep, rev, subName, dev@dev)
    }
    }  
}
  
applyAnimation.polyline <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
    if (doNotAnimate(animSet, animation))
        return()
    
  # If we only have one line
    if (is.null(x$id) && is.null(x$id.lengths)) {
        x$id <- rep(1L, length(x$x))
    }

  # Multiple lines exist
    if (is.null(x$id)) {
        n <- length(x$id.lengths)
        id <- rep(1L:n, x$id.lengths)
    } else {
        n <- length(unique(x$id))
        id <- x$id
    }

  # Repeating animation parameters so that each element can have
  # distinct values
    begin <- rep(animSet$begin, length.out = n)
    interp <- rep(animSet$interp, length.out = n)
    dur <- rep(animSet$duration, length.out = n)
    rep <- rep(animSet$rep, length.out = n)
    rev <- rep(animSet$revert, length.out = n)

    for (i in 1:n) {
        subName <- subGrobName(x$name, i)
        
        if ("x" %in% names(animSet$animations)) {
            au <- ithAnimUnit(animSet$animations$x, x$x, i)
            xx <- au$values
            timeid <- au$timeid
        } else {
            xx <- x$x[x$id == i]
        }
        if ("y" %in% names(animSet$animations)) {
            au <- ithAnimUnit(animSet$animations$y, x$y, i)
            yy <- au$values
            timeid <- au$timeid
        } else {
            yy <- x$y[x$id == i]
        }

        if (any(c("x", "y") %in% names(animSet$animations))) {
            loc <- locToInches(xx, yy, dev)
            svgAnimatePoints(cx(loc$x, dev), cy(loc$y, dev), timeid,
                             begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
        }
        # Any other attribute
        if (!(animation %in% c("x", "y"))) {
            svgAnimate(animation,
                       paste(ithValue(animSet$animations[[animation]], i),
                             collapse=";"),
                       begin[i], interp[i], dur[i], rep[i], rev[i], subName, dev@dev)
        }
    }
    }
}

# Possibly multiple line segments, each of which become <polyline> elements
applyAnimation.segments <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {

        if (doNotAnimate(animSet, animation))
            return()
        
        # We may be dealing with multiple rects that need animating
        n <- max(length(x$x0), length(x$y0),
                 length(x$x1), length(x$y1))

        # Rep the original x/y/width/height out to be the same length
        x$x0 <- rep(x$x0, length.out=n)
        x$y0 <- rep(x$y0, length.out=n)
        x$x1 <- rep(x$x1, length.out=n)
        x$y1 <- rep(x$y1, length.out=n)
  
        # Repeating animation parameters so that each element can have
        # distinct values
        begin <- rep(animSet$begin, length.out = n)
        interp <- rep(animSet$interp, length.out = n)
        dur <- rep(animSet$duration, length.out = n)
        rep <- rep(animSet$rep, length.out = n)
        rev <- rep(animSet$revert, length.out = n)
        
        for (i in 1:n) {
            subName <- subGrobName(x$name, i)

            if ("x0" %in% names(animSet$animations))
                x0i <- ithUnit(animSet$animations$x0, x$x0, i)
            else
                x0i <- x$x0[i]
            if ("y0" %in% names(animSet$animations))
                y0i <- ithUnit(animSet$animations$y0, x$y0, i)
            else
                y0i <- x$y0[i]
            if ("x1" %in% names(animSet$animations))
                x1i <- ithUnit(animSet$animations$x1, x$x1, i)
            else
                x1i <- x$x1[i]
            if ("y1" %in% names(animSet$animations))
                y1i <- ithUnit(animSet$animations$y1, x$y1, i)
            else
                y1i <- x$y1[i]
    
            if (any(c("x0", "y0", "x1", "y1") %in%
                    names(animSet$animations))) {
                nvals <- max(length(x0i), length(y0i),
                             length(x1i), length(y1i))
                loc0 <- locToInches(x0i, y0i, dev)
                loc1 <- locToInches(x1i, y1i, dev)
                xx <- rbind(convertX(loc0$x, "inches", valueOnly=TRUE),
                            convertX(loc1$x, "inches", valueOnly=TRUE))
                yy <- rbind(convertY(loc0$y, "inches", valueOnly=TRUE),
                            convertY(loc1$y, "inches", valueOnly=TRUE))
                svgAnimatePoints(cx(unit(xx, "inches"), dev),
                                 cy(unit(yy, "inches"), dev),
                                 rep(1:nvals, each=2), # timeid
                                 begin[i], interp[i], dur[i],
                                 rep[i], rev[i], subName, dev@dev)
            }
            # Any other attribute
            if (!(animation %in% c("x0", "y0", "x1", "y1"))) {
                svgAnimate(animation,
                           paste(ithValue(animSet$animations[[animation]], i),
                                 collapse=";"),
                           begin[i], interp[i], dur[i],
                           rep[i], rev[i], subName, dev@dev)
            }
            
        }
    }
}

applyAnimation.polygon <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur,
                   animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
        
        if (doNotAnimate(animSet, animation))
            return()
    
        # If we only have one polygon
        if (is.null(x$id) && is.null(x$id.lengths)) {
            x$id <- rep(1L, length(x$x))
        }

        # Multiple polygons exist
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }

        # Repeating animation parameters so that each element can have
        # distinct values
        begin <- rep(animSet$begin, length.out = n)
        interp <- rep(animSet$interp, length.out = n)
        dur <- rep(animSet$duration, length.out = n)
        rep <- rep(animSet$rep, length.out = n)
        rev <- rep(animSet$revert, length.out = n)

        for (i in 1:n) {
            subName <- subGrobName(x$name, i)
        
            if ("x" %in% names(animSet$animations)) {
                au <- ithAnimUnit(animSet$animations$x, x$x, i)
                xx <- au$values
                timeid <- au$timeid
            } else {
                xx <- x$x[x$id == i]
            }
            if ("y" %in% names(animSet$animations)) {
                au <- ithAnimUnit(animSet$animations$y, x$y, i)
                yy <- au$values
                timeid <- au$timeid
            } else {
                yy <- x$y[x$id == i]
            }

            if (any(c("x", "y") %in% names(animSet$animations))) {
                loc <- locToInches(xx, yy, dev)
                svgAnimatePoints(cx(loc$x, dev), cy(loc$y, dev), timeid,
                                 begin[i], interp[i], dur[i], rep[i], rev[i],
                                 subName, dev@dev)
            }
            # Any other attribute
            if (!(animation %in% c("x", "y"))) {
                svgAnimate(animation,
                           paste(ithValue(animSet$animations[[animation]], i),
                                 collapse=";"),
                           begin[i], interp[i], dur[i], rep[i], rev[i],
                           subName, dev@dev)
            }
        }
    }
}

applyAnimation.pathgrob <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur,
                   animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
        
        if (doNotAnimate(animSet, animation))
            return()

        # NOTE:  only ever drawing ONE line
        begin <- animSet$begin
        interp <- animSet$interp
        dur <- animSet$duration
        rep <- animSet$rep
        rev <- animSet$revert

        subName <- subGrobName(x$name, 1)

        # Rather than looping through 'n' different shapes
        # need to generate a set of animation values for a
        # single shape consisting of multiple sub-paths
        # at multiple time points

        # If we only have one sub-path
        if (is.null(x$id) && is.null(x$id.lengths)) {
            x$id <- rep(1L, length(x$x))
        }

        # Multiple sub-paths
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }
            
        # NOTE:  according to the SVG spec, I think the animated values
        #        HAVE to follow the original series of M, L, Z to be
        #        valid;  otherwise behaviour of browser is undefined?
        if ("x" %in% names(animSet$animations)) {
            au <- as.animUnit(animSet$animations$x, attr(x$x, "unit"))
            xx <- au$values
            pathid <- au$id
            timeid <- au$timeid
        } else {
            xx <- x$x
        }
        if ("y" %in% names(animSet$animations)) {
            au <- as.animUnit(animSet$animations$y, attr(x$y, "unit"))
            yy <- au$values
            pathid <- au$id
            timeid <- au$timeid
        } else {
            yy <- x$y
        }

        # Only one sub-path
        if (is.null(pathid)) {
            pathid <- rep(id, length.out=length(timeid))
        }
    
        if (any(c("x", "y") %in% names(animSet$animations))) {
            loc <- locToInches(xx, yy, dev)
            svgAnimatePath(cx(loc$x, dev), cy(loc$y, dev), pathid, timeid,
                           begin, interp, dur, rep, rev, subName, dev@dev)
        }
        # Any other attribute
        if (!(animation %in% c("x", "y"))) {
            svgAnimate(animation,
                       paste(ithValue(animSet$animations[[animation]], 1),
                             collapse=";"),
                       begin, interp, dur, rep, rev, subName, dev@dev)
        }
    }  
}

applyAnimation.rastergrob <- function(x, animSet, animation, group, dev) {

    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur,
                   animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {

        # Raster may have NULL width or height
        x$width <- widthDetails(x)
        x$height <- heightDetails(x)

        # We may be dealing with multiple rasters that need animating
        n <- max(length(x$x), length(x$y), length(x$width), length(x$height))

        # Rep the original x/y/width/height out to be the same length
        x$x <- rep(x$x, length.out=n)
        x$y <- rep(x$y, length.out=n)
        x$width <- rep(x$width, length.out=n)
        x$height <- rep(x$height, length.out=n)
  
        # Repeating animation parameters so that each element can have
        # distinct values
        begin <- rep(animSet$begin, length.out = n)
        interp <- rep(animSet$interp, length.out = n)
        dur <- rep(animSet$duration, length.out = n)
        rep <- rep(animSet$rep, length.out = n)
        rev <- rep(animSet$revert, length.out = n)

        angle <- current.rotation()
        
        for (i in 1:n) {
            subName <- subGrobName(x$name, i)

            # If x AND y change, need to transform together
            # If width/height changes, have to animate x/y as well
            # because SVG <rect> does not have justification
            if ("x" %in% names(animSet$animations))
                xi <- ithUnit(animSet$animations$x, x$x, i)
            else
                xi <- x$x[i]
            if ("y" %in% names(animSet$animations))
                yi <- ithUnit(animSet$animations$y, x$y, i)
            else
                yi <- x$y[i]
            if ("width" %in% names(animSet$animations))
                widthi <- ithUnit(animSet$animations$width, x$width, i)
            else
                widthi <- x$width[i]
            if ("height" %in% names(animSet$animations))
                heighti <- ithUnit(animSet$animations$height, x$height, i)
            else
                heighti <- x$height[i]
            lb <- leftbottom(xi, yi, widthi, heighti,
                             x$just, x$hjust, x$vjust, dev)
    
            switch(animation,
                   x={
                       dim <- dimToInches(widthi, heighti, dev)
                       additive <- "replace"
                       if (angle != 0) {
                           svgAnimateRotation(angle,
                                              cx(lb$x, dev), cy(lb$y, dev),
                                              begin[i], interp[i], dur[i],
                                              rep[i], rev[i],
                                              id=subName, svgdev=dev@dev)
                           additive <- "sum"
                       }
                       svgAnimateTranslation(cx(lb$x, dev),
                                             ch(dim$h, dev) + cy(lb$y, dev),
                                             begin[i], interp[i], dur[i],
                                             rep[i], rev[i],
                                             additive,
                                             subName, dev@dev)
                   },
                   y={
                       # If we are also animating "x" then this has
                       # already been done
                       if (!"x" %in% names(animSet$animations)) {
                           dim <- dimToInches(widthi, heighti, dev)
                           additive <- "replace"
                           if (angle != 0) {
                               svgAnimateRotation(angle,
                                                  cx(lb$x, dev), cy(lb$y, dev),
                                                  begin[i], interp[i], dur[i],
                                                  rep[i], rev[i],
                                                  id=subName, svgdev=dev@dev)
                               additive <- "sum"
                           }
                           svgAnimateTranslation(cx(lb$x, dev),
                                                 ch(dim$h, dev) +
                                                 cy(lb$y, dev),
                                                 begin[i], interp[i], dur[i],
                                                 rep[i], rev[i],
                                                 additive,
                                                 subName, dev@dev)
                       }
                   },
                   width={
                       dim <- dimToInches(widthi, heighti, dev)
                       # If x is also animated,
                       # this has already been handled above
                       if (!any(c("x", "y") %in% names(animSet$animations))) {
                           if (angle != 0) {
                               svgAnimateRotation(angle,
                                                  cx(lb$x, dev), cy(lb$y, dev),
                                                  begin[i], interp[i], dur[i],
                                                  rep[i], rev[i],
                                                  id=subName, svgdev=dev@dev)
                               additive <- "sum"
                           }
                           svgAnimateTranslation(cx(lb$x, dev),
                                                 ch(dim$h, dev) +
                                                 cy(lb$y, dev),
                                                 begin[i], interp[i], dur[i],
                                                 rep[i], rev[i],
                                                 additive,
                                                 subName, dev@dev)
                       }
                       svgAnimateScale(cw(dim$w, dev), -ch(dim$h, dev), 
                                       begin[i], interp[i], dur[i],
                                       rep[i], rev[i],
                                       id=paste(subName, "scale",
                                           sep=getSVGoption("id.sep")),
                                       svgdev=dev@dev)
                   },
                   height={
                       # If "width" is also animated,
                       # this has already been done
                       if (!"width" %in% names(animSet$animations)) {
                           dim <- dimToInches(widthi, heighti, dev)
                           if (!any(c("x", "y") %in%
                                    names(animSet$animations))) {
                               if (angle != 0) {
                                   svgAnimateRotation(angle,
                                                      cx(lb$x, dev),
                                                      cy(lb$y, dev),
                                                      begin[i], interp[i],
                                                      dur[i], rep[i], rev[i],
                                                      id=subName,
                                                      svgdev=dev@dev)
                                   additive <- "sum"
                               }
                               svgAnimateTranslation(cx(lb$x, dev),
                                                     ch(dim$h, dev) +
                                                     cy(lb$y, dev),
                                                     begin[i], interp[i],
                                                     dur[i], rep[i], rev[i],
                                                     additive,
                                                     subName, dev@dev)
                           }
                           svgAnimateScale(cw(dim$w, dev), -ch(dim$h, dev), 
                                           begin[i], interp[i], dur[i],
                                           rep[i], rev[i],
                                           id=paste(subName, "scale",
                                                 sep=getSVGoption("id.sep")),
                                           svgdev=dev@dev)
                       }
                   },
                   # Any other attribute
                   {
                       svgAnimate(animation,
                                  paste(ithValue(animSet$animations[[animation]], i),
                                        collapse=";"),
                                  begin[i], interp[i], dur[i],
                                  rep[i], rev[i], subName, dev@dev)
                   })
        }
    }
}

applyAnimation.xspline <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur,
                   animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } else {
        
        if (doNotAnimate(animSet, animation))
            return()

        # If we only have one xspline
        if (is.null(x$id) && is.null(x$id.lengths)) {
            x$id <- rep(1L, length(x$x))
        }

        # Multiple xsplines exist
        if (is.null(x$id)) {
            n <- length(x$id.lengths)
            id <- rep(1L:n, x$id.lengths)
        } else {
            n <- length(unique(x$id))
            id <- x$id
        }

        # Repeating animation parameters so that each element can have
        # distinct values
        begin <- rep(animSet$begin, length.out = n)
        interp <- rep(animSet$interp, length.out = n)
        dur <- rep(animSet$duration, length.out = n)
        rep <- rep(animSet$rep, length.out = n)
        rev <- rep(animSet$revert, length.out = n)

        for (i in 1:n) {
            subName <- subGrobName(x$name, i)
        
            if ("x" %in% names(animSet$animations)) {
                au <- ithAnimUnit(animSet$animations$x, x$x, i)
                xx <- au$values
                timeid <- au$timeid
            } else {
                xx <- x$x[x$id == i]
            }
            if ("y" %in% names(animSet$animations)) {
                au <- ithAnimUnit(animSet$animations$y, x$y, i)
                yy <- au$values
                timeid <- au$timeid
            } else {
                yy <- x$y[x$id == i]
            }
            
            if (any(c("x", "y") %in% names(animSet$animations))) {

                getSplinePoints <- function(x, y, grob) {
                    tempSpline <- grob
                    tempSpline$x <- x
                    tempSpline$y <- y
                    tempSpline$id <- NULL
                    tempSpline$id.lengths <- NULL
                    xsplinePoints(tempSpline)                    
                }
                
                # for each time period, need to convert control
                # points into (x, y) and then generate new timeid
                nval <- length(timeid)
                xx <- rep(xx, length.out=nval)
                yy <- rep(yy, length.out=nval)
                xxx <- split(xx, timeid)
                yyy <- split(yy, timeid)
                points <- mapply(getSplinePoints, xxx, yyy,
                                 MoreArgs=list(grob=x), SIMPLIFY=FALSE)
                timeid <- rep(1:length(points),
                              sapply(points, function(p) length(p$x)))
                xpoints <- do.call("unit.c", lapply(points, function(p) p$x))
                ypoints <- do.call("unit.c", lapply(points, function(p) p$y))
                if (x$open) {
                    # animating a polyline element
                    loc <- locToInches(xpoints, ypoints, dev)
                    svgAnimatePoints(cx(loc$x, dev), cy(loc$y, dev), timeid,
                                     begin[i], interp[i], dur[i],
                                     rep[i], rev[i], subName, dev@dev)
                } else {
                    # animating a path element
                    loc <- locToInches(xpoints, ypoints, dev)
                    pathid <- rep(id, length.out=length(timeid))
                    svgAnimatePath(cx(loc$x, dev), cy(loc$y, dev),
                                   pathid, timeid,
                                   begin[i], interp[i], dur[i],
                                   rep[i], rev[i], subName, dev@dev)
                } 
            }
            # Any other attribute
            if (!(animation %in% c("x", "y"))) {
                svgAnimate(animation,
                           paste(ithValue(animSet$animations[[animation]], i),
                                 collapse=";"),
                           begin[i], interp[i], dur[i], rep[i], rev[i],
                           subName, dev@dev)
            }
        }
    }  
}


applyAnimation.grob <- function(x, ...) {
    # If we got here, then we've hit something that is not yet implemented
    stop(paste("Animation of ", paste(class(x), collapse=":"),
               " objects is not yet implemented",
               sep=""))
}

applyAnimation.gTree <- function(x, animSet, animation, group, dev) {
    if (group) {
        svgAnimate(animation,
                   paste(ithValue(animSet$animations[[animation]], 1),
                         collapse=";"),
                   animSet$begin, animSet$interp, animSet$dur, animSet$rep, animSet$rev,
                   x$name, dev@dev)
    } 
}

applyAnimationSet <- function(x, animationSet, group, dev) {
    x$name <- getID(x$name, "grob", FALSE)
    animations <- animationSet$animations
    for (i in names(animations)) 
        applyAnimation(x, animationSet, i, group, dev)
}

animate <- function(x, dev) {
    UseMethod("animate")
}

animate.grob <- function(x, dev) {
    lapply(x$animationSets,
           function(as) {
               applyAnimationSet(x, as, FALSE, dev)
           })
    lapply(x$groupAnimationSets,
           function(as) {
               applyAnimationSet(x, as, TRUE, dev)
           })
}

animate.gTree <- function(x, dev) {
    lapply(x$groupAnimationSets,
           function(as) {
               applyAnimationSet(x, as, TRUE, dev)
           })
    # If you want to do something with the 'animationSets'
    # for your gTree then you have to write your own
    # animate() method
}

# NOTE that this has to be a primToDev() method
# NOT a grobToDev() method
# OTHERWISE, viewports will not be set up correctly
primToDev.animated.grob <- function(x, dev) {
    animate(x, dev)
    NextMethod()    
}

# Ensure the animation is retained on a forced grob
forceGrob.animated.grob <- function(x) {
    y <- NextMethod()
    if (inherits(y, "forcedgrob")) {
        y$animationSets <- x$animationSets
        y$groupAnimationSets <- x$groupAnimationSets
        class(y) <- unique(c("animated.grob", class(y)))
    }
    y
}
