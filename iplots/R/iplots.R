#==========================================================================
# iplots - interactive plots for R
# Package version: 1.1-7
#
# $Id: iplots.R 173 2013-11-19 18:21:53Z urbanek $
# (C)Copyright 2003-13 Simon Urbanek, 2006 Tobias Wichtrey
# Authors: Simon Urbanek, Tobias Wichtrey, Alex Gouberman
#
# This copy of iplots is licensed under GPL v2.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#==========================================================================
# -- global variables (in the package environment)
#    (note: since 1.1-4 we hold them above imports in an unlocked env)
# .iplots.fw     - framework object
# .iplots        - list of iplot objects
# .iplot.curid   - current plot ID
# .iplot.current - .iplots[[.iplot.curid]]
# .iplots.env is the environment of this package
#==========================================================================

#==========================================================================
# used objects:
#==========================================================================
#
# ivar
#   vid (int), name (char), obj (jobjRef [SVar])
#
# iplot
#   id (int), iname (char), obj (jobjRef)
#
# from rJava:
#
# jobjRef
#   jobj (int), jclass (char)
#

#==========================================================================
# initialization
#==========================================================================

# we have parts that require user interaction and so
# must be disabled during make check
.inside.make.check <- function() ("CheckExEnv" %in% search())

setClass("iset", representation(obj="jobjRef", name="character"))
setClass("ivar", representation(obj="jobjRef", vid="integer", name="character", iset="iset"))

.onAttach <- function(...) {
  if (length(.issue.warning)) packageStartupMessage(.issue.warning)
}

# library initialization: Add "<iplots>/java/iplots.jar" to classpath,
# initialize Java and create an instance on the Framework "glue" class
.onLoad <- function(lib, pkg) {
  .jpackage(pkg, lib.loc = lib)

  ## variables from iplots 0.x-x
  ipv <- c(".iplots",".iplots.fw",".iset.selection",".isets",".iplot.curid",".iplot.current")
  mipv <- match(ipv, ls(envir=.GlobalEnv, all.names=TRUE))
  if (any(!is.na(mipv))) {
    rm(list=ipv[!is.na(mipv)], envir=.GlobalEnv) 
    warning("iPlots currently don't support saving of sessions. Data belonging to iPlots from your previous session will be discarded.")
  }

  ## insert an environment above imports that will hold variables we can modify
  ## this allows pre-namespace semantics under namespaces
  .imports <- parent.env(topenv())
  .i.par <- parent.env(.imports)
  ipe <- new.env(parent=.i.par)
  attr(ipe, "name") <- "volatiles:iplots"
  parent.env(.imports) <- ipe

  ipe$.restricted.el <- FALSE # restricted event loop use (on OS X in the GUI and shell)
  ipe$.issue.warning <- NULL

  headless <- identical(.jcall("java/lang/System","S","getProperty","java.awt.headless"), "true")
  # disable compatibility mode on Macs (experimental!)
  if (length(grep("^darwin",R.version$os))) {
      # this is a JGR 1.5 hack - it allows us to find out whether we are running
      # inside JGR without touching any classes
      if (!any(.jcall("java/lang/System","S","getProperty","main.class")=="org.rosuda.JGR.JGR")) {
          #.jcall("java/lang/System","S","setProperty","com.apple.eawt.CocoaComponent.CompatibilityMode","false")
          .restricted.el <<- TRUE
          if (!nchar(Sys.getenv("R_GUI_APP_VERSION"))) {
            ## start Quartz only if we are not in headless mode
            if (!headless) {
              # fire up event loop by simply starting a Quartz device
              grDevices::quartz("dummy", 2, 2)
              grDevices::dev.off()
              # improve response time
              # would need QuartzCocoa_SetLatency(10) call
            }
          } else {
              ## disable all handlers as they conflict with the GUI
              .jcall("java/lang/System","S","setProperty","register.about","false")
              .jcall("java/lang/System","S","setProperty","register.open","false")
              .jcall("java/lang/System","S","setProperty","register.preferences","false")
              .jcall("java/lang/System","S","setProperty","register.quit","false")
          }
          ipe$.issue.warning <- "Note: On Mac OS X we strongly recommend using iplots from within JGR.\nProceed at your own risk as iplots cannot resolve potential ev.loop deadlocks.\n'Yes' is assumed for all dialogs as they cannot be shown without a deadlock,\nalso ievent.wait() is disabled.\nMore recent OS X version do not allow signle-threaded GUIs and will fail.\n"
      } else {
          # don't mess JGR up
          .jcall("java/lang/System","S","setProperty","register.about","false")
          .jcall("java/lang/System","S","setProperty","register.open","false")
          .jcall("java/lang/System","S","setProperty","register.preferences","false")
          .jcall("java/lang/System","S","setProperty","register.quit","false")
      }
  }
  
  assign(".iplots.fw", if (headless || nchar(Sys.getenv("NOAWT"))) NULL else .jnew("org/rosuda/iplots/Framework"), ipe)
  
   # we need to reset everything for sanity reasons
  assign(".iset.selection", vector(), ipe)
  assign(".isets", list(), ipe)
  .isets[[1]]<<-list()
  ipe$.iplots<-list()
  ipe$.iplot.curid<-1
  ipe$.iplot.current<-NULL

  if (!is.jnull(.iplots.fw) && (.inside.make.check() || .restricted.el || nchar(Sys.getenv("NOAWT")))) .jcall(.iplots.fw,, "setNoInteractionFlag", TRUE)
}

# helper function to identify a class in a strstr manner (not nice)
.class.strstr <- function(o, class) {
  if (!inherits(o, "jobjRef")) return(FALSE);
  if (length(grep(class, .jclass(o)))) TRUE else FALSE
}

#==========================================================================
# iSet API
#==========================================================================

# select a current dataset
iset.set <- function(which = iset.next()) {
  if (length(which)!=1)
    stop("You must specify exactly one dataset.")

  ci<-.jcall(.iplots.fw,"I","curSetId")+1
  .iset.save(ci)

  if (is.character(which)) {
    nso<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","selectSet",which)
    if (is.null(nso))
      stop("There is no iSet with the name `",which,"`.")
    ci<-.jcall(.iplots.fw,"I","curSetId")+1
    .iplots<<-.isets[[ci]]$iplots
    .iplot.curid<<-.isets[[ci]]$iplot.cid
    .iplot.current<<-.isets[[ci]]$iplot.cur
    ci
  } else {
    if (!is.numeric(which))
      stop("The `which' parameter must be a name or an ID.")

    nso<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","selectSet",as.integer(which-1:1))
    if (is.null(nso))
      stop("Invalid iSet ID.")
    ci<-.jcall(.iplots.fw,"I","curSetId")+1
    .iplots<<-.isets[[ci]]$iplots
    .iplot.curid<<-.isets[[ci]]$iplot.cid
    .iplot.current<<-.isets[[ci]]$iplot.cur
    .jstrVal(.jcall(.iplots.fw,"S","getSetName"))
  }
}

iset.cur <- function() {
  .jcall(.iplots.fw,"I","curSetId")+1
}

iset.next <- function(which=iset.cur()) {
  if (!is.numeric(which))
    stop("'which' must be a number.")
  w<-as.integer(which)-1:1
  t<-.jcall(.iplots.fw,"I","countSets")
  ((w+1) %% t)+1
}

iset.prev <- function(which=iset.cur()) {
  if (!is.numeric(which))
    stop("'which' must be a number.")
  w<-as.integer(which)-1:1
  t<-.jcall(.iplots.fw,"I","countSets")
  ((w-1:1) %% t)+1
}

iset.list <- function() {
  tot<-.jcall(.iplots.fw,"I","countSets")
  l<-vector()
  for (i in 1:tot)
    l[.jstrVal(.jcall(.iplots.fw,"S","getSetName",as.integer(i-1:1)))]<-i
  l
}

.iset.save <-function (ci=iset.cur()) {
  .isets[[ci]]<<-list(iplots=.iplots,iplot.cid=.iplot.curid,iplot.cur=.iplot.current)
}

iset.new <- function(name=NULL, payload=NULL) {
  if (!is.null(payload) && !is.list(payload)) stop("payload must be a list")
  .iset.save()
  if (is.null(name)) name<-.jnull("java/lang/String")
  ci<-.jcall(.iplots.fw,"I","newSet",name)+1
  .iplots<<-list()
  .iplot.current<<-NULL
  .iplot.curid<<-0
  iset.set(ci)
  if (!is.null(payload)) {
      n <- names(payload)
      if (is.null(n)) n <- paste("V",1:length(payload),sep='')
      for (i in 1:length(payload)) ivar.new(n[i], payload[[i]])
  }
  iset(ci)
}

iset.rm <- function(which = iset.cur()) {
    if (inherits(which, "iset")) {
        which <- .jcall(.iplots.fw, "I", "indexOfSet", which@obj)
        if (which < 0)
            stop("invalid set (possibly already deleted")
        which <- which + 1
    }
    if (is.character(which)) {
        which <- .jcall(.iplots.fw, "I", "getSetIdByName", as.character(which)[1])
        if (which < 0) stop("iSet of that name doesn't exist")
        which <- which + 1
    }
    ocs <- iset.cur()
    if (ocs == which) {
        iset.set()
        if (iset.cur() == ocs)
            stop("cannot delete the last iSet, there must be at least one iSet at all times")
    } else .iset.save(ocs)
    if (!.jcall(.iplots.fw, "Z", "removeSetById", as.integer(which-1)))
        stop("cannot remove iSet, invalid id")

    l <- .isets[[which]]
    lapply(l$iplots, .iplot.dispose)
    .isets[[which]]<<-NULL
    ci <- iset.cur()
    .iplots<<-.isets[[ci]]$iplots
    .iplot.curid<<-.isets[[ci]]$iplot.cid
    .iplot.current<<-.isets[[ci]]$iplot.cur
    ci
}

iplot.location <- function(x, y, relative=FALSE, plot=iplot.cur()) {
    if (is.numeric(plot)) plot<-.iplots[[plot]]
    if (inherits(plot, "iplot")) {
        f<-.jcall(plot$obj,"Ljava/awt/Frame;","getFrame")
        if (!is.jnull(f)) {
            if (missing(x) && missing(y)) {
                p <- .jcall(f,"Ljava/awt/Rectangle;","getBounds")
                return (c(x=.jcall(p,"D","getX"), y=.jcall(p,"D","getY"), width=.jcall(p,"D","getWidth"), height=.jcall(p,"D","getHeight")))
            }
            if (missing(x))
                { x <- if (is.logical(relative) && isTRUE(!relative)) iplot.location(plot=plot)[1] else 0 }
            if (missing(y))
                { y <- if (is.logical(relative) && isTRUE(!relative)) iplot.location(plot=plot)[2] else 0 }
            loc <- c(0,0)
            if (inherits(relative, "iplot"))
                loc <- iplot.location(plot=relative)
            else if (isTRUE(relative))
                loc <- iplot.location(plot=plot)
            .jcall(f,"V","setLocation",as.integer(x[1]+loc[1]), as.integer(y[1]+loc[2]))
        }
    } else stop("invalid plot")
}

iplot.size<-function(width, height, plot=iplot.cur()) {
    if (is.numeric(plot))
        plot <- .iplots[[plot]]
    query <- (missing(width) && missing(height))
    if (inherits(plot, "iplot")) {
        width <- if (missing(width)) .jcall(plot$obj,"I","getWidth") else as.integer(width)[1]
        height <- if (missing(height)) .jcall(plot$obj,"I","getHeight") else as.integer(height)[1]
        if (!query) {
            .jcall(plot$obj,,"setSize", width, height)
            .jcall(.jcall(plot$obj,"Ljava/awt/Frame;","getFrame"),,"pack")
        }
        c(width=width, height=height)
    } else stop("invalid plot")
}

iset <- function(which=iset.cur()) {
  if (!is.numeric(which)) {
    which <- .jcall(.iplots.fw, "I", "getSetIdByName", as.character(which)[1])
    if (which<0) stop("iSet of that name doesn't exist")
  } else which <- which - 1:1
  so <- .jcall(.iplots.fw, "Lorg/rosuda/ibase/SVarSet;", "getSet", as.integer(which))
  if (is.jnull(so)) stop("iSet not found")
  new("iset", obj=so, name=.jcall(so, "S", "getName"))
}

isets <- function() {
  tot<-.jcall(.iplots.fw,"I","countSets")
  l <- lapply(1:tot, iset)
  names(l) <- unlist(lapply(l, function(x) x$name))
  l
}

# ivar object constructor
.new.ivar <- function(vid, name, obj, iset=iplots::iset()) new("ivar", vid=vid, name=name, obj=obj, iset=iset)

# iset accessor
"$.iset" <- function(x, name) {
  vid <- .jcall(x@obj, "I", "indexOf", name)
  if (vid < 0) return(NULL)
  vobj <- .jcall(x@obj, "Lorg/rosuda/ibase/SVar;", "at", vid)
  if (is.jnull(vobj)) return (NULL)
  .new.ivar(vid=vid, name=.jcall(vobj, "S", "getName"), obj=vobj, iset=x)
}

"$<-.iset" <- function(x, name, value) {
    vid <- .jcall(x@obj, "I", "indexOf", as.character(name))
    if (vid < 0) { # new name
        ivar.new(name, value)
        return(x)
    }
    # existing name = replace
    vobj <- .jcall(x@obj, "Lorg/rosuda/ibase/SVar;", "at", vid)
    x[[name]] <- value
}

"[[.iset" <- function(x, i) {
  vid <- if (is.numeric(i)) as.integer(i-0.75) else .jcall(x@obj, "I", "indexOf", as.character(i))
  if (vid < 0) stop("subscript out of bounds")
  vobj <- .jcall(x@obj, "Lorg/rosuda/ibase/SVar;", "at", vid)
  if (is.jnull(vobj)) stop("subscript out of bounds")
  .new.ivar(vid=vid, name=.jcall(vobj, "S", "getName"), obj=vobj, iset=x)
}

"[.iset" <- function(x, i=1:(dim(x)[1]), j=1:length(x)) {
    # trick to get negative indexing right
    if (is.numeric(j)) j <- (1:length(x))[j]
    l <- lapply(j, function(p) x[[p]][i])
    names(l) <- unlist(lapply(j, function(p) x[[p]]@name))
    as.data.frame(l,row.names=1:length(l[[1]]))
}

"[<-.iset" <- function(x, i=1:(dim(x)[1]), j=1:length(x), value) {
  j <- (1:length(x))[j]
  if (length(j)==0 || length(i)==0) return(x)
  for (vi in 1:length(j)) {
    v <- x[[j[vi]]]
    if (!is.null(v)) {
      d<-dim(value)[2]
      if (is.null(d)||is.na(d)) {
        va <- if (is.list(value)) value[[((vi-1) %% length(vi))+1]] else value
        if (length(va)<length(i)) va <- rep(va, length.out=length(i))
        v[i] <- va
      } else {
        di <- ((vi-1) %% d)+1
        va <- value[,di]
        if (length(va)<length(i)) va <- rep(va, length.out=length(i))
        v[i] <- va        
      }
    }
  }
  x
}

"[[<-.iset" <- function (x, i, value) {
  if (is.numeric(i)) {
    if (i<1) stop("index out of bounds")
    lx <- length(x)
    if (i > lx+1) stop("index out of bounds (assignment would create uncontiguous iSet)")
    if (i == lx+1) {
      ivar.new(deparse(substitute(value)),value)
      return(x)
    }
  }
  iv <- x[[i]]
  if (!is.null(iv))
      ivar.update(iv, value)
  else
      ivar.new(i, value)
  x
}

names.iset <- function(x)
  unlist(lapply(1:length(x), function(p) { v<-.jcall(x@obj,"Lorg/rosuda/ibase/SVar;","at",as.integer(p-1)); if (is.jnull(v)) NA else .jcall(v,"S","getName") }))

"names<-.iset" <- function(x, value) {
  lx <- length(x)
  v <- as.character(value)[1:lx]
  if (length(v) < lx) v <- c(v, paste("NA:",1:(lx-length(v)),sep=''))
  lapply(1:lx, function(p) { var<-.jcall(x@obj,"Lorg/rosuda/ibase/SVar;","at",as.integer(p-1)); if (!is.jnull(var)) .jcall(var,"V","setName", v[p]) })
  x
}

dim.iset <- function(x)
  c(.jcall(x@obj, "I", "length"), .jcall(x@obj, "I", "count"))

length.iset <- function(x)
  .jcall(x@obj, "I", "count")

length.ivar <- function(x)
  .jcall(x@obj, "I", "size")

print.iset <- function (x, ...) {
  xo <- x@obj
  vc <- .jcall(xo, "I", "count")
  vl <- vector()
  if (vc>0) vl <- names(x)
  print(paste("iSet[",.jcall(xo,"S","getName"),"]{",paste(vl,collapse=','),"}",sep=''),...)
}

print.ivar <- function(x, ...) {
  print(paste("iVar{",.jcall(x@obj,"S","getName")," from ",x@iset@name,"}",sep=''),...)
}

# S4 objects are sometimes printed using `show' so we have to map show to print
setMethod("show", "ivar", function(object) print.ivar(object))
setMethod("show", "iset", function(object) print.iset(object))

#==========================================================================
# variables handling (iVar)
#==========================================================================

.ivar.check.length <- function(cont) {
  dsl <- .jcall(.iplots.fw,"I","getCurVarSetLength")
  if (dsl>0 && length(cont)!=dsl) {
    cat("iSet and data length differ. Please observe the dialog box (it may be hidden by the R window).\n")
    flush.console()
  }
}

# re-format the variable name to prevent collision with existing variables
#.ivar.valid.name <- function(name) {
#  as.character(.jcall(.iplots.fw,"S","getNewTmpVar",as.character(name)))
#}


# create a new variable (we should make this a method to handle various types)
ivar.new <- function (name=deparse(substitute(cont)), cont) {
  if (!is.character(name) || length(name)>1)
    stop("variable name must be a single string")
  if(is.factor(cont) || is.character(cont)) {
    cont <- as.factor(cont)
    .ivar.check.length(cont)
    id<-.jcall(.iplots.fw,"I","newVar",name,as.integer(cont),as.character(levels(cont)))
    if (id==-2) stop("Operation canceled by user.")
    if (id==-3) {
      iset.new()
      id<-.jcall(.iplots.fw,"I","newVar",name,as.integer(cont),as.character(levels(cont)))
      if (id<0)
        stop("Unable to create an iVariable");
    }
    if (id>=0) {
      vobj<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVar;","getVar",id)
      .jcall(.jcast(vobj,"org/rosuda/ibase/SVarFixFact"),"V","createMissingsCat")
      .jcall(vobj,"V","categorize",TRUE)
      return(.new.ivar(vid=id,name=.jcall(vobj,"S","getName"),obj=vobj))
    }
    return (NULL)
  }
  if (!is.numeric(cont)) stop("iSet can handle only factors, numeric, integer and character vectors.")
  .ivar.check.length(cont)
  id<-.jcall(.iplots.fw,"I","newVar",name,.jarray(cont))
  if (id==-2) stop("Operation canceled by user.")
  if (id==-3) {
    iset.new()
    id<-.jcall(.iplots.fw,"I","newVar",name,.jarray(cont))
    if (id<0)
      stop("Unable to create an iVariable");
  }
  if (id>=0) {
    o <- .jcall(.iplots.fw,"Lorg/rosuda/ibase/SVar;","getVar",id)
    .new.ivar(vid=id,name=.jcall(o,"S","getName"),obj=o)
  } else NULL
}

ivar.new.map <- function (name, x, y) {
  id<-.jcall(.iplots.fw,"I","newVar",name,.jarray(x),.jarray(y))
  if (id==-2) stop("Operation canceled by user.")
  if (id==-3) {
    iset.new()
    id<-.jcall(.iplots.fw,"I","newVar",name,.jarray(x),.jarray(y))
    if (id<0)
      stop("Unable to create an map iVariable");
  }
  if (id>=0) {
    o <- .jcall(.iplots.fw,"Lorg/rosuda/ibase/SVar;","getVar",id)
    .new.ivar(vid=id, name=.jcall(o,"S","getName"), obj=o)
  } else NULL
}

# retrieve data back to R
ivar.data <- function(var) {
  if (!inherits(var,"ivar"))
    stop("parameter is not an iVariable.")
  vid<-as.integer(var@vid)
  if(.jcall(.iplots.fw,"I","varIsNum",vid)!=0) .jcall(.iplots.fw,"[D","getDoubleContent",vid) else as.factor(.jcall(.iplots.fw,"[S","getStringContent",vid))
}

# update contents of an existing variable (undocumented!)
ivar.update <- function (var, cont, batch = FALSE) {
  r <- -1
  if (!inherits(var, "ivar"))
    stop("invalid variable")
  if (is.factor(cont))
    r <- .jcall(.iplots.fw, "I", "replaceVar",
                var@vid, as.integer(unclass(cont)-1:1), levels(cont))
  else if (is.numeric(cont) || is.character(cont))
    r <- .jcall(.iplots.fw, "I", "replaceVar", var@vid, .jarray(cont))
  else {
    vf <- factor(as.character(cont))
    r <- .jcall(.iplots.fw, "I", "replaceVar",
                var@vid, as.integer(unclass(vf)-1:1), levels(vf))
  }
  if (r<0) stop("Unable to replace iSet variable contents (either type or length doesn't match)")
  if (!batch)
    iset.updateVars()
}

"[.ivar" <- function(x,i,...) if (missing(i)) ivar.data(x) else ivar.data(x)[i]
"[<-.ivar"<- function(x,...,value) { a <- ivar.data(x); a[...] <- value; ivar.update(x, a) }


#`<-.ivar` <- ivar.update

#==========================================================================
# iplot management functions
#==========================================================================

iplot.set <- function (which=iplot.next()) {
  a<-try(.iplots[[which<-as.integer(which)]])
  if (!is.list(a)) {
    warning("There is no such plot")
  } else {
    .iplot.current<<-a
    .iplot.curid<<-which
  }
}

iplot.list <- function () {
  if (exists(".iplots")) .iplots else list()
}

# close doesn't dispose of the plot back-end, use .iplot.dispose instead
.iplot.close <- function(plot) {
    .jcall(.jcall(plot$obj,"Ljava/awt/Frame;","getFrame"),"V","dispose")
}

# dispose implies close i.e. the window is closed but also other structures are released
.iplot.dispose <- function(plot) {
    .jcall(plot$obj,"V","dispose")
}

iplot.off <- function(plot=iplot.cur()) {
  if (length(.iplots)==0) {
    if (inherits(plot, "iplot")) {
      warning("Plot is not on the plot list (which is empty), attempting to close anyway.")
      .iplot.dispose(plot)
    }
    return()
  }

  if (inherits(plot,"iplot")) {
    po <- .jcast(plot$obj, "java/lang/Object")
    for (i in 1:length(.iplots))
      if (.jcall(.iplots[[i]]$obj, "Z", "equals", po)) {
        plot <- i
        break
      }
    if (!is.integer(plot)) {
      warning("Plot is not on the plot list, attempting to close anyway.")
      .iplot.dispose(plot)
      return()
    }
  }
  
  j<-list()
  k <- 1
  for (i in 1:length(.iplots))
    if (i!=plot) {
      j[[k]]<-.iplots[[i]]
      k<-k+1
    } else
      .iplot.dispose(.iplots[[i]])
  .iplots <<- j
  if (length(j)==0) {
    .iplot.current<<-NULL
    .iplot.curid<<-0
  } else
    iplot.set(iplot.prev())
  invisible()
}

# this function should not be called directly by the user
# is creates a new R iplot-object for an existing Java plot object
.iplot.new <- function (plotObj, class) {
  if (!exists(".iplots")) {
    .iplots <<- list()
  }
  a<-list(id=(length(.iplots)+1),obj=plotObj)
  class(a)<-c("iplot", class)
  attr(a,"iname")<-.jstrVal(.jcall(plotObj,"S","getTitle"))
  .iplot.current <<- .iplots[[.iplot.curid <<- (length(.iplots)+1)]] <<- a
}

iplot.cur <- function() { .iplot.curid }
iplot.next <- function(which = iplot.cur()) {
 (which%%length(.iplots))+1
}
iplot.prev <- function(which = iplot.cur()) {
 ((which-2)%%length(.iplots))+1
}

# generic print for i[lots
print.iplot <- function(x, ...) { cat("ID:",x$id," Name: \"",attr(x,"iname"),"\"\n",sep="") }

# low-level plot calls
.iplot.setXaxis <- function(ipl,x1,x2) { .jcall(.jcall(ipl,"Lorg/rosuda/ibase/toolkit/Axis;","getXAxis"),"Z","setValueRange",x1,x2-x1); }
.iplot.setYaxis <- function(ipl,x1,x2) { .jcall(.jcall(ipl,"Lorg/rosuda/ibase/toolkit/Axis;","getYAxis"),"Z","setValueRange",x1,x2-x1); }

.iplot.iPlot <- function (x,y,...) {
  a<-.iplot.new(.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ScatterCanvas;","newScatterplot",x@vid,y@vid), "ixyplot")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iHist <- function (var, ...) {
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/HistCanvas;","newHistogram",var@vid),"ihist")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iBar  <- function (var, ...) {
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/BarCanvas;","newBarchart",var@vid),"ibar")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iMap  <- function (var, ...) {
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/MapCanvas;","newMap",var@vid),"imap")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iBox <- function (x, y, ...) {
  if (!inherits(x,"ivar")) {
    vv <- integer(length(x))
    for (i in seq.int(length(x))) {
      v <- x[[i]]
      vv[i] <- if (inherits(v, "ivar")) v@vid else v
    }
    if (length(x) == 1L && !is.null(y))
      a <- .iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ParallelAxesCanvas;","newBoxplot", as.integer(vv), y@vid), "iboxplot")
    else if (is.null(y))
      a <- .iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ParallelAxesCanvas;","newBoxplot",as.integer(vv)),"iboxplot")
    else stop("Cannot have multivariate `x' and a grouping variable")
  } else {
    if (is.null(y))
      a <- .iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ParallelAxesCanvas;","newBoxplot",x@vid),"iboxplot")
    else
      a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ParallelAxesCanvas;","newBoxplot",x@vid,y@vid),"iboxplot")
  }
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iHammock <- function (vars, ...) {
  vv<-vector()
  for (v in vars) {
    if (inherits(v, "ivar"))
      vv<-c(vv,v@vid)
    else
      vv<-c(vv,v)
  }
  if (length(vv)<2)
    stop("At least 2 valid variables are necessary for a hammock plot")
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/HamCanvas;","newHammock",as.integer(vv)),"ihammock")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iPCP  <- function (vars, ...) {
  vv<-vector()
  for (v in vars) {
    if (inherits(v, "ivar"))
      vv<-c(vv,v@vid)
    else
      vv<-c(vv,v)
  }
  if (length(vv)<2)
    stop("At least 2 valid variables are necessary for a pcp plot")
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/ParallelAxesCanvas;","newPCP",as.integer(vv)),"ipcp")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

.iplot.iMosaic <- function (vars, ...) {
  if (length(vars) < 2L)
    stop("At least 2 valid variables are necessary for a mosaic plot")
  a<-.iplot.new(lastPlot<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/plots/MosaicCanvas;","newMosaic",as.integer(vars)),"imosaic")
  if (length(list(...))>0) iplot.opt(...,plot=a)
  a
}

# processes arguments and returns two components:
# - vars : IDs of variables in the call
# - opts : any other (named) options to the call
#
# variables are either
#  * unnamed arguments to the call
#  * list passed as "vars" argument
# they are converted from their raw form to var IDs pushing them to iset
# if necessary
.var.list <- function (..., .apply.coersion.fn) {
  l <- list(...)
  vars <- NULL
  opts <- l
  if (!length(names(l))) {
    vars <- l
    sl <- substitute(list(...))
    sl[[1]] <- NULL
    names(vars) <- unlist(lapply(sl,deparse))
    if (length(vars)==1 && is.list(vars))
      vars <- vars[[1]]
    opts <- NULL
  } else if (any(names(l)=="")) {
    vars <- l
    sl <- substitute(list(...))
    sl[[1]] <- NULL
    sl[names(l) != ""] <- NULL
    vars[names(l) != ""] <- NULL
    opts[names(opts) == ""] <- NULL
    if (length(vars) == 1L && is.list(vars[[1]]))
      vars <- vars[[1]]
    else
      names(vars) <- unlist(lapply(sl, deparse))
  } else {
    m <- pmatch(names(l), "vars")
    if (!is.na(m)) {
      vars <- l[[which(m==1)]]
      opts[[which(m==1)]] <- NULL
    }
  }
  
  if (is.null(vars)) stop("Missing data")

  len <- length(vars)
  if (inherits(vars, "iset"))
    vars <- lapply(seq.int(len), function(i) vars[[i]])
  if (inherits(vars, "ivar"))
    vv <- vars@vid
  else {
    if (!is.list(vars)) vars <- list(V=vars)
    vv <- integer(length(vars))
    if (length(vars) != length(names(vars)))
      names(vars) <- rep("V",length(vars))
    for (v in 1:length(vars))
      vv[v] <- if (inherits(vars[[v]], "ivar")) vars[[v]]@vid else ivar.new(names(vars)[v], if (missing(.apply.coersion.fn)) vars[[v]] else .apply.coersion.fn(vars[[v]]))@vid
  }
  if (!length(vv)) stop("Missing data")
  #do.call(".iplot.iMosaic",c(list(vars=vv), opts))
  list(vars=vv, opts=opts)
}

#==========================================================================
# user-level plot calls
#==========================================================================


imosaic <- function(...) {
  l <- .var.list(..., .apply.coersion.fn=as.factor)
  do.call(".iplot.iMosaic",c(list(vars=l$vars), l$opts))
}

iplot <- function(x, y=NULL, xlab=NULL, ylab=NULL, ...) {
  if (inherits(x,"ivar") || inherits(y,"ivar")) {
    lx <- if (inherits(x,"ivar")) .jcall(x@obj,"I","size") else length(x)
    ly <- if (inherits(y,"ivar")) .jcall(y@obj,"I","size") else length(y)
  } else {
    xlabel <- if (!missing(x))
        deparse(substitute(x))
    ylabel <- if (!missing(y))
        deparse(substitute(y))
    xy <- xy.coords(x, y, xlabel, ylabel, log="")
    xlab <- if (is.null(xlab))
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab))
        xy$ylab
    else ylab
    nx<-xlab
    ny<-ylab
    if (!is.factor(x)) x<-xy$x
    if (!is.factor(y)) y<-xy$y
    lx<-length(x)
    ly<-length(y)
  }
  if (lx<2)
    stop("iplot coordinates must specify at least two points")
  if (!inherits(y,"ivar") && length(y)==1) ry<-rep(y,lx)
  if (lx!=ly)
    stop(paste("Incompatible vector lengths (",lx," and ",ly,").\nBoth vectors x and y must be of the same length.",sep=''))

  if (!inherits(x,"ivar"))
    x<-ivar.new(nx, x)
  if (!inherits(y,"ivar"))
    y<-ivar.new(ny, y)
  if (is.null(x)) stop("Invalid X variable")
  if (is.null(y)) stop("Invalid Y variable")
  .iplot.iPlot(x,y,...)
}

ibar <- function(var, ...) {
  len<-length(var)
  if (inherits(var,"ivar")) len<-.jcall(var@obj,"I","size")
  if (len<2)
    stop("ibar requires at least two data points")
  if (!inherits(var, "ivar") && (is.vector(var) || is.factor(var)) && length(var)>1)
     var<-ivar.new(deparse(substitute(var))[1], as.factor(var))
  .iplot.iBar(var, ...)
}

imap <- function(x, y=NULL, ...) {
  if (inherits(x, "map")) {
    y <- x$y
    x <- x$x
  }
  if (!inherits(x, "ivar")) {
    if (is.null(y))
      stop("Missing y coordinates")
    var<-ivar.new.map(deparse(substitute(var))[1], as.numeric(x), as.numeric(y))
  }
  .iplot.iMap(var, ...)
}

ihist <- function(var, ...) {
  varn <- deparse(substitute(var))[1]
  len <- length(var)
  isNum <- is.numeric(var)
  if (inherits(var,"ivar")) {
    len<-.jcall(var@obj,"I","size")
    isNum <- .jcall(var@obj, "Z", "isNum")
  }
  if (len<2)
    stop("ihist requires at least two data points")
  if (!isNum) {
    if (inherits(var, "ivar"))
      stop("Variable must be numeric.")
    else {
      var <- as.numeric(var)
      varn <- paste("as.numeric(",varn,")",sep='')
    }
  }
  if (!inherits(var, "ivar"))
    var<-ivar.new(varn, var);
  .iplot.iHist(var, ...)
}

ibox <- function(x, y=NULL, ...) {
  if (inherits(x, "iset"))
    x <- lapply(seq.int(length(x)), function(i) x[[i]])
  if (!is.null(y) && !inherits(y, "ivar"))
    y <- ivar.new(deparse(substitute(y))[1], as.factor(y));
  if (is.list(x)) {
    vv <- integer()
    for (i in seq.int(length(x))) {
      var <- x[[i]]
      if (!inherits(var, "ivar")) {
        if (!is.numeric(var))
          var <- as.numeric(var)
        varname <- names(x)[[i]]
        if (!is.null(varname))
          var <- ivar.new(varname, var)
        else
          var <- ivar.new("V", var)
      }
      if (inherits(var, "ivar")) vv <- c(vv, var@vid)
    }
    .iplot.iBox(vv,y, ...)
  } else {
    len <- length(x)
    if (inherits(x,"ivar")) len<-.jcall(x@obj,"I","size") else if (!is.numeric(x)) x <- as.numeric(x)
    if (len < 2L)
      stop("ibox requires at least two data points")
    if (!inherits(x,"ivar")) x<-ivar.new(deparse(substitute(x))[1], x);
    .iplot.iBox(x, y, ...)
  }
}

ihammock <- function(...) {
  l <- .var.list(..., .apply.coersion.fn=as.factor)
  if (length(l$vars) < 2) stop("At least two variables are required for hammock plots")
  do.call(".iplot.iHammock",c(list(vars=l$vars), l$opts))
}

ipcp <- function(...) {
  l <- .var.list(...)
  if (length(l$vars)<2) stop("At least two variables are required for PCP")
  do.call(".iplot.iPCP",c(list(vars=l$vars), l$opts))
}

iplot.opt <- function(..., plot=iplot.cur()) {
  if (is.numeric(plot)) plot<-.iplots[[as.integer(plot)]]
  if (length(list(...))==0)
    .iplot.getPar(plot)
  else
    .iplot.opt(...,plot=plot)
}

.iplot.opt <- function(xlim=NULL, ylim=NULL, col=NULL, main=NULL, ..., plot=iplot.cur()) {
  repaint = FALSE

  optlist <- list(...)

  # now, this is ugly, we should use proper class dispatch
  if (inherits(plot, "imosaic")) {
    if (!is.null(optlist$type)) {
      t <- pmatch(optlist$type,c("observed","expected","fluctuation","same.bin.size","multiple.barcharts"))
      if (any(is.na(t))) stop("Invalid type for imosaic")
      .jcall(plot$obj, "Ljava/lang/Object;","run",.jcast(.iplots.fw,"java.lang.Object"), c("observed","expected","fluctuation","samebinsize","multiplebarcharts")[t])
      optlist$type <- NULL
    }
  }
  
  if(length(optlist)) {
    for(i in 1:length(optlist)) .jcall(plot$obj,"V","setOption",names(optlist)[i],optlist[[i]])
    repaint = TRUE
  }

  if (is.numeric(plot)) plot<-.iplots[[as.integer(plot)]]
  if (!is.null(xlim)) .iplot.setXaxis(plot$obj,xlim[1],xlim[2])
  if (!is.null(ylim)) .iplot.setYaxis(plot$obj,ylim[1],ylim[2])
  if (!is.null(main))
    .jcall(.jcall(plot$obj,"Ljava/awt/Frame;","getFrame"),"V","setTitle",as.character(main))
  if (!is.null(col)) iset.brush(col)
  if (!(is.null(xlim) && is.null(ylim))) {
    repaint = TRUE
  }

  if(repaint) {
    .jcall(plot$obj,"V","updateObjects")
    .jcall(plot$obj,"V","setUpdateRoot",as.integer(0))
    .jcall(plot$obj,"V","repaint")
  }

  # dirty temporary fix
  .jcall(plot$obj,"V","forcedFlush")
  invisible()

}

.iplot.getPar <- function(plot=.iplot.current) {
  p <- list()
#  if (.class.strstr(plot$obj,"ScatterCanvas")) {
#    p$xlim<-.Java(.Java(plot$obj,"getXAxis"),"getValueRange")
#    p$ylim<-.Java(.Java(plot$obj,"getYAxis"),"getValueRange")
#  }
  p
}

iplot.data <- function(id=NULL) {
  if (!is.null(id)) {
    v<-.jcall(.iplot.current$obj,"Lorg/rosuda/ibase/SVar;","getData",as.integer(id-1:1))
    if (!is.null(v)) {
      if(.jcall(.iplots.fw,"I","varIsNum",v)!=0)
        .jcall(.iplots.fw,"[D","getDoubleContent",v,evalArray=TRUE)
      else
        as.factor(.jcall(.iplots.fw,"[S","getStringContent",v))
    } else {
      NULL
    }
  } else {
    i<-1
    l<-list()
    while (!is.null(a<-iplot.data(i))) {
      l[[i]]<-a
      if (i==1) names(l)<-"x"
      if (i==2) names(l)<-c("x","y")
      i<-i+1
    }
    l
  }
}

# BaseCanvas methods

iplot.zoomIn<-function(x1,y1,x2,y2) {.jcall(.iplot.current$obj,,"performZoomIn",as.integer(x1),as.integer(y1),as.integer(x2),as.integer(y2));invisible();}
iplot.zoomOut<-function(x,y) {.jcall(.iplot.current$obj,,"performZoomOut",as.integer(x),as.integer(y));invisible();}
iplot.resetZoom<-function() {.jcall(.iplot.current$obj,"V","resetZoom");invisible();}
iplot.rotate<-function(i){.jcall(.iplot.current$obj,"V","rotate",as.integer(i));invisible();}
#//iplot.setMargins<-function(left,right,top,bottom) {
#//	v<-c(left,right,top,bottom);
#//	.jcall(.iplot.current$obj,"V","setDefaultMargins",.jarray(v));invisible();}

iplot.backend <- function(type = NULL) {
  getypes <- c("awt","swing","opengl")
  if (!is.null(type)) {
    ge <- pmatch(type[1],getypes)
    if (any(is.na(ge))) stop("invalid backend type")
    .jcall(.iplots.fw,"V","setGraphicsEngine",as.integer(ge-1:1))
  }
  getypes[.jcall(.iplots.fw,"I","getGraphicsEngine")+1]
}

iplot.setExtendedQuery <- function(str, plotID=.iplot.curid) {
    if(is.null(str) || str==FALSE)
        .jcall(.iplots.fw,,"useExtQueryString",as.integer(plotID),FALSE)
    else
        .jcall(.iplots.fw,,"setExtQueryString",as.integer(plotID),as.character(str)[1])
    invisible()
}

#==========================================================================
# selection/highlighting API
#==========================================================================

iset.select <- function(what, mode="replace", mark=TRUE, batch=FALSE) {
  m<-.jcall(.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","getCurrentSet"),"Lorg/rosuda/ibase/SMarker;","getMarker")
  if (mode=="replace") iset.selectNone(batch=batch)
  if (mode=="intersect")
      stop("I'm sorry, mode='intersect' is not yet supported in this version.")
  if (is.logical(what)) {
    for(i in 1:length(what)) { if (what[i]) .jcall(m,"V","set",as.integer(i-1:1),mark) }
  } else {
    for(i in what) { .jcall(m,"V","set",as.integer(i-1:1),mark) }
  }
  if (!batch) .jcall(.iplots.fw,"V","updateMarker")
  invisible()
}

iset.selected <- function() {
  m<-.jcall(.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","getCurrentSet"),"Lorg/rosuda/ibase/SMarker;","getMarker")
  .jcall(m,"[I","getSelectedIDs",evalArray=TRUE)+1
}

iset.selectAll <- function(batch=FALSE) { .jcall(.jcall(.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","getCurrentSet"),"Lorg/rosuda/ibase/SMarker;","getMarker"),"V","selectAll",as.integer(1)); if (!batch) .jcall(.iplots.fw,"V","updateMarker"); }
iset.selectNone <- function(batch=FALSE) { .jcall(.jcall(.jcall(.iplots.fw,"Lorg/rosuda/ibase/SVarSet;","getCurrentSet"),"Lorg/rosuda/ibase/SMarker;","getMarker"),"V","selectNone"); if (!batch) .jcall(.iplots.fw,"V","updateMarker"); }


#==========================================================================
# brushing API
#==========================================================================

# in paper: iset.color(color, what=iset.selected())
iset.col <- function(col=NULL) { iset.brush(col) }
iset.brush <- function(col=NULL) {
  if (is.null(col) || is.na(col)) col<-as.integer(c(0,0))
  if (is.numeric(col) && !is.integer(col)) col<-as.integer(col)
  if (is.factor(col)) col<-as.integer(as.integer(col)+1)
  .jcall(.iplots.fw,"V","setSecMark",col);
  # dirty temporary fix
  .jcall(.iplot.current$obj,"V","forcedFlush");
  invisible()
}

#** iset.cols(): return colors

iset.updateVars <- function() { .jcall(.iplots.fw,"V","updateVars"); }

#==========================================================================
# iobj API
#==========================================================================

# internal function - creates a new object of the Java-class <type>
.iobj.new <- function(plot, type) {
  pm<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/PlotManager;","getPlotManager")
  a<-list(obj=.jnew(paste("org/rosuda/ibase/toolkit",type,sep="/"),pm),pm=pm,plot=plot)
  class(a)<-"iobj"
  plot$curobj<-a
  # dirty temporary fix
  .jcall(plot$obj,"V","forcedFlush");
  a
}

.get.plot.obj <- function(plot) {
  if (is.numeric(plot))
    plot <- .iplots[[plot]]
  if (is.null(plot))
    stop("There is such plot")
  if (!inherits(plot, "iplot"))
    stop("Invalid plot object")
  plot
}

iobj.list <- function(plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  pm<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/PlotManager;","getPlotManager")
  i<-.jcall(pm,"I","count")
  l<-list()
  if (i>0) {
    for(j in 1:i) {
      a<-list(obj=.jcall(pm,"Lorg/rosuda/ibase/toolkit/PlotObject;","get",as.integer(j-1:1)),pm=pm,plot=plot)
      class(a)<-"iobj"
      l[[j]]<-a
    }
  }
  l
}

iobj.get <- function(pos, plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  if (!inherits(plot,"iplot"))
    stop("The specified plot is no iplot")
  pm<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/PlotManager;","getPlotManager")
  a<-list(obj=.jcall(pm,"Lorg/rosuda/ibase/toolkit/PlotObject;","get",as.integer(pos-1:1)),pm=pm,plot=plot)
  class(a)<-"iobj"
  a
}

iobj.cur <- function(plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  pm<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/PlotManager;","getPlotManager")
  oo<-.jcall(pm,"Lorg/rosuda/ibase/toolkit/PlotObject;","getCurrentObject");
  if (is.null(oo)) {
    NULL
  } else {
    a<-list(obj=oo,pm=pm,plot=plot)
    class(a)<-"iobj"
    a
  }
}

.iobj.equal <- function(a,b) {
  if (!inherits(a,"iobj") || !inherits(b,"iobj")) stop("wrong types to compare")
  (!is.null(a) && !is.null(b) && .jcall(a$obj,"Z","equals",.jcast(b$obj,"java/lang/Object")))
}

`==.iobj` <- .iobj.equal
`!=.iobj` <- function(a,b) !.iobj.equal(a,b)

iobj.next <- function(which=iobj.cur(), plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  l <- iobj.list(plot)
  if (length(l) == 0) return(numeric())
  if (length(l) == 1) return(1)
  if (inherits(which, "iobj")) which <- base::which(unlist(lapply(l,.iobj.equal,which)))
  which <- if (length(which) != 1) 1 else which + 1
  which <- ((which-1:1) %% length(l)) + 1
  which
}

iobj.prev <- function(which=iobj.cur(), plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  l <- iobj.list(plot)
  if (length(l) == 0) return(numeric())
  if (length(l) == 1) return(1)
  if (inherits(which, "iobj")) which <- base::which(unlist(lapply(l,.iobj.equal,which)))
  which <- if (length(which) != 1) 1 else which - 1
  which <- ((which-1:1) %% length(l)) + 1
  which
}

.iplot.get.by.pm <-function (which) { # get plot from the list by pm entry
  for (i in .iplots)
    if (i$pm == which) return(i)
  NULL;
}

iobj.rm <- function(which=iobj.cur(), plot = iplot.cur()) {
  plot <- .get.plot.obj(plot)
  if (is.vector(which) && (inherits(which[1],"iobj") || is.numeric(which))) {
    for (i in which) iobj.rm(i, plot=plot)
    return()
  }
  if (is.numeric(which)) which<-iobj.get(which, plot=plot)
  .jcall(which$pm,"V","rm",.jcast(which$obj,"org/rosuda/ibase/toolkit/PlotObject"))
  which$plot$curobj<-.jcall(which$pm,"I","count")
  .jcall(which$pm,"V","update")
  rm(which)
}

.iobj.opt.PlotText <- function(o,x=NULL,y=NULL,txt=NULL,ax=NULL,ay=NULL) {
  if (!is.null(x) || !is.null(y)) {
    if (is.null(x)) x<-.jcall(o$obj,"[D","getX",evalArray=TRUE)
    if (is.null(y)) y<-.jcall(o$obj,"[D","getY",evalArray=TRUE)
    .jcall(o$obj,"V","set",as.double(x),as.double(y))
  }
  if (!is.null(ax))
    .jcall(o$obj,"V","setAX",as.double(ax))
  if (!is.null(ay))
    .jcall(o$obj,"V","setAY",as.double(ay))
  if (!is.null(txt)) {
    .jcall(o$obj,"V","set",as.character(txt))
  }
}

itext <- function(x, y=NULL, labels=seq(along=x), ax=NULL, ay=NULL, ..., plot=iplot.cur()) {
  plot <- .get.plot.obj(plot)

  pt<-.iobj.new(plot,"PlotText")
  c<-xy.coords(x,y,recycle=TRUE)
  iobj.opt(pt,x=c$x,y=c$y,ax=ax,ay=ay,txt=labels)
  if (length(list(...))>0)
    iobj.opt(pt,...)
  pt
}

.raw.read <- function(con) {
  x <- 20e6
  y <- raw()
  while(TRUE) {
    z <- readBin(con, raw(), x)
    if (length(y)) {
      if (length(z) < x)
        return(c(y, z))
      y <- c(y, z)
    } else {
      if (length(z) < x)
        return(z)
      y <- z
    }
  }
}

iraster <- function(x1, y1, x2, y2, img, ..., plot=iplot.cur()) {
  if (inherits(img, "raster") || inherits(img, "nativeRaster") || is.matrix(img) || is.array(img)) img <- writePNG(img)
  else if (inherits(img, "connection")) img <- .raw.read(img)
  else if (!is.raw(img)) {
    img <- as.character(img)[1L]
    sz <- file.info(img)$size
    if (any(is.na(sz))) stop("cannot read file `", img, "'")
    img <- readBin(img, raw(), sz)
  }
  plot <- .get.plot.obj(plot)
  if (length(x1) == 4L && missing(y1)) {
    y1 <- x1[2L]
    x2 <- x1[3L]
    y2 <- x1[4L]
    x1 <- x1[1L]
  } else if (length(x1) == 2L && missing(y1) && !missing(x2) && length(x2) == 2L) {
    y1 <- x1[2L]
    x1 <- x1[1L]
    y2 <- x2[2L]
    x2 <- x2[1L]
  }
  if (length(x1) != 1L || length(x2) != 1L || length(y1) != 1L || length(y2) != 1L)
    stop("Invalid image coordinates")
  
  pt <- .iobj.new(plot, "PlotImage")
  .jcall(pt$obj, "V", "set", as.double(x1), as.double(y1), as.double(x2), as.double(y2))
  .jcall(pt$obj, "V", "setImage", img)
  if (length(list(...)))
    .iobj.opt(pt, ...)
  else {
    .jcall(plot$obj, "V","forcedFlush")
    .jcall(plot$obj, "V","repaint")
  }
  pt
}

.iobj.opt.get.PlotText <- function(o)
  list(x=.jcall(o$obj,"[D","getX",evalArray=TRUE),
       y=.jcall(o$obj,"[D","getY",evalArray=TRUE),
       ax=.jcall(o$obj,"[D","getAX",evalArray=TRUE),
       ay=.jcall(o$obj,"[D","getAY",evalArray=TRUE),
       txt=.jstrVal(.jcall(o$obj,"S","getText")))


iobj.opt <- function(o=iobj.cur(),...) {
  if (is.list(o) && inherits(o[[1]],"iobj")) {
    r<-list()
    for (i in o) r<-c(r, iobj.opt(o=i, ...))
    r
  } else {
    if (length(list(...))==0)
      .iobj.opt.get(o)
    else
      .iobj.opt(o,...)
  }
}

.iobj.opt.get <- function(o) {
  if (!is.null(o)) {
    if (.class.strstr(o$obj, "PlotText"))
      .iobj.opt.get.PlotText(o)
  }
}

.iobj.opt.PlotPolygon <- function(o,x,y=NULL) {
  .co<-xy.coords(x,y)
  x<-.co$x
  y<-.co$y
  .jcall(o$obj,"V","set",as.numeric(x),as.numeric(y))
}

.iobj.opt <- function(o=iobj.cur(),...,col=NULL, fill=NULL, layer=NULL, reg=NULL, visible=NULL, coord=NULL, update=TRUE, a=NULL, b=NULL) {
  if (is.numeric(o)) o<-iobj.get(o)
  if (!is.null(layer)) .jcall(o$obj,"V","setLayer",as.integer(layer))
  if (length(list(...))>0) {
    if (.class.strstr(o$obj,"PlotText"))
      .iobj.opt.PlotText(o,...)
    else if (.class.strstr(o$obj,"PlotPolygon"))
      .iobj.opt.PlotPolygon(o,...)
    else
      .jcall(o$obj,"V","set",...)
  }
  if (!is.null(col)|| !is.null(fill)) .iobj.color(o,col,fill)
  if (!is.null(reg)||!is.null(a)||!is.null(b)) .iabline.set(a=a,b=b,reg=reg,obj=o)
  if (!is.null(visible))
    .jcall(o$obj,"V","setVisible",as.logical(visible))
  if (!is.null(coord) && length(coord)>0) {
    if (length(coord)==1)
      .jcall(o$obj,"V","setCoordinates",as.integer(coord))
    else
      .jcall(o$obj,"V","setCoordinates",as.integer(coord[1]),as.integer(coord[2]))
  }
  if (update) {
    .jcall(o$obj,"V","update")
    # dirty temporary fix
    .jcall(.iplot.current$obj,"V","forcedFlush")
    if (.class.strstr(o$obj, "PlotImage")) ## not sure why, but PlotImage requires explicit repaint
      .jcall(o$plot$obj, "V", "repaint")
  }
}

iobj.set <- function(which=iobj.next()) {
  if (is.numeric(which)) which<-iobj.get(which)
  if (is.null(which)) stop("opject doesn't exist")
  .jcall(which$pm,"V","setCurrentObject",which$obj);
}

.iobj.color <- function(obj=iobj.cur(), col=NULL, fill=NULL) {
  if (is.numeric(obj)) obj<-iobj.get(as.integer(obj))
  if (!is.null(col))
    if (is.na(col))
      .jcall(obj$obj,"V","setDrawColor",NULL)
    else
      .jcall(obj$obj,"V","setDrawColor",.jnew("org/rosuda/ibase/toolkit/PlotColor",col))

  if (!is.null(fill))
    if (is.na(fill))
      .jcall(obj$obj,"V","setFillColor",NULL)
    else
      .jcall(obj$obj,"V","setFillColor",.jnew("org/rosuda/ibase/toolkit/PlotColor",fill))

  .jcall(obj$obj,"V","update")
  # dirty temporary fix
  .jcall(.iplot.current$obj,"V","forcedFlush");
}

print.iobj <- function(x, ...) {
  cat(.jstrVal(x$obj),"\n")
}

ilines <- function(x,y=NULL,col=NULL,fill=NULL,visible=NULL,plot=iplot.cur()) {
  plot <- .get.plot.obj(plot)
  .co<-xy.coords(x,y)
  x<-.co$x
  y<-.co$y
  if (length(x)==1 && length(y)==1) {
    x<-c(x,x)
    y<-c(y,y)
  }
  pp<-.iobj.new(plot,"PlotPolygon")
  .jcall(pp$obj,"V","set",as.numeric(x),as.numeric(y))
  if (!is.null(col) || !is.null(fill))
    .iobj.color(pp,col,fill) # includes "update"
  else
    .jcall(pp$obj,"V","update")
    if (!is.null(visible))
      .jcall(pp$obj,"V","setVisible",visible)
  pp
}

.iabline.AB <- function(a=NULL, b=NULL, reg=NULL, coef=NULL, ...) {
  if (!is.null(reg)) a<-reg
  if (!is.null(a) && is.list(a)) {
    t <- as.vector(coefficients(a))
    if (length(t) == 1) {
      a<-0
      b<-t
    } else {
      a<-t[1]
      b<-t[2]
    }
  }
  if (!is.null(coef)) {
    a <- coef[1]
    b <- coef[2]
  }
  if (is.null(b)) b<-0
  list(a=a,b=b)
}

iabline <- function(a=NULL, b=NULL, reg=NULL, coef=NULL, ..., plot=iplot.cur()) {
  plot <- .get.plot.obj(plot)
  l<-.iabline.AB(a,b,reg,coef)
  a<-l$a
  b<-l$b
  ax<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/Axis;","getXAxis")
  if (is.null(ax)) {
    stop("The plot has no X axis")
  } else {
    r<-.jcall(ax,"[D","getValueRange",evalArray=TRUE)
    mi<-min(r)
    mx<-max(r)
    ilines(c(mi,mx),c(a+b*mi,a+b*mx),...,plot=plot)
  }
}

.iabline.set <- function(...,obj=iobj.cur()) {
  l<-.iabline.AB(...)
  a<-l$a
  b<-l$b
  plot<-obj$plot;
  ax<-.jcall(plot$obj,"Lorg/rosuda/ibase/toolkit/Axis;","getXAxis")
  if (is.null(ax)) {
    stop("The plot has no X axis")
  } else {
    r<-.jcall(ax,"[D","getValueRange",evalArray=TRUE)
    mi<-min(r)
    mx<-max(r)
    iobj.opt(o=obj,c(mi,mx),c(a+b*mi,a+b*mx))
  }
  invisible()
}

ievent.wait <- function() {
  if (.inside.make.check()) {
    cat("NOTE: ievent.wait is likely being run from within `make check', returning NULL to prevent `make check' from stalling.\n")
    return(NULL)
  }
  if (.restricted.el) {
      warning("ievent.wait cannot be used on Mac OS X without JGR, because waiting will deadlock the R, Java and system event loops. Returning NULL.")
      return(NULL)
  }
  msg<-.jcall(.iplots.fw,"Lorg/rosuda/ibase/NotifyMsg;","eventWait")
  if (!is.null(msg)) {
    o<-list(obj=msg,msg=.jcall(msg,"I","getMessageID"),cmd=.jcall(msg,"S","getCommand"),pars=.jcall(msg,"I","parCount"))
    class(o)<-"ievent"
    return(o)
  }
  return(NULL)
}

iset.sel.changed <- function (iset=iset.cur()) {
  a <- iset.selected()
  if (length(a)!=length(.iset.selection))
    b <- TRUE
  else
    b <- !(sum(a==.iset.selection)==length(a)) # this is stupid, but works :P
  if (b) .iset.selection <<- a
  b
}

.iDebug <- function(level=1) {
  .jcall(.iplots.fw,"V","setDebugLevel",as.integer(level))
}
