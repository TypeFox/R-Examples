#-- ref.r ------------------------------------
# Jens Oehlschlaegel
# created:            14.08.97
# simplified:         30.08.03
# performance tested: 30.08.03
# documented:         30.08.03
# (gpl) 2003
#---------------------------------------------


if (!exists("is.R")){
  is.R <- function()exists("version") && !is.null(vl <- version$language) && vl == "R"
}

#! \name{ref}
#! \alias{ref}
#! \alias{print.ref}
#! \title{ creating references }
#! \description{
#!   Package \code{ref} implements references for S (R/S+).
#!   Function \code{ref} creates references.
#!   For a memory efficient wrapper to matrixes and data.frames which allows nested subsetting see \code{\link{refdata}}
#! }
#! \usage{
#! ref(name, loc = parent.frame())
#! }
#! \arguments{
#!   \item{name}{ name of an (existing) object to be referenced }
#!   \item{loc}{ location of the referenced object, i.e. an environment in R or a frame in S+ }
#! }
#! \details{
#!   In S (R/S+) paramters are passed by value and not by reference.
#!   When passing big objects, e.g. in recursive algorithms, this can quickly eat up memory.
#!   The functions of package \code{ref} allow to pass references in function calls.
#!   The implementation is purely S and should work in R and S+.
#!   Existence of the referenced object is not checked by function \code{ref}.
#!   Usually \code{\link{as.ref}} is more convenient and secure to use.
#!   There is also a print method for references.
#! }
#! \value{
#!   a list with
#!   \item{ name }{  name of the referenced object }
#!   \item{ loc }{ location of the referenced object, i.e. an environment in R or a frame in S+ }
#!   and class "ref"
#! }
#! \note{
#!  Using this type of references is fine for prototyping in a non-objectoriented programming style.
#!  For bigger projects and safer programming you should consider the approach suggested by Henrik Bengtsson
#!  at \url{http://www.maths.lth.se/help/R/ImplementingReferences} (announced to be released as package "oo" or "classes")
#! }
#! \section{WARNING}{
#!  Usually functions in S have no side-effects except for the main effect of returning something.
#!  Working with references circumvents this programming style and can have considerable side-effects.
#!  You are using it at your own risk.
#! }
#! \section{R 1.8 WARNING}{
#!  Changing parts of referenced objects has been slowed down by order of magnitudes since R version 1.8, see performance test examples on the help page for \code{\link{deref}}.
#!  Hopefully the old performance can be restored in future versions.
#! }
#! \section{S+ WARNING}{
#!  Package ref should generally work under R and S+. However, when changing very small parts of referenced objects, using references under S+ might be inefficient (very slow with high temporary memory requirements).
#! }
#! \section{Historical remarks}{
#!   This package goes back to an idea submitted April 9th 1997 and code offered on August 17th 1997 on s-news.
#!   The idea of implementing references in S triggered an intense discussion on s-news. The status reached in 1997 can be summarized as follows:\cr
#!   \enumerate{
#!     \item{\bold{advantage}}{passing by reference can save memory compared to passing by value}
#!     \item{\bold{disadvantage}}{passing by reference is more dangerous than passing by value}
#!     \item{\bold{however}}{the implementation is purely in S, thus rather channels existing danger than adding new danger}
#!     \item{\bold{restriction}}{assigning to a subsetted part of a referenced object was inefficient in S+ (was S+ version 3)}
#!   }
#!   Due to the last restriction the code was never submitted as a mature library.
#!   Now in 2003 we have a stable version of R and astonishingly assigning to a subsetted part of a referenced object \emph{can} be implemented efficient.
#!   This shows what a great job the R core developers have done. In the current version the set of functions for references was dramatically simplified, the main differences to 1997 beeing the following:
#!   \enumerate{
#!     \item{\bold{no idempotence}}{ \code{\link{deref}} and  \code{\link{deref<-}} now are a simple function and no longer are methods. This decision was made due top performance reasons. As a consequence, \code{deref()} no longer is idempotent: one has to know whether an object is a reference. Function \code{\link{is.ref}} provides a test. }
#!     \item{\bold{no write protection}}{ The 1997 suggestion included a write protection attribute of references, allowing for read only references and allowing for references that could only be changed by functions that know the access code. Reasons for this: there is no need for readonly references (due to copy on modify) and oop provides better mechanisms for security. }
#!     \item{\bold{no static variables}}{ The suggestion made in 1997 did include an implementation of static variables realized as special cases of references with a naming convention which reduced the risc of name collisions in the 1997 practice of assigning to frame 0. Now R has namespaces and the oop approach of Henrik Bengtsson using environments is to be prefered over relatively global static objects. }
#!   }
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{
#!  \code{\link{as.ref}}, \code{\link{deref}}, \code{\link{deref<-}}, \code{\link{exists.ref}}, \code{\link{is.ref}}, \code{\link{print.ref}}, \code{\link{HanoiTower}}
#! }
#! \examples{
#!   v <- 1
#!   r <- ref("v")
#!   r
#!   deref(r)
#!   cat("For more examples see ?deref\n")
#! }
#! \keyword{ programming }

ref <-
if(is.R()){
  function(name, loc=parent.frame())
    # creator for ref objects
  {
    temp <- list(
      name  = name
    , loc = loc
    )
    class(temp) <- "ref"
    temp
  }
}else{
  function(name, loc=sys.parent())
    # creator for ref objects
  {
    temp <- list(
      name  = name
    , loc = loc
    )
    class(temp) <- "ref"
    temp
  }
}


print.ref <-
if(is.R()){
  function(x, ...)
    # prints reference, optionally referenced object
  {
    cat("reference to", x$name, "in ")
    print(x$loc)
    if (exists.ref(x)){
      cat("(object exists)\n")
    }else{
      cat("(no object)\n")
    }
    invisible()
  }
}else{
  function(x, ...)
    # prints reference, optionally referenced object
  {
    if (x$loc==1)
      cat("reference to", x$name, "in where=", x$loc)
    else
      cat("reference to", x$name, "in frame=", x$loc)
    if (exists.ref(x)){
      cat(" (object exists)\n")
    }else{
      cat(" (no object)\n")
    }
    invisible()
  }
}



#! \name{as.ref}
#! \alias{as.ref}
#! \title{ coercing to reference }
#! \description{
#!   This function RETURNs a reference to its argument.
#! }
#! \usage{
#! as.ref(obj)
#! }
#! \arguments{
#!   \item{obj}{ an object existing in the current environment/frame }
#! }
#! \value{
#!   an object of class "ref"
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ref}}, \code{\link{deref}} }
#! \examples{
#!   v <- 1
#!   r <- as.ref(v)
#!   r
#!   deref(r)
#! }
#! \keyword{ programming }

as.ref <-
if(is.R()){
  function(obj)
    # return reference to obj
    # if obj is already a reference return the reference directly
  {
    obj.name <- substitute(obj)
    obj.loc <- parent.frame()
    if (!is.name(obj.name))
      stop("obj must be a named object")
    obj.name <- deparse(obj.name)
    if (is.ref(obj)){
      obj
    }else{
      ref(obj.name, obj.loc)
    }
  }
}else{
  function(obj)
    # return reference to obj
    # if obj is already a reference return the reference directly
  {
    obj.name <- substitute(obj)
    obj.loc <- sys.parent()
    if (!is.name(obj.name))
      stop("obj must be a named object")
    obj.name <- deparse(obj.name)
    if (is.ref(obj)){
      obj
    }else{
      ref(obj.name, obj.loc)
    }
  }
}


#! \name{deref}
#! \alias{deref}
#! \alias{deref<-}
#! \title{ dereferencing references }
#! \description{
#!   This functions allow to access a referenced object. \code{deref(ref)} returns the object, and \code{deref(ref) <- value} assigns to the referenced object.
#! }
#! \usage{
#! deref(ref)
#! deref(ref) <- value
#! #the following does not pass R CMD CHECK
#! #deref<-(ref, value)
#! #deref(ref)[1] <- value  # subsetted assignment appears to be inefficent in S+.
#! }
#! \arguments{
#!   \item{ref}{ a reference as returned by \code{\link{ref}} or \code{\link{as.ref}} }
#!   \item{value}{ a value to be assigned to the reference }
#! }
#! \details{
#!   \code{deref} and \code{deref<-} provide convenient access to objects in other environments/frames.
#!   In fact they are wrappers to \code{\link{get}} and \code{\link{assign}}.
#!   However, convenient does not neccessarily means efficient.
#!   If performance is an issue, the direct use of \code{\link{new.env}}, \code{\link{substitute}} and \code{\link{eval}} may give better results.
#!   See the examples below.
#! }
#! \value{
#!   \code{deref} returns the referenced object.
#!   \cr \code{"deref<-"} returns a reference to the modified object, see \code{\link{ref}}.
#! }
#! \references{ Writing R Extensions }
#! \author{ Jens Oehlschlägel }
#! \note{ Subsetted assignment appears to be inefficent in S+. Note the use of \code{\link{substitute}} in the examples. }
#! \seealso{ \code{\link{ref}}, \code{\link{as.ref}},  \code{\link[base]{get}},  \code{\link[base]{assign}},  \code{\link[base]{substitute}},  \code{\link[base]{eval}} }
#! \examples{
#!   # Simple usage example
#!   x <- cbind(1:5, 1:5)          # take some object
#!   rx <- as.ref(x)               # wrap it into a reference
#!   deref(rx)                     # read it through the reference
#!   deref(rx) <- rbind(1:5, 1:5)  # replace the object in the reference by another one
#!   deref(rx)[1, ]                # read part of the object
#!   deref(rx)[1, ] <- 5:1         # replace part of the object
#!   deref(rx)                     # see the change
#!   cat("For performance test examples see ?deref\n")
#!
#!  \dontrun{
#!   ## Performance test examples showing actually passing by reference
#!   # define test matrix size of nmatrix by nmatrix
#!   nmatrix <- 1000
#!   # you might want to use less loops in S+
#!   # you might want more in R versions before 1.8
#!   nloop   <- 10     
#!
#!   # Performance test using ref
#!   t1 <- function(){ # outer function
#!     m <- matrix(nrow=nmatrix, ncol=nmatrix)
#!     a <- as.ref(m)
#!       t2(a)
#!     m[1,1]
#!   }
#!   # subsetting deref is slower (by factor 75 slower since R 1.8 compared to previous versions
#!   # , and much, much slower in S+) ...
#!   t2 <- function(ref){
#!     cat("timing", timing.wrapper(
#!       for(i in 1:nloop)
#!         deref(ref)[1,1] <- i
#!     ), "\n")
#!   }
#!   if (is.R())gc()
#!   t1()
#!   # ... than using substitute
#!   t2 <- function(ref){
#!     obj <- as.name(ref$name)
#!     loc <- ref$loc
#!     cat("timing", timing.wrapper(
#!       for(i in 1:nloop)
#!         eval(substitute(x[1,1] <- i, list(x=obj, i=i)), loc)
#!     ), "\n")
#!   }
#!   if (is.R())gc()
#!   t1()
#!
#!
#!   # Performance test using Object (R only)
#!   # see Henrik Bengtsson package(oo)
#!   Object <- function(){
#!     this <- list(env.=new.env());
#!     class(this) <- "Object";
#!     this;
#!   }
#!   "$.Object" <- function(this, name){
#!     get(name, envir=unclass(this)$env.);
#!   }
#!   "$<-.Object" <- function(this, name, value){
#!     assign(name, value, envir=unclass(this)$env.);
#!     this;
#!   }
#!   # outer function
#!   t1 <- function(){
#!     o <- Object()
#!     o$m <- matrix(nrow=nmatrix, ncol=nmatrix)
#!       t2(o)
#!     o$m[1,1]
#!   }
#!   # subsetting o$m is slower ...
#!   t2 <- function(o){
#!     cat("timing", timing.wrapper(
#!       for(i in 1:nloop)
#!         o$m[1,1] <- i
#!     ), "\n")
#!   }
#!   if (is.R())gc()
#!   t1()
#!   # ... than using substitute
#!   t2 <- function(o){
#!     env <- unclass(o)$env.
#!     cat("timing", timing.wrapper(
#!       for(i in 1:nloop)
#!         eval(substitute(m[1,1] <- i, list(i=i)), env)
#!     ), "\n")
#!   }
#!   if (is.R())gc()
#!   t1()
#!
#!   }
#! }
#! \keyword{ programming }


deref <-
if(is.R()){
  function(ref)
    # returns referenced object
  {
    get(ref$name, envir=ref$loc)
  }
}else{
  function(ref)
    # returns referenced object
  {
    if (ref$loc==1)
      get(ref$name, where=ref$loc)
    else
      get(ref$name, frame=ref$loc)
  }
}

"deref<-" <-
if(is.R()){
  function(ref, value)
    # assigns value to referenced object
  {
    assign(ref$name, value, envir=ref$loc)
    ref
  }
}else{
  function(ref, value)
    # assigns value to referenced object
  {
    if (ref$loc==1)
      assign(ref$name, value, where=ref$loc)
    else
      assign(ref$name, value, frame=ref$loc)
    ref
  }
}


#! \name{is.ref}
#! \alias{is.ref}
#! \alias{exists.ref}
#! \title{ checking (for) references }
#! \description{
#!   \code{is.ref} checks whether an object inherits from class "ref". \cr
#!   \code{exists.ref} checks whether an referenced object exists.
#! }
#! \usage{
#!   is.ref(x)
#!   exists.ref(ref)
#! }
#! \arguments{
#!   \item{x}{ an object that might be a reference }
#!   \item{ref}{ a reference as returned from \code{\link{ref}} or \code{\link{as.ref}} }
#! }
#! \value{
#!   logical scalar
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ref}}, \code{\link[base]{exists}}, \code{\link[base]{inherits}}, \code{\link[base]{class}} }
#! \examples{
#!   v <- 1
#!   good.r <- as.ref(v)
#!   bad.r <- ref("NonExistingObject")
#!   is.ref(v)
#!   is.ref(good.r)
#!   is.ref(bad.r)
#!   exists.ref(good.r)
#!   exists.ref(bad.r)
#! }
#! \keyword{ programming }

exists.ref <-
if(is.R()){
  function(ref)
  {
    exists(ref$name, envir=ref$loc)
  }
}else{
  function(ref)
  {
    if (ref$loc==1)
      exists(ref$name, where=ref$loc)
    else
      exists(ref$name, frame=ref$loc)
  }
}


is.ref <- function(x)
{
  inherits(x, "ref")
}


#! \name{sleep.wrapper}
#! \alias{sleep.wrapper}
#! \alias{memsize.wrapper}
#! \alias{timing.wrapper}
#! \title{ wrapper to get some measures for all platforms }
#! \description{
#!   interrupts execution for specified no. of seconds
#! }
#! \usage{
#! sleep.wrapper(time)
#! }
#! \arguments{
#!   \item{time}{ no. of seconds }
#! }
#! \value{
#!   NULL
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link[base]{Sys.sleep}} }
#! \keyword{ internal }

sleep.wrapper <- if(is.R()){
  function(time)Sys.sleep(time)
}else{
  function(time)sleep(time)
}

memsize.wrapper <- if(is.R()){
  function()sum(gc()[,2])
}else{
  function()memory.size()/ 1048576
}

timing.wrapper <- if(is.R()){
  function(expr)system.time(expr)[1]
}else{
  function(expr)sys.time(expr)[1]
}


#! \name{HanoiTower}
#! \alias{HanoiTower}
#! \alias{move.HanoiTower}
#! \alias{print.HanoiTower}
#! \alias{plot.HanoiTower}
#! \title{ application example for references }
#! \description{
#!   This is an example for using references in S (R/S+) with package \code{ref}.
#!   \code{HanoiTower} implements a recursive algorithm solving the Hanoi tower problem.
#!   It is implemented such that the recursion can be done either by passing the HanoiTower \emph{by reference} or \emph{by value} to the workhorse function \code{move.HanoiTower}.
#!   Furthermore you can choose whether recursion should use \code{\link{Recall}} or should directly call \code{move.HanoiTower}.
#!   As the HanoiTower object is not too big, it can be extended by some garbage MBytes, that will demonstrate the advantage of passing references instead of values.
#!   The deeper we recurse, the more memory we waist by passing values (and the more memory we save by passing references).
#!   Functions \code{move.HanoiTower} and \code{print.HanoiTower} are internal (not intended to be called by the user directly).
#! }
#! \usage{
#!   HanoiTower(n = 5
#!   , parameter.mode = c("reference", "value")[1]
#!   , recursion.mode = c("recall", "direct")[1]
#!   , garbage = 0
#!   , print = FALSE
#!   , plot = TRUE
#!   , sleep = 0
#!   )
#! }
#! \arguments{
#!   \item{n}{ number of slices }
#!   \item{parameter.mode}{ one of "reference" or "value" deciding how to pass the HanoiTower object }
#!   \item{recursion.mode}{ one of "recall" or "direct" deciding how to call recursively }
#!   \item{garbage}{ no. of bytes to add to the HanoiTower size }
#!   \item{print}{ TRUE print the HanoiTower changes }
#!   \item{plot}{ FALSE not to plot the HanoiTower changes }
#!   \item{sleep}{ no. of seconds to wait between HanoiTower changes for better monitoring of progress }
#! }
#! \details{
#!   The Hanoi Tower problem can be described as follows: you have n slices of increasing size placed on one of three locations a,b,c such that the biggest slice is at the bottom, the next biggest slice on top of it and so forth with the smallest slice as the top of the tower.
#!   Your task is to move all slices from one stick to the other, but you are only allowed to move one slice at a time and you may never put a bigger slice on top of a smaller one.
#!   The recursive solution is: to move n slices from a to c you just need to do three steps: move n-1 slices to b, move the biggest slice to c and move n-1 slices from b to c. If n equals 1, just move from a to c.
#! }
#! \value{
#!   invisible()
#! }
#! \author{ Jens Oehlschlägel }
#! \seealso{ \code{\link{ref}}, \code{\link[base]{Recall}} }
#!
#! \examples{
#!     HanoiTower(n=2)
#!
#!  \dontrun{
#!     # small memory examples
#!     HanoiTowerDemoBytes <- 0
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "reference"
#!     , recursion.mode  = "direct"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "reference"
#!     , recursion.mode  = "recall"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "value"
#!     , recursion.mode  = "direct"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "value"
#!     , recursion.mode  = "recall"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     rm(HanoiTowerDemoBytes)
#!
#!     # big memory examples
#!     HanoiTowerDemoBytes <- 100000
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "reference"
#!     , recursion.mode  = "direct"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "reference"
#!     , recursion.mode  = "recall"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "value"
#!     , recursion.mode  = "direct"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     if (is.R())
#!       gc()
#!     HanoiTower(
#!       parameter.mode  = "value"
#!     , recursion.mode  = "recall"
#!     , garbage = HanoiTowerDemoBytes
#!     )
#!     rm(HanoiTowerDemoBytes)
#!   }
#! }
#! \keyword{ programming }


HanoiTower <- function(
  n = 5
, parameter.mode  = c("reference", "value")[1]
, recursion.mode  = c("recall", "direct")[1]
, garbage = 0     # bytes added to emulate bigger object size
, print   = FALSE
, plot    = TRUE
, sleep   = 0
)
  ## outer function
{
  cat("\nHanoiTower() start\n")

  parameter.mode <- match.arg(parameter.mode, c("reference", "value"))
  recursion.mode <- match.arg(recursion.mode, c("recall", "direct"))

  # create a local object which shall be modified by a recursive function, either by ref or by val
  Tower <- list(
    Tower = list(
      a=n:1
    , b=numeric()
    , c=numeric()
    )
  , garbage = 1:(garbage %/% 4)
  )
  class(Tower) <- "HanoiTower"

  # create a local object that shall always be manipulated by ref
  Info <- list(
    print = print
  , plot  = plot
  , sleep = sleep
  , frames = 0
  , MBytes = 0
  )
  Info.ref <- as.ref(Info)

  # show initial object size
  cat("object.size(HanoiTower)=", object.size(Tower), "\n")

  # display first Tower
  if (print)
    print.HanoiTower(Tower)
  if (plot)
    plot.HanoiTower(Tower)

  # inform user about recursion mode
  cat(ifelse(recursion.mode=="direct","recursion.mode='direct'  [uses move.HanoiTower()]","recursion.mode='recall'  [uses Recall()]"),"\n")

  # do the work according to parameter mode
  if(parameter.mode=="value") {
    cat("parameter.mode='value'     [each recursive call copies Tower]\n")
    # Tower is a local object
    seconds <- timing.wrapper(move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, n, 1, 3))
  }else{
    cat("parameter.mode='reference' [reference to TowerObject in local frame of HanoiTower()]\n")
    # Tower parameter to move.HanoiTower is a reference into this frame
    seconds <- timing.wrapper(move.HanoiTower(as.ref(Tower), Info.ref, parameter.mode, recursion.mode, n, 1, 3))
  }
  # here's the work done
  # get Info
  cat("\nHanoiTower() done\n")
  cat("seconds     ", round(seconds, 2), "\n")
  cat("max(frames) ", Info$frames, "\n")
  cat("max(MBytes) ", round(Info$MBytes, 2), "\n")
  invisible()
}


move.HanoiTower <- function(Tower, Info.ref, parameter.mode, recursion.mode, n=1, from=1, to=1)
  ## (inner) recalling function
{

  # emulate a manipulation to the big object
  # (mainly because otherwise S+ is so clever to recognize that the big object is nver changed)
  if (is.ref(Tower)){
    if (length(deref(Tower)$garbage))
      deref(Tower)$garbage[1] <- 1
  }else{
    if (length(Tower$garbage))
      Tower$garbage[1] <- 1
  }

  if (n==1){

    # parameter Tower is either value or reference
    if (is.ref(Tower)){

      # write actions in the final recursion branches
      nfrom <- length(deref(Tower)$Tower[[from]])
      deref(Tower)$Tower[[to]][length(deref(Tower)$Tower[[to]])+1] <- deref(Tower)$Tower[[from]][nfrom]
      length(deref(Tower)$Tower[[from]]) <- nfrom-1

    }else{

      # write actions in the final recursion branches
      nfrom <- length(Tower$Tower[[from]])
      Tower$Tower[[to]][length(Tower$Tower[[to]])+1] <- Tower$Tower[[from]][nfrom]
      length(Tower$Tower[[from]]) <- nfrom-1

    }

    # show progress
    sleep.wrapper(deref(Info.ref)$sleep)
    if (deref(Info.ref)$print)
      print.HanoiTower(Tower)
    if (deref(Info.ref)$plot)
      plot.HanoiTower(Tower)

    # update performance Info
    deref(Info.ref)$frames <- max( deref(Info.ref)$frames, sys.nframe() )
    deref(Info.ref)$MBytes <- max( deref(Info.ref)$MBytes, memsize.wrapper() )

  }else{

  # recall actions
    free <- (1:3)[-c(from, to)]

    if (parameter.mode=="reference"){
      if (recursion.mode=="direct"){
        move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, n-1, from, free)
        move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, 1  , from, to  )
        move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, n-1, free, to  )
      }else{
        Recall(Tower, Info.ref, parameter.mode, recursion.mode, n-1, from, free)
        Recall(Tower, Info.ref, parameter.mode, recursion.mode, 1  , from, to  )
        Recall(Tower, Info.ref, parameter.mode, recursion.mode, n-1, free, to  )
      }
    }else{
      if (recursion.mode=="direct"){
        Tower <- move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, n-1, from, free)
        Tower <- move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, 1  , from, to  )
        Tower <- move.HanoiTower(Tower, Info.ref, parameter.mode, recursion.mode, n-1, free, to  )
      }else{
        Tower <- Recall(Tower, Info.ref, parameter.mode, recursion.mode, n-1, from, free)
        Tower <- Recall(Tower, Info.ref, parameter.mode, recursion.mode, 1  , from, to  )
        Tower <- Recall(Tower, Info.ref, parameter.mode, recursion.mode, n-1, free, to  )
      }
    }
  }
  Tower
}


print.HanoiTower<- function(x, ...){
  # parameter x is either deref or reference
  if (is.ref(x))
    Tower <- deref(x)$Tower
  else
    Tower <- x$Tower
  nmax <- max(unlist(lapply(Tower, length)))
  outp <- matrix("", nrow=nmax, ncol=3)

  if (length(Tower[[1]]))
    outp[1:length(Tower[[1]]), 1] <- as.character(Tower[[1]])
  if (length(Tower[[2]]))
    outp[1:length(Tower[[2]]), 2] <- as.character(Tower[[2]])
  if (length(Tower[[3]]))
    outp[1:length(Tower[[3]]), 3] <- as.character(Tower[[3]])

  outp <- outp[nmax:1,]
  outp <- rbind(outp,rep("-",3),c("A","B","C"))
  outp <- data.frame(outp)
  names(outp) <- c("         .","    .","    .")
  print(outp)
  invisible()
}

plot.HanoiTower <- function(x, ...){
  # parameter x is either deref or reference
  if (is.ref(x))
    Tower <- deref(x)$Tower
  else
    Tower <- x$Tower
  n <- sum(unlist(lapply(Tower,length)))
  plot(1, 1, xlim=c(0,3*n), ylim=c(0,n), type="n")
  for (i in 1:3) if (length(Tower[[i]])>0)
    for (j in 1:length(Tower[[i]])){
      x <- (i-0.5)*n
      y <- j-0.5
      polygon(x+Tower[[i]][j]*c(0.5,-0.5,-0.5,0.5),y+c(0.5,0.5,-0.5,-0.5), density=0)
      text(x,y,Tower[[i]][j], adj=0.5)
    }
  invisible()
}

