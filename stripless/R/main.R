#' @name strucplot
#' @title Structured Trellis Plots
#' @description Structured Trellis Plots Without Strip Labels
#' @details  A structured display is a trellis display without strip
#'  labels. The factors and levels are
#'  encoded by plot layout and spacing and decoded by a separate legend.
#'  See the \strong{Overview} below for a detailed explanation.
#'
#' @concept structured plots,structured trellis plots, factorial designs,plot
#'  factorial designs, fractional factorial designs, display factorial designs
#'
#' @section Overview: The trellis display paradigm breaks down when there are
#'  more than 2 or 3 conditioning variables, because the plotting area becomes
#'  cluttered with multiple layers of strip labels identifying the panel
#'  settings. This is especially a problem if one wants to show the structure or
#'  results of studies that are (fractions of) factorial designs with the design
#'  factors as conditioning variables. For example, in industrial type
#'  experiments, it is not uncommon to have 5 or more experimental factors.
#'
#'  It is also often the case that there are multiple responses -- e.g.
#'  several different characteristics of a product or a functional response like
#'  an IR spectrum, MRI scan, or surface plot. It can be desirable to display
#'  such results so that direct visual comparison of the complex responses can be made.
#'
#'  The \code{strucplot} function enables such functionality by omitting the
#'  strip labels and laying out the panels in a regular array, a 'xyLayout' in
#'  which the position of the panels rather than their strip labels identifies
#'  the variable settings.
#'
#'  Because the lattice package already has these capabilities, \code{strucplot}
#'  is implemented as a wrapper to lattice's \code{xyplot} function. This
#'  provides a convenient interface, simplifies the code and, most
#'  importantly, allows the user to access all of lattice's (and the underlying
#'  grid graphics') functionality in the usual way. Only two extra arguments --
#'  'xyLayout' and 'spacings' -- are used to do this, although the default
#'  'spacings' argument normally need not be changed. (There is also a third
#'  'center' argument for 2 level factorial designs, which are commonly used in
#'  industrial experiments, that is explained below.)
#'
#'@section How does strucplot() work?:
#'
#'  Suppose that the data consist of a numeric vector response y for 4
#'  conditioning factors, f1, f2, f3, and f4, where f1 and f2 each have 2
#'  levels, f3 has 3 levels, and f4 has 4 levels. (Because the conditioning
#'  variables are coerced to factors by \code{factor}, if level orderings other
#'  than that given by these coercions are wanted, the user should do them
#'  explicitly before calling \code{strucplot}).
#'
#'  Then the call, \code{strucplot(~y|f1*f2*f3*f4)} would produce a trellis plot
#'  without strip labels with 12 (= 3 x 4) rows and 4 (= 2 x 2) columns of
#'  panels, some of which may be empty if the corresponding factor settings are
#'  missing. The default xyLayout argument that produces this is: \code{xyLayout
#'  = list(x = 1:2, y = 3:4)}. It splits the conditioning variables as evenly as
#'  possible into 2 groups with the x component getting the first 2 variables,
#'  f1 and f2, and the y component getting the second 2 variables, f3 and f4.
#'  (if there are an odd number of variables, the x component gets one more
#'  variable than the y). This means the levels of the x variables, f1 and f2,
#'  vary across each row; and the levels of the y variables, f3 and f4, vary
#'  down each column.  Since there are 4 combinations of levels for the x
#'  variables and 12 for y, this gives a 12 row by 4 column display.
#'
#'  The panels are displayed in each direction in reverse lexicographic order,
#'  where the 'alphabets' are the factor levels. This means that the first
#'  variable changes the fastest; the second the next fastest, and so on. Using
#'  (i,j) to denote setting in which the first factor is at the ith level and
#'  the second is at the jth, this translates to:
#'
#'  \describe{ \item{Row ordering}{ (f1,f2): (1,1), (2,1), (1,2), (2,2)}
#'  \item{Column ordering from top down}{ (f3, f4): ((1,1), (2,1), (3,1), (1,2),
#'  (2,2), ... , (1,4), (2,4), (3,4)} }
#'
#'  If one component is missing, if it is x, there will only be 1 column; if it
#'  is y, only 1 row. The nonmissing component must still be correctly specified
#'  to provide the panel ordering.
#'
#'@section Panel spacing:
#'
#'  Variable spacing between the panels hierarchically groups them to identify
#'  their settings. The default spacing for both x and y directions is 0:9. This
#'  means that panels corresponding to the first, fastest changing, variable are
#'  separated by 0 units (= character heights); groups of panels at each fixed
#'  level of the second next fastest changing variable are separated by 1 unit;
#'  groups of groups of panels at fixed levels of the third are separated by 2
#'  units; and so forth.
#'
#'  For the example, this means that the row spacing would look like ('X'
#'  indicates a panel): XX  XX . And for columns it would be going down: XXX XXX
#'  XXX  XXX . The spacings can be different for x and y, but this is usually
#'  unnecessary.
#'
#'@section Effective xyLayout specification:
#'
#'  The default layout is often enhanced by changing the order of the factors
#'  and the xyLayout; for example, ordering the factors from those with the
#'  least change among levels to the most, or vice-versa; or setting the "most
#'  important" factors along rows to facilitate visual comparison.
#'
#'  The order of the variables -- and hence which vary in the x or y direction
#'  -- is given both by the left to right order of the conditioning in the
#'  formula and the xyLayout argument. Thus, in the example, conditioning with
#'  \code{~|f3*f1*f2*f4} and setting the layout with \code{xyLayout = list(x=4,
#'  y=3:1)} is equivalent to \code{~|f4*f2*f1*f3} and \code{xyLayout =
#'  list(x=1,y=2:4)} and produces a display with 12 rows and 4 columns in which
#'  the row panels now correspond to the 4 f4 levels and the column panels to
#'  the levels of f2, f1, and f3. This redundancy is deliberate: it allows
#'  changing layouts via the xyLayout argument to avoid rewriting a formula with
#'  long names in a different order.
#'
#'@section 2 level designs with a center point -- the 'center' argument:
#'
#'  Finally, the 'center' argument, a logical with default = \code{TRUE},
#'  controls the display when the conditioning factors are arranged as a 2 level
#'  factorial design with a single "pseudo"-center point. "pseudo" here means
#'  that the settings of the factors at the center need not be exactly in the
#'  middle for numeric factors. If the design is not of this form, the 'center'
#'  argument is ignored.
#'
#'  For such designs, when \code{center = TRUE},  a more compact display will be
#'  drawn in which a center panel corresponding to the center point is shown as
#'  the single panel in its row and column, but all other empty panels
#'  corresponding to settings where some of the conditioning variables are at
#'  their mid levels and some are not, are omitted. Examples are given below.
#'
#'  If it is desired to show all the empty panels, which can be useful to
#'  informatively represent the actual design sparsity, set  \code{center = FALSE}.
#'
#'@param obj Argument determining method dispatch. For the formula method, a
#'  \code{xyplot} type formula of the form \code{~y|f1*f2*...*fn} or
#'  \code{y~x|f1*f2*...*fn} (where "..." means the actual variable names).
#'  Instead of explicitly specifying the conditioning variables,
#'  i.e. the variables after the |, you can instead use a "." after |. This is
#'  interpreted to mean "all variables in the data argument \emph{except}
#'  those to the left of the |". For example, the second formula above could be
#'  written as \code{y~x|.} where f1, f2, ..., fn and possibly x and y are the
#'  variables in data (typically columns of a data frame).
#'
#'  Note 1: Extended formulas and 3-d formulas (see
#'  \code{\link[lattice]{xyplot}}) are not implemented.
#'
#'  Note 2: For the lm method, the model object should contain a model
#'  component. See \code{\link[stats]{predict.lm}} and Help pages for predict
#'  methods for classes inheriting from "lm" for other arguments and for
#'  exactly what is plotted.
#'
#'  It is preferable that this be the the first argument of the call.
#'
#'@param data For the formula method, a data frame (or more precisely, anything
#'  that is a valid envir argument in eval, e.g., a list or an environment)
#'  containing values for any variables in the formula, as well as groups and
#'  subset if applicable. If not found in data, or if data is unspecified, the
#'  variables are looked for in the environment of the formula.
#'
#'@param xyLayout In its most general form, a list with named "x" and "y"
#'  components, either or both of which can be missing (or NULL). If there are
#'  \emph{n} conditioning factors and both x and y are given, then
#'  the combined components (e.g. via \code{unlist}) must be a permutation of
#'  the sequence 1, 2,\ldots, n (with no duplicates).
#'
#'  The integers in x specify the indices of the conditioning factors and their
#'  levels in the x direction hierarchy. Correspondingly for y. For a fuller
#'  explanation of how this controls the plot layout, see the \strong{Overview}
#'  section and the examples.
#'
#'  For convenience, the xyLayout argument can be specified in several other
#'  ways. The basic idea is that the full list will be \emph{filled in}
#'  "appropriately" if possible. Specifically, this means:
#'
#'  \itemize{
#'  \item If both components are missing, NULL, or of zero length the list is
#'  constructed by splitting the conditioning factors equally in the horizontal
#'  and vertical directions, with horizontal gettting one more if the number of
#'  conditioning variables is odd;
#'   \item If only one component of the list is given, the missing component gets
#'  the remaining factors, if any, in order(note that this means if both
#'  directions are not in order low to high, then both must be explicitly
#'  specified in the xyLayout list);
#'  \item If component names are missing or empty(i.e. ""), the first one is
#'  "x" and the second is "y". Nonmissing component names must be "x" and/or "y"
#'  and must be unique;
#'  \item A 1 or 2 column matrix can be used instead of a list with an even
#'  number of conditioning variables. If the column names are "x" or "y", they
#'  will be used. Otherwise, the x component is first and y second;
#'  \item An unnamed vector can be given in place of a list with a single unnamed
#'  component (which would be assumed to be the "x" component)
#'  }
#'
#'  As a result, all the following would produce exactly the same layout for
#'  n = 5 conditioning factors when used as the xyLayout argument:
#'  \itemize{
#'    \item \code{list(x = 1:3, y = 4:5)}
#'    \item \code{list (x = 1:3)}
#'    \item \code{list (y = 4:5)}
#'    \item argument missing or NULL
#'    \item \code{1:3}
#'    \item \code{cbind(y = 4:5)}
#'    \item \code{matrix(1:3,ncol=1)}
#'  }
#'  Note that arbitrary integer vectors can be used, not just
#'  (ordered) sequences; e.g. the following are correct and equivalent when
#'  there are 6 conditioning factors:
#'
#'  \itemize{
#'    \item \code{list(y = c(1,5,3), x = c(2,4,6))}
#'    \item \code{cbind(y=c(1,5,3), x=c(2,4,6))}
#'    \item \code{cbind(y = c(1,5,3))}
#'  }
#'  \emph{But be careful}: \code{xyLayout = c(2,4,6)} is equivalent to
#'  \code{xyLayout = list(x = c(2,4,6), y = c(1,3,5))}, which is different than the
#'  previous, since order matters!
#'
#'
#'@note Because xyLayout and the number of levels in the conditioning variables
#'  determine the plot structure, a 'layout' argument in the call will be
#'  ignored.
#'
#'@seealso \code{\link[lattice]{lattice}}
#'
#'@param spacings A \strong{list} with x and y components that are nondecreasing
#'  integer sequences. The ith value in the sequence gives the spacing in
#'  character heights between the sets of panels at each level for the ith
#'  conditioning variable in the component's direction. See the
#'  \strong{Overview} for further explanation.
#'
#'@param center Logical, default TRUE. If the conditioning factors constitute a
#'  2-level design with a center point and center is TRUE, a compact trellis
#'  plot that omits panels for the missing settings of the factors is drawn. If
#'  FALSE or the design is not of this form, all panels are shown, which may be
#'  informative in any case to better visualize the design sparsity.
#'
#'@param newdata For the \code{lm} method (and objects inheriting from the
#'  \code{lm} class). The 'newdata' argument for the associated \code{predict}
#'  method. Consult \code{\link[stats]{predict.lm}} for details.
#'
#'@param ylab For the \code{lm} method. Optional y axis label for predicted values.
#'
#'@param predictArgs For the \code{lm} method. A \strong{list} of optional named
#'  arguments for the relevant \code{predict} method.
#'
#'@param col Fill color for the data frame method.
#'
#'@param \dots Further arguments to panel functions, methods, or
#'\code{\link[lattice]{xyplot}}.
#'
#'@return The formula method returns an object of type c("structured",
#'  "trellis); or c("doe","structured","trellis") for 2-level designs with a
#'  center point. This is a trellis object with additional \code{structure} and
#'  \code{formula} attributes. Other methods may extend the class to allow
#'  special print/plot methods.
#'
#'  The \code{structure} attribute provides xyLayout information. It is a list
#'  with \code{x} and \code{y} components, each of which in turn is a list, one
#'  of which may be empty (i.e. of length 0). The components of the x list have
#'  names of the conditioning variables in the horizontal direction of the
#'  layout with values the variable values; analogously for the y list.
#'
#'  The \code{formula} attribute gives the actual formula used in the trellis
#'  call, i.e. with the actual conditioning variable names substituted for ".".
#'
#' @example man/Examples/strucplot_exmpl.R
#
strucplot<- function(obj,...)UseMethod("strucplot")
#
#' @describeIn strucplot Formula Method
strucplot.formula <- function(
  obj
  ## A formula of the kind used in xyplot
  ,data = list()
  ## as in \code{xyplot}
  ,xyLayout = list()
  ## A list with named x and/or y integer vector components that when
  ## concatenated give the order in which the conditioning variables are used in
  ## the display. That is c(x,y) will be used as the perm.cond argument for the
  ## xyplot call. The x and y components indicate which variables change in the
  ## x (horizontal) and y (vertical) directions, respectively. A missing x
  ## component means only 1 row; a missing y component, only 1 column.
  ,spacings = list(x= 0:9,y= 0:9)
  ##  list with x and y components that are nondecreasing integer sequences
##  giving x and y spacings vectors between the panels at the various hierachies
  ,center=FALSE
 ## if TRUE,only for 2^n factorial designs with a center point. The display will
 ## omit all the empty panels corresponding to the middle levels of the variables
 ## except for the center point panel in which all conditioning variables are at
 ## their middle levels. Otherwise will be ignored.
  ,...
  ## further arguments
)
{
  ## create data argument containing the variables to evaluate obj
  spf <- strucParseFormula(obj,data)
  obj <- attr(spf,"form")
   dots <- list(...)
   if(identical(dots[["outer"]],TRUE))stop("'outer = TRUE' option not permitted for strucplot plots")
   d <- spf$condition
   lend <- length(d)
   k <- ceiling(lend/2)
   xy <- c(x= "x", y = "y")
  if(is.null(dots[["horizontal"]])||
     !is.logical(dots[["horizontal"]])||
     length(dots[["horizontal"]])>1) dots[["horizontal"]] <- FALSE
  if(is.null(spf[["left"]])) {
    ## No LHS, so must modify formula and create one to satisfy xyplot
    nm <-paste(sample(c(LETTERS,letters),12),collapse="")
    data[[nm]]<- rep(0,length(spf$right))
    if(dots[["horizontal"]]){
      obj[[3]] <- obj[[2]]
      obj[[2]] <- obj[[3:2]]
      obj[[2]]<- as.name(nm)
      dots[["ylab"]] <- ""
      dots[["scales"]][["y"]] <- list(draw = FALSE)
    } else {
      obj[[3]] <- obj[[2]]
      obj[[c(3,2)]]<- as.name(nm)
      obj[[2]] <- obj[[c(2,2)]]
      dots[["xlab"]] <- ""
      dots[["scales"]][["x"]] <- list(draw=FALSE)
    }
  }
  xyLayout <- xyLayout(xyLayout, n = lend)
  levLen <-tryCatch(lapply(xyLayout,
       function(w){
        if(identical(w,integer(0)))1
        else sapply(d[w],function(x)length(levels(x)))
      }),
    error=function(e)stop("Can't compute number of levels in conditioning factors"))

  ## Check that spacings are nondecreasing
  if(!all(sapply(spacings,function(x)all(diff(x)>=0))))
    stop("spacings must be nondecreasing",call.=FALSE)
  if(!center || !check2lvl(d)){
    center <- FALSE ## to set class attribute of return
    dots[["layout"]] <- sapply(levLen,prod)
    dots[["perm.cond"]] <- unlist(xyLayout)
    dots[["skip"]] <- FALSE
    dots[["between"]] <- lapply(c(x="x",y="y"),function(i){
      if(identical(xyLayout[[i]],integer(0)))0
      else makeSpacings(spacings[[i]],levLen[[i]])
    }
    )
  } else  tryCatch(
      {  ## 2^n design with center point
   ## must create an artificial factor of right form to index the panels
    nm <-paste(sample(c(LETTERS,letters),12),collapse="")
    obj[[c(3,3)]] <- as.name(nm)
    dord <- d[unlist(xyLayout)]
    data[[nm]] <- factor(interaction(dord,drop=TRUE),
      levels = make_2level_with_center(lapply(dord,levels)))
    nfac <- sapply(xyLayout,length)
    dots[["layout"]] <- 2^nfac + (nfac>0)
    dots[["perm.cond"]] <- 1
    dots[["skip"]] <-if(any(nfac <1))FALSE else {
      half <- 2^(nfac-1)
      rowskip<- rep(c(FALSE,TRUE,FALSE),times=c(half[1],1,half[1]))
      unlist(rep(list(rowskip,!rowskip,rowskip),times=c(half[2],1,half[2])))
    }
    dots[["between"]] <-lapply(xy,function(i){
        k <- nfac[i]
        if(!k)0 else {
          if(k==1)rep(spacings[[i]][1],2) else {
            space <- rep(spacings[[i]],length.out=k)
            sp <- makeSpacings(space,rep(2,k-1))
            c(sp,rep(space[k],2),sp)
          }
        }
      })
  }, ##tryCatch expression end
   error=function(e)stop("Cannot display 2 level factorial design with center point"))
if(is.null(dots[["as.table"]])) dots[["as.table"]]  <- TRUE
dots[["strip"]] <- FALSE
dots[["drop.unused.levels"]] <- FALSE
## (To override possible user provided values)
  ## return trellis object for plotting
structure(do.call(xyplot,c(list(obj,data=data), dots)),
  class=c(if(center)"doe","structured","trellis","list"),
  formula = obj,
  structure = lapply(xyLayout,function(i)lapply(d[i],levels)) )
}
#
#' @describeIn strucplot Default method prints an error message
strucplot.default <- function(...)
stop(paste("No strucplot method",
           "Use 'methods(strucplot)' to find available strucplot methods",
           sep = "\n"))

#' @describeIn strucplot Data frame method with no response to show design structure
strucplot.data.frame <- function(
  obj ## design data frame
  ,col = "darkblue" ## color for filling the panels
  ,...
  ## as in strucplot.formula + arguments for panel.fill
  )

{
  y <- rep(0,nrow(obj))
  dat <- obj
  obj <- formula(paste0("~y|",paste0(names(obj),collapse="*")))
  dots <- list(...)
#   strucplot.formula(obj,data=dat, ...,
#    scales=list(draw=FALSE),
#    ylab= "" ,
#    col=col,
#    panel=function(y,col,...)if(length(y))lattice::panel.fill(col=col,...) else NA
#  )
  dots$panel <- function(y,col,...)if(length(y))lattice::panel.fill(col=col,...)
              else NA
  dots$scales <- list(draw = FALSE)
  dots$ylab <- ""
  dots$col <- col
  do.call(strucplot,c(list(obj=obj,data=dat),dots))
}
#
#'@describeIn strucplot Converts to a data frame and calls the data.frame method
strucplot.list <- function(obj,...)
{
  strucplot(
    tryCatch(data.frame(obj)), error=function(e)
      stop("Unable to convert to data frame",call.= FALSE), ...)
}
#
#' @describeIn strucplot Converts to a data frame and calls the data.frame method
strucplot.matrix <- strucplot.list
#
#' @describeIn strucplot Plots predicted values for fitted models inheriting
#'   from class "lm"
strucplot.lm <- function(obj
  ,newdata = model.frame(obj)[-1]
  ## newdata argument (a data frame) for predict method
  ,ylab = "Predicted Response"
  ,predictArgs = NULL
  ## a **list** of optional named arguments for the predict method.
  ## Note that the panel function must be able to accommodate them
  ,...
  )
{
  nm <- deparse(substitute(obj))
  predicted <-
    tryCatch(do.call(predict,c(list(obj,newdata=newdata),predictArgs)),
    error = function(e)
      stop(paste("'predict' method failed. Was model frame saved in",nm,"?"),call. = FALSE))
  obj <- formula(paste("~predicted",paste(names(newdata),collapse="*"),sep="|"))
  out <- strucplot.formula(obj,data=newdata,ylab=ylab,...)
  class(out)<- c("modelFit",class(out))
  out
}

