#' @title Easy Barcharts

#' @description Wrapper for \code{\link{barchart}} in package \code{lattice}.  Creates a 
#' barchart from raw data using formula-data syntax similar to that of \code{\link{xtabs}},
#' or from a table.  Defaults to a "standard"
#' barchart in which the bars are vertical and unstacked.  Supports percentage barcharts.
#' 
#' @rdname barchartGC
#' @usage barchartGC(x,data=parent.frame(),type="frequency",flat=FALSE,auto.key=TRUE,
#'                        horizontal=FALSE,stack=FALSE,...)
#' @param x Either a formula or an object that can be coerced to a table.  If formula, it must be 
#' of the form ~var or ~var1+var2.
#' @param data Usually a data frame that supplies the variables in \code{x}.  Variables not in the data
#' argument are searched for in the parent environment.
#' @param type Possible values are "frequency" and "percent".
#' @param flat If set to TRUE, will produce barchart that resembles the layout of \code{xtabs}
#' @param auto.key Provides a simple key
#' @param horizontal Determines orientation of the bars (overrriden by flat)
#' @param stack Determines whether bars for tallies are stacked on eac other or placed
#' next to one another (overrriden by flat)
#' @param ... other arguments passed to \code{barchart}:  these include main, sub, and
#' xlab, which are likely to be familiar to students from other \code{lattice} graphical
#' functions.  An error is possible if other arguments
#' pertaining to legends are passed (hopefully anyone interested in such will have moved on
#' to \code{barchart}).
#' @return A trellis object describing the barchart.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #barchart of counts for one factor variable:
#' barchartGC(~sex,data=m111survey)
#' 
#' #barchart with percentages and title:
#' barchartGC(~sex,data=m111survey,
#'    type="percent",
#'    main="Distribution of Sex")
#' 
#' #barchart of counts, to study the relationship between
#' #two factor variables:
#' barchartGC(~sex+seat,data=m111survey)
#' 
#' #percentage barchart, two factor variables:
#' barchartGC(~sex+seat,data=m111survey,type="percent")
#' 
#' #From tabulated data:
#' sexseat <- xtabs(~sex+seat,data=m111survey)
#' barchartGC(sexseat,type="percent",main="Sex and Seating Preference")
#' 
#' #from tabulated data:
#' dieTosses <- c(one=8,two=18,three=11,four=7,five=9,six=7)
#' barchartGC(dieTosses,main="60 Rolls of a Die")
#' 
#' # a "flat" barchart, pictorial version of xtabs() 
#' barchartGC(~sex+seat,data=m111survey,flat=TRUE,ylab="Sex")
#' 
#' # a "flat" barchart, pictorial version of xtabs() 
#' barchartGC(~sex+seat,data=m111survey,type="percent",flat=TRUE,ylab="Sex")
barchartGC <-
  function(x,data=parent.frame(),
           type="frequency",flat=FALSE,auto.key=TRUE,horizontal=FALSE,stack=FALSE,...)  {
    
    
    levelTol <- 5 #above this number, set legend vertically to the right
    
    #handle incorrect type specifications
    if (type %in% c("frequency","count","counts",
                    "freq","Freq","Counts","Count","Frequency",
                    "fre","fr","f","Fre","Fr","F",
                    "Coun","Cou","Co","C",
                    "coun","cou","co","c")) {
      type <- "frequency"
    }
    
    if (type %in% c("percentage","percent","perc","%",
                    "Percentage","Percent","Perc",
                    "per","per","p","Per","Pe","P")) {
      type <- "percent"
    }
    
    # handle arugments when user wants a barchart that looks like xtabs()
    if (flat==TRUE) {
      stack <- TRUE
      horizontal <- TRUE
    }
    
    
    #handle user-specified x or y axis labels
    
    ellipses <- list(...)
    
    if (type=="frequency" && horizontal==FALSE) {
      if (is.null(ellipses$ylab)) {
        ellipses$ylab <- "frequency"
      }
    }
      
      if (type=="percent" && horizontal==FALSE) {
        if (is.null(ellipses$ylab)) {
          ellipses$ylab <- "percent"
        }
      }
    
    if (type=="frequency" && horizontal==TRUE) {
      if (is.null(ellipses$xlab)) {
        ellipses$xlab <- "frequency"
      }
    }
    
    if (type=="percent" && horizontal==TRUE) {
      if (is.null(ellipses$xlab)) {
        ellipses$xlab <- "percent"
      }
    }
    
    
    if (is(x,"formula"))  { #we have a formula
      prsd <- ParseFormula(x)
      pullout <- as.character(prsd$rhs)
      
      if (length(pullout) == 1) {  #one variable
        varname <- pullout[1]
        variable <- simpleFind(varName=varname,data=data)
        tab <- xtabs(~variable) 
        
        #set up the "other arguments" to be passed to barchart
        otherArgs <- c(list(stack=stack,horizontal=horizontal,
                            auto.key=auto.key),
                       ellipses)
        
        
        if (type=="frequency") {
        args <- c(list(x=tab),otherArgs)
        return(do.call(lattice::barchart,args))
        }
      
        if (type=="percent") {
        perctab <- 100*tab/sum(tab)
        perctab[is.nan(perctab)] <- 0
        args <- c(list(x=perctab),otherArgs)
        return(do.call(lattice::barchart,args))
          }
      
      } # end one variable
      
      if (length(pullout)==3) { #two variables
        expname <- pullout[2]
        respname <- pullout[3]
        explanatory <- simpleFind(varName=expname,data=data)
        response <- simpleFind(varName=respname,data=data)
        tab <- table(explanatory,response)
        
        #make groups appear in order like xtabs
        if (horizontal==TRUE) {
          tab <- tab[nrow(tab):1,]
        }
        
        #make legend work well with orientation of chart, if possible.  If
        #there are lots of levels in exp variable, though, then put legend to the right
        
        
        if (stack==FALSE && horizontal==TRUE && auto.key==TRUE) {
          
            space <- ifelse(ncol(tab) <= levelTol,"top","right")
            key=simpleKeyRev(columns=1,space=space,
                  text=colnames(tab),
                  rectangles=TRUE,
                  points=FALSE)
            
            otherArgs <- c(list(stack=stack,horizontal=horizontal,
                                key=key),
                           ellipses) 
            
            tab <- tab[,rev(colnames(tab))]
            
            if (type=="frequency") {
              args <- c(list(x=tab),otherArgs)
              return(do.call(lattice::barchart,args))
            }
            if (type=="percent") {
              perctab <- 100*prop.table(tab,margin=1)
              perctab[is.nan(perctab)] <- 0
              args <- c(list(x=perctab),otherArgs)
              return(do.call(lattice::barchart,args))
            }   

        }
        
        
        if (stack==TRUE && horizontal==FALSE && auto.key==TRUE) {
          
          space <- ifelse(ncol(tab) <= levelTol,"top","right")
          key=simpleKeyRev(columns=1,space=space,
                           text=colnames(tab),
                           rectangles=TRUE,
                           points=FALSE)
          
          otherArgs <- c(list(stack=stack,horizontal=horizontal,
                              key=key),
                         ellipses) 
          
          tab <- tab[,rev(colnames(tab))]
          
          if (type=="frequency") {
            args <- c(list(x=tab),otherArgs)
            return(do.call(lattice::barchart,args))
          }
          if (type=="percent") {
            perctab <- 100*prop.table(tab,margin=1)
            perctab[is.nan(perctab)] <- 0
            args <- c(list(x=perctab),otherArgs)
            return(do.call(lattice::barchart,args))
          }   
          
        }
        
        
        if (stack==FALSE && horizontal==FALSE && auto.key==TRUE) {
          if (ncol(tab) <= levelTol) {
            auto.key=list(space="top",columns=ncol(tab))
          } else auto.key=list(space="right",columns=1)
        }
        
        if (stack==TRUE && horizontal==TRUE && auto.key==TRUE) {
          if (ncol(tab) <= levelTol) {
            auto.key=list(space="top",columns=ncol(tab))
          } else auto.key=list(space="right",columns=1)
        }

        
        #set up the "other arguments" to be passed to barchart
        otherArgs <- c(list(stack=stack,horizontal=horizontal,
                            auto.key=auto.key),
                       ellipses)
        
        if (type=="frequency") {
        args <- c(list(x=tab),otherArgs)
        return(do.call(lattice::barchart,args))
        }
        if (type=="percent") {
          perctab <- 100*prop.table(tab,margin=1)
          perctab[is.nan(perctab)] <- 0
          args <- c(list(x=perctab),otherArgs)
          return(do.call(lattice::barchart,args))
        }   
      }
      
      
    } #end check for formula
      
    
    if (!is(x,"formula")) {  #we have tabular data
      x <- as.table(x)
      if (length(dim(x))==1) {#one variable
        
        #set up the "other arguments" to be passed to barchart
        otherArgs <- c(list(stack=stack,horizontal=horizontal,
                            auto.key=auto.key),
                       ellipses)
        
        if (type=="frequency") {
          args <- c(list(x=x),otherArgs)
          return(do.call(lattice::barchart,args))
        }
        if (type=="percent") {
          perctab <- 100*x/sum(x)
          perctab[is.nan(perctab)] <- 0
          args <- c(list(x=perctab),otherArgs)
          return(do.call(lattice::barchart,args))
        }     
      }
      if (length(dim(x))>1) {#two variables
        
        #make groups appear in order like xtabs
        if (horizontal==TRUE) {
          x <- x[nrow(x):1,]
        }
        
        #make legend work well with orientation of chart, if possible.  If
        #there are lots of levels in exp variable, though, then put legend to the right
        
        if (stack==FALSE && horizontal==TRUE && auto.key==TRUE) {
          
          space <- ifelse(ncol(x) <= levelTol,"top","right")
          key=simpleKeyRev(columns=1,space=space,
                           text=colnames(x),
                           rectangles=TRUE,
                           points=FALSE)
          
          otherArgs <- c(list(stack=stack,horizontal=horizontal,
                              key=key),
                         ellipses) 
          
          x <- x[,rev(colnames(x))]
          
          if (type=="frequency") {
            args <- c(list(x=x),otherArgs)
            return(do.call(lattice::barchart,args))
          }
          if (type=="percent") {
            perctab <- 100*prop.table(x,margin=1)
            perctab[is.nan(perctab)] <- 0
            args <- c(list(x=perctab),otherArgs)
            return(do.call(lattice::barchart,args))
          }   
          
        }
        
        
        if (stack==TRUE && horizontal==FALSE && auto.key==TRUE) {
          
          space <- ifelse(ncol(x) <= levelTol,"top","right")
          key=simpleKeyRev(columns=1,space=space,
                           text=colnames(x),
                           rectangles=TRUE,
                           points=FALSE)
          
          otherArgs <- c(list(stack=stack,horizontal=horizontal,
                              key=key),
                         ellipses) 
          
          x <- x[,rev(colnames(x))]
          
          if (type=="frequency") {
            args <- c(list(x=x),otherArgs)
            return(do.call(lattice::barchart,args))
          }
          if (type=="percent") {
            perctab <- 100*prop.table(x,margin=1)
            perctab[is.nan(perctab)] <- 0
            args <- c(list(x=perctab),otherArgs)
            return(do.call(lattice::barchart,args))
          }   
          
        }
        
        
        if (stack==FALSE && horizontal==FALSE && auto.key==TRUE) {
          if (ncol(x) <= levelTol) {
            auto.key=list(space="top",columns=ncol(x))
          } else auto.key=list(space="right",columns=1)
        }
        
        if (stack==TRUE && horizontal==TRUE && auto.key==TRUE) {
          if (ncol(x) <= levelTol) {
            auto.key=list(space="top",columns=ncol(x))
          } else auto.key=list(space="right",columns=1)
        }
        
        #set up the "other arguments" to be passed to barchart
        otherArgs <- c(list(stack=stack,horizontal=horizontal,
                            auto.key=auto.key),
                       ellipses)
        
        if (type=="frequency") {
          args <- c(list(x=x),otherArgs)
          return(do.call(lattice::barchart,args))
        }
        if (type=="percent") {
          perctab <- 100*prop.table(x,margin=1)
          perctab[is.nan(perctab)] <- 0
          args <- c(list(x=perctab),otherArgs)
          return(do.call(lattice::barchart,args))
      }
      
      }
        
    } # end tabular processing
    
  } #end barchartGC

#' @title Reversed Simple Key Function
#' 
#' Utility function for barchartGC.
#' 
#' @rdname simpleKeyRev
#' @usage simpleKeyRev(text, points = TRUE,
#' rectangles = FALSE,
#' lines = FALSE,
#' col = add.text$col,
#' cex = add.text$cex,
#' alpha = add.text$alpha,
#' font = add.text$font,
#' fontface = add.text$fontface,
#' fontfamily = add.text$fontfamily,
#' lineheight = add.text$lineheight,
#' ...)
#' @export
#' @keywords internal
simpleKeyRev <-
  function(text, points = TRUE,
           rectangles = FALSE,
           lines = FALSE,
           col = add.text$col,
           cex = add.text$cex,
           alpha = add.text$alpha,
           font = add.text$font,
           fontface = add.text$fontface,
           fontfamily = add.text$fontfamily,
           lineheight = add.text$lineheight,
           ...)
  {
    add.text <- lattice::trellis.par.get("add.text")
    foo <- seq_along(text)
    ans <-
      list(text = list(lab = text),
           col = col, cex = cex, alpha = alpha,
           font = font,
           fontface = fontface,
           fontfamily = fontfamily,
           ...)
    if (points) ans$points <-
      lattice::Rows(lattice::trellis.par.get("superpose.symbol"), foo)
    if (rectangles) {                               #modification is in here
      temp <-lattice::Rows(lattice::trellis.par.get("superpose.polygon"), foo)
      temp$col <- rev(temp$col)
      ans$rectangle <- temp
    }
    
# this bit not needed: (problem with updateList (where is it from?))
#     if (lines) ans$lines <-
#       updateList(Rows(trellis.par.get("superpose.symbol"), foo), ## for pch
#                  Rows(trellis.par.get("superpose.line"), foo))

    ans
  } #end SimpleKeyRev

# for easy sourcing during development process
# simpleFind <- function(varName,data) {
#   tryCatch({get(varName,envir=as.environment(data))},
#            error=function(e) {
#              get(varName,inherits=T)
#            }
#   ) 
#}
