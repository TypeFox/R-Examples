setClass("gTreetcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )

## map a list to a tree
setMethod(".gtree",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   offspring = NULL,
                   hasOffspring = NULL,                 # for defining offspring. FUN
                                        # of children only.
                   offspring.data = NULL,
                   col.types = NULL, # data frame with logical
                   icon.FUN = NULL,                      # return stock name --called
                                        # on offspring, returns a
                                        # vector of length nrow
                   chosencol = 1,
                   multiple = FALSE,
                   handler = NULL,
                   action=NULL,
                   container=NULL,
                   ...
                   ) {
            
            force(toolkit)
            theArgs <- list(...)
            
            if(is.null(offspring)) {
              message(gettext("Need to have specified an offspring function\n"))
              return(NA)
            }
            
            ## allow multiple if asked
            if(as.logical(multiple))
              selectmode <- "extended"
            else
              selectmode <- "browse"

            
            ## get base offspring
            os <- offspring(c(), offspring.data)
            if(!inherits(os,"data.frame"))
              os <- as.data.frame(os)
            
            ## icons
            icons <- rep("", nrow(os))
            if(!is.null(icon.FUN)) 
              icons <- icon.FUN(os)
            ## fix icons - allow for stock or file or "" or null or NA
            ## are icons "", NA, filename or stockname?
            icons <- sapply(icons,function(i) {
              findTkIcon(i)
            })


            
            ## need os, hasChild, icons, size: d, n, m
            l <- .treeGetOffspring(os, hasOffspring)
            os <- l$children
            whichHaveOffspring <- l$offspring
            d <- dim(os); m <- d[1]; n <- d[2]
            
            
            ## make tree view and 
            ## pack into scrolled window
            if(is.logical(container) && container)
              container <- gwindow(visible=TRUE)
            
            tt <- getWidget(container)
            gp <- ttkframe(tt, width=20*n) # default width

            xscr <- ttkscrollbar(gp, orient="horizontal",
                                 command=function(...)tkxview(tr,...))
            yscr <- ttkscrollbar(gp, 
                                 command=function(...)tkyview(tr,...))

            if(n >= 2) 
              tr <- ttktreeview(gp, columns = 2:n,
                                ## this works, but the above is
                                ##cleaner. It gives one extra column
                                ##when n = 1
                                ## columns = as.tclObj(columns, drop=FALSE),
                                ## but the following fails -- extra columns
                                ## columns = columns,
                                displaycolumns="#all",
                                selectmode=selectmode,
                                xscrollcommand=function(...)tkset(xscr,...),
                                yscrollcommand=function(...)tkset(yscr,...))
            else                        # no columns argument
              tr <- ttktreeview(gp,   
                                displaycolumns="#all",
                                selectmode=selectmode,
                                xscrollcommand=function(...)tkset(xscr,...),
                                yscrollcommand=function(...)tkset(yscr,...))

            tkgrid(tr,row=0,column=0, sticky="news")
            tkgrid(yscr,row=0,column=1, sticky="ns")
            tkgrid(xscr, row=1, column=0, sticky="ew")
            ## see tkFAQ 10.1 -- makes for automatic resizing
            tkgrid.columnconfigure(gp, 0, weight=1)
            tkgrid.rowconfigure(gp, 0, weight=1)
            tkpack(gp, expand=TRUE, fill="both")

            
            ## call in autoscroll if requested -- has issues with sizing
            if(getWithDefault(theArgs$do.autoscroll, TRUE) &&
               windowingsystem() != "aqua") {
              tcl("autoscroll::autoscroll", xscr)
              tcl("autoscroll::autoscroll", yscr)
            }

            ## turn on alternating shading if more than 1 column
            ## XXX 

            obj <- new("gTreetcltk", block=gp, widget=tr,
              ID=getNewID(),
              e = new.env(),
              toolkit=toolkit)

            tag(obj,"offspring") <- offspring
            tag(obj,"offspring.data") <- offspring.data
            tag(obj,"icon.FUN")  <-  icon.FUN

            ## warn if chosencol not = 1
            if(chosencol != 1) {
              gwCat(gettext("In gWidgetstcltk, chosencol is always the first column\n"))
            }
            tag(obj,"chosencol") <- chosencol
            tag(obj,"multiple") <- multiple
          
            ## put in children, 
            .treeAddOffspring(tr, parent="", os, whichHaveOffspring, icons)
            
            ## now add a handler to row-exapnd
            ## Tree view open handler. No need to delete on close, as we
            ## delete children on open.
            tkbind(tr, "<<TreeviewOpen>>",function(W, x,y) {
                 sel <- unlist(strsplit(tclvalue(tcl(W,"selection"))," "))[1]
                 ## check if  children, if not return
                 children <- unlist(strsplit(tclvalue(tcl(W,"children",sel))," "))
                 if(length(children) == 0)
                   return()
                 lapply(children, function(i) tcl(W,"delete",i))
                 ## add in children
                 path <- .treeGetPath(W)
                 
                 os <- offspring(path, tag(obj, "offspring.data"))

                 ## icons
                 icons <- rep("", nrow(os))
                 if(!is.null(icon.FUN)) 
                   icons <- icon.FUN(os)
                 ## fix icons - allow for stock or file or "" or null or NA
                 ## are icons "", NA, filename or stockname?
                 icons <- sapply(icons,function(i) {
                   findTkIcon(i)
                 })
              
                l <- .treeGetOffspring(os, hasOffspring)
                 os <- l$children
                 whichHaveOffspring <- l$offspring
                 
                 .treeAddOffspring(W, parent=sel, os, whichHaveOffspring,
                                   icons=icons)
               })
            
            
            if(!is.null(handler)) {
              id <- addhandlerdoubleclick(obj,handler,action)
              tag(obj, "handler.id") <- id
            }
            
            ## attach to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj,...)
            }
            
            return(obj)
          })

## Some helper functions
.treeAddOffspring <- function(tr, parent="",os, hasChild, icons) {
  d <- dim(os); m <- d[1]; n <- d[2]
  if(m == 0) return()


  ## coerce to a character matrix, but worry about single row data frames
  os <- as.data.frame(lapply(os,as.character),stringsAsFactors=FALSE)
  os <- as.matrix(os)
  nms <- colnames(os)                   # os is a matrix now!


  ## widths and column headings
  if(parent == "" && n > 1) {
    ## column headings
    tcl(tr,"heading","#0", text="")
    lapply(2:n, function(j) {
      width <- max(nchar(c(nms[j],os[,j,drop=TRUE]))) * 8 + 15
      tcl(tr,"column",j, "-width",width) # was -2??
      tcl(tr,"heading",j,text=nms[j])
    })
  } else if(parent == "" && n == 1) {
    tcl(tr,"heading","#0", text=nms[1])
    tcl(tr,"column","#0", "-width", max(nchar(c(nms[1],os[,1,drop=TRUE]))) * 8 + 15)
#    tcl(tr,"column",0,"-width",0)
  }
  
  lapply(1:m, function(i) {
    if(n > 1) {
      values <- os[i,2:n,drop=TRUE]
    } else {
      values<- NA
    }
    icon <- tcl("image","create","photo",file=icons[i])
    ## need to use 16 by 16 icons
    
    if(length(values) > 1 || !is.na(values)) {
      tpath <- tcl(tr,"insert",parent,"end",text=os[i,1],image=icon,
                 values = values)
    } else {                            # no values
      tpath <- tcl(tr,"insert",parent,"end",text=os[i,1],image=icon)
    }
    ## add offspring
    if(!is.na(hasChild[i]) && as.logical(hasChild[i])) {
      tcl(tr,"insert",tclvalue(tpath),"end", text="")
    }
  })
  
  
}

## for svalue
.treeGetPath <- function(tr) {
  ## return path from selection
  sel <- unlist(strsplit(tclvalue(tcl(tr,"selection"))," "))

  if(is.null(sel)) return(c())
  path <- tclvalue(tcl(tr,"item",sel,"-text"))
  parent <- tclvalue(tcl(tr,"parent",sel))
  while(parent != "") {
    path <- c(tclvalue(tcl(tr,"item",parent,"-text")),path)
    parent <- tclvalue(tcl(tr,"parent",parent))
  }

  return(path)
}

.treeGetSelectedValue <- function(tr) {
  ## return path from selection
  sel <- unlist(strsplit(tclvalue(tcl(tr,"selection"))," "))
  if(length(sel) == 0) return(NA)
  sapply(sel, function(i) tclvalue(tcl(tr,"item",i,"-text")))
}

## Take the data frame and massage it to return
## icons if asked, and figure out offspring
.treeGetOffspring <- function(children, hasOffspring) {

  ## do we expand?
  ## how to determine if offspring are needed?
  ## default to hasOffspring, then second column, then default to FALSE
  if(!is.null(hasOffspring)) {
    offspring <- hasOffspring(children)
  } else {
    ## if second column is logical, we use that
    if(dim(children)[2] > 2 && is.logical(children[,2])) {
      offspring <- children[,2]
      children <- children[,-2, drop=FALSE]
    } else {
      offspring <- rep(FALSE, nrow(children))
    }
  }
  
  return(list(children=children, offspring=offspring))
}


## number of columns 
.treeNoColumns <- function(tr) {
  colnames <- tclvalue(tcl(tr,"cget","-columns"))
  colnames <- unlist(strsplit(colnames," "))
  length(colnames) + 1 # 1 for key
}

## width of columns
.treeColWidths <- function(tr) {
  n <- .treeNoColumns(tr)
  widths <- sapply(2:n , function(j) {
    tclvalue(tcl(tr, "column", j - 2, "-width"))
  })
  widths <- c(tclvalue(tcl(tr, "column", "#0", "-width")), widths)
  return(as.numeric(widths))
}


## help out with gtree
.treeByReturnVector = function(df, FUN,...) {
  tmp = by(df, factor(1:nrow(df)), FUN)
  sapply(tmp, function(x) x)
}

## has different arguments, but we mask this with ...
## this has offspringdata as first argument
setMethod(".update",
          signature(toolkit="guiWidgetsToolkittcltk",object="gTreetcltk"),
          function(object,toolkit,...) {
            obj <- object     # rename, object from update generic
            tr <- getWidget(obj)
            
            theArgs <- list(...)
            offspring <- tag(obj,"offspring")
            hasOffspring <- tag(obj,"hasOffspring")
            icon.FUN <- tag(obj,"icon.FUN")
            
            offspring.data <- theArgs$offspring.data
            if(is.null(offspring.data)) {
              if(length(theArgs) > 1)
                offspring.data <- theArgs[[1]]
              else
                offspring.data <- NULL
            }
            ## what should now be in this part of the tree
            os <- offspring(c(), offspring.data)

            ## icons
            icons <- rep("", nrow(os))
            if(!is.null(icon.FUN)) 
              icons <- icon.FUN(os)
            ## fix icons - allow for stock or file or "" or null or NA
            ## are icons "", NA, filename or stockname?
            icons <- sapply(icons,function(i) {
              findTkIcon(i)
            })

            
            l <- .treeGetOffspring(os, hasOffspring)
            os <- l$children
            whichHaveOffspring <- l$offspring
            d <- dim(os); m <- d[1]; n <- d[2]
            

            ## delete what is there
            children <- unlist(strsplit(tclvalue(tcl(tr,"children",""))," "))
            if(length(children) > 0)
              lapply(children, function(i) tcl(tr,"delete",i))

            ## add children
            .treeAddOffspring(tr, parent="", os, whichHaveOffspring,
                              icons=icons)
            invisible()
          })

## index returns the indices
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTreetcltk"),
          function(obj, toolkit, index=NULL, drop=NULL,...) {

            tr <- getWidget(obj)
            
            index <- getWithDefault(index, FALSE)


            if(index) {
              sel <- unlist(strsplit(tclvalue(tcl(tr,"selection"))," "))
              if(length(sel) == 0)
                return(NULL)
              vals <- lapply(sel, function(i) {
                ind <- numeric(0)
                parent <- i
                while(parent != "") {
                  ind <- c(as.numeric(tcl(tr, "index", parent)) + 1, ind)
                  parent <- tclvalue(tcl(tr, "parent", parent))
                }
                ind
              })
              if(length(vals) == 1)
                vals <- vals[[1]]
              return(vals)
            } else {
              ## give key
              whichCol <- tag(obj,"chosencol")
              if(whichCol != 1)
                gwCat(gettext(" svalue only returns first column with gWidgetstcltk\n"))
              
              
              ## possible multiple selection
              ## return path from selection
              sel <- unlist(strsplit(tclvalue(tcl(tr,"selection"))," "))
              if(length(sel) == 0) return(NULL)
              vals <- sapply(sel, function(i) tclvalue(tcl(tr,"item",i,"-text")))
              ## strip off names
              attr(vals,"names") <- rep(NULL,length(vals))
              
              return(vals)
            }
          })

##' set selection by index
##' 
##' @param value a numeric path, or list of numeric paths specifying the selection
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gTreetcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {

                   index <- getWithDefault(index, TRUE)
                   if(!index) {
                     gwCat(gettext("Need to have index=TRUE (or NULL)"))
                     return(obj)
                   }

                   if(is.atomic(value))
                     value <- list(value)
                   ## 0-based
                   value <- lapply(value, function(i) i)

                   tr <- getWidget(obj)
                   selected <- as.character(tcl(tr, "selection"))
                   lapply(selected, function(sel) tcl(tr, "selection", "toggle", sel))
                   lapply(value, function(path) {
                     parent <- ""
                     for(i in path) {
                       parent <- as.character(tcl(tr, "children", parent))[i]
                       tcl(tr, "see", parent)
                     }
                     tcl(tr, "selection", "add", parent)
                   })

                   return(obj)

                                 
                 })
### need to figure this out
## return the path in values
setMethod("[",
          signature(x="gTreetcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, guiToolkit("tcltk"), i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gTreetcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            tr <- getWidget(x)

            ## return path from selection
            sel <- unlist(strsplit(tclvalue(tcl(tr,"selection"))," "))

            if(is.null(sel)) return(c())
            path <- tclvalue(tcl(tr,"item",sel,"-text"))
            parent <- tclvalue(tcl(tr,"parent",sel))
            while(parent != "") {
              path <- c(tclvalue(tcl(tr,"item",parent,"-text")),path)
              parent <- tclvalue(tcl(tr,"parent",parent))
            }
            

            if(missing(i))
              return(path)
            else
              return(path[i])
          })

### methods
## row-activated in gtable gives double click
setMethod(".addhandlerdoubleclick",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gTreetcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            addhandler(obj, "<Double-1>",handler,action)
           })


setReplaceMethod(".size", 
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gTreetcltk"),
                 function(obj, toolkit, ..., value) {
                   tr <- getWidget(obj)

                   n <- .treeNoColumns(tr)
                   widths <- .treeColWidths(tr)
                   curWidth <- sum(widths)
                   widths <- floor((1+widths) * value[1]/curWidth)
                   
                   ## set width
                   sapply(2:n, function(j) {
                     tcl(tr, "column", j - 2, width=widths[j], stretch=TRUE, anchor="w")
                   })
                   tcl(tr,"column","#0","-width", widths[1])
                   
                   ## set height
                   height <- value[2]
                   tcl(tr,"configure", height = floor(height/16))
                   
                   return(obj)
                 })

