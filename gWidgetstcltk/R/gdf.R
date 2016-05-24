##################################################
### Gdf
## now we use the tablelist code




## gGrid cover gDf and gTable
setClass("gGridtcltk",
         contains="gComponenttcltk",
         prototype=prototype(new("gComponenttcltk"))
         )
setClass("gDftcltk",
         contains="gGridtcltk",
         prototype=prototype(new("gComponenttcltk"))
         )


## constructor for editing a data frame
setMethod(".gdf",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   items = NULL,
                   name = deparse(substitute(items)),
                   do.subset = FALSE,
                   container=NULL,...)  {

            force(toolkit)

            tt  <- getWidget(container)
            block <- ttkframe(tt)

            ##
            xscr <- ttkscrollbar(block, orient="horizontal",
                                 command=function(...) tkxview(widget,...))
            yscr <- ttkscrollbar(block, orient="vertical",
                                 command=function(...) tkyview(widget,...))
            

            widget <- tkwidget(block, "tablelist::tablelist",
                                resizablecolumns=1,
                                xscrollcommand=function(...) tkset(xscr,...),
                                yscrollcommand=function(...) tkset(yscr,...))
            
            tcl(widget, "configure", selecttype="cell")
            
            
            tkgrid(widget, row=0, column=0, sticky="news")
            tkgrid(yscr, row=0, column=1, sticky="ns")
                        tkgrid(xscr, row=1, column=0, sticky="ew")
            tkgrid.columnconfigure(block, 0, weight=1)
            tkgrid.rowconfigure(block, 0, weight=1)
            
            ##  tcl("autoscroll::autoscroll", xscr)
            ##  tcl("autoscroll::autoscroll", yscr)
            
            
            ## new object
            obj <- new("gDftcltk",block=block, widget=widget,
              toolkit=toolkit, ID=getNewID(), e = new.env())


            items <- as.data.frame(items)
            tl_configure_columns(widget, names(items))

            obj[] <- items
            tag(obj, "head") <- head(items, n=1)
            
            ## add to container
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj, ...)
            }

            return(obj)
            
          })


##
####################################################



## gWidget methods
setReplaceMethod(".size",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gGridtcltk"),
          function(obj, toolkit,  ..., value) {
            width = as.integer(value[1])
            height = as.integer(value[2])
            ## size
            tkconfigure(getWidget(obj),
                        maxwidth=width,
                        maxheight=height)
            
            return(obj)
          })



## data frame methods
## get selected value
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gGridtcltk"),
          function(obj, toolkit, index=NULL, drop=NULL,...) {
            message("svalue not implemented")
          })
          
          
## set by index value selected value
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",obj="gGridtcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                  message("svalue<- not implemented")
                 })


## refers to the entire data frame
## index returned by svalue(index=T) works here
setMethod("[",
          signature(x="gGridtcltk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j,..., drop=drop)
          })

setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {

            widget <- x@widget
            
            opar <- options("warn"); on.exit(options(opar))
            options(list(warn=-1)) # quiet for coerce_raw
            d <- dim(x)

            head <- tag(x, "head")
            l <- lapply(seq_len(d[2]), function(j) {
              coerce_raw(head[[j]], tl_get_column_raw(widget, j))
            })
            
            m <- structure(l,
                           .Names=tl_get_column_names(widget),
                           row.names=seq_len(d[1]),
                           class="data.frame")
            
            m[i,j, ...]

          })

## [<-
setReplaceMethod("[",
                 signature(x="gGridtcltk"),
                 function(x, i, j,..., value) {
                   .leftBracket(x, x@toolkit, i, j,...) <- value
                   return(x)
                 })


setReplaceMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x, toolkit, i, j, ..., value) {
            widget <- x@widget
            
            if(!missing(i) || !missing(j)) {
              tmp <- x[]
              tmp[i,j] <- value
              value <- tmp
            }

            tl_load_data(widget, value)
            
            return(x)
            
          })
                 
## data frame like
setMethod(".dim", 
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x,toolkit) {
            widget <- x@widget

            c(tl_no_rows(widget), tl_no_cols(widget)) 
          })

setMethod(".length",
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x,toolkit) return(dim(x)[2]))


## no dimnames for gGrid, only names
setMethod(".dimnames",
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x,toolkit) {
            tktable <- getWidget(x)

            toVector <- function(i) sapply(i, function(j) paste(j, collapse=" "))
            
            d <- dim(x)
            dimnames <- list(rownames=NULL,
                             colnames=names(x))
            dimnames
          })
          

setReplaceMethod(".dimnames",
                 signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
                 function(x, toolkit,  value) {
                   
                   if(!is.list(value))
                     stop("value is a list with first element the row names, and second the column names")
                   rnames = make.row.names(value[[1]])
                   cnames = value[[2]]
                   d = dim(x)
                   if(is.null(rnames) || length(rnames) != d[1])
                     stop("Row names are the wrong size")
                   if(is.null(cnames) || length(cnames) != (d[2]))
                     stop("Column names are the wrong size")

                   ## set column names
                   names(x) <- cnames

                   ## set row names
                   ## ignore
                          
                   return(x)
                 })


setMethod(".names",
          signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
          function(x, toolkit) {
            widget <- x@widget
            tl_get_column_names(widget)
          })


setReplaceMethod(".names",
                 signature(toolkit="guiWidgetsToolkittcltk",x="gGridtcltk"),
                 function(x, toolkit, value) {
                   tl_set_column_names(x@widget,value)
                   return(x)
                 })

## ##################################################

## .gDfaddPopupMenu <- function(obj) {
##   ## global variables to record row, column of menu popup
##   x0 <- NA; y0 <- NA
##   tktable <- getWidget(obj)
  
##   menu <- tkmenu(tktable)

##   insert <- function(x,ind, y) {
##     if(ind == 0) {
##       c(y,x)
##     } else if(ind >= length(x)) {
##       c(x,y)
##     } else {
##       c(x[1:(ind-1)], y, x[ind:length(x)])
##     }
##   }
  
##   ## return row, column of popup area
##   getWhere <-  function() {
##     where <- paste("@",x0,",",y0, sep="")
##     ind <- tcl(tktable,"index",where)
##     ind <- as.numeric(unlist(strsplit(as.character(ind),",")))
##     ind
##   }

##   ## GUI to write expression to evaluate to fill in column
##   transformVariable <- function(col) {
##     ## obj is main object
##     w <- gwindow(gettext("Transform variable"), parent=obj,width=300, height=400)
##     g <- ggroup(horizontal=FALSE, container = w)
##     glabel("To transform a variable you define a function body.", container = g)
##     glabel("You can use 'x' for the data frame and the column names.", container = g)
##     glabel("", container = g)
##     glabel("function(x) {", container=g)
##     glabel("\twith(x, {", container = g)
##     out <- gtext("", container = g); size(out) <- c(300,100)
##     glabel("\t})", container = g)
##     glabel("}", container = g)
## ##    gseparator(container =g, expand=TRUE)
##     bg <- ggroup(container = g)
##     cancelButton <- gbutton("cancel", container = bg, handler = function(h,...) dispose(w))
##     okButton <- gbutton("ok", container = bg, handler = function(h,...) {
##       str <- paste("x <- obj[,]",
##                    "f <- function(x) { with(x,{",
##                    svalue(out),
##                    "})}",
##                    "f(x)",
##                    sep="\n", collapse="\n")
##       val <- try(eval(parse(text=str)))
##       if(!inherits(val,"try-error")) {
##         obj[,col] <- val
##         dispose(w)
##       } else {
##         galert(gettext("Error in function body"), parent = w)
##       }
##     })
##     size(w) <- c(300,200)
##   }

##   columnEmpty <- function(col) {
##     val <- obj[,col] ## XXX write me
##     return(FALSE)
##   }
##   rowEmpty <- function(row) {
##     val <- obj[,row] ## XXX write me
##     return(FALSE)
##   }

##   ## confirm a delete
##   confirmDelete <- function(msg="Really delete? There is non empty data") {
##     out <- tkmessageBox(icon="question",
##                         message=gettext(msg),
##                         type="yesno",
##                         parent=tktable)
##     ifelse(as.character(out) == "yes",TRUE, FALSE)
##   }
  
##   formatColumn <- function(col, type) {
##     ## use tktable tag to format column to type.
##   }

##   ## make the menu
##   tkadd(menu,"command",label=gettext("Transform Variable"), command = function() {
##     ind <- getWhere()
##     transformVariable(ind[2])
##   })
##   tkadd(menu,"separator")
##   ##
##   tkadd(menu,"command",label=gettext("Insert Variable"), command = function() {
##     ind <- getWhere()
##     tcl(tktable,"insert", "cols", ind[2])
##     classes <- tag(obj, "classes")
##     tag(obj,"classes") <- insert(classes, ind[2]+1, "character")

##     val <- ginput("New variable name:", parent=obj)
##     if(!is.na(val))
##       names(obj)[ind[2] + 1] <- val
##   })
##   tkadd(menu,"command",label=gettext("Delete Variable"), command = function() {
##     ind <- getWhere()
##     if(columnEmpty(ind[2]) || confirmDelete())
##       tcl(tktable,"delete","cols",ind[2])
##     tag(obj, "classes") <- tag(obj, "classes")[-ind[2]]
##   })

##   tkadd(menu,"command",label=gettext("Rename Variable"), command = function() {
##     ind <- getWhere()
##     j <- ind[2]
##     oldName <- names(obj)[j]
##     val <- ginput("New variable name:", oldName, icon="question", parent=obj)
##     if(!is.na(val))
##       names(obj)[j] <- val
##   })

  
##   tkadd(menu,"command",label=gettext("Insert Case"), command = function() {
##     ind <- getWhere()
##     tcl(tktable,"insert","rows",ind[1])

##     val <- ginput("New case name:", parent=obj)
##     if(is.na(val)) 
##       val <- "NA"                       # fill in
##     rownames(obj)[ind[1] + 1] <- val

##   })
##   tkadd(menu,"command",label=gettext("Delete Case"), command = function() {
##     ind <- getWhere()
##     if(rowEmpty(ind[1]) || confirmDelete())
##       tcl(tktable,"delete","rows",ind[1])
##   })
##   tkadd(menu,"command",label=gettext("Rename case"), command = function() {
##     ind <- getWhere()
##     i <- ind[1] 
##     oldName <- rownames(obj)[i]
##     val <- ginput("New case name:", oldName, icon="question", parent=obj)
##     if(!is.na(val))
##       rownames(obj)[i] <- val
##   })

##   tkadd(menu,"separator")

##   setClass <- function(type) {
##     ind <- getWhere()
##     tclvalue(typeVar) <- type
##     classes <- tag(obj,"classes")
##     classes[ind[2]] <- type
##     tag(obj,"classes") <- classes
##     formatColumn(col=ind[2], type=type)
##   }
  
##   typeVar <- tclVar("numeric")          # for selecting type via radiobutton
##   tkadd(menu, "radiobutton", label="numeric", variable=typeVar, command=function() setClass("numeric"))
##   tkadd(menu, "radiobutton", label="integer", variable=typeVar, command=function() setClass("integer"))
##   tkadd(menu, "radiobutton", label="factor", variable=typeVar, command=function() setClass("factor"))
##   tkadd(menu, "radiobutton", label="character", variable=typeVar, command=function() setClass("character"))
##   tkadd(menu, "radiobutton", label="logical", variable=typeVar, command=function() setClass("logical"))
##   tkadd(menu, "radiobutton", label="other", variable=typeVar, command=function() {
##     ## need to popup dialog to get function name for other.
##     galert("other is not written", parent=obj)
##     setClass("character")
##   })
  
##   popupCommand <- function(x,y,X,Y) {
##     ## before popping up we have some work to do
##     x0 <<- x; y0 <<- y;
##     classMenuItems <- 7:12 + 2
##     ind <- getWhere() ## row, column
##     ## fix menu basd on where
##     tkentryconfigure(menu, 0, state=ifelse(ind[2]==0,"disabled","normal"))
##     tkentryconfigure(menu, 2, state=ifelse(ind[2]==0,"disabled","normal"))
##     tkentryconfigure(menu, 3, state=ifelse(ind[2]==0,"disabled","normal"))
##     tkentryconfigure(menu, 4, state=ifelse(ind[2]==0,"disabled","normal"))
##     tkentryconfigure(menu, 5, state=ifelse(ind[1]==0,"disabled","normal"))
##     tkentryconfigure(menu, 6, state=ifelse(ind[1]==0,"disabled","normal"))
##     tkentryconfigure(menu, 7, state=ifelse(ind[1]==0,"disabled","normal"))

##     for(i in classMenuItems)
##       tkentryconfigure(menu, i, state=ifelse(ind[2]==0,"disabled","normal"))

##     ## adjust class depends on which column
##     if(ind[2] == 0) {
##       tclvalue(typeVar) <- FALSE
##     } else {
##       theClass <- tag(obj,"classes")[ind[2]]
##       if(theClass %in% c("numeric","integer","character","factor","logical"))
##         tclvalue(typeVar) <- theClass
##       else
##         tclvalue(typeVar) <- "other"
##     }
##     ## popup
##     tkpopup(menu,X,Y)
##   }
##   ## mac binding, just 3 for all
##   if( as.character(tcl("tk","windowingsystem")) == "aqua" ) {
##     tkbind(tktable, "<2>", popupCommand)
##     tkbind(tktable, "<Control-1>", popupCommand)
##   }
##   tkbind(tktable, "<3>", popupCommand)
## }


## ## getFromIndex -- not using tcl array variable
## tktable.get <- function(tktable, i, j) {
##   val <- tkget(tktable, paste(i,j, sep=","))
##   as.character(val)
## }

## ## set From Index -- not using tcl array variable
## tktable.set <- function(tktable, i, j, value) 
##   tkset(tktable, paste(i, j, sep=","), as.character(value))



## ## take a data frame or matrix make a character matrix
## ## basically sapply(mat,format) but also has dimnames
## toCharacterMatrix <- function(x, rNames, cNames) {
##   mat <- as.data.frame(x, stringsAsFactors=FALSE)
##   mat <- as.data.frame(lapply(mat, format), stringsAsFactors=FALSE)
##   if(!missing(rNames)) 
##     mat <- cbind(rNames,mat)
##   mat[,1] <- as.character(mat[,1])
  
##   if(!missing(cNames)) 
##     mat <- rbind(c(rep("", !missing(rNames)), cNames), mat)
##   return(mat)
## }

## ## fill in a tclArray object from character matrix
## ## modifies ta in place -- passed through environment
## fillTclArrayFromCharMat <- function(ta, cm) {
##   ## cm[,1] contains column names, while cm[1,] has rownames
##   lapply(2:ncol(cm), function(j)
##          ta[[0, j - 1]] <- as.tclObj(cm[1, j], drop = TRUE))
##   for(j in 1:ncol(cm)) 
##     lapply(2:nrow(cm), function(i) 
##       ta[[i - 1, j - 1]] <- as.tclObj(cm[i, j], drop = TRUE))
## }

## ## tclArray -> DataFrame
## tclArrayToDataFrame <- function(ta, tktable, classes) {
##   d <- tkindex(tktable, "end") # get size from tktable
##   d <- as.numeric(unlist(strsplit(as.character(d), ",")))
##   l <- list()
##   for (j in 1:d[2]) {
##     vals <- sapply(1:d[1], function(i) {
##       val <- ta[[i,j]]
##       ifelse(is.null(val), NA, tclvalue(val))
##     })
##     l[[j]] <- try(switch(classes[j],
##                      factor=factor(vals),
##                      as(vals, classes[j])),
##                   silent=TRUE)
##     if(inherits(l[[j]], "try-error")) l[[j]] <- vals ## character
##   }
##   ind <- which(classes == "character")
##   if(length(ind)) {
##     ## convert NA to ""
##     for(i in ind) {
##       tmp <- l[[i]]
##       tmp[is.na(tmp)] <- ""
##       l[[i]] <- tmp
##     }
##   }
  
##   df <- as.data.frame(l)
##   ## fix character -- turned to factor above through as.data.frame
##   if(length(ind)) {
##     df[,ind] <- as.character(df[,ind])
##   }
##   ## dimnames
##   getTclValueWithDefault <- function(val, default) {
##     if(is.null(val))
##       default
##     else
##       tclvalue(val)
##   }
##   colnames(df) <- sapply(1:d[2], function(j) getTclValueWithDefault(ta[[0,j]], sprintf("X%s",j)))
##   rownames(df) <- make.row.names(sapply(1:d[1], function(i) getTclValueWithDefault(ta[[i,0]], as.character(i))))
##   return(df)
## }






## helper function here
## unlike make.names this doesn't put "X" prefix
make.row.names <- function(x) {
  dups = duplicated(x)
  if(any(dups))
    x[dups] <- make.unique(x)[dups]
  return(unlist(x))
}



###

## Code for interfacing with tablelist5.6 which is loaded in
## zzz.R


## Events are:  <<TablelistCellUpdated>> <<TablelistSelect>> 


## Configure tbl
tl_configure_columns <- function(tbl, nms) {
  .Tcl(sprintf("%s configure -columns {%s}",
               tbl$ID,
               paste(sprintf("0 {%s} left", nms), collapse="\n")
      ))
  sapply(seq_along(nms), function(j) tl_set_column_editable(tbl, j))
}

## Load Data
## helper to load a row
tl_insert_row <- function(tbl, row) {
  if(length(row) == 1 && grepl(" ", row))
    row <- paste("{", row, "}", sep="")
  tcl(tbl, "insert", "end", unlist(lapply(row, as.character)))
}

tl_clear_data <- function(tbl) {
  tcl(tbl, "delete", "0", "end")
}

tl_load_data <- function(tbl, items) {
  ## need to clear old first!
  tl_clear_data(tbl)
  sapply(seq_len(nrow(items)), function(i)
         tl_insert_row(tbl, items[i,,drop=TRUE]))
}

## return tcl cell index
tl_get_cellindex <- function(tbl, i, j) {
  tcl(tbl, "cellindex", sprintf("%s, %s", i-1, j-1))
}


## Get Data
## get cell infor -- raw = text
tl_get_cell_raw <- function(tbl, i, j) {
  raw <- tcl(tbl, "cellcget", tl_get_cellindex(tbl, i, j), "-text")
  tclvalue(raw)
}

## returns text value for column -- must coerce to ...
tl_get_column_raw1 <- function(tbl, j) {
  m <- tl_no_rows(tbl)
  sapply(seq_len(m), function(i) tl_get_cell_raw(tbl, i, j))
}

##helper
parse_tcl <- function(x) {

  ctr <- 0
  y <- strsplit(x, "")[[1]]
  tmp <- character(0)
  cur <- ""
  
  push_chr <- function(cur, i) {
    if(cur == "") i else paste(cur, i, sep="")
  }
  commit_cur <- function() {
    if(nchar(cur) > 0)
      tmp <<- c(tmp, cur)
    cur <<- ""
  }
  for(i in y) {
    if(i == "{") {
      if(ctr == 1) {
        commit_cur()
      }
      ctr <- ctr + 1
    } else if(i == "}") {
      if(ctr == 2) {
        commit_cur()
      }
      ctr <- ctr - 1
    } else if(i == " ") {
      if(ctr == 1) {
        commit_cur()
      } else {
        cur <- push_chr(cur, i)
      }
    } else {
      cur <- push_chr(cur, i)
    }
  }
  commit_cur()
  tmp
}

tl_get_column_raw <- function(tbl, j) {
  raw <- tcl(tbl, "getcolumns", j-1, j-1)
  parse_tcl(tclvalue(raw))
}


## return character matrix
tl_get_raw <- function(tbl) {
  do.call(cbind, lapply(seq_len(tl_no_cols(tbl)), function(j) tl_get_column_raw(tbl, j)))
}

## coerce
coerce_raw <- function(x, values) UseMethod("coerce_raw")
coerce_raw.default <- function(x, values) as.character(values)
coerce_raw.integer <- function(x, values) as.integer(values)
coerce_raw.numeric <- function(x, values) as.numeric(values)
coerce_raw.logical <- function(x, values) as.logical(values)
coerce_raw.factor <- function(x, values) factor(values)


## names
tl_set_column_name <- function(tbl, j, nm) {
  tcl(tbl, "columnconfigure", j-1, title=nm)
}

tl_set_column_names <- function(tbl, nms) {
  for(j in seq_along(nms)) tl_set_column_name(tbl, j, nms[j])
}


tl_get_column_name <- function(tbl, j) {
  tail(as.character(tcl(tbl, "columnconfigure", j-1, title=NULL)), n=1)
}

tl_get_column_names <- function(tbl) {
  sapply(seq_len(tl_no_cols(tbl)), function(j) tl_get_column_name(tbl, j))
}

## remove column
tl_remove_column <- function(tbl, j) {
  tcl(tbl, "deletecolumns", j-1, j-1)
}


## sort by column
tl_sort_bycolumn <- function(tbl, j, decreasing=FALSE) {
  dir <- if(decreasing) "decreasing" else "increasing"
  tcl(tbl, "sortbycolumn", j-1, sprintf("-%s", dir))
}


## size
tl_no_rows <- function(tbl) as.numeric(tcl(tbl, "childcount", "root"))
tl_no_cols <- function(tbl) as.numeric(tcl(tbl, "columncount"))


##
tl_set_focus_on_cell <- function(tbl, i, j) {
  tcl(tbl, "see", sprintf("%s, %s", i-1, j-1))
}


## show/hide column
tl_hide_row <- function(tbl, i, hide=TRUE) {
  hide <- if(hide) 1 else 0
  tcl(tbl, "rowconfigure", i-1, hide=hide)
}

tl_hide_column <- function(tbl, j, hide=TRUE) {
  hide <- if(hide) 1 else 0
  tcl(tbl, "columnconfigure", j-1, hide=hide)
}

## toggle editabbility of column
tl_set_column_editable <- function(tbl, j, editable=TRUE) {
  editable <- if(editable) "yes" else "no"
  tcl(tbl, "columnconfigure", j-1, editable=editable)
}


