# TODO

# 2. Make selection caption work
# 4. Add menu
# 5. Autosize window
# 7. Fix t.to error
# 8. Remove cursor when not edit mode
# 9. Doesn't select on factor matches properly
# 10. Row key press paint update doesn't work
# 11. Page-down etc doesn't call selection

# Written by Tom Taverner <t.taverner@gmail.com>
# for the U.S. Department of Energy (PNNL, Richland, WA, USA)
# Website: http://omics.pnl.gov/software
#
# Notice: This computer software was prepared by Battelle Memorial Institute,
# hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
# Department of Energy (DOE).  All rights in the computer software are reserved
# by DOE on behalf of the United States Government and the Contractor as
# provided in the Contract.
#
# NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR
# IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
#
# This notice including this sentence must appear on any copies of this computer
# software.

#' A package for an editing data frames for RGtk2. ImIproves on base edit.data.frame function found in utils
#' @name RGtk2DfEdit-package
#' @docType package

# Bugs:
# Prints Logical NA as TRUE
# The cell render function is a kludge
# Some platforms don't paint column selection, I can't reproduce this

require(RGtk2)

"\\001\\011\\000" # wide dots
"\\001\\011\\000" # dashes

############################################################
# Style RC
############################################################

# don't make the focus line width bigger than the v-separator
style <- 'style "treeview-style" { 
GtkTreeView::horizontal-separator = 0
GtkTreeView::vertical-separator = 0
GtkTreeView::focus-line-width = 0
GtkTreeView::focus-line-pattern = "\\000\\000\\000"
GtkTreeView::grid-line-pattern = "\\000\\000\\000"
GtkTreeView::grid-line-width = 0
bg_pixmap[NORMAL]		= "<none>"
bg_pixmap[INSENSITIVE]	= "<none>"
bg_pixmap[SELECTED]		= "<none>"
bg_pixmap[ACTIVE]		= "<none>"
bg_pixmap[PRELIGHT]		= "<none>" 
}
class "GtkTreeView" style "treeview-style" 
style "paned-style" { 
GtkPaned::handle-size = 5
}
class "GtkPaned" style "paned-style"'

style <- 'style "treeview-style" { 
GtkTreeView::focus-line-width = 0
}
class "GtkTreeView" style "treeview-style" 
style "paned-style" { 
GtkPaned::handle-size = 5
}
class "GtkPaned" style "paned-style"'

gtkRcParseString(style)

DATA_OBJECTS = c("data.frame", "matrix", "array")

# 
#old_warning <- warning

#12345.67 is "12,345.67" in the US, "12 345,67" in France and "12.345,67" 

STACK_LIMIT <- 1E8
# Don't bother drawing row selection rectangles if there are more than this
MAX_ROWS_TO_DRAW_SELECTION <- 1000
VERSION_STRING <- "version 0.6.1"

COLUMN_OFFSET <- 1

# Optionally, we include a blank row at the bottom
# 0 if it isn't, 1 if it is
EXTRA_ROW <- 0

ToInternalColIdx <- function(x) x+COLUMN_OFFSET    
ToExternalColIdx <- function(x) x-COLUMN_OFFSET

# replace these with "" in pretty_print


if(.Platform$OS.type == "windows") {
  PLATFORM_OS_TYPE <- "windows"
} else if (.Platform$OS.type == "unix"){
  if (Sys.info()["sysname"] == "Darwin")
    PLATFORM_OS_TYPE <- "mac"
  else{ 
    PLATFORM_OS_TYPE <- "unix"    }
}

# for pasting
if(PLATFORM_OS_TYPE == "windows"){
  NEWLINE_CHAR <- "\r\n"
  HEADER_BOX_MARGIN <- 2 # how to get this?
} else {
  NEWLINE_CHAR <- "\n"
  HEADER_BOX_MARGIN <- 4
}

DO_CAIRO <- TRUE
if (Sys.info()["sysname"] == "Darwin")
  DO_CAIRO <- FALSE
  
DEFAULT_COLNAMES <-  c(LETTERS, sort(as.vector(outer(LETTERS, LETTERS, FUN=paste, sep=""))))

MergeAdjacent <- function(v, warn=T){
      dv <- diff(v)
      if(warn && sum(dv < 1)) 
        warning("Vector is not increasing")
    
      if(length(v)){
        w1 <- which(dv != 1)
        cbind(start=v[c(1, w1+1)],
          end=v[c(w1, length(v))])
          #, 
          #length=1+length(w1))
      } else {
        cbind(start=integer(0), end=integer(0))#, length=0)
      }
    }  
  
  # http://www.mail-archive.com/r-help@r-project.org/msg38009.html  
make.call = function(func, my.arg){
  Call <- match.call()
  fn <- deparse(substitute(func))
  Call[[1]] <- as.name(fn)
  Call$func <- NULL
  Call$my.arg <- NULL
  my.formals <- formals(func)
  Call[[2]] <- as.name(deparse(substitute(my.arg)))
  return(Call)
}

findxInParseTree <- function(txt){
  .env00 <- new.env() 
  .env00$flag <- FALSE  
  dfs <- function(ll){
    if(!is.call(ll)) {
      if (length(ll) && as.character(ll) == "x") .env00$flag <- TRUE
      return()
    } else {
      for(ii in 1:length(ll)) dfs(ll[[ii]])
    }
  }
  pt <- as.list(parse(text=txt))
  dfs(pt[[1]])
  return(.env00$flag)
}

  # Remove leading and trailing white space from text and replace all blanks with "NA"
StripBlanks <- function(x){
  x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
  x[x=="NA"|nchar(x) == 0] <- NA 
  return(x)
}

# Keep all our coercions in here
# to.levels means you're turning a factor into its levels, rather than
#   as.integer()
# Coercing a factor to a factor will keep its level ordering.
myDataTypeCoercions <- function(typ, x, to.levels=FALSE){
  rv <- NA
  tryCatch(
    if(to.levels){
      stopifnot(is.factor(x))
      rv <- levels(x)[as.integer(x)]
      rv <- myDataTypeCoercions("character", rv)
    } else if(typ=="integer"){
      rv <- as.integer(x)
    } else if (typ== "logical"){
      rv <- as.logical(x)
    } else if (typ=="numeric") {
      rv <- as.numeric(x)
    } else if ("factor"%in%typ){
      sbx <- StripBlanks(x)    
      if(inherits(x, "factor")){
        theLevels = unique(StripBlanks(levels(x)))
      } else {
        theLevels = unique(sbx, na.last = TRUE)
      }
      rv <- factor(sbx, levels=theLevels)
    } else {
      x <- as.character(x) 
      x[is.na(x)] <- ""
      rv <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)  
      rv[nchar(rv) == 0] <- ""
    }, warning = function(w) return(rv))
  return(rv)
}

GetClasses <- function(x){
 nc <- ncol(x)
 rv <- rep(NA, nc)
 for(jj in 1:nc){
   cc <- x[,jj]
   cx <- class(cc)
     # class might be "ordered factor"   
   if("factor"%in%cx) {   
     cx <- "factor"
   } else if(!is.atomic(cc) || length(cx) != 1) {
     cx <- "character"
   }
   rv[jj] <- cx
 }
 return(rv)
}
  
# Coerce frame2 to theClasses1 column classes 
CoerceDataTypes <- function(frame2, theClasses1, to.levels=FALSE){ 
  stopifnot(is.data.frame(frame2))
  theClasses2 <- GetClasses(frame2)
  if(NA%in%theClasses1) stop("Trying to coerce to NA class") 
  stopifnot(length(theClasses1) == length(theClasses2))
  xx <- theClasses1 != theClasses2  
  xx[is.na(xx)] <- TRUE 
  for(jj in which(xx)){ 
    frame2[,jj] <- myDataTypeCoercions(theClasses1[jj], frame2[,jj], to.levels)
  } # for jj  
  return(frame2) 
}

# We normally coerce when we change cells, but sometimes we might want to paste
# and accept the default data frame coercion
ChangeCells <- function(df, nf, row.idx, col.idx, do.coercion=T){ 
  dmm <- dim(df)
  if(!is.data.frame(nf)){
    nf <- data.frame(nf, check.names=FALSE)
  }                                 
  if(missing(row.idx)) row.idx <- 1:dmm[1] 
  if(missing(col.idx)) col.idx <- 1:dmm[2]
  stopifnot(ncol(nf) == length(col.idx) && nrow(nf) == length(row.idx)) 

  oldf = df[row.idx, col.idx, drop=F]
 
  theClasses <- GetClasses(df)
  newClasses <- GetClasses(nf)

  if(do.coercion){
    nf <- CoerceDataTypes(nf, theClasses[col.idx])
  }

  idxf <- which(theClasses == "factor")
  idxf <- idxf[idxf%in%col.idx]
  
  if(length(idxf)){
    for(jj in idxf){
      nf.col = which(col.idx==jj)
      xx <- df[[jj]] 
      lvls <- levels(xx) 
      to.model <- as(nf[,nf.col], class(lvls)) 
      to.model[!nchar(to.model)] <- NA  # 3-16-10
        # we're changing the levels
      if((!all(to.model%in%lvls))||(!all(xx[row.idx]%in%xx[-row.idx]))){  
        x <- as.vector(xx) 
        x[row.idx] <- to.model 
        df[,jj] <- myDataTypeCoercions("factor", x)
      } else {
        df[row.idx,jj] <- to.model
      }
    }
    cc <- !col.idx%in%idxf 
    df[row.idx, col.idx[cc] ] <- nf[,cc,drop=F]
  } else { # no factors to change
    if(!do.coercion){
        # if we paste into a blank or NA column,
        # want to make sure the new frame takes attributes
		for(jj in seq(length=length(newClasses))){
		   df.jj = col.idx[jj]

		   if (newClasses[jj] != theClasses[df.jj] && all(df[,df.jj] == ""))
		       class(df[,df.jj]) <- newClasses[jj]
		}
    }
    df[row.idx, col.idx] <- nf
  }

  return(list(df = df,  
    undo = list( 
      func = "ChangeCells", 
      args = list(nf=oldf, row.idx=row.idx, col.idx=col.idx, do.coercion=T) 
    ) 
  )) 
}


SetFactorAttributes <- function(df, col.idx, info){ 
  idx <- col.idx
  theCol <- df[[idx]]
  stopifnot("factor"%in%class(theCol))
  lvls <- levels(theCol)
  old <- list(levels=lvls)
  if(length(lvls) > 1){
    old$contrasts <- contrasts(theCol)
    old$contrast.name <- attr(theCol, "contrast.name")
  }

  if(identical(info$ordered, TRUE)) 
    theCol <- as.ordered(theCol)    
  if(!is.null(info$levels))
    theCol <- factor(theCol, levels=info$levels)
  if(!is.null(info$contrasts))
    contrasts(theCol) <- info$contrasts
  if(!is.null(info$contrast.name)) 
    attr(theCol, "contrast.name") <- info$contrast.name
  
  df[[idx]] <- theCol

  return(list(df = df,  

    undo = list( 
      func = "SetFactorAttributes", 
      args = list(col.idx=idx, info=old) 
    ) 
  )) 
}

  # to.levels: when coercing a factor, use levels(x)[as.numeric(x)]
CoerceColumns <- function(df, theClasses, col.idx, to.levels=FALSE){
  idx <- col.idx
  if(length(theClasses) == 1 && length(idx) > 1) 
    theClasses <- rep(theClasses, length(idx))
  stopifnot(length(idx) == length(theClasses)) 
  stopifnot(max(idx) <= ncol(df))

  old.c <- GetClasses(df)[idx] # the previous class of the column
  df[,idx] <- CoerceDataTypes(df[,idx,drop=F], theClasses, to.levels)

  return(list(df = df,
    undo=list(  
      func = "CoerceColumns", 
      args = list(theClasses = old.c, col.idx=idx)
    ) 
  )) 
} 

ChangeColumnNames <- function(df, theNames, col.idx){
  idx <- col.idx
  stopifnot(length(idx) == length(theNames)) 
  stopifnot(max(idx) <= ncol(df))

  oldNames <- colnames(df)[idx] 
  colnames(df)[idx] <- theNames 
 
  return(list(df = df,
    undo=list(  
      func = "ChangeColumnNames", 
      args = list(theNames = oldNames, col.idx=idx)
    ) 
  )) 
} 
 
ChangeRowNames <- function(df, theNames, row.idx){
  idx <- row.idx
  stopifnot(length(idx) == length(theNames))   
  stopifnot(max(idx) <= nrow(df))
  rdf <- rownames(df)
  oldNames <- rdf[idx] 
  theNames <- make.unique(c(rdf[-idx], theNames))[length(rdf)+1:length(idx)-length(idx)] 
  rownames(df)[idx] <- theNames
  df[idx, 1] <- theNames
  return(list(df = df,
    undo=list(  
      func = "ChangeRowNames", 
      args = list(theNames = oldNames, row.idx=idx)
    ) 
  )) 
} 

# Makes indexes 
InsertIndex <- function(orig, insertions){ 
  all.idx <- seq(length=orig)
  li <- seq(length=length(insertions)) 
  del.idx <- insertions + li - 1 
  ins.idx <- orig + li 
  for(jj in li){ 
    after.idx = max(del.idx[jj]-1, 0)
    all.idx <- append(all.idx, ins.idx[jj], after=after.idx)
  }
  return(all.idx) 
} 
# deleting indexed rows and columns 
DeleteRows <- function(df, row.idx){
  idx <- row.idx   
  new_df <- df[-idx,,drop=F]
  new_df = maintain_null_rownames(df, new_df)
  list(df = new_df,
    undo = list( 
      func = "InsertRows", 
      args = list(
        nf = df[idx,,drop=F],
        row.idx=idx-1:length(idx)+1)))
}

InsertRows <- function(df, nf, row.idx){
  idx <- row.idx
  ddf <- dim(df)

  stopifnot(ddf[2] == dim(nf)[2] && dim(nf)[1] == length(idx))

  rv <- InsertNARows(df, row.idx)
  new_df <- rv$df
  nf <- CoerceDataTypes(nf, GetClasses(df))  
  new_df[row.idx,] <- nf
  new_df = maintain_null_rownames(df, new_df)
  rv$df <- new_df
  return(rv)
}

RowNamesAreNull = function(df)
  identical(rownames(df), as.character(1:dim(df)[1]))

maintain_null_rownames = function(df, new_df){
    # set this flag to TRUE if the row names are 1:dim(df)[1]
    # and remember that the last rowname is actually " ", so omit it
  if(EXTRA_ROW) { 
    row_names_null_flag <- identical(rev(rownames(df))[-1], as.character((dim(df)[1]-1):1))
  }  else {
    row_names_null_flag <- RowNamesAreNull(df)
  }
    # insert the new frame into the old frame at the insertion position
  if(row_names_null_flag && dim(df)[1]>0){
    #cat("row names are null\n")
    rownames(new_df) <- new_df[,1] <- seq(length=nrow(new_df))
  }
  return(new_df)
}
    
InsertNARows <- function(df, row.idx){
  idx <- row.idx
  ddf <- dim(df)
  xx <- InsertIndex(ddf[1], idx)  
  lidx <- length(idx)

  nf <- data.frame(rbind(rep(NA, ddf[2])))[rep(1, lidx),]
  colnames(nf) <- colnames(df)
  rownames(nf) <- make.unique(c(rownames(df), idx+1:lidx-1))[ddf[1]+1:lidx]
  nf[,ddf[2]] <- ""
  nf[,1] <- rownames(nf)
  nf <- CoerceDataTypes(nf, GetClasses(df))
  
  new_df = rbind(df, nf)[xx,,drop=F]
  new_df = maintain_null_rownames(df, new_df)

  list(df = new_df,  
    undo = list( 
      func = "DeleteRows", 
      args = list(
        row.idx=idx+1:length(idx)-1)))
}

# Check if we have null type column names
ColNamesAreNull = function(df)
identical(
    colnames(df)[seq(from=2, length.out=(ncol(df)-2))], 
    DEFAULT_COLNAMES[seq(from=1, length.out=(ncol(df)-2))]
  )

# if the old df (df) has default colnames, make sure the new_df does too
maintain_null_colnames = function(df, new_df){
  stopifnot(ncol(df)>=2)
  col_names_null_flag <- ColNamesAreNull(df)
  if(col_names_null_flag && ncol(new_df) > 2)
    colnames(new_df) <- c("rows", DEFAULT_COLNAMES[seq(from=1, length.out=(ncol(new_df)-2))], "")
  return(new_df)
}

DeleteColumns <- function(df, col.idx) {
  idx <- col.idx
  new_df <- df[,-idx,drop=F]
  new_df = maintain_null_colnames(df, new_df)

  list(df = new_df,  
    undo = list( 
      func = "InsertColumns", 
      args = list(
        nf = df[,idx,drop=F],
        col.idx=idx-1:length(idx)+1))) 
}

InsertColumns <- function(df, nf, col.idx){
  idx <- col.idx
  ddf <- dim(df)

  stopifnot(ddf[1] == dim(nf)[1] && dim(nf)[2] == length(idx))
  xx <- InsertIndex(ddf[2], idx) 

  new_df = cbind(df, nf)[,xx,drop=F]
  new_df = maintain_null_rownames(df, new_df)

  list(df = new_df,  
    undo = list( 
      func = "DeleteColumns", 
      args = list(
        col.idx=idx+1:length(idx)-1)))
}

InsertNAColumns <- function(df, col.idx, NA.opt=""){
  idx <- col.idx
  ddf <- dim(df)
  xx <- InsertIndex(ddf[2], idx) 
  lidx <- length(idx)
  nf <- data.frame(X=cbind(rep(NA.opt, ddf[1])),stringsAsFactors=FALSE, check.names=FALSE)[,rep(1, length(idx)),drop=F]

  colnames(nf) <- make.unique(c(colnames(df), DEFAULT_COLNAMES[idx+1:lidx-2]))[ddf[2]+1:lidx]

  new_df <- cbind(df, nf)[,xx,drop=F]

  new_df = maintain_null_colnames(df, new_df)

  list(df = new_df,  
    undo = list( 
      func = "DeleteColumns", 
      args = list(
        col.idx=idx+1:length(idx)-1)))
}

# insert new frame dat
GetTaskPasteIn <- function(theFrame, dat, insert.row, insert.col, 
  do.rownames = F, do.colnames = F, do.coercion=T){ 
  stopifnot(class(dat) == "data.frame")

  theRownames <- NULL 
  theColnames <- NULL 
  if(do.rownames)
    theRownames <- rownames(dat)
  if (do.colnames)
    theColnames <- colnames(dat) 
   
  dd <- dim(dat)
  dm <- dim(theFrame)
  if (is.null(insert.row)) insert.row <- 1
  if (is.null(insert.col)) insert.col <- 1
  ins <- c(insert.row, insert.col)
  
  ins.end <- c(ins[1]+dd[1]-1, ins[2]+dd[2]-1) # end of insertion range  
  tasks <- list()

    # need to expand the frame by adding NA rows/cols
  if(dm[1] < ins[1]+dd[1]-1){ # Rows
    #Indexing: 1 means insert so it forms new row 1, so there are N+1 insert
    # positions with an N row data frame.
    #idx <- (dm[1]+1):(ins[1]+dd[1])
    #idx <- rep(dm[1], length(idx))
    idx <- rep(dm[1]+1, ins[1]+dd[1]-dm[1]-1+EXTRA_ROW)
    tasks[[length(tasks)+1]] <- list(func="InsertNARows", args=list(row.idx=idx)) 
  }
  if(dm[2] < ins[2]+dd[2]){ # Columns
    #idx <- (dm[2]+1):(ins[2]+dd[2])
    idx <- rep(dm[2], dd[2]+ins[2]-dm[2])
    #idx <- rep(dm[2], length(idx))
    tasks[[length(tasks)+1]] <- list(func="InsertNAColumns", args=list(col.idx=idx))
  }

  row.idx <- ins[1]:ins.end[1]
  col.idx <- ins[2]:ins.end[2]
  
  if(do.rownames && !is.null(theRownames))
    tasks[[length(tasks)+1]] <- list(func="ChangeRowNames", 
        arg = list(theNames = theRownames, row.idx=row.idx))

  if(do.colnames && !is.null(theColnames))
    tasks[[length(tasks)+1]] <- list(func="ChangeColumnNames", 
        arg = list(theNames = theColnames, col.idx=col.idx))  

  if(nrow(dat) && ncol(dat))
    tasks[[length(tasks)+1]] <- list(func="ChangeCells", 
      args=list(nf=dat, row.idx=row.idx, col.idx=col.idx, do.coercion=do.coercion))

  return(tasks)
} 

# Do a user call
DoUserCall <- function(call_name, fargs, .local){
  handler <- .local$changed.handler[[call_name]]  
  if(!is.null(handler$func) && is.function(handler$func) ){
    fargs <- append(list(obj = .local$group.main), fargs)
    if(!is.null(handler$data))
      fargs <- append(fargs, list(data=handler$data))

    tryCatch({
      do.call(handler$func, fargs)
    }, error = function(e) warning(e))
  }
}

# Do a task list
# returns list with df=new frame, undo = list
DoTask <- function(.local, df, task, handler=NULL){
#print("DoTask called")
  w <- TransientWindow("Updating...", .local)
  on.exit(w$destroy())
  
  undo <- list()
  for(taskItem in task){
    arg <- taskItem$arg
    func.name <- taskItem$func
    arg$df <- df

    # This can fail when it's too big
#print(arg)
    rv <- do.call(taskItem$func, arg) 
    df <- rv$df
    undo[[length(undo)+1]] <- rv$undo
        
    if("col.idx"%in%names(taskItem$arg)) taskItem$arg$col.idx <- ToExternalColIdx(taskItem$arg$col.idx)
    DoUserCall(func.name, taskItem$arg, .local)
    if("col.idx"%in%names(taskItem$arg)) taskItem$arg$col.idx <- ToInternalColIdx(taskItem$arg$col.idx)
  }
  return(list(df=df, undo=undo))
}
 

  bgColor <- list(as.GdkColor(c(237,236,235)*256),as.GdkColor(c(235,234,219)*256))[[(PLATFORM_OS_TYPE == "windows")+1]]
  selectedColor <- as.GdkColor(c(198, 213, 253)*256) # Linux
  whiteColor <- as.GdkColor(c(255, 255, 255)*256)
  #selectedColor <- as.GdkColor(c(255,255,255)*256) # Linux
  selectedColumnColor <- as.GdkColor(c(198, 213, 253)*256) # Linux  
  selectedTextColor <- as.GdkColor("black")
  
  myValidInputKeys <- c(GDK_space:GDK_asciitilde, GDK_Delete, GDK_BackSpace)

  myMetaKeys <- c(GDK_Shift_L, GDK_Shift_R, GDK_Control_L, GDK_Control_R)
  myShiftKeys <- c(GDK_Shift_L, GDK_Shift_R)
  myValidNavigationKeys <- c(GDK_Down, GDK_Left, GDK_Up, GDK_Right,
  GDK_Page_Up, GDK_Page_Down, GDK_Return, GDK_ISO_Left_Tab, GDK_Tab, GDK_Home, GDK_End, GDK_KP_Enter, GDK_KP_Left:GDK_KP_Down)
  myValidKeys <- c(myValidInputKeys, myMetaKeys, myValidNavigationKeys, GDK_Insert)

   # make our own key-letter list
  myGDKKeys <- sapply(apropos("GDK_"), get)
  myGDKKeys <- sort(unlist(myGDKKeys[sapply(myGDKKeys, class)=='numeric']))


  exclude.factor <- NA

#' Convenience function to call data frame editor in its own window
#'
#' @param items A data frame (?) to display graphically
#' @param dataset.name Name for data set
#' @param size (height, width)
#' @param col.width (width)
#' @return An object of class GtkDfEdit for which a few RGtk2-style methods are defined
#' @export
# Edit object in place argument
dfedit <- function(items, dataset.name = deparse(substitute(items)), 
  size=c(600, 300), col.width = 64, editable=TRUE, 
  autosize = is.null(dim(items))||ncol(items)<25,
update=TRUE, modal=TRUE){
  
  stopifnot(modal%in%c(TRUE, FALSE))
  if(missing(items)||is.null(items)) {
    items <- data.frame(NULL)
    #items <- data.frame(matrix("", 1E2, 26), stringsAsFactors=FALSE, check.names=FALSE)
    dataset.name <- "Untitled"
  }

  obj <- gtkDfEdit(items, dataset.name, size.request=size, col.width = col.width, editable=editable, autosize=autosize, update=update)
  #dialog <- gtkDialog(dataset.name, NULL, "modal", "Close", 1,show = FALSE)
  if(identical(editable, TRUE) && identical(modal, TRUE)) dialog <- gtkDialog("Data Frame Editor", NULL, c("destroy-with-parent"), "gtk-ok", 1, "gtk-cancel", 0, show = FALSE)    
  else dialog <- gtkDialog("Data Frame Viewer", NULL, c("destroy-with-parent"), "gtk-close", 0, show = FALSE)
   
  #dialog$setModal(FALSE)
  #dialog$setPosition(GtkWindowPosition["center-on-parent"])

  infoBox <- gtkHBoxNew()
  selLabel = gtkLabelNew("")
  selLabel$setAlignment(0, 0.5)
  selLabel$setSizeRequest(150, -1)
  infoBox$packStart(selLabel, padding=2, expand=FALSE)

  .local <- obj$getData(".local")

  dialog[["vbox"]]$packStart(infoBox, padding=2, expand=FALSE)

  dialog[["vbox"]]$add(obj)
  dialog[["vbox"]]$setFocusChild(obj)

  obj$setActionHandler("OnLoad", function(obj, sourc, typ) {
   .local <- obj$getData(".local")
   dialog$setTitle(paste(.local$dataset.name, "- dfedit", VERSION_STRING))
  })

     obj$setActionHandler("Selection", function(obj, selections, ...) {
      selLabel$setText("")
      dstr <- ""
      if(length(selections)==1){
        sel.col <- c(selections[[1]]$start['col.idx'], selections[[1]]$end['col.idx'])
        sel.col.let <- DEFAULT_COLNAMES[c(selections[[1]]$start['col.idx'], selections[[1]]$end['col.idx'])]
        sel.row <- c(selections[[1]]$start['row.idx'], selections[[1]]$end['row.idx'])
        if(length(sel.row)==2 && length(sel.col)==2){
          sel.col <- sort(sel.col)
          sel.row <- sort(sel.row)
          if(sel.row[1]==sel.row[2] && sel.col[1]==sel.col[2]){
            range.txt <- paste(sel.col.let[1], sel.row[2], sep="")
            dstr <- as.character(.local$theFrame[sel.row[1], sel.col[1]+COLUMN_OFFSET])
          } else {
            range.txt <- 
               paste(sel.row[2]-sel.row[1]+1, 
               " R x ", sel.col[2] - sel.col[1]+1, " C", sep="")
  #             paste(sel.col.let[1], sel.row[1], 
  #             ":", sel.col.let[2], sel.row[2], sep="")
          }
        }
      } else if(length(selections)>1){
        range.txt <- "Multiple ranges"
      } else {
        range.txt <- ""
      }
      selLabel$setText(range.txt)
    })

  view <- .local$view
  #dialog$showAll()
  #win <- gtkWindowNew()  
  #win$setTitle(dataset.name)
  if(modal){
    dialog$showAll()
    if (dialog$run() == 1){
      rv <- invisible(obj$getDataFrame())
      dialog$destroy()
      return(invisible(rv))
    } else {
	   dialog$destroy()
	   return(invisible(items))
    }
  } else {
    gSignalConnect(dialog, "response", gtkWidgetDestroy)
    dialog$showAll()
    return(invisible(obj))
  }
}




dfview <- function(items, dataset.name = deparse(substitute(items)), 
  size=c(600, 300), col.width = 64, editable=FALSE, 
  autosize = is.null(dim(items))||ncol(items)<25,
  update=FALSE){
  rv <- dfedit(items, dataset.name, 
  size=size, col.width=col.width, editable=editable, autosize=autosize, update=update)
  invisible(rv)
}
   
# Our handler:
#      isText <- theClass == "factor" || theClass == "character"
#	    renderer <- .local$allColumns[[kk-1]]$renderer
#      renderer.set <- renderer['ellipsize-set']
#      if( isText && !renderer.set) {
#	        renderer['ellipsize'] <- PangoEllipsizeMode['end']
#	        renderer['ellipsize-set'] <- TRUE
#      } else if (!isText && renderer.set) {
#	        renderer['ellipsize'] <- PangoEllipsizeMode['none']
#	        renderer['ellipsize-set'] <- FALSE
#      }

DoUndo <- function(.local){  
  if(!length(.local$undoStack)) return(TRUE)
  undo <- .local$undoStack[[length(.local$undoStack)]]
  rv <- DoTask(.local, .local$theFrame, rev(undo), .local$changed.handler)
  UpdateDfEditor(.local, rv$df)
   if(length(.local$undoStack)){
     #.local$redoStack[[length(.local$redoStack)+1]] <- rev(undo)
     .local$undoStack[[length(.local$undoStack)]] <- NULL
   }
}

DoRedo <- function(.local){
  if(!length(.local$redoStack)) return(TRUE)
  redo <- .local$redoStack[[length(.local$redoStack)]]
  rv <- DoTask(.local, .local$theFrame, redo, .local$changed.handler)
  UpdateDfEditor(.local, rv$df)
   if(length(.local$undoStack))
     .local$undoStack[[length(.local$undoStack)]] <- NULL
}

TransientWindow <- function(txt, .local=NULL){
    w <- gtkWindowNew("", show=F) #GTK_WINDOW_POPUP
    gtkWidgetSetCanFocus(w, TRUE)
    gtkWindowSetSkipTaskbarHint(w, TRUE)
    gtkWidgetSetSizeRequest(w, 100, 30)
    gtkWindowSetDecorated(w, FALSE)
    lt <- gtkLabelNew(txt)
    gtkAdd(w, lt)
    gtkWindowSetPosition(w, GtkWindowPosition["center-on-parent"])              
    gtkWindowSetTransientFor(w, .local$toplevel)    
    gtkWindowSetModal(w, TRUE)
    gtkWidgetShowAll(w)
    #gtkWidgetGrabFocus(w)
    return(w)
}


# restrict means you can't select the last index
GetSelectedRows <- function(tv, .local, ignore.last.row=FALSE){
   # block transient interactions with a popup
#  w <- TransientWindow("Getting row selection...", .local)  
#  on.exit(w$destroy())  
  sr <- integer(0)
   #tryCatch({    
    selection <- gtkTreeViewGetSelection(tv)
    rv <- gtkTreeSelectionGetSelectedRows(selection)$retval
    lrv <- length(rv)
    if(lrv){
      #sr <- vector("integer", lrv)
      #for(ii in seq(length=lrv))
      #  sr[ii] <- gtkTreePathGetIndices(rv[[ii]])
      sr <- sapply(rv, gtkTreePathGetIndices)
      if(!ignore.last.row){
        if(is.numeric(sr)) {
          sr <- sr + 1
          } else {
          sr <- integer(0)        
        }
      } else { # remove last row
        if(is.numeric(sr)) {
          sr <- sr + 1
          sr <- sr[!sr%in%dim(.local$theFrame)[1]]
        } else {
          sr <- integer(0)        
        }
      }
    }
  #},
  #error=function(e) {
  #  warning(e)
  #  integer(0)
  #})
  return(sr)
}
    
# Make a data frame to put into the model.
# We have an extra row for names and a blank at the end.
# Also, data frames with all-NAs turn into logical, which displays badly.
MakeInternalDataFrame <- function(dataset, add.rows=EXTRA_ROW, add.columns=T, NA_replace = NA_character_){
  
  if(is.null(dim(dataset))) dataset <- cbind(dataset)
  if(!length(rownames(dataset)) && dim(dataset)[1])
      rownames(dataset) <- seq(length=dim(dataset)[1])
  if(!length(colnames(dataset)) && dim(dataset)[2])
      colnames(dataset) <- DEFAULT_COLNAMES[seq(length=dim(dataset)[2])]
      
   # Fix bug with NA row name, 09-15-2010
  row_names_dataset <- row.names(dataset)
  row.problem.idx <- is.na(row_names_dataset) | duplicated(row_names_dataset)
  if(any(row.problem.idx)){
    dataset <- dataset[!row.problem.idx,,drop=F]
    message(paste("Removed missing or duplicated row names at positions:", paste(which(row.problem.idx), collapse=", ")))
  }  

  if(add.columns){
   theFrame <- data.frame(rows = row.names(dataset), dataset, 
     " " = vector("character", nrow(dataset)), 
     check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    theFrame <- data.frame(dataset, stringsAsFactors=FALSE)
  }
  
    # turn all-NA columns to whatever default we want
#  if(dim(theFrame)[2] > 0) {
#    which.allNA <- which(colSums(is.na(theFrame)) == dim(theFrame)[1])
#    for(jj in which.allNA) theFrame[,jj] <- NA_replace
#  }
#
  theClasses <- GetClasses(theFrame)  
  for(jj in which(theClasses == "character")) 
    theFrame[,jj] <- myDataTypeCoercions("character", theFrame[,jj])

  if(add.rows){ 
    if (add.columns){
		  blankRow <- data.frame(rows=" ", rbind(rep(NA, dim(dataset)[2])), " " = "", row.names = " ", stringsAsFactors=F)
   } else {
		blankRow <- data.frame(rbind(rep(NA, dim(dataset)[2])), stringsAsFactors=F)
   }
		blankRow <- CoerceDataTypes(blankRow, theClasses)
		for(jj in which(theClasses == "factor")){
  		theFrame[,jj] <- myDataTypeCoercions("factor", theFrame[,jj])
		  blankRow <- SetFactorAttributes(df=blankRow, col.idx=jj, 
        info=list(levels=levels(theFrame[,jj]), ordered=is.ordered(theFrame[,jj])))$df  
      }
		names(blankRow) <- names(theFrame)
		theFrame <- rbind(theFrame, blankRow)
  }
  return(theFrame)
}

# Make a data frame to extract from the model.
MakeExternalDataFrame <- function(theFrame, .local=NULL){

  if(EXTRA_ROW) new.frame <- theFrame[-nrow(theFrame),-c(1, ncol(theFrame)), drop=F]
  else 
    new.frame <- theFrame[,-c(1, ncol(theFrame)), drop=F]

  if(!is.null(.local)){
	  dataset.class <- .local$dataset.class
	  tryCatch({
		if(length(dataset.class)==0 || "data.frame"%in%dataset.class){ 
		  class(new.frame) <- dataset.class
		} else if(length(dataset.class)==1) {
		  new.frame <- as(new.frame, dataset.class)
		}}, error = function(e) warning(paste("Couldn't create data of class", paste(dataset.class, collapse = ", "))))

		# Copy any other attributes over
	  dataset.attributes <- .local$dataset.attributes
	  copy.attr <- dataset.attributes[!names(dataset.attributes)%in%c('class', 'dim', 'dimnames', 'names', 'row.names')]    
	  if(length(copy.attr)) for(ii in 1:length(copy.attr)) attr(new.frame, names(copy.attr)[ii]) <- copy.attr[[ii]]
  }

  return(new.frame)
}

# return last usable row path
MakeLastPath <- function(df){
  if(!EXTRA_ROW){
    if(dim(df)[1] == 0) return(gtkTreePathNewFromString("0"))
  } else { # adding an extra row
    if(dim(df)[1] == 1) return(gtkTreePathNewFromString("0"))
  } 
  return(gtkTreePathNewFromString(as.character(dim(df)[1]-1-EXTRA_ROW)))
}

quick_entry <- function(msg, handler, data=NULL, win=NULL) {
  dialog <- gtkDialog(msg, NULL, c("destroy-with-parent"), "gtk-ok", 1, "gtk-cancel", 0, show = F)
  #dialog$setDecorated(FALSE)
  entry = gtkEntryNew()        
  entry$setSizeRequest(300, -1)  
  gSignalConnect(entry, "key-press-event", function(obj, evt){
    if(evt[["keyval"]] == GDK_Return) {
      handler(entry$getText(), data)
      dialog$destroy()     
    }
    if(evt[["keyval"]] == GDK_Escape) {
      dialog$destroy()     
    }
    return(FALSE)
    })
   
  if(!is.null(win)) {
     checkPtrType(win, "GtkWindow")
     dialog$setPosition(GtkWindowPosition["center-on-parent"])
     dialog$setTransientFor(win)
  }
  gSignalConnect(dialog, "response", function(dlg, arg1, user.data) {
     if(arg1==1)
       handler(entry$getText(), user.data)
     dialog$destroy()
  }, data=data)                                    
  dialog[["vbox"]]$packStart(entry, TRUE, TRUE, 10)
  dialog$showAll()
}

quick_query <- function(message, handler, data, win=NULL) {
  dialog <- gtkDialog("Query", NULL, c("modal", "destroy-with-parent"), "gtk-ok", 1, "gtk-cancel", 0, show = F)
  label <- gtkLabel(message)     
  if(!is.null(win)) {
     checkPtrType(win, "GtkWindow")
     dialog$setPosition(GtkWindowPosition["center-on-parent"])
     dialog$setTransientFor(win)
  }
  gSignalConnect(dialog, "response", function(dlg, arg1, user.data) {
    if(arg1) handler(data)
  dialog$destroy()
  })
  dialog[["vbox"]]$packStart(label, TRUE, TRUE, 10)
  dialog$showAll()
}

quick_message <- function(message, caption="Warning", win=NULL) {
  dialog <- gtkDialog("Message", NULL, c("modal","destroy-with-parent"), "gtk-ok", 1,
                      show = FALSE)
  label <- gtkLabel(message)     
  if(!is.null(win)) {
     checkPtrType(win, "GtkWindow")
     dialog$setPosition(GtkWindowPosition["center-on-parent"])
     dialog$setTransientFor(win)
  }
  gSignalConnect(dialog, "response", gtkWidgetDestroy)  
  dialog[["vbox"]]$packStart(label, TRUE, TRUE, 10)
  dialog$showAll()
}
    
#warning <- quick_message

IsCtrl <- function(stat) (as.flag(stat) & GdkModifierType['control-mask'])

IsShift <- function(stat) (as.flag(stat) & GdkModifierType['shift-mask'])

IsCtrlShift <- function(stat) (as.flag(stat) & GdkModifierType['control-mask'] & GdkModifierType['shift-mask'])

CtrlLetter <- function(keyval, stat, let)
  keyval == let && (as.flag(stat) & GdkModifierType['control-mask'])
ShiftLetter <- function(keyval, stat, let)
  keyval == let && (as.flag(stat) & GdkModifierType['shift-mask'])  
CtrlShiftLetter <- function(keyval, stat, let)
  keyval == let && (as.flag(stat) & GdkModifierType['control-mask'] & GdkModifierType['shift-mask'])

   # this is all horribly platform dependent
CopyToClipboard <- function(dat, do.rownames=F, do.colnames=F){
   
  write.function <- function(dat, p)
    write.table(dat, p, row.names=F, col.names=F, quote=F, sep="\t")

  if(do.rownames) {
    t.dat <- t(cbind(rownames(dat), dat))
  } else {
    t.dat <- t(dat)      
  }

  dat2 <- paste( apply(t.dat, 2, function(ll) {
      ll[is.na(ll)] <- ""
      paste(ll, collapse="\t")
    }), collapse=NEWLINE_CHAR)
        
  if(do.colnames) {
    dat.cn <- ""
    if(do.rownames) dat.cn <- "\t" 
    dat.cn <- paste(dat.cn, paste(colnames(dat), collapse="\t"), sep="")
    dat2 <- paste(dat.cn, dat2, sep=NEWLINE_CHAR)
  }
  
  if(.Platform$OS.type == "windows") {
    #write.function(dat, "clipboard")
    get("writeClipboard", envir=.GlobalEnv)(dat2)
  } else if (.Platform$OS.type == "unix"){
    if (Sys.info()["sysname"] == "Darwin")
      a <- pipe("pbcopy", "w")
    else {
      if(!length(system('which xclip', intern=T))){
        quick_message("xclip must be installed to copy")
        return(FALSE)
      }
      a <- pipe("xclip -selection c", open="w")
    }
  }

  if(.Platform$OS.type != "windows" ) 
    tryCatch(
      write.function(dat2, a),
      error = stop,
      finally = close(a)
    )
}  

# Where to read from
GetPasteClipboardPipe <- function(){
  if(PLATFORM_OS_TYPE == "windows") {
    p <- "clipboard"
  } else if (PLATFORM_OS_TYPE == "unix"){
    if(!length(system('which xsel', intern=T))){
      quick_message("xsel must be installed to paste")
      return(FALSE)
    }
    p <- pipe("xsel -o -b", open="r")
  } else if (PLATFORM_OS_TYPE == "mac"){
    p <- pipe("pbpaste")
  } else{ 
    stop("Unrecognized platform type")
  }
}

CloseClipboardPipe <- function(p)
  if (PLATFORM_OS_TYPE == "unix" || PLATFORM_OS_TYPE=="mac" )
    close(p)

# Dialog for reading in tables
#

# Options dialog menu
# setwd("~/rgtk2extras/pkg/RGtk2Extras/R/")

# This wraps read.table and tries to handle a few common errors
# First common problem. Bug in read.table: if nrows <= 5 it gives an error
# So if it throws an error, just try scanning 1 row at a time.
# Error handling: first 5 rows can lack EOF and give an error

# Second common problem: numbers with commas or full stop delimiters.
# Use 'num.string.ignore' argument to filter out "%" "," etc
 
# quote="", na,strings="", comment.char=""
# Get pasted data from clipboard. 
# Arguments are for read.table
ReadTableFromClipboard = function(...){
  CHAR_MAX = 10000 # max chars in first 5 lines if they fail
  rv <- tryCatch({
  b = GetPasteClipboardPipe()
  rv = read.table(b, ...)
  },
  error = function(e) { 
        # if the read table fails,
        # we might need to append EOF
    if(identical(as.numeric(gregexpr("incomplete final line", e$message)), 1)){
      tryCatch({
        # Do this by calling pushBack on another paste pipe
        b2 = GetPasteClipboardPipe()
#print(readLines(b2))
        b2_lines = readChar(b2, CHAR_MAX)
#        tfn = tempfile(fileext=".txt")

#        cat(c(b2_lines, "\n"), file=tfn)
#print(readLines(tfn))
#        rt2 = read.table(tfn, ...)
#        unlink(tfn)
#print(b2_lines)
        if(nchar(b2_lines)==CHAR_MAX)
          warning(paste("Maximum of", 
          CHAR_MAX, "characters in first 5 lines exceeded."))
        pushBack(c(b2_lines, ""), b2)
        rv = read.table(b2, ...)
      }, error = stop, 
      finally = CloseClipboardPipe(b2))
    } else {
      # something else went wrong
      stop(e$message)
    }
  }, 
  finally= CloseClipboardPipe(b))
  return(rv)
}

CheckStringsForNumeric = function(rv, num.string.ignore = "%"){
  # Go through columns. If the first 5 lines give intelligible numbers on removing
  # the ignored characters, read as numbers.

  suppressWarnings({
  if(!is.null(num.string.ignore) && nrow(rv)) {
    num.string.ignore = paste("[", num.string.ignore, "]", sep="")
    test.idx = 1:min(5, nrow(rv))
    for(jj in seq(length=ncol(rv))){
      first.entries = rv[test.idx,jj]
      where.na = is.na(first.entries)
        # if removing the separator/% char doesn't introduce new NA's and is OK
      if(!all(where.na) && is.character(first.entries[!where.na])){
         read.nums = as.numeric(gsub(num.string.ignore, "", first.entries))
         if(any(is.numeric(read.nums)) && sum(is.na(read.nums)) <= sum(where.na)){
           rv[,jj] <- as.numeric(gsub(num.string.ignore, "", rv[,jj]))
         }
      }
    }
  }
  })
  #gsub(",", "", "1,234,567", fixed = TRUE)
  return(rv)
}

# Returns true only if the argument is 1x1, and 
# its whitespace-trimmed content doesn't contain the delimiter
LooksLikeASingleItem = function(dat, 
  delimit=paste(" \t;", ifelse(identical(options("OutDec")[[1]], "."), ",", ""), sep="")
  ){
dat[1,1] <- StripBlanks(dat[1,1])
!is.null(dat) && identical(as.numeric(dim(dat)), c(1,1)) && 
!length(grep(paste("[", delimit, "]", sep=""), dat[1,1]))
}

ReadFromClipboardWrapper <- function(.local, ...){
      # Check to see if the user's pasting in something 
      # that looks like a regular old number, otherwise,
      # go to the ReadFromClipboard dialog
  dat <- NULL
  tryCatch({
	  dat <- ReadTableFromClipboard(nrow=2, stringsAsFactors=FALSE)
	  if(LooksLikeASingleItem(dat)){
	   dat <- CheckStringsForNumeric(dat, num.string.ignore="[$%]")
	   rfc <- list(xx=dat, 
		 other.opts = list(paste.row.labels=F, paste.col.labels=F))
	   } else {
		 dat <- NULL
	   }
   }, error = function(e) {})

  if(is.null(dat)){
    rfc <- ReadFromClipboard(...) # character vector
	dat <- rfc$xx
  }
  if(!is.null(rfc)){
      #cmd.list <- rfc$cmd.list
      other.opts = rfc$other.opts
	  row.idx <- GetSelectedRows(.local$view, .local, ignore.last.row=EXTRA_ROW)
	  sc <- GetSelectedColumns(.local, restrict=FALSE)
	    # if there are no rows, we want to paste in
	  if(!length(row.idx)){
	    row.idx <- 1
	  }
	    # if there are no columns, we want to paste in
	  if(length(sc) == 0 && length(.local$allColumns) == 1) 
	    sc <- 1 # 
	  if(length(sc) && length(row.idx)){
	    col.idx <- min(sc)
	    row.idx <- row.idx[1]
	    task <- GetTaskPasteIn(.local$theFrame, dat, row.idx, col.idx+COLUMN_OFFSET, do.colnames=other.opts$paste.col.labels, do.rownames = other.opts$paste.row.labels, do.coercion=F)
	    DoTaskWrapper(.local, task)
	  }
  }
}

# ... args are arguments to read.table, which we might override or ignore (currently ignoring)
# returns the data frame and options
ReadFromClipboard <- function(handler=NULL, user.data=NULL, parent.window=NULL, pipe=NULL, fromFile=FALSE, fileName=NULL, sep=NULL, ...){
## read table dialog

  if(fromFile)
   stopifnot(!is.null(fileName) && !is.na(fileName) && nchar(fileName))# && file.exists(fileName)) # might be URL
       
  function.opts = list(...)
  dialog <- gtkDialog("Read Table", NULL, c("destroy-with-parent"), "gtk-ok", 1, "gtk-cancel", 0, show = F)
  dialog$setPosition(GtkWindowPosition["center-on-parent"])              
  dialog$setTransientFor(parent.window)          
  
  #dialog$setDecorated(FALSE)
  box0 <- gtkVBoxNew(FALSE, 5)  

  if (0){ 
  fr1 <- gtkFrameNew(label="Choose File")
  box1 <- gtkHBoxNew(FALSE, 5)
  fr1$add(box1)
  
  entry1 <- gtkEntryNew()
  entry1$setAlignment(0)
  button1 <- gtkButtonNewWithLabel("Select File...")

  box1$packStart(entry1, TRUE, TRUE, 0)
  box1$packStart(button1, FALSE, FALSE, 0)

  nlines1_adj <- gtkAdjustment(5, 1, 1000, 1, 1)
  nlines1 <- gtkSpinButton(nlines1_adj, 1.0, 0)
  nlines1$setValue(5)
  label1 <- gtkLabelNew("Lines To Preview:")

  box1$packEnd(nlines1, FALSE, FALSE, 0)
  box1$packEnd(label1, FALSE, FALSE, 0)

  gSignalConnect(button1, "clicked", data=entry1, function(obj, label){
    tryCatch({
      f <- my_choose_files(file.path(getwd(), "*.*"), multi=F)
      if(nchar(f) > 0){
        setwd(dirname(f))
        # modify settings
        ext <- strsplit(f, "[.]")[[1]]
        ext <- ext[length(ext)]
        if(ext == "csv"){
          sep1$setActive(2)
        } else if (ext == "txt") {
          sep1$setActive(3)
        } else if (ext == "prn") {
          sep1$setActive(1)
        } else {
          sep1$setActive(0)
        }
        
        entry1$setText(f)
        preview.txt <- paste(readLines(con = f, n = nlines1$getValue()), collapse="\n")
        tb$getBuffer()$setText(preview.txt)
      }
    }, error = function(e){print(e)})
  })
    
  box1 <- gtkHBoxNew(FALSE, 5)

  fr0 <- gtkFrameNew(label="File Preview")
  scroll <- gtkScrolledWindowNew()
  scroll$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
  
  tb <- gtkTextViewNew()
  tb$setEditable(FALSE)
  tb$setSizeRequest(400, 100)
  scroll$add(tb)
  fr0$add(scroll)

  tb$getBuffer()$setText("No File Loaded")
  box0$add(fr1)
  #box0$add(fr0)

  }

  fr2 <- gtkFrameNew(label="Table Read Options")
  
  box2 <- gtkVBoxNew(FALSE, 5)

  sep.types = list("Tab (\"\t\")" = "\t", "Whitespace" = "", "Space (\" \")" = " ", "Comma (\",\")" = ",", "Semicolon (\";\")" = ";")

  sep1 <- gtkComboBoxNewText()
  for (ll in names(sep.types)) sep1$appendText(ll)
  sep1$setTooltipText("The field separator character. Values on each line of the file are separated by this character. White space is one or more spaces, tabs, newlines or carriage returns.")
  # if (identical(.Options$dec[[1]], ".")) sep1$setActive(2)
  # else 
  if(is.null(sep)) {
    sep1$setActive(0)
  } else {
    stopifnot(sep%in%sep.types)
    sep1$setActive(which(sep==sep.types)-1)
  }
  
  hbox <- gtkHBoxNew(FALSE, 0)
  hbox$add(gtkLabelNew("Separator"))
  hbox$add(sep1)
#  box2$packStart(hbox, FALSE, FALSE, 0)

NUMBER_FORMATS = data.frame(row.names= c("thousands", "decimal"),
   "US/UK (12,345.67)" = c(",", "."),
   "France (12 345,67)" = c(" ", ","), 
   "Germany (12.345,67)" = c(".", ","), stringsAsFactors=F, check.names=FALSE
)

  dec1 <- gtkComboBoxNewText()
  dec1$setTooltipText("The characters used for decimal points and thousands separators")
  dec1$setSizeRequest(120, -1)
  #for (ll in c(".", ",")) dec1$appendText(ll)
  for (ll in colnames(NUMBER_FORMATS)) dec1$appendText(ll)
  # Default to France if we aren't in US/UK
  dec1$setActive(ifelse(identical(options("OutDec")[[1]], "."), "0", "1"))
#  hbox <- gtkHBoxNew(FALSE, 0)
  hbox$add(gtkLabelNew("Format"))
  hbox$add(dec1)
  box2$packStart(hbox, FALSE, FALSE, 0)

  opt.box <- gtkHBoxNew(FALSE, 0)
  opt.box1 <- gtkVBoxNew(FALSE, 0)
  opt.box2 <- gtkVBoxNew(FALSE, 0)
  opt.box$add(opt.box1)
  #opt.box$add(opt.box2)
  box2$add(opt.box)
  
  opt.table <- gtkTableNew(6, 2, homogeneous=FALSE)
  opt.table$setRowSpacings(5)
  opt.box$packEnd(opt.table, TRUE, TRUE, 5)

  row.names0 <- gtkCheckButtonNewWithLabel(label="Use Row Labels")

#  row.names0$setActive(TRUE)
  opt.box1$packStart(row.names0, FALSE, FALSE, 0)

  paste.row.labels <- gtkCheckButtonNewWithLabel(label="Row Labels")
  paste.col.labels <- gtkCheckButtonNewWithLabel(label="Column Labels")
  pasteoptsbox <- gtkVBoxNew()
  pasteoptsbox$add(paste.row.labels)
  pasteoptsbox$add(paste.col.labels)
 pasteoptsbox$setTooltipText("Whether to paste row and column labels along with the data")

  check.names1 <- gtkCheckButtonNewWithLabel(label="Check Names")
  check.names1$setTooltipText(" If TRUE then the names of the variables in the data frame are checked to ensure that they are syntactically valid variable names.  If necessary they are adjusted (by make.names) so that they are, and also to ensure that there are no duplicates.")
  pasteoptsbox$add(check.names1)

  fr.pasteopts <- gtkFrameNew(label="Paste Options")
  fr.pasteopts$add(pasteoptsbox)

  header1 <- gtkCheckButtonNewWithLabel(label="Use Header Row")
  header1$setTooltipText("Whether the file contains the names of the variables as its first line.")
  header1$setActive(FALSE)

  #check.opts = function(theArg, checkButton)
  #  if(theArg%in%names(function.opts) && function.opts[[theArg]]%in%c(TRUE, FALSE)) gtkToggleButtonSetActive(checkButton, (function.opts[[theArg]]))
  #check.opts("header", header1)

  opt.box1$packStart(header1, FALSE, FALSE, 0)
  
  # Stolen from Peter Dalgaard
  ASCII <- c("Any White Space", sapply(1:255, function(i) parse(text=paste("\"\\", structure(i,class="octmode"), "\"", sep=""))[[1]]))


  row.names1 <- gtkSpinButton(gtkAdjustment(1, 1, 1000, 1, 1), 1.0, 0)
  rowTooltip <- "Switching this on allows the column of the table giving row names to be specified. If it is off, and the first row contains one fewer field than the number of columns, the first column in the input is used for the row names."
  row.names0$setTooltipText(rowTooltip)
  row.names1$setTooltipText(rowTooltip)
  row.names1$setValue(1)

  row.names1$setSensitive(row.names0$getActive())
  gSignalConnect(row.names0, "toggled", function(widget){
    row.names1$setSensitive(widget$getActive())
  })
    
  #hbox <- gtkHBoxNew(FALSE, 0)
  label = gtkLabelNew("Use Column")
  label$setAlignment(1, 0.5)                
  opt.table$attach(label, 0, 1, 0, 1, xpad=10)
  opt.table$attach(row.names1, 1, 2, 0, 1)

  skip1 <- gtkSpinButton(gtkAdjustment(0, 0, 1000, 1, 1), 1.0, 0)
  skip1$setValue(0)
skip1$setTooltipText("the number of lines of the data file to skip before beginning to read data.")
  label = gtkLabelNew("Rows To Skip")
  label$setAlignment(1, 0.5)
  opt.table$attach(label, 0, 1, 1, 2, xpad=10)
  opt.table$attach(skip1, 1, 2, 1, 2)
  
  colClasses1 <- gtkComboBoxNewText()
  for (ll in c("Default", "logical", "integer", "numeric", "complex", "character", "raw", "factor", "Date", "POSIXct")) colClasses1$appendText(ll)
  colClasses1$setActive(0)
  #hbox <- gtkHBoxNew(FALSE, 0)
  label = gtkLabelNew("Column Classes")
  label$setAlignment(1, 0.5)  
  opt.table$attach(label, 0, 1, 2, 3, xpad=10)
  opt.table$attach(colClasses1, 1, 2, 2, 3)
  
  na.strings1 <- gtkEntryNew()
  na.strings1$setTooltipText("strings which are to be interpreted as NA values.  Blank fields are also considered to be missing values.")
  na.strings1$setText("NA")
  hbox <- gtkHBoxNew(FALSE, 0)
  label = gtkLabelNew("NA String")
  label$setAlignment(1, 0.5)  
  opt.table$attach(label, 0, 1, 3, 4, xpad=10)
  opt.table$attach(na.strings1, 1, 2, 3, 4)

  quote1 <- gtkComboBoxNewText()
  for (ll in c("", "'", "\"", "\"'")) quote1$appendText(ll)
quote1$setTooltipText("the set of quoting characters. To disable quoting altogether, use \"\" ")
  quote1$setActive(0)
  #hbox <- gtkHBoxNew(FALSE, 0)
  label = gtkLabelNew("Quotes")
  label$setAlignment(1, 0.5)  
  opt.table$attach(label, 0, 1, 4, 5, xpad=10)
  opt.table$attach(quote1, 1, 2, 4, 5)
  
  fill1 <- gtkCheckButtonNewWithLabel(label="Fill Unequal Length Rows")
  fill1$setTooltipText("If TRUE then in case the rows have unequal length, blank fields are implicitly added.")
  fill1$setActive(TRUE)
  opt.box1$packStart(fill1, FALSE, FALSE, 0)

  stringsAsFactors1 <- gtkCheckButtonNewWithLabel(label="Read Strings As Factors")
stringsAsFactors1$setTooltipText("Convert character variables to factors.")
  stringsAsFactors1$setActive(FALSE)
  opt.box1$packStart(stringsAsFactors1, FALSE, FALSE, 0)

  num.ignore.string1 <- gtkEntryNew()
num.ignore.string1$setTooltipText("Ignore these characters if they appear inside or around numeric variables. For example, reading '12.3%', '12 345.67' and '$512.00' as numbers should ignore '%', ' ' and '$' respectively.")

  num.ignore.string1$setText(paste("%$", NUMBER_FORMATS["thousands", dec1$getActiveText()], sep=""))
  hbox <- gtkHBoxNew(FALSE, 0)
  label = gtkLabelNew("Numbers Ignore")
  label$setAlignment(1, 0.5) 
  opt.table$attach(label, 0, 1, 5, 6, xpad=10)
  opt.table$attach(num.ignore.string1, 1, 2, 5, 6)

  opt.box1$add(fr.pasteopts)
  
  fr2$add(box2)
  
  #box0$packStart(fr1, FALSE, FALSE, 0)
  #box0$packStart(fr0, FALSE, FALSE, 0)
  
  fr3 <- gtkFrameNew(label="Table Preview")
  box3 <- gtkVBoxNew(FALSE, 5)
  #preview.button <- gtkButtonNewWithLabel("  Preview Table  ")
  #hbox <- gtkHBoxNew(FALSE, 0)
  #hbox$packEnd(preview.button, FALSE, FALSE, 5)
  #box3$packStart(hbox, FALSE, FALSE, 0)
  
  tb2 <- gtkTextViewNew()
  tb2$setSizeRequest(400, 200)
  tb2$getBuffer()$setText("No Table loaded")
  box3$packStart(tb2)

  read_table <- function(nrows=-1, doMsg=FALSE){
#    if(nchar(entry1$getText()) == 0) stop("No file selected")

    row.names.arg <- NULL
    if(row.names0$getActive()) row.names.arg <- as.integer(row.names1$getText())

    col.classes.arg <- NA

    cmd.list = list(header= header1$getActive(),
      sep = sep.types[[sep1$getActiveText()]],
      skip= as.integer(skip1$getText()),
      #row.names = row.names.arg,
      nrows = nrows,
      colClasses= col.classes.arg,
      stringsAsFactors = stringsAsFactors1$getActive(),
      fill = fill1$getActive(),
      dec = NUMBER_FORMATS["decimal", dec1$getActiveText()],
      quote = quote1$getActiveText(),
      strip.white=TRUE,
      comment.char = "",
      check.names=check.names1$getActive()
    )
    
    if(fromFile){
      stopifnot(!is.null(fileName) && !is.na(fileName) && nchar(fileName))# && file.exists(fileName)) # could be url
      cmd.list$file <- fileName
      GetTable <- read.table
    } else { 
      GetTable <- ReadTableFromClipboard
    }
    
    if(colClasses1$getActiveText() != "Default")  {
      #if(col.classes.arg != "character" && !is.null(row.names.arg)){
      cmd.list$nrows=5
      xx <- do.call(GetTable, cmd.list) # read it with row names NULL 
      NN <- dim(xx)[2]
      if(NN){
        cmd.list$colClasses <- c("character", rep(colClasses1$getActiveText(), NN-1))
      }
      cmd.list$nrows=nrows
   }

    xx <- do.call(GetTable, cmd.list) # read it with row names NULL 

    # ignore a percentage sign too
    #num.string.ignore = paste("[", NUMBER_FORMATS["thousands", dec1$getActiveText()], num.ignore.string1$getText(), "]", sep="")

    
    num.string.ignore = num.ignore.string1$getText()
    xx <- CheckStringsForNumeric(xx, num.string.ignore)
     # options for pasting into the table after reading
    other.opts = list(
      paste.row.labels = paste.row.labels$getActive(), 
      paste.col.labels = paste.col.labels$getActive()
    )

    if(!is.null(row.names.arg) && length(dim(xx))==2 && 0 < row.names.arg && row.names.arg <= ncol(xx)){
	  #dupes <- duplicated(xx[,row.names.arg])
      #if(any(dupes)) {
	  #	  xx <- xx[!dupes,,drop=F]
      #  #if(doMsg) message("Dropping ", sum(dupes), " row(s) from table to create unique row names")
      #}
      xr <- as.character(xx[,row.names.arg])
      dupes <- duplicated(xr)
      nas <- is.na(xr)
      if(any(nas)) xr[nas] <- ""
      if(any(dupes)) xr <- make.unique(xr)
      rownames(xx) <- xr
      xx <- xx[,-row.names.arg,drop=F]
    }
    return(list(xx=xx, cmd.list=cmd.list, other.opts=other.opts))
  } # end read_table

update_view = function(obj=NULL, label=NULL){
    tryCatch({    
      preview <- read_table(500)$xx
#fr3$setLabel(paste("Table Size:", nrow(preview), "R x", ncol(preview), "C"))
      obj <- gtkDfEdit(preview, size.request = c(400, 200), col.width=100, editable=FALSE)
      theObj = gtkContainerGetChildren(box3)[[1]]
      if(!is.null(theObj)) 
         gtkWidgetDestroy(theObj)
      gtkBoxPackEnd(box3, obj)
    }, error = function(e) {
      tb3 <- gtkTextViewNew()
      tb3$setWrapMode(GtkWrapMode['word_char'])
      tb3$setSizeRequest(400, 200)
      tb3$getBuffer()$setText(as.character(e$message))
      theObj = gtkContainerGetChildren(box3)[[1]]       
      if(!is.null(theObj)) 
         gtkWidgetDestroy(theObj)
      box3$packEnd(tb3)
    })
  }

  #gSignalConnect(preview.button, "clicked", update_view)

sapply(list(row.names0, header1, fill1, stringsAsFactors1, check.names1),
   function(x) tryCatch(gSignalConnect(x, "toggled", update_view), error=function(e) print(e)))

  sapply(list(num.ignore.string1, sep1, skip1, colClasses1, quote1, na.strings1),
   function(x) tryCatch(gSignalConnect(x, "changed", update_view), error=function(e) print(e)))

	gSignalConnect(dec1, "changed", 
	  function(...){
	    all.ign <- paste(c("[", NUMBER_FORMATS["thousands",], "]"), collapse="")
		txt <- gsub(all.ign, "", num.ignore.string1$getText())
        new.ign <- NUMBER_FORMATS["thousands", dec1$getActiveText()]
	    num.ignore.string1$setText(paste(txt, new.ign, sep=""))
	  } 
	)

	gSignalConnect(row.names0, "toggled", 
	  function(...){
		paste.row.labels$setActive(row.names0$getActive())
	  }
	)

	gSignalConnect(header1, "toggled", 
	  function(...){
		paste.col.labels$setActive(header1$getActive())
	  } 
	)

  fr3$add(box3)

  box0$packStart(fr3)
  box0$packStart(fr2, FALSE, FALSE, 0)
  
  dialog[["vbox"]]$packStart(box0, TRUE, TRUE, 10)
  dialog$showAll()
  update_view()

  xx <- NULL
  rv <- NULL
  if(dialog$run() == 1){# && nchar(entry1$getText())){
      w <- TransientWindow("Copying...")
      on.exit(w$destroy())      
  tryCatch({               
   rv <- read_table(doMsg=TRUE)
   xx <- rv$xx
   cmd.list = rv$cmd.list
          
   #aa <- basename(entry1$getText())
   #newname <- make.names(strsplit(aa, "[.]")[[1]][1])
   #assign(newname, xx, envir=.GlobalEnv)         
   
  }, error = function(e) {print(e)})
   dialog$destroy()
 } else {
   dialog$destroy()
 }
 return(rv)
  #return(xx)
}

    
# from gWidgets
gtkMenuPopupHack <- function (object, parent.menu.shell = NULL, parent.menu.item = NULL, 
  func = NULL, data = NULL, button, activate.time) {
  checkPtrType(object, "GtkMenu")
  if (!is.null(parent.menu.shell)) 
      checkPtrType(parent.menu.shell, "GtkWidget")
  if (!is.null(parent.menu.item)) 
      checkPtrType(parent.menu.item, "GtkWidget")
  if (!is.null(func)) 
     func <- as.function(func)
  button <- as.numeric(button)
  activate.time <- as.numeric(activate.time)
  w <- .RGtkCall("S_gtk_menu_popup", object, parent.menu.shell, 
      parent.menu.item, func, data, button, activate.time, 
      PACKAGE = "RGtk2")
  return(invisible(w))
}


###############################################################################
# Factor editor
###############################################################################
##########################################
# Blocking editor

BlockSizeHandler <- function(the.column, data){
  .local <- data$.local
  row.idx = data$row.idx
  col.idx = data$col.idx
  entry.frame <- data.frame(X=the.column)
  task <- list(
   list(func="ChangeCells", 
    arg = list(nf=entry.frame, row.idx=row.idx, col.idx=col.idx))
   )
  DoTaskWrapper(.local, task)
}

 # start.lvl = level name to start cycling at
DoBlockSize <- function(the.column, loc.window, handler, data, start.lvl=NULL){

  .BlockEnv <- new.env()

  .BlockEnv$the.column <- the.column 
  .BlockEnv$data <- data

  UpdateColumn <- function(block.size, dummy){
    total.len <- length(.BlockEnv$the.column)
    .BlockEnv$the.column <-  gl(length(lvls), block.size, total.len, labels=lvls)
    #handler(the.column, data)
  }
  
  if(!is.factor(the.column)) stop("Can't block over non-factors")  
  lvls <- unique(levels(the.column))
  if(!is.null(start.lvl)){
    if(!start.lvl%in%lvls) stop("start.level must be in levels")
    ww <- which(lvls == start.lvl)[1]
    if(ww > 1) lvls <- c(lvls[ww:length(lvls)], lvls[1:(ww-1)])
  }
  UpdateColumn(1)
  
  ok.handler <- function(block.size, dummy){
    handler(.BlockEnv$the.column, .BlockEnv$data) 
  }
 
  MakeSpinDialog(loc.window, title = "Blocking", label="Select Block Size", spin.handler=UpdateColumn, ok.handler=ok.handler, data=NULL)
}

MakeSpinDialog <- function(loc.window, title, label, spin.handler, ok.handler, data=NULL, maximum=100.0, minimum=1.0){
  window2 <- gtkWindowNew(show=F)
  window2$setTitle(title)
  window2$setPosition(GtkWindowPosition["center-on-parent"])  
  window2$setTransientFor(loc.window)
  window2$setModal(TRUE)
  box0 <- gtkVBoxNew(FALSE, 5)
  
  window2$add(box0)
  window2$setResizable(FALSE)
  window2$setPosition(GtkWindowPosition["center-on-parent"])  
  if(!is.null(loc.window)) window2$setTransientFor(loc.window)
  window2$setModal(TRUE)

  box1 <- gtkHBoxNew(FALSE, 5)
  box0$packStart(box1, FALSE, FALSE, 10)
  box1$packStart(gtkLabelNew(label), FALSE, FALSE, 10)
  spinner_adj <- gtkAdjustment(1, 1, maximum, minimum, 5.0)
  spinner <- gtkSpinButton(spinner_adj, 1.0, 0)
  spinner$setValue(1)
  box1$add(spinner)
  spinner$grabFocus()

  box2 <- gtkHBoxNew(FALSE, 5)
  box0$packEnd(box2, FALSE, FALSE, 0)
  button <- gtkButtonNewWithLabel("    OK    ")
  gSignalConnect(spinner, "value-changed", function(...){
    spin.handler(spinner$getValue(), data)
  })
  gSignalConnect(spinner, "key-press-event", function(obj, evt){
    if(evt[["keyval"]] == GDK_Return) {
      spinner$update()
      ok.handler(spinner$getValue(), data) 
      window2$destroy()
    }
    if(evt[["keyval"]] == GDK_Escape) {
      window2$destroy()
    }    
    FALSE
  })
  gSignalConnect(button, "clicked", function(handler, data){
    spinner$update()  
    ok.handler(spinner$getValue(), data)
    window2$destroy()
  })
  
  box2$packEnd(button, FALSE, FALSE, 0)
  button <- gtkButtonNewWithLabel("   Cancel   ")
  gSignalConnect(button, "clicked", function(...){
    window2$destroy()
    })
  box2$packEnd(button, FALSE, FALSE, 5)
  window2$show()
}

FactorEditorHandler <- function(the.column, col.idx, data){
  .local <- data
  #print(.FactorEnv$the.column)
#  the.column <- myDataTypeCoercions("factor", the.column)
#  x <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)
#  x[x=="NA"|nchar(x) == 0] <- NA
  entry.frame <- data.frame(X=the.column) 
  the.contrasts <- NULL
  if(length(levels(the.column)) > 1) the.contrasts <- contrasts(the.column)
  task <- list(list(func="ChangeCells", 
     arg = list(nf=entry.frame, col.idx=col.idx)),
    list(func="SetFactorAttributes", 
    arg = list(col.idx=col.idx, info=list(levels=levels(the.column), 
     contrasts=the.contrasts, contrast.name=attr(the.column, "contrast.name"), ordered=is.ordered(the.column))))
     )
  DoTaskWrapper(.local, task)
}

DoFactorEditor <- function(theFrame, toplevel, col.idx=integer(0), 
  handler=NULL, data=NULL){

  .FactorEnv <- new.env()
  .FactorEnv$col.idx <- col.idx

  UpdateView <- function(){
    contr <- contrast.coding[[which(sapply(rev(grb$getGroup()), gtkToggleButtonGetActive))]]
    x <- .FactorEnv$the.column
    if(ordered.checkbox$getActive() != is.ordered(x))
      x <- factor(x, levels = levels(x), ordered = ordered.checkbox$getActive())
    if(length(levels(x)) > 1) {
      contrasts(x) <- contr
      attr(x, "contrast.name") <- contr
    }
    handler(x, .FactorEnv$col.idx , data)
    .FactorEnv$the.column <- x
  }

  cell.edited <- function(cell, path.string, new.text, data){
    xx <- .FactorEnv$xx      
    if(!nchar(new.text))# || new.text%in%xx) 
      stop("New name must exist")# and be unique")
    checkPtrType(data, "GtkListStore")
    model <- data
    path <- gtkTreePathNewFromString(path.string)
    iter <- model$getIter(path)$iter

    i <- path$getIndices()[[1]]+1

      # editing the level names # 101010
    #zz <- theFrame[,.FactorEnv$col.idx]
    zz <- .FactorEnv$the.column# <- factor(.FactorEnv$the.column, levels=unique(xx))    
    lzz <- levels(zz)    
    ww <- which(xx[i]==lzz)
    
    if(length(ww)){
      lzz[ww] <- new.text
      levels(zz) <- lzz
      .FactorEnv$the.column <- zz
      UpdateView()            
    }
    xx[i] <- new.text
    model$set(iter, 0, new.text)
    .FactorEnv$xx <- xx
    UpdateLabel()    
  }
    
  add.item <- function(button, data) {
    xx <- .FactorEnv$xx      
     for(k in 1:(length(xx)[1]+1)){
       nl <- paste("Level", k, sep="_")
       if(!nl%in%xx) break
     }
     xx <- c(xx, nl)
     .FactorEnv$xx <- xx     
     iter <- model$append()$iter
     model$set(iter, 0, xx[length(xx)])
            
    .FactorEnv$the.column <- factor(.FactorEnv$the.column, levels=unique(xx))
     UpdateView()
     UpdateLabel()
   }


  remove.item <- function(widget, data)
  {
     xx <- .FactorEnv$xx   
     checkPtrType(data, "GtkTreeView")
     treeview <- data
     model <- treeview$getModel()
     selection <- treeview$getSelection()
   
     selected <- selection$getSelected()
     if (selected[[1]]){
        iter <- selected$iter
         
        path <- model$getPath(iter)
        i <- path$getIndices()[[1]]+1
        model$remove(iter)
        
        xx <- xx[-i]
        .FactorEnv$xx <- xx           
        path$prev()
        selection$selectPath(path)        
        .FactorEnv$the.column <- factor(.FactorEnv$the.column, levels=unique(xx))
        UpdateView()            
        UpdateLabel()        
      }
  }

  move.item.up <- function(widget, data)
  {
     xx <- .FactorEnv$xx
     checkPtrType(data, "GtkTreeView")
     treeview <- data

     model <- treeview$getModel()
     selection <- treeview$getSelection()
   
     selected <- selection$getSelected()
     if (selected[[1]])
       {
         iter <- selected$iter
           
         path <- model$getPath(iter)
         i <- path$getIndices()[[1]]+1
         if(i == 1) return()
         model$set(iter, 0, xx[i-1])
         path$prev()
         selection$selectPath(path)
         selected <- selection$getSelected()
         iter <- selected$iter
         model$set(iter, 0, xx[i])
         tmp <- xx[i-1]
         xx[i-1] <- xx[i]
         xx[i] <- tmp
         .FactorEnv$xx <- xx
        .FactorEnv$the.column <- factor(.FactorEnv$the.column, levels=unique(xx))
         UpdateView()
         UpdateLabel()                 
       }
  }

  move.item.down <- function(widget, data)
  {
     xx <- .FactorEnv$xx
     checkPtrType(data, "GtkTreeView")
     treeview <- data

     model <- treeview$getModel()
     selection <- treeview$getSelection()
   
     selected <- selection$getSelected()
     if (selected[[1]])
       {
         iter <- selected$iter
           
         path <- model$getPath(iter)
         i <- path$getIndices()[[1]]+1
         if(i == length(xx)) return()
         model$set(iter, 0, xx[i+1])
         gtkTreePathNext(path)
         selection$selectPath(path)
         selected <- selection$getSelected()
         iter <- selected$iter
         model$set(iter, 0, xx[i])
         tmp <- xx[i+1]
         xx[i+1] <- xx[i]
         xx[i] <- tmp
         .FactorEnv$xx <- xx         
        .FactorEnv$the.column <- factor(.FactorEnv$the.column, levels=unique(xx))
         UpdateView()                
         UpdateLabel()              
       }
  }

  edit.item <- function(widget, data) {
     checkPtrType(data, "GtkTreeView")
     treeview <- data

     model <- treeview$getModel()
     selection <- treeview$getSelection()
   
     selected <- selection$getSelected()
     if (selected[[1]])
       {
         iter <- selected$iter         
         path <- model$getPath(iter)         
         treeview$setCursorOnCell(path,
           treeview$getColumns()[[1]],
           treeview$getColumns()[[1]]$getCellRenderers()[[1]],
           TRUE)
         UpdateLabel()           
       }
  }
   
  contrast.coding <- list(
  "Treatment (default) - contrasts each level with the first.\nThe first level is omitted." = "contr.treatment", 
  "Helmert - contrast the second level with the first, the \nthird with the average of the first two, and so on." = "contr.helmert",  
  "Polynomial - contrasts based on orthogonal polynomials." = "contr.poly", 
  "Sum - sum to zero contrasts." = "contr.sum", 
  "SAS - treatment contrasts with base level set to be the\nlast level of the factor." = "contr.SAS")
  
  UpdateLabel <- function(){
    xx <- .FactorEnv$xx  
    contr <- contrast.coding[[which(sapply(rev(grb$getGroup()), gtkToggleButtonGetActive))]]

    contrF <- paste("Control (first level) is: <b>", xx[1], "</b>", sep="")
    contrL <- paste("Control (last level) is: <b>", rev(xx)[1], "</b>", sep="")
    
    contr.msg <- list(contr.treatment = contrF, contr.helmert = contrL, 
      contr.poly = "L and Q terms (see MASS, p156)", 
      contr.sum = contrL, contr.SAS = contrL)[[contr]]
    
    lab1$setMarkup(contr.msg)
    
  }  

  .FactorEnv$the.column <- theFrame[,.FactorEnv$col.idx]
  .FactorEnv$col.original <- .FactorEnv$the.column  
  contr <- attr(.FactorEnv$the.column, "contrast.name")    

  box4 <- gtkVBoxNew(FALSE, 5)  
  grb <- gtkRadioButtonNewWithLabel(NULL, label=names(contrast.coding)[1])
  grb$setActive(TRUE)
  box4$packStart(grb, FALSE, FALSE, 0)
  for(ii in 2:length(names(contrast.coding))){
    grb <- gtkRadioButtonNewWithLabel(group=grb$getGroup(), 
      label=names(contrast.coding)[ii])
    if(!is.null(contr) && contr==unlist(contrast.coding)[ii]) grb$setActive(TRUE)
    gSignalConnect(grb, "toggled", function(grb, ...) if(grb$getActive()) UpdateLabel())
    box4$packStart(grb, FALSE, FALSE, 0)

  }

  data.factors <- which(GetClasses(theFrame) == "factor")
  if(!length(colnames(theFrame)[data.factors])) stop("No data columns are of type \"factor\"")

  # called with no selection  
  if(!length(col.idx)) {
    col.idx <- data.factors[1]
    .FactorEnv$col.idx <- col.idx
  } else if(length(col.idx==1) && !is.factor(theFrame[,col.idx])) {
    stop(paste("Data column:", col.idx, "is not of type \"factor\""))
  }
    
#  window <- gtkWindowNew(show=F)
  dialog <- gtkDialog("Factor Editor", NULL, "modal", "gtk-ok", 1, "gtk-cancel", 0, show = T)    
  dialog$setPosition(GtkWindowPosition["center-on-parent"])
  dialog$setTransientFor(toplevel) 
      
#  if(!is.null(toplevel)) {
#    window$setTransientFor(toplevel)  
#    window$setModal(TRUE)
#  }  

  dialog$setTitle("Factor Editor")
#  window$setBorderWidth(5)
  box0 <- gtkVBoxNew(FALSE, 5)
  dialog[["vbox"]]$add(box0)  
#  window$add(box0)

  fr0 <- gtkFrameNew(label="Column Selection")

  cb1 <- gtkComboBoxNewText()  
  cb1$show()
  cb1["width-request"] <- 75
  
  for (item in colnames(theFrame)[data.factors]) # omit "rows" heading
    cb1$appendText(item)  

  theIdx <- which(data.factors %in% col.idx)-1
  cb1$setActive(theIdx)
  fr0$add(cb1)  
  box0$packStart(fr0, FALSE, FALSE, 5)    
  fr1 <- gtkFrameNew(label="Factor Level Order")
  box0$packStart(fr1, TRUE, TRUE, 5)
  box1 <- gtkHBoxNew(FALSE, 5)
  fr1$add(box1)
  box1.5 <- gtkVBoxNew(FALSE, 0)  
  box1$packStart(box1.5, TRUE, TRUE, 5)

  box1.6 <- gtkHBoxNew(FALSE, 0)  
  box1.5$packStart(box1.6, FALSE, TRUE, 5)  
  lab1 <- gtkLabelNew("")  
  box1.6$packStart(lab1, FALSE, FALSE, 5)
                                            
  MakeModel <- function(){
    col.original <- .FactorEnv$col.original
    if(!is.factor(col.original)) {
      col.original <- as.factor(integer(0))
      #stop("Can't edit non-factors")
      warning("Can't edit non-factors")    
      }
    .FactorEnv$xx <- na.omit(cbind(levels(col.original)))
#print(col.original)
#print(levels(col.original))
    #.FactorEnv$xx <- cbind(levels(col.original))
    model <- gtkListStoreNew("gchararray")
    xx <- .FactorEnv$xx
    sapply(xx, function(x) model$set(model$append()$iter, 0, x))
    UpdateLabel()  
    return(model)
  }
    
  sw <- gtkScrolledWindowNew(NULL, NULL)
  sw$setPolicy("automatic", "automatic")
  box1.5$packStart(sw, TRUE, TRUE, 0)
  treeview <- gtkTreeViewNew()  
  
  model <- MakeModel()  
  treeview$setModel(model)
  
  gSignalConnect(cb1, "changed", function(widget){
      new.idx <- which(widget$getActiveText()==colnames(theFrame))
      .FactorEnv$col.idx <- new.idx
      .FactorEnv$the.column <- theFrame[,new.idx]
      .FactorEnv$col.original <- .FactorEnv$the.column
      model <- MakeModel()
      treeview$setModel(model)
  })   
  
  treeview$setRulesHint(TRUE)
  treeview$setHeadersVisible(FALSE)
  treeview$setEnableSearch(FALSE)  
  treeview$getSelection()$setMode("single")
  renderer <- gtkCellRendererTextNew()
  renderer$setData("column", 0)

  renderer['editable-set'] <- TRUE
  renderer['editable'] <- TRUE
  treeview$insertColumnWithAttributes(-1, "Name", renderer, text = 0)

  sw$setShadowType(as.integer(1))
  sw$add(treeview)
  sw$setSizeRequest(300, -1)
#  window$setResizable(FALSE)
  dialog$setResizable(FALSE)

  gSignalConnect(renderer, "edited", cell.edited, model)

  box2 <- gtkVBoxNew(FALSE, 5)
  box1$packStart(box2, FALSE, FALSE, 5)

  button1 <- gtkButtonNewWithLabel("Move Up")
  box2$packStart(button1, FALSE, FALSE, 0)

  button2 <- gtkButtonNewWithLabel("Move Down")
  box2$packStart(button2, FALSE, FALSE, 0)

  button3 <- gtkButtonNewWithLabel("Edit Name")
  box2$packStart(button3, FALSE, FALSE, 0)

  button4 <- gtkButtonNewWithLabel("Add Level")
  box2$packStart(button4, FALSE, FALSE, 0)

  button5 <- gtkButtonNewWithLabel("Remove Level")
  box2$packStart(button5, FALSE, FALSE, 0)
  
  gSignalConnect(button1, "clicked", move.item.up, data=treeview)
  gSignalConnect(button2, "clicked", move.item.down, data=treeview)
  gSignalConnect(button3, "clicked", edit.item, data=treeview)
  gSignalConnect(button4, "clicked", add.item, model)
  gSignalConnect(button5, "clicked", remove.item, treeview)   
  
  ordered.checkbox <- gtkCheckButtonNewWithLabel(label="Levels Are Ordered")
  ordered.checkbox$setActive(is.ordered(.FactorEnv$the.column))
  box0$packStart(ordered.checkbox)

  expander <- gtkExpanderNew(label="Factor Contrasts (For Experts)")
  fr2 <- gtkFrameNew(label="Contrast Coding")
  #box0$packStart(fr2, FALSE, FALSE, 5)
  box0$packStart(expander, FALSE, FALSE, 5)
  expander$add(fr2)
  fr2$add(box4)

#  box6 <- gtkHBoxNew(FALSE, 5)
   
  UpdateLabel()    
  if(dialog$run() == 1){
    UpdateView()
    dialog$destroy()
  } else {
    .FactorEnv$the.column <- .FactorEnv$col.original  
    UpdateView()
    dialog$destroy()    
  }
  
#  button <- gtkButtonNewWithLabel("  Cancel  ")
#  gSignalConnect(button, "clicked", function(...) {
#   .FactorEnv$the.column <- .FactorEnv$col.original
#    UpdateView()
#    window$destroy()
#    })
#  box6$packEnd(button, FALSE, FALSE, 5)
#  
#  box0$packEnd(box6, FALSE, FALSE, 0)
#  button <- gtkButtonNewWithLabel("    OK    ")
#  gSignalConnect(button, "clicked", function(...){    	
#    UpdateView()
#    window$destroy()
#  })
#  box6$packEnd(button, FALSE, FALSE, 0)
#  
#  window$show()
}

###############################################################################
# End factor editor
###############################################################################

#TableEntryDialog <- function(.localenv){

MakeRadioDialog <- function(title, items, handler=NULL, data=NULL){
  window <- gtkWindowNew(show=F) 
  if(!is.null(data$.local)) {
    window$setPosition(GtkWindowPosition["center-on-parent"])  
    window$setTransientFor(data$.local$toplevel)  
    window$setModal(TRUE)
  }  
  window$setTitle(title)
  window$setBorderWidth(5)
  box0 <- gtkVBoxNew(FALSE, 5)
  window$add(box0)
  window$setResizable(FALSE)
  window$resize(window$getSize()$width, 1)

  xx <- list()
  for(jj in 1:length(items)){
    fr1 <- gtkFrameNew(label=names(items[[jj]]))
    xx[[jj]] <- MakeRadiobuttonGroup(items[[jj]][[1]])
    fr1$add(xx[[jj]]$box)
    box0$add(fr1)
  }

  fr2 <- gtkVBoxNew()
  box18000 <- gtkHBoxNew(FALSE, 5)
  fr2$add(box18000)
  button <- gtkButtonNewWithLabel(" Cancel ")
  gSignalConnect(button, "clicked", function(...){
    window$destroy()
  })

  box18000$packEnd(button, FALSE, FALSE, 0)
  button <- gtkButtonNewWithLabel("   OK   ")
  box18000$packEnd(button, FALSE, FALSE, 0)
  gSignalConnect(button, "clicked", function(...){
    results<- sapply(xx, function(item) 
        which(rev(sapply(item$grb$getGroup(), gtkToggleButtonGetActive))))
    res <- list()
    for(jj in 1:length(items)) 
      res[[names(items[[jj]])]] <- items[[jj]][[1]][results[jj]]
    if(!is.null(handler)) handler(res, data)
    window$destroy()
  })
  box0$packStart(fr2, FALSE, FALSE, 0)  
  window$show()
}

TABLE_IMPORT_OPTIONS <- list(
    list("Data Type" = c("Numeric Data", "Character Data")), 
    list("Name Options" = c("Row And Column Names", "Row Names Only", "Column Names Only", "No Row Or Column Names"))
    )
 
TABLE_IMPORT_OPTIONS_FUNC <- function(results){
  do.colnames <- F; do.rownames <- F
  if(results$"Name Options"%in%c("Row And Column Names", "Row Names Only"))
    do.rownames <- T
  if(results$"Name Options"%in%c("Row And Column Names", "Column Names Only"))
    do.colnames <- T
  colClasses <- NA
  if(results$"Data Type"%in%c("Numeric Data"))
    colClasses <- "numeric"
  if(results$"Data Type"%in%c("Character Data"))
    colClasses <- "character"
  list(colClasses=colClasses, do.rownames=do.rownames, do.colnames=do.colnames)
}

###############################################################################
# Sort dialog
###############################################################################

MakeRadiobuttonGroup <- function(labels, theChoice=1){
  box3 <- gtkVBoxNew()
  grb1 <- gtkRadioButtonNewWithLabel(NULL, label=labels[1])
  if(theChoice == 1) grb1$setActive(TRUE)
  box3$add(grb1)
  for(jj in 2:length(labels)){
    grb1 <- gtkRadioButtonNewWithLabel(group=grb1$getGroup(),label=labels[jj])
    if(jj == theChoice) grb1$setActive(TRUE)
    box3$add(grb1)
  }
  return(list(box=box3, grb=grb1))
}

MakeComboEntry <- function(items, box1){
  box2 <- gtkHBoxNew()
	cb1 <- gtkComboBoxNewText()  


	cb1$show()
	cb1["width-request"] <- 100		
	for (item in items) # omit "rows" heading
	  cb1$appendText(item)  
	cb1$setActive(0)
  box5 <- gtkVBoxNew()		
  box5$packStart(cb1, TRUE, FALSE, 0)    
  box2$packStart(box5, TRUE, TRUE, 10)
  xx1 <- MakeRadiobuttonGroup(c("Low to high", "High to low"))
  box3 <- xx1$box
  grb1 <- xx1$grb
  xx2 <- MakeRadiobuttonGroup(c("Default", "Character", "Numeric"))
  box4 <- xx2$box
  grb2 <- xx2$grb
	box2$packStart(box4, FALSE, FALSE, 0)    
	box2$packStart(box3, FALSE, FALSE, 0)
  box1$packStart(box2, FALSE, FALSE, 0)
  return(list(col=cb1, ord=grb1, typ=grb2, box=box2))
}
  
# Returns the ordering of the table
DoSortDialog <- function(theFrame, handler, .localenv){
  window <- gtkWindowNew(show=F) 
  if(!is.null(.localenv$toplevel)) {
    window$setPosition(GtkWindowPosition["center-on-parent"])    
    window$setModal(TRUE)
    window$setTransientFor(.localenv$toplevel)      
  }
  window$setTitle("Sort Options")
  window$setBorderWidth(5)
  box0 <- gtkVBoxNew(FALSE, 5)
  window$add(box0)
  items <- colnames(theFrame)
  fr0 <- gtkFrameNew(label="Order By Column")
  box1 <- gtkVBoxNew(FALSE, 5)
  box0$add(fr0)
  fr0$add(box1)

  .sl <- new.env()
  .sl$theList <- list()
  .sl$theList[[length(.sl$theList)+1]] <- MakeComboEntry(items, box1)
  fr1 <- gtkFrameNew(label="Add/Remove Columns")
  box9000 <- gtkHBoxNew(FALSE, 5)
  fr1$add(box9000)
  button <- gtkButtonNewWithLabel("- Column")
  gSignalConnect(button, "clicked", function(obj, ...){
    if(length(.sl$theList)<2) return(FALSE)
    .sl$theList[[length(.sl$theList)]]$box$destroy()
    .sl$theList[[length(.sl$theList)]] <- NULL
    window$resize(window$getSize()$width, 1)
  })
  box9000$packEnd(button, FALSE, FALSE, 0)
  button <- gtkButtonNewWithLabel("+ Column")
  gSignalConnect(button, "clicked", function(obj, data=.sl){
    .sl$theList[[length(.sl$theList)+1]] <- MakeComboEntry(items, box1)
  })
  box9000$packEnd(button, FALSE, FALSE, 0)
  box0$packStart(fr1, FALSE, FALSE, 0)

  fr2 <- gtkVBoxNew()
  box18000 <- gtkHBoxNew(FALSE, 5)
  fr2$add(box18000)
  button <- gtkButtonNewWithLabel(" Cancel ")
  gSignalConnect(button, "clicked", function(obj, data=.sl){
    window$destroy()
  })
  box18000$packEnd(button, FALSE, FALSE, 0)
  button <- gtkButtonNewWithLabel("   OK   ")
  box18000$packEnd(button, FALSE, FALSE, 0)
  gSignalConnect(button, "clicked", function(obj, data=.sl){
		opts <- lapply(.sl$theList, function(item){
				list(col = item$col$getActiveText(),
  				ord = which(rev(sapply(item$ord$getGroup(), gtkToggleButtonGetActive))),
  				typ = which(rev(sapply(item$typ$getGroup(), gtkToggleButtonGetActive))))
		})
		dataset.order <- do.call("order", lapply(opts, function(item){
		  xx <- theFrame[[item$col]]
		  if(item$typ == 1) {
        xrank <- xtfrm(xx)
      } else if (item$typ == 2) {
		    xrank <- xtfrm(as.character(xx))
      } else if (item$typ == 3) {
		    xrank <- xtfrm(as.numeric(xx))
      } else {
        stop("Sort error")
      }
			(c(-1, 1)[(item$ord==1)+1])*xrank
		}))
    handler(dataset.order, .localenv)
    window$destroy()    
  })
  box0$packStart(fr2, FALSE, FALSE, 0)  
  window$setResizable(FALSE)
  window$show()
}
###############################################################################
# End sort dialog
###############################################################################
                

DoTaskWrapper <- function(.local, task, do.undo = TRUE){
  #tryCatch({
#print("DoTaskWrapper")
    rv <- DoTask(.local, .local$theFrame, task, .local$changed.handler)
    rows.changed <- NULL
      # Special for changecells to minimize edit time
    if(length(task) == 1 && task[[1]]$func == "ChangeCells")
      rows.changed <- task[[1]]$arg$row.idx
    #print("Updating here...")
    #print(dim(.local$model))
    UpdateDfEditor(.local, rv$df, rows.changed)    
    #print("Updating finished")
    #print(dim(.local$model))
      
    if(do.undo){
      .local$undoStack[[length(.local$undoStack)+1]] <- rv$undo         
      
      if(object.size(.local$undoStack) > STACK_LIMIT){
        warning("Stack full")    
        jj <- 0

        while(object.size(.local$undoStack) > STACK_LIMIT && length(.local$undoStack))
          .local$undoStack[[jj <- jj + 1]] <- NULL
        if(object.size(.local$undoStack) > STACK_LIMIT){ 
          warning("This edit is too large to support undo")
          .local$undoStack <- list()
        }
      }
    }
  #}, error = function(e) {
  #  warning("An error has occurred in performing the task")
  #  old_warning(e)
  #})      
}


  # update the scroll in one treeview from another
ViewOnScrollChanged <- function(obj, data){
  sw2 <- data$sw2
  .local <- data$.local
  # Take over the event loop! 
  # See http://wiki.laptop.org/go/PyGTK/Smooth_Animation_with_PyGTK
  while(gtkEventsPending())
    gtkMainIteration()
  gtkAdjustmentSetValue(sw2, gtkAdjustmentGetValue(obj))
  .local$allow.key.flag <- TRUE  # For paging, however, this has a problem...
}


############################################################
# Turns out that scroll_row_timeout has no horizontal autoscroll.
# I based this on gtktreeview::gtk_tree_view_vertical_autoscroll, no offset
HScroll <- function(data){ 
  .local <- data
    # which side of the middle are we?
  view <- .local$view
  sw.ha <- .local$sw.ha
  ptr <- view$getBinWindow()$getPointer()
  vr <- view$getVisibleRect()$visible.rect
  sw.ha.value <- sw.ha$value
  direction <- ifelse (ptr$x - sw.ha.value <= vr$width/2, -1, 1) 

  val <- sw.ha.value + direction*sw.ha[["step_increment"]]
  if(0 <= val && val <= sw.ha[["upper"]] - sw.ha[["page_size"]]) {
    sw.ha$setValue(val)
  } else if (val < 0) {
    sw.ha$setValue(0)
  } else {
    sw.ha$setValue(sw.ha[["upper"]] - sw.ha[["page_size"]])
  }
  TRUE
}

# Bind this to button-release, focus-out and enter-notify
RemoveHScrollTimeout <- function(obj, event, data){
  .local <- data
  try({
    gSourceRemove(.local$hScrollTimeout) 
    .local$doingHScroll <- FALSE
  }, silent=T) 
  FALSE
}

AddHScrollTimeout <- function(obj, event, data){
  .local <- data
  sw.ha <- .local$sw.ha
  view <- .local$view
  if (!.local$doingHScroll){
    ptr <- obj$getBinWindow()$getPointer()
    if (as.flag(ptr$mask) & GdkEventMask['button-press-mask']){
    x <- ptr$x - sw.ha$getValue()
    y <- ptr$y
    vr <- view$getVisibleRect()$visible.rect
    h <- vr$height
    w <- vr$width
    z1 <- y > h/w*x 
    z2 <- y > h - h/w*x
    if((z1 && !z2) || (!z1 && z2)){ 
      .local$hScrollTimeout <- gTimeoutAdd(.local$SCROLL_ROW_TIMEOUT, HScroll, data=.local)
      .local$doingHScroll <- TRUE
    }
    }} #if, if
    TRUE
}
############################################################
# End autoscroll

PaintSelectionOnTimeout <- function(.local, selection=.local$ss, widget=.local$view){
  # We want to add a timeout to the LAST key pressed.
  try({
    gSourceRemove(.local$do.paint.timeout) 
    .local$do.paint.timeout <- NULL
  }, silent=T)

  .local$do.paint.timeout <- gTimeoutAdd(100, function(){
  	.local$do.paint <- TRUE
    UpdateSelectionRectangle(.local, selection, widget)
    return(FALSE)
	})
}
  
####
############################################################
MoveCursor <- function(widget, direction, .local, stat=as.integer(0)){
  cursor.info <- gtkTreeViewGetCursor(widget)
  path <- cursor.info$path
  column <- cursor.info$focus.column
  allColumns <- .local$allColumns
  page.flag <- FALSE
  if(is.null(path)){
    row.idx = 1
    new.idx = 2
  } else {
    row.idx <- as.integer(gtkTreePathGetIndices(path))+1
    new.idx <- GetColIdx(cursor.info$focus.column)  
  }

  col.idx <- new.idx
  if(IsShift(stat)){
    if(is.null(.local$selections))
      .local$selections <- list(list(start=c(row.idx=row.idx, col.idx=col.idx)))
  }

if (direction %in% c("page_down", "page_up")){
    page.flag <- TRUE
    visible_row_range <-  as.integer(diff(sapply(gtkTreeViewGetVisibleRange(widget)[-1], gtkTreePathGetIndices)))
    row.idx <- row.idx + 
       ifelse(direction=="page_down", 1, -1)*visible_row_range 
    if(direction=="page_down") row.idx <- row.idx-1
    else if(direction=="page_up") row.idx <- row.idx
    if(row.idx < 1) 
      row.idx <- 1
    else if (row.idx > nrow(.local$theFrame)) 
      row.idx <- nrow(.local$theFrame)
    path <- gtkTreePathNewFromString(row.idx-1)
    gtkTreeViewScrollToCell(widget, path, column=NULL, use.align=TRUE, row.align=0)
    gtkTreeViewSetCursorOnCell(widget, path, NULL, NULL, FALSE)
  } else if (direction == "right"){
     if( new.idx < length(allColumns)){
       new.idx <- new.idx+1
     } else {
       if(.local$editable) insert.dialog(.local, row.idx=row.idx-1, col.idx=new.idx, insert.type ="Columns", win=.local$toplevel)
     } 
  } else if (direction == "left") {
    if(1 < new.idx)
      new.idx <- new.idx-1
  } else if (direction == "down") {
    if (row.idx < nrow(.local$theFrame)){
      gtkTreePathNext(path)
      row.idx <- row.idx+1
    } else {
       if(.local$editable) insert.dialog(.local, row.idx=row.idx-EXTRA_ROW, col.idx=new.idx, insert.type ="Rows", win=.local$toplevel)
    }            
  } else if (direction == "up") {
    gtkTreePathPrev(path)
    if(row.idx > 1) row.idx <- row.idx-1
  } else if (direction%in%c("home", "end")) {
    row.idx <- ifelse(direction=="home", 1, nrow(.local$theFrame))
    path <- gtkTreePathNewFromString(row.idx-1)
    gtkTreeViewScrollToCell(widget, path, column=NULL, use.align=TRUE, row.align=0)
    gtkTreeViewSetCursorOnCell(widget, path, NULL, NULL, FALSE)
  }

  if(as.flag(stat) & GdkModifierType['shift-mask']){
    stopifnot(length(.local$selections)>0)
    .local$selections[[1]]$end <- c(row.idx=row.idx, col.idx=new.idx)
    gdkWindowInvalidateRect(.local$viewGetBinWindow, NULL, FALSE)
    PaintSelectionOnTimeout(.local)
  } else {
      # this won't get called, we caught these in ViewKeyPress
    if (as.flag(stat) & GdkModifierType['control-mask']) {
    } else {
      if(!is.null(.local$selections)){
        .local$selections <- NULL
        UpdateSelectionRectangle(.local)
      }
      UpdateColumnSelection(.local, new.idx)          
    }
  }
      
  if(!page.flag){
    new.col <- allColumns[[new.idx]]$column
    renderer <- allColumns[[new.idx]]$renderer
    gtkTreeViewSetCursorOnCell(widget, path, new.col, renderer, FALSE)
  }

    # ScrollToCell seems to stop working after we add a new column.
    # This routine is from gtktreeview::gtk_tree_view_scroll_to_cell
  if(direction%in%c("left", "right")){
    cell_rect <- gtkTreeViewGetBackgroundArea(widget, path, new.col)$rect
    vis_rect <- gtkTreeViewGetVisibleRect(widget)$visible.rect
    dest_x <- vis_rect$x
    dest_y <- vis_rect$y
    if (cell_rect$x < vis_rect$x)
      dest_x <- cell_rect$x
    if (cell_rect$x+cell_rect$width > vis_rect$x + vis_rect$width)
      dest_x <- cell_rect$x + cell_rect$width - vis_rect$width
    gtkTreeViewScrollToPoint(widget, dest_x, dest_y)
  }

selections <- list(
    list(start=c(row.idx=row.idx, col.idx=new.idx),
         end=c(row.idx=row.idx, col.idx=new.idx)))
   # update the selection only if we didn't select anything else
  if(!length(.local$selections)) DoUserCall("Selection", list(selections=selections), .local)

  #} # if 0
  return(TRUE)
}  

RowNamesClearContents <- function(row.idx, .local){
  if(length(row.idx)==0) return(FALSE)
  nf <- .local$theFrame[row.idx,, drop=F]
  nf[,-1] <- ""
  task <- list(list(func="ChangeCells", 
    arg = list(nf=nf, row.idx=row.idx)))
  DoTaskWrapper(.local, task)
}  

GetColIdx <- function(column){
  #tryCatch(
    as.integer(column["title"])#,
  #  error = function(e) integer(0))
}

RowNamesKeyPress <- function(widget, event, data) {
  .local <- data
  keyval <- event[["keyval"]] 
  stat <- as.flag(event[["state"]])    
  #print(keyval)
  if (CtrlLetter(keyval, stat, GDK_z)){
    #print("Ctrl-z")
    DoUndo(.local)
    return(TRUE)
  } else if(keyval == GDK_Delete && .local$editable) {
    #row.idx <- GetSelectedRows(.local$view.rn, .local)
    #RowNamesClearContents(row.idx, .local)
    DeleteSelection(.local)
    #.local$allow.key.flag <- TRUE                  
    return(TRUE)
  } else if (CtrlLetter(keyval, stat, GDK_c)) {
     # block transient interactions with a popup
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())   
      nf <- GetSelectedCells(.local)
     
      if(!is.null(nf)) CopyToClipboard(nf, 
        do.rownames=!RowNamesAreNull(.local$theFrame), do.colnames=FALSE)
      return(TRUE)
      #return(TRUE)
  } else if (CtrlLetter(keyval, stat, GDK_v) && .local$editable){
	   dat <- ReadFromClipboard(.local, row.names=1, sep="\t", stringsAsFactors=F)$xx # character vector
  cursor.info <- gtkTreeViewGetCursor(widget)  
  row.idx <- as.integer(gtkTreePathGetIndices(cursor.info$path))+1
	  if(!is.null(dat)){
		if(!length(row.idx)) row.idx <- 1
		col.idx <- 1
		row.idx <- row.idx[1]
		task <- GetTaskPasteIn(.local$theFrame, dat, row.idx, col.idx+COLUMN_OFFSET, do.colnames=F, do.rownames = F)
		DoTaskWrapper(.local, task)
   }
  }
           
    # ignore last row
  cursor.info <- gtkTreeViewGetCursor(widget)  

  if(keyval%in%myValidNavigationKeys){
    #cat("*")
    if( .local$scrollRowNames) { 
      .local$scrollRowNames <- FALSE
     while(gtkEventsPending())
            gtkMainIteration()
      gtkPropagateEvent(widget, event)
      #gtkTreeSelectionUnselectAll(.local$ss)
      gtkAdjustmentSetValue(.local$sw.view.va, gtkAdjustmentGetValue(.local$sw.rn.va))
      .local$scrollRowNames <- TRUE 
      return(TRUE)
   }
  }


  return(FALSE)
}

  # Apply command to range [sr, sc]
CommandData <- function(.local, sr, sc) {
  command.dialog <- list(
    title = "Apply Command",
    label = "Apply a command or function to cell selection.\nSelection is stored as variable x.\n",
    txt.stringItem = "", label = "Enter Code (Ctrl-Enter To Run)", multi=T, signal.on.startup = F,
      signal = c("default", "run.it"),    
    do.apply.trueFalseItem = FALSE, label = "Apply Function",
      signal = c("default", "toggle.sensitive", "apply.over"),
      apply.over.radiobuttonItem = c(1,2), item.labels = c("Rows", "Columns"), label = "Over", indent=10,    
      signal = c("default", function(apply.over, position){
        if(get.value(apply.over) == 1) set.value(position, "Right")
        if(get.value(apply.over) == 2) set.value(position, "Bottom")        


      }, "position"),
    insert.trueFalseItem = TRUE, label = "Put Output In Table",
      signal = c("default", "toggle.sensitive", "position"),    
    position.radiobuttonItem = c("Replace", "Bottom", value="Right"), label = "Position to Insert", indent=10,
    otherTable.stringItem = "", label = "Put Output In Another Table", sensitive=F
  )
  
  command <- function(txt, insert, position, do.apply, apply.over, otherTable){
   txt <- sub("^=(.*?)$", "\\1", txt, perl=TRUE)
   try.short.fn <- paste("function(x) {", txt, "}")             
   dat <- NULL
   .e = new.env()   
   .e$output <- NULL
   x <- as.matrix(.local$theFrame[sr, sc+1, drop=F])   
   .e$xx <- NULL
   tryCatch({
     .e$output <- eval(parse(text=txt), envir=.GlobalEnv)
     if(!is.function(.e$output)){ # not a function, maybe something like 1:10                        
        .e$xx <- .e$output
        if(findxInParseTree(txt) && exists("x", envir = .GlobalEnv)) 
          quick_message("  Warning! A variable called 'x' exists already and will be used.  ", win=.local$toplevel)
     }
   }, error = function(e) { # try putting "function(x)" on the front and see if that works
     tryCatch({    
       .e$output <- eval(parse(text=try.short.fn), envir=.GlobalEnv)
     }, error = function(e) quick_message("  Sorry, that didn't make sense to R  ^_^;;  ", win=.local$toplevel))
  })  
  output <- .e$output    
  xx <- .e$xx  
  if(is.function(output)) {
    if(do.apply){
      xx <- apply(x, apply.over, output) 
    } else {
      xx <- eval(make.call(output, x))  
    }
  }
  dat <- NULL
  tryCatch({  
    if(!is.null(xx) && length(xx) && !(is.list(xx) && !is.data.frame(xx))){
      if (position == "Replace"){
          dat <- array(xx, c(length(sr), length(sc)))
          if(is.atomic(dat)) DoTaskWrapper(.local,list(list(func="ChangeCells", args=list(nf=dat, row.idx=sr, col.idx=sc+1))))
      } else if (position == "Right"){
         dat <- t(ragged.cbind(xx))
         new.dat <- array(NA, c(dim(.local$theFrame)[1], dim(dat)[2]))
         new.dat[sr, 1:dim(dat)[2]] <- dat
         #new.dat <- data.frame(dat, stringsAsFactors=F)         
         task <- list(list(func="InsertColumns",          
                 arg = list(nf = new.dat, col.idx=max(sc)+1+rep(1, dim(new.dat)[2]))))
        DoTaskWrapper(.local, task)                       
      } else if (position == "Bottom") {
         dat <- ragged.cbind(xx)
         new.dat <- array(NA, c(dim(dat)[1], dim(.local$theFrame)[2]))
         new.dat[1:dim(dat)[1], sc+1] <- dat
         new.dat <- data.frame(new.dat, stringsAsFactors=FALSE, check.names=FALSE)
         task <- list(list(func="InsertRows",                   
                 arg = list(nf = new.dat, row.idx=max(sr)+rep(1, dim(new.dat)[1]))))
        DoTaskWrapper(.local, task)      
      }
    }
    #}, error = function(e) quick_message("  Output couldn't be put into cells  ", win=.local$toplevel))
    }, error = function(e) cat("Output couldn't be put into cells\n"))    
  }  

  run.dialog(func=command, dlg.list=command.dialog, win=.local$toplevel)
  
 }
    
insert.dialog <- function(.local, row.idx=NULL, col.idx=NULL, insert.type ="Rows", win=NULL) {
  
  choice.list <- list("Columns" = "InsertNAColumns", "Rows" = "InsertNARows")

  dialog <- gtkDialog("Insert Rows/Columns", NULL, c("modal", "destroy-with-parent"), "gtk-ok", 1, "gtk-cancel", 0, show = F)
  #dialog$setDecorated(FALSE)

  do.insert <- function(){                
    choice.idx <- which(sapply(rev(grb$getGroup()), gtkToggleButtonGetActive))
    choice <- names(choice.list)[choice.idx]
    func <- choice.list[[choice.idx]]
    n <- as.integer(spinner$getValue())
    stopifnot(is.integer(n) && n > 0)
    arg <- list()
    if(choice == "Rows"){
      idx <- row.idx
      stopifnot(!is.null(idx))
      theIdx <- rep(idx+1, n)            
      arg$row.idx <- theIdx      
    } else {
      idx <- col.idx                     
      stopifnot(!is.null(idx))      
      theIdx <- rep(idx+1, n)                  
      arg$col.idx <- theIdx

    }
    DoTaskWrapper(.local,list(list(func=func,arg=arg)))
  }

  box0 <- gtkVBoxNew(FALSE, 5)
  
  grb.keypress <- function(obj, evt){
    keyval <- evt[["keyval"]]
    if(keyval == GDK_Return){
      do.insert() 
      dialog$destroy()
      .local$view$grabFocus()
      return(TRUE)
    }
    if(keyval == GDK_Right){
      spinner$spin(GtkSpinType["step-forward"], 1)
      return(TRUE)
    }    
    if(keyval == GDK_Left){
      spinner$spin(GtkSpinType["step-backward"], 1)
      return(TRUE)
    }        
    FALSE
  }

  fr2 <- gtkFrameNew(label="Insert")   
  box4 <- gtkVBoxNew(FALSE, 5)
  grb <- gtkRadioButtonNewWithLabel(NULL, label=names(choice.list)[1])
  gSignalConnect(grb, "key-press-event", grb.keypress)  
  grb$setActive(TRUE)
  box4$packStart(grb, FALSE, FALSE, 0)
  grb <- gtkRadioButtonNewWithLabel(group=grb$getGroup(), label=names(choice.list)[2])
  gSignalConnect(grb, "key-press-event", grb.keypress)
  if(!is.null(insert.type) && insert.type==names(choice.list)[2]) grb$setActive(TRUE)
  box4$packStart(grb, FALSE, FALSE, 0)
  fr2$add(box4)
  box0$add(fr2)


  box1 <- gtkHBoxNew(FALSE, 5)
  box1$packStart(gtkLabelNew("Number To Insert"), FALSE, FALSE, 10)
  spinner_adj <- gtkAdjustment(1, 1, .Machine$integer.max, 1, 10.0)
  spinner <- gtkSpinButton(spinner_adj, 1.0, 0)
  spinner$setValue(1)
  box1$add(spinner)
  spinner$grabFocus()
  box0$add(box1)

  gSignalConnect(spinner, "key-press-event", function(obj, evt){
    if(evt[["keyval"]] == GDK_Return) {
      spinner$update()
      do.insert() 
      dialog$destroy()  
    }
    if(evt[["keyval"]] == GDK_Escape) {
      dialog$destroy()
    }    
    .local$view$grabFocus()    
    FALSE
  })
  
  if(!is.null(win)) {
     checkPtrType(win, "GtkWindow")
     dialog$setPosition(GtkWindowPosition["center-on-parent"])
     dialog$setTransientFor(win)
  }
  gSignalConnect(dialog, "response", function(dlg, arg1, ...) {
     if(arg1==1) do.insert() 
     dialog$destroy()
     .local$view$grabFocus()     
  })                                    
  dialog[["vbox"]]$packStart(box0, TRUE, TRUE, 10)
  dialog$showAll()
}    

DeleteSelection <- function(.local){
  selections <- .local$selections
  if(length(selections)){
    xors <- get.xor.selection(.local$selections)
    # drop rows and cols that don't contain cells
    xbb <- as.integer(c( 
             rownames(xors)[c(1, nrow(xors))], 
             colnames(xors)[c(1, ncol(xors))] ))
    row.idx <- xbb[1]:xbb[2]
    col.idx <- xbb[3]:xbb[4]+COLUMN_OFFSET
    theSelection <- .local$theFrame[row.idx, 
      col.idx, drop=FALSE]
    theSelection[xors] <- NA
    nf <- theSelection
    task <- list(list(func="ChangeCells", 
      arg = list(nf=nf, row.idx=row.idx, col.idx=col.idx)))
    DoTaskWrapper(.local, task)
  }
}

GetCursorRowColIdx <- function(widget){
  cursor.info <- gtkTreeViewGetCursor(widget)
  path <- cursor.info$path
  if(is.null(path)) return(NULL)       
  col.idx <- GetColIdx(cursor.info$focus.column)  
  row.idx <- as.integer(gtkTreePathGetIndices(path))+1                
  return(c(row.idx, col.idx))
}

ViewKeyPress <- function(widget, event, data) {
  .local <- data
  retval <- TRUE       
  event.time <- event[["time"]]
  keyval <- event[["keyval"]]
  .local$last.view.event <- event
  .local$last.key.pressed <- keyval
  
      # remove click info
  .local$last_click_info <- NULL
  #cat(c(".", "*")[.local$allow.key.flag+1])
  if(0 && !is.null(.local$entry)) {
   #print("ViewKeyPress unrealizing entry")
  .local$entry$unrealize()
  .local$entry <- NULL
  }     

  if(.local$allow.key.flag && event.time > .local$last.time && is.null(.local$entry)){ #4-21-10
  #if( allow.key.flag && event.time > .local$last.time){  

#  if(.local$allow.key.flag) if(!gtkWidgetIsFocus(.local$view)) { #added 4-21-2010
#      print("ViewKeyPress: focus is FALSE")
#     gtkWidgetGrabFocus(.local$toplevel)      
#     gtkWidgetGrabFocus(.local$view)
#  }  

  
    #tryCatch({
      .local$last.time <- event.time
      .local$allow.key.flag <- FALSE                            
      stat <- event[["state"]] 
  
      allColumns <- .local$allColumns
      view <- .local$view
      model <- .local$model
      
      # Paging events
    if (CtrlLetter(keyval, stat, GDK_c) || CtrlLetter(keyval, stat, GDK_C)) {
      namesflag <- CtrlLetter(keyval, stat, GDK_C)

      gsr <- GetSelectedRows(.local$view, .local)
     # block transient interactions with a popup
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())   
      nf <- GetSelectedCells(.local)
      if(!is.null(nf))
        CopyToClipboard(nf, 
        do.rownames=namesflag && !RowNamesAreNull(.local$theFrame), do.colnames=namesflag && !ColNamesAreNull(.local$theFrame))
      retval <- TRUE
      .local$allow.key.flag <- TRUE                          
      #return(TRUE)
    } else if (CtrlLetter(keyval, stat, GDK_v) && .local$editable){
      # print("Ctrl-v")
      ReadFromClipboardWrapper(.local)
      retval <- TRUE      
      .local$allow.key.flag <- TRUE                          
      #return(TRUE)
    } else if (CtrlLetter(keyval, stat, GDK_z)){
      #print("Ctrl-z")
      DoUndo(.local)
      .local$allow.key.flag <- TRUE
      retval <- TRUE
      #return(TRUE)
    } else if (CtrlLetter(keyval, stat, GDK_a)){
      #print("Ctrl-a")
      SelectAll(.local)
      .local$allow.key.flag <- TRUE        
      retval <- TRUE
      #return(TRUE)
    } else if(keyval == GDK_Delete  && .local$editable) {
       # replace set of selection rectangles with disjoint
       # rectangles in order to minimize paint/display properly
      DeleteSelection(.local)
      .local$allow.key.flag <- TRUE                  
      retval <- TRUE      
      #return(TRUE)   
    } else if(keyval%in%myMetaKeys ||  (stat == GdkModifierType['control-mask'])){
       # start shift-selection if user pressed shift
      if(keyval%in%myShiftKeys){
        rc <- GetCursorRowColIdx(widget)
        if(!is.null(rc)){
        .local$selections <- list(list(start=c(row.idx=rc[1], col.idx=rc[2])))
#print("Shift pressed")
#print(unlist(.local$selections))
        }
      }
      .local$allow.key.flag <- TRUE                
      retval <- TRUE        
      #return(TRUE)
    } else if (keyval%in%myValidNavigationKeys){
      #.local$rectangles <- list()      
      if (keyval == GDK_Up || ShiftLetter(keyval, stat, GDK_Return))
        moveDirection <- "up"
      else if (keyval %in% c(GDK_Down, GDK_Return))
        moveDirection <- "down"
      else if(keyval %in% c(GDK_Right, GDK_Tab))
        moveDirection <- "right"
      else if(keyval %in% c(GDK_Left, GDK_ISO_Left_Tab))
        moveDirection <- "left"
      else if(keyval %in% c(GDK_Page_Down))
        moveDirection <- "page_down"
      else if(keyval %in% c(GDK_Page_Up))
        moveDirection <- "page_up"
      else if(keyval %in% c(GDK_Home))
        moveDirection <- "home"
      else if(keyval %in% c(GDK_End))
        moveDirection <- "end"
      else
        moveDirection <- NULL

      if(!is.null(moveDirection)){
        MoveCursor(widget, moveDirection, .local, stat)
        retval <- TRUE        
      } else {
        retval <- FALSE
      }

      .local$allow.key.flag <- TRUE                  
    } else if (keyval==GDK_Insert){
        if(.local$editable){
		    #print("Starting insert")
		    cursor.info <- gtkTreeViewGetCursor(widget)
		    path <- cursor.info$path          
		    col.idx <- GetColIdx(cursor.info$focus.column)
		    row.idx <- as.integer(gtkTreePathGetIndices(path))                
		    if(!is.null(path)){
		      insert.dialog(.local, row.idx, col.idx, insert.type="Rows", win=.local$toplevel)
		    }
        }
        .local$allow.key.flag <- TRUE           
        retval <- TRUE                   
    } else if (keyval == GDK_equal && length(sc <- GetSelectedColumns(.local))) {
      if(length(sr <- GetSelectedRows(.local$view, .local))==0) 
          sr <- 1:(dim(.local$theFrame)[1]-EXTRA_ROW)
      CommandData(.local, sr, sc)
      .local$allow.key.flag <- TRUE        
      retvale <- TRUE
    } else if (keyval%in%myValidInputKeys){ # valid input keys
      cursor.info <- gtkTreeViewGetCursor(widget)
      col.idx <- GetColIdx(cursor.info$focus.column)
      renderer <- allColumns[[col.idx]]$renderer        
      path <- cursor.info$path 
      row.idx <- as.integer(gtkTreePathGetIndices(path))+1
      ddf <- dim(.local$theFrame)
        # ignore last row
      if(row.idx > ddf[1]-EXTRA_ROW || col.idx > ddf[2]-2) {
        .local$allow.key.flag <- TRUE        
        retval <- TRUE          
      } else {             
        if(.local$editable){
           gtkTreeViewSetCursorOnCell(widget, cursor.info$path,
             cursor.info$focus.column, renderer, TRUE)
           if(keyval != GDK_space){
             gtkPropagateEvent(view, event)
           }
         } else {
           .local$allow.key.flag <- TRUE
         }
         retval <- FALSE
           # Bug spotted by Liviu 28/9/12
         if(keyval == GDK_space){
           retval <- TRUE
         }
      } # if row.idx
      #.local$allow.key.flag <- TRUE           
    }
    #}, error = function(e){
    #    warning(e)
    #    #cat("E")
    #    .local$allow.key.flag <- TRUE
    #   }) # TryCatch       
  } else {# if allow.key.flag is false or event.time
#print("ViewKeyPress Blocked")
  }
  if(!keyval%in%myValidKeys) { # otherwise bad things
    .local$allow.key.flag <- TRUE
  }
  
  return(retval)
}

GetSelectedColumns <- function(.local, restrict=TRUE) {
 rv <- .local$selectedColumns
 if(restrict) rv <- rv[!rv%in%(dim(.local$theFrame)[2]-1)]
 return(rv)
}

RendererEditingStarted <- function(renderer, entry, path, data) {  
  .local <- data$.local
  #print(paste("RendererEditingStarted called, local key flag is", .local$allow.key.flag))

  if(!is.null(.local$entry)) {
    #cat ("Already using entry!!!\n")
  }
  .local$entry <- entry
  checkPtrType(entry, "GtkEntry")
  col.idx <- data$col.idx
  #print(entry)

  theFrame <- .local$theFrame
  view <- .local$view
    # don't paint rectangles
  .local$do.paint <- FALSE
  myBlack = as.GdkColor("black")
  myWhite = as.GdkColor("white")
  gtkEntrySetHasFrame(entry, TRUE)
  gtkWidgetModifyBase(entry, as.integer(1), myBlack)  
  gtkWidgetModifyBase(entry, as.integer(3), myBlack)      
  gtkWidgetModifyText(entry, as.integer(1), myWhite)              
  gtkWidgetModifyText(entry, as.integer(3), myWhite)   
            
  #entry$setAlignment(1)
  
  gtkWidgetSetEvents(entry, GdkEventMask["all-events-mask"])	  	  
  gObjectSetData(entry, "index", col.idx)
  gObjectSetData(entry, "renderer", renderer)
  gObjectSetData(entry, "path.str", path)
  gObjectSetData(entry, "path", gtkTreePathNewFromString(path))
  gObjectSetData(entry, "text", gtkEntryGetText(entry))    

   # if the used pressed space, we want to insert it the normal way
   # this is surprisingly hard to do with a treeview...
   # However we could probably avoid the propagate call
   # this way!
  if(identical(.local$last.key.pressed, GDK_space))
   gSignalConnect(entry, "map-event", 
     function(entry, event, data=.local) {
       entry$setText(" ")
       #entry$setText(.local$last.key.pressed)
       gtkEditableSelectRegion(entry, 1, 1)
       })

  # Insert cursor at point where we clicked
  gSignalConnect(entry, "map-event", function(entry, event, data=.local){
#print("RendererEditingStarted: map-event insert") 
   # This denies focus to the entry
if(0){
    insertion_point <- 0
    pixel_size <- entry$getLayout()$getPixelSize()
    click_info <- .local$last_click_info
    col_width <- click_info$column$getWidth()
    text_width <- pixel_size$width
    text_pos <- (click_info$cell.x - col_width + text_width)/(text_width)
    if(length(text_pos) != 1 || text_pos < 0 || text_pos > 1) text_pos <- 0
    insertion_point <- round(text_pos*nchar(entry$getText()))   
    gtkEditableSelectRegion(entry, insertion_point, insertion_point)
}

#print("RendererEditingStarted: map-event insert returning")
    return(FALSE)  
 })
  # 7-4-10 - disable right-click menu, seems to confuse people
  IgnoreRClick <- function(entry, event, data){
    if(event[["button"]] == as.integer(3)) return(TRUE)
    return(FALSE)
  }

  if(1) gSignalConnect(entry, "focus-in-event", function(entry, event) {
    #print("RendererEditingStarted: entry focused in")
    return(FALSE)
  })

  gSignalConnect(entry, "button-press-event", IgnoreRClick)
  gSignalConnect(entry, "button-release-event", IgnoreRClick)  
    # We've trapped the temporary GtkEntry the renderer creates.
  gSignalConnect(entry, "key-press-event", RendererEntryKeyPress, data=.local)
    # you can leave the entry in 2 ways, focusing out or pressing a key from within
  .local$entry.focus.out <- gSignalConnect(entry, "unrealize", function(entry, ...){
    #print("RendererEditingStarted: entry focused out")
    #if(!identical(.local$ignore_entry_focus_out, NULL)) {
    #  .local$ignore_entry_focus_out <- NULL
    #  return(FALSE)
    #}
    #.local$ignore_entry_focus_out <- TRUE
    #gtkPropagateEvent(entry, event)           
    EnterTextIntoCell(entry, .local)    
    #print(.local$entry)
    .local$allow.key.flag <- TRUE # unlock
    .local$entry <- NULL                 
      # changed 082210
    .local$entry.focus.out <- NULL
    #print("RendererEditingStarted: unlocked from focus-out")

    return(FALSE)
  }) # focus out event
  
 .local$doingEntryIntoCompletion <- FALSE  
 typ <- GetClasses(.local$theFrame)[col.idx]
 if(typ == "factor"){
   create.completion.model <- function(factor.levels) {
     store <- gtkListStoreNew("gchararray")
     sapply(factor.levels, function(x) store$set(store$append()$iter, 0, x))
     return(store)
   }
   factor.levels <- levels(theFrame[,col.idx])
   completion <- gtkEntryCompletionNew()
     # Assign the completion to the entry
 
     # Create a tree model and use it as the completion model
   completion.model <- create.completion.model(factor.levels)
   entry$setData("levels", factor.levels)    
   completion$setModel(completion.model)
   completion$setTextColumn(0)
   completion$setMinimumKeyLength(0)
   completion$setInlineCompletion(TRUE)
   completion$setPopupSingleMatch(FALSE)
   .local$doingEntryIntoCompletion <- TRUE
     # move cursor down if you hit return after you select a match
        
   gSignalConnect(completion, "match-selected", after=F, function(widget, completion.model, iter){
     #entry$setText(gtkTreeModelGetValue(completion.model, iter, 0)$value)
#print("RendererEditingStarted: Focused out via match-selection")
     EnterTextIntoCell(entry, .local)
     MoveCursor(.local$view, "down", .local)
     return(FALSE)
   })
  entry$setCompletion(completion)
  entry$setData("completion", completion)

  entry$setData("pos", 0)  # match position

  entry$setData("current.match", character(0))  
  }
  #print("RendererEditingStarted returning")

  return(FALSE)
}

EnterTextIntoCell <- function(entry, .local){
  #w <- TransientWindow("Updating...", .local)  
  #on.exit(w$destroy())
#print("EnterTextIntoCell entered")
  col.idx <- gObjectGetData(entry,"index")  
  row.idx <- as.integer(gObjectGetData(entry, "path.str"))+1      
  keyval <- gObjectGetData(entry,"keyval")
  #} else {
    theText <- gtkEntryGetText(entry)          
  #}
#print(c("theText", theText))
  nf <- data.frame(theText, stringsAsFactors=F, check.names=F)
  task <- list(list(func="ChangeCells", 
    arg = list(nf=nf, row.idx=row.idx, col.idx=col.idx)))
  DoTaskWrapper(.local, task)
#print("EnterTextIntoCell returning")
}

RendererEntryKeyPress <- function(entry, event, data) {
  .local <- data
  retval <- TRUE    
  event.time <- event[["time"]]
  keyval <- event[["keyval"]]
  #cat(c("A", "B")[.local$allow.renderer.key.flag+1])
  if(1){# .local$allow.renderer.key.flag && event.time >= .local$last.time){
    .local$last.time <- event.time
    retval <- FALSE
    keyval <- event[["keyval"]]
    stat <- as.flag(event[["state"]])
    gObjectSetData(entry, "keyval", keyval)
  	gObjectSetData(entry,"stat", stat)
  	
    view <- .local$view
      # control keys for popup selection
    if(.local$doingEntryIntoCompletion) { # we're in the popup. or something.
      retval <- FALSE
              
      if (keyval %in% c(GDK_Return)){        
          # this is to handle returns on the entry completion
        current.match <- entry$getData("current.match") 
        if(length(current.match)) 
          entry$setText(current.match)
          
        if(nchar(entry$getText()) > 0) {  
          #EnterTextIntoCell(entry, .local)  
          #.local$allow.key.flag <- TRUE # unlock
          #.local$entry <- NULL
#print("RendererEntryKeyPress, returning on completion")                
          MoveCursor(view, "down", .local)
        } 
        
       #.local$allow.key.flag <- TRUE # unlock
       #.local$entry <- NULL              
        #
        retval <- FALSE
      } else if(0 && keyval == GDK_Escape){
        #entry$setText(entry$getData("text"))
        #MoveCursor(view, "down", .local)
        retval <- TRUE
      }
        # Keep track of entry completion
        # this is really hacky. first, get the matches in the entry
      ll <-  entry$getData("levels")
        # this may not be perfect...
      matches <- ll[regexpr(tolower(entry$getText()),tolower(ll)) == 1]
      if( length(matches) > 1) {  # because the popup is only for these matches
         # now keep track of the up/down arrow selection...
        if (keyval %in% c(GDK_Up, GDK_Down)){
          xx <- -1 
          if(keyval%in%c(GDK_Down)) 
            xx <- 1          
          xx <- (entry$getData("pos") + xx)%%(length(matches)+1)
        } else {
          xx <- 0        
        }
        #print(xx)
        #print(matches[xx])
        entry$setData("pos", xx)
        entry$setData("current.match", matches[xx])
      }
      # just a regular text cell
      # calling MoveCursor ends the focus of the entry
    }  else if (!keyval%in%myValidInputKeys && !keyval%in%myMetaKeys){  
      # EnterTextIntoCell(entry, .local)  # update the cell
      # .local$allow.key.flag <- TRUE # unlock
      # .local$entry <- NULL
 
      if (keyval == GDK_Escape){
        gtkEntrySetText(entry, gObjectGetData(entry, "text"))   
      } else if (keyval == GDK_Up || ShiftLetter(keyval, stat, GDK_Return)){
        MoveCursor(view, "up", .local)
      } else if ( keyval == GDK_Return || keyval == GDK_Down ) {
        MoveCursor(view, "down", .local)
      }  else if (keyval == GDK_Left || keyval == GDK_ISO_Left_Tab) {
        MoveCursor(view, "left", .local)
      } else if (keyval == GDK_Right || keyval == GDK_Tab) {
        MoveCursor(view, "right", .local)
      } 
      retval <- FALSE
    } else { # not entry into completion, or invalid keys
    }

    if(keyval%in%myMetaKeys) # ignore these
      retval <- TRUE
  } 

  # new code here  
  #print("RendererEntryKeyPress setting flag to TRUE")
  #.local$allow.key.flag <- TRUE
  #print(.local$view$isFocus()) # Here's the problem....
  return(retval)
}

ColumnSetState <- function(allColumns, ii, state, .local){    
  if(length(ii) != 1 || ii < 1 || length(allColumns) < ii)
	  stop(paste("trying to set column index outside bounds:", ii))	
  #if(ii != length(allColumns))
  gtkWidgetSetState(allColumns[[ii]]$button,
    ifelse(state, as.integer(1), as.integer(0)))
  .local$allColumns <- allColumns
  if(state) {
    .local$selectedColumns <- sort(union(ii, .local$selectedColumns)) 
  } else {
    .local$selectedColumns <- sort(setdiff(.local$selectedColumns, ii)) 
  }
}  

UpdateColumnSelection <- function(.local, selectedColumns.new) {
    allColumns <- .local$allColumns
    selectedColumns <- GetSelectedColumns(.local, restrict=FALSE)
    ll <- length(allColumns)
    if(ll > 1){ # for zero columns, don't select anything
        # Don't activate the last column - it's left blank
      for (ii in setdiff(selectedColumns.new, selectedColumns))
        ColumnSetState(allColumns, ii, TRUE, .local)
      for (ii in setdiff(selectedColumns, selectedColumns.new))
        ColumnSetState(allColumns, ii, FALSE, .local)
    }
}

SelectAll <- function(.local){
  allColumns <- .local$allColumns
  UpdateColumnSelection(.local, 1:length(allColumns))
  if(length(.local$rectangles)){
    .local$rectangles <- list()
    .local$viewGetBinWindow$invalidateRect(NULL, FALSE)
  }
  .local$do.paint <- TRUE
  
   # block transient interactions with a popup
  w <- TransientWindow("Selecting...", .local)
  on.exit(w$destroy())
  
  #gtkTreeSelectionSelectAll(gtkTreeViewGetSelection(.local$view))    
  .local$selections <- list(list(
     start=c(row.idx=1, col.idx=1),
     end=c(row.idx=nrow(.local$theFrame), col.idx=ncol(.local$theFrame)-2)
     )) 
  UpdateSelectionRectangle(.local)      
}

RowNamesButtonPress <- function(widget, event, data) {
  .local <- data
   # kill any open entries   
  typ <- event[['type']]
  stat <- event[['state']]
  info <- widget$getPathAtPos(event[["x"]], event[["y"]])
  if(is.null(info$path)) return(TRUE)         

  row.idx <- info$path$getIndices()+1
  col.idx <- 1
        
  allColumns <- .local$allColumns

  #UpdateColumnSelection(.local, 1:(length(allColumns)))
  .local$do.paint <- TRUE
    
  #UpdateColumnSelection(.local, integer(0))
  #if(length(.local$rectangles)){
  #  .local$rectangles <- list()
  #  .local$viewGetBinWindow$invalidateRect(NULL, FALSE)             
  #}

  if (.local$scrollRowNames) # i.e., we're not propagating event
  if (event[["button"]] == as.integer(1)){
    if(row.idx > dim(.local$theFrame)[1]-EXTRA_ROW) return(TRUE)
    if(typ == as.integer(4)){ # single clicked     
    selections <- .local$selections

  #if(0){
   theFrameNCols <- ncol(.local$theFrame)-1-COLUMN_OFFSET
   rowSel <- list(start=c(row.idx=row.idx, col.idx=1),
        end=c(row.idx=row.idx, col.idx=theFrameNCols))

   # selectedRows <- GetSelectedRows(.local)  	
   if (as.flag(stat) & GdkModifierType['shift-mask']) { # range
    if(is.null(selections) || !length(selections))
      selections <- list(list(start=rowSel$start))
    # For the tricky case where we've selected a cell range
    # and shift-click - unsure of best behavior here
    selections[[length(selections)]]$start['col.idx'] <- 1
    selections[[length(selections)]]$end <- rowSel$end

  } else if (as.flag(stat) & GdkModifierType['control-mask']) { # range
    # Simpler: Just add the column to whatever you clicked.
    selections <- .local$selections

    if(length(selections)){
      xors <- get.xor.selection(selections)
      airx <- as.integer(rownames(xors))
      selectedRows <- airx[rowSums(xors)>0]
        # if the column index is not in the range 
        # of selected columns,
        # we want to add it
      whichidx <- which(airx==row.idx)
      # if the entire column is already selected and TRUE
      # then set it to FALSE 
      if(row.idx%in%selectedRows){
        #dont.add.flag <- TRUE
        if(all(xors[whichidx,])
          && ncol(xors)==theFrameNCols){
          xors[whichidx,] <- FALSE
          selections <- get.xor.rectangles(xors)
        } else {
          xors[whichidx,] <- FALSE
          selections <- get.xor.rectangles(xors)
          selections[[length(selections)+1]] <- rowSel
        }
      # otherwise set the whole column to TRUE
      } else {
         selections[[length(selections)+1]] <- rowSel
      }
    # if length(selections)
    } else {
       selections[[length(selections)+1]] <- rowSel
    }

  } else {
    gtkTreeSelectionUnselectAll(.local$ss)    
    selections <- list()
    selections[[length(selections)+1]] <- rowSel
 }
 .local$selections <- selections
  #} # if 0
  .local$do.paint <- TRUE

    UpdateSelectionRectangle(.local)
    gdkWindowInvalidateRect(.local$viewGetBinWindow, NULL, FALSE)
    #gdkWindowInvalidateRect(.local$view.rnGetBinWindow, NULL, FALSE)


    # prevent a repeated click from editing
	  cursor.info <- gtkTreeViewGetCursor(widget)
      if(is.null(cursor.info$path))
        return(FALSE)

	  cursor.row.idx <- as.integer(gtkTreePathGetIndices(cursor.info$path))+1
      if(identical(cursor.row.idx, row.idx))
        return(TRUE)

    DoUserCall("RowClicked", list(idx = info$path$getIndices()+1), .local)

     # scroll main view in sync
    if( .local$scrollRowNames) { 
      .local$scrollRowNames <- FALSE
      gtkPropagateEvent(widget, event)
      gtkTreeSelectionUnselectAll(.local$ss)
      gtkAdjustmentSetValue(.local$sw.view.va, gtkAdjustmentGetValue(.local$sw.rn.va))
      .local$scrollRowNames <- TRUE 
      return(TRUE)
   }
     
   } else if(typ == as.integer(5)){ # ignore double click
     return(TRUE)
   }      
  } else if (event[["button"]] == as.integer(3)){ # our popup menu
    m <- Row3rdButtonMenu(.local, row.idx)
    gtkMenuPopupHack(m, button = event$GetButton(), 
      activate.time = gdkEventGetTime(event))
    return(TRUE)
  } 
  return(FALSE)
}


# given the selection, creates a bounding box mask 041210
get.xor.selection <- function(selections){
	mtx <- apply(as.data.frame(selections), 1, range)
  stopifnot(nrow(mtx)==2 && ncol(mtx)> 0 && ncol(mtx)%%2==0)
	ary <- array(FALSE, apply(mtx, 2, diff)+1)
	msq <- function(xc) {xc <- unlist(xc); seq.int(from=xc[1], to=xc[2])}
	dimnames(ary) <- list(msq(mtx[,1]), msq(mtx[,2]))
	for(ii in seq(length=length(selections))){
		x1 = as.data.frame(selections[[ii]])
		x1.idx <- x1 - mtx[1,] + 1
    if(!(nrow(x1.idx)==2 && ncol(x1.idx)==2)) break;
		sii <- msq(x1.idx[1,])
		sjj <- msq(x1.idx[2,])
		ary[sii, sjj] <- !ary[sii, sjj]
	}
	return(ary)
}

get.xor.rectangles <- function(xors){
	rects = list()

    row.offset <- as.integer(rownames(xors)[1])-1
    names(row.offset)<-NULL
    col.offset <- as.integer(colnames(xors)[1])-1
    names(col.offset)<-NULL
	while(any(xors)){
	wxors <- which(xors, arr.ind=T)
	  # get first element which isn't in first row containing true
	top.row = wxors[1,1]
	left.col = wxors[1,2]

	first.col.containing.true <-  xors[top.row:nrow(xors),left.col]

	row.height <- which(diff(first.col.containing.true)!=0)[1]
    if(is.na(row.height)) row.height <- length(first.col.containing.true)
    bottom.row <- top.row+row.height-1

	xx <- colSums(!xors[top.row:bottom.row,,drop=F])

	right.col <- which(xx[-(1:left.col)]!=0)[1] + left.col-1
	if(is.na(right.col)) right.col <- length(xx)
	stopifnot(all(xors[top.row:bottom.row, left.col:right.col]))
	xors[top.row:bottom.row, left.col:right.col] <- FALSE

    names(top.row) <- names(left.col) <- names(bottom.row) <- names(right.col) <- NULL

	rects[[length(rects)+1]] <- list(
		 start = c(row.idx=top.row+row.offset, 
                  col.idx=left.col+col.offset),
		 end = c(row.idx=bottom.row+row.offset, 
                col.idx=right.col+col.offset)
	  )
#print("A")
#print(rects[[length(rects)]])
#print("B")
	}

return(rects)
}

# Takes selections as xor-ing one another
# Return selection list, i.e. list of start, end row.idx col.idx
GetSelectedCells <- function(.local){
   # replace set of selection rectangles with disjoint
   # rectangles in order to minimize paint/display properly
  selections <- .local$selections
  if(length(selections)){
    xors <- get.xor.selection(.local$selections)
    # drop rows and cols that don't contain cells
    xbb <- as.integer(c( 
             rownames(xors)[c(1, nrow(xors))], 
             colnames(xors)[c(1, ncol(xors))] ))
    wxx <- xbb[3:4] == ncol(.local$theFrame)-1
    theSelection <- NULL
    if(!all(wxx)){
      xbb[3:4][wxx] <- ncol(.local$theFrame)-2
      theSelection <- .local$theFrame[xbb[1]:xbb[2], 
        xbb[3]:xbb[4]+COLUMN_OFFSET, drop=FALSE]
        # remove entirely unselected rows and columns
      if(any(rowSums(xors)==0)){
        theSelection <- theSelection[rowSums(xors)!=0,,drop=F]
        xors <- xors[rowSums(xors)!=0,,drop=F]
      }
      if(any(colSums(xors)==0)){
        theSelection <- theSelection[,colSums(xors)!=0,drop=F]
        xors <- xors[,colSums(xors)!=0,drop=F]
      }
      if(!all(xors)){
        theSelection[!xors] <- NA
      }
    }
  }
  return(theSelection)
}

ViewButtonPress <- function(widget, event, data) {

  .local <- data
  .local$last.view.event <- event
  model <- .local$model

  allColumns <- .local$allColumns

  info <- widget$getPathAtPos(event[["x"]], event[["y"]])

  if (is.null(info$column)) return(FALSE)

  renderer <- info$column$getCellRenderers()[[1]]

    # store last click info here
  .local$last_click_info <- info
  
  if(is.null(info$path)) return(TRUE)
  row.idx <- info$path$getIndices()+1
  col.idx <- GetColIdx(info$column)
  
  typ <- event[['type']]
  stat <- event[['state']]

  if (event[["button"]] == as.integer(3)){ # our popup menu
  	m <- Cell3rdButtonMenu(.local, row.idx, col.idx)
    gtkMenuPopupHack(m, button = event$GetButton(),
      activate.time = gdkEventGetTime(event))
    return(TRUE)
  } else if (event[["button"]] == as.integer(1)){    
    if(typ == as.integer(4)){ # single clicked
        # we want to keep going if user clicks on column to resize - hard to tell           
      if(is.null(info$path) || (identical(row.idx, 1) && info$cell.x < 5)) {   
        .local$rectangles <- list()
        return(FALSE)
      }

        # if it isn't editable, 
        # make sure we don't drop into edit mode
	  if(!.local$editable) {
		  cursor.info <- gtkTreeViewGetCursor(widget)
		  path <- cursor.info$path
		  if(is.null(path)){
			cursor.row.idx = 1
			cursor.col.idx = 2
		  } else {
			cursor.row.idx <- as.integer(gtkTreePathGetIndices(path))+1
			cursor.col.idx <- GetColIdx(cursor.info$focus.column)  
		  }
	      if(identical(cursor.row.idx, row.idx) && identical(col.idx, cursor.col.idx)){
  	        return(TRUE)
          }
	  }
        
      selectedColumns <- GetSelectedColumns(.local)
      if (as.flag(stat) & GdkModifierType['shift-mask']) { # range
        if(is.null(.local$start.column.select)) .local$start.column.select <- col.idx
        selectedColumns <- .local$start.column.select:col.idx
    
       if(is.null(.local$selections)) return(TRUE)
       stopifnot(length(.local$selections[[length(.local$selections)]])>0)
       .local$selections[[length(.local$selections)]]$end <- 
          c(row.idx=row.idx, col.idx=col.idx)

          # If SHIFT is held down and no start is selected, do nothing
        if(is.null(.local$start.column.select)) {
          return(TRUE)
#          .local$selections <- list()
#          .local$selections[[length(.local$selections)+1]] <- 
#            list(start=c(row=row.idx, col=col.idx))
       }


	   } else if (as.flag(stat) & GdkModifierType['control-mask']) { # range
         #cursor.row.idx <- as.integer(gtkTreePathGetIndices(path))+1
		 selectedColumns <- union(selectedColumns, col.idx)
		 #selectedColumns <-
		 # c(union, setdiff)[[(col.idx%in%selectedColumns+COLUMN_OFFSET)]](selectedColumns, col.idx)
         if(is.null(.local$start.column.select)) 
           .local$start.column.select <- col.idx
          # Treat ctrl-button down as starting a new selection range
        .local$selections[[length(.local$selections)+1]] <- list(start=
           c(row.idx=row.idx, col.idx=col.idx))

        return(TRUE)
      } else { # no control or shift
         selectedColumns <- col.idx
        .local$start.column.select <- col.idx

        .local$selections <- list()
        .local$selections[[1]] <- list()
        .local$selections[[1]]$start <- c(row.idx=row.idx, col.idx=col.idx)

        #gsr <- GetSelectedRows(.local$view, .local)      
       gsc <- GetSelectedColumns(.local)

         # This is to prevent clicking twice on the same row from opening
         # an edit entry in a different cell
#        if(length(gsr) && length(gsc) && gsr[1] == row.idx && (length(gsr) > 1 || length(gsc) > 1 || col.idx != gsc) ) {

        # if you've clicked twice on a selected cell, 
        # ignore it
        # Test whether the first selected row is equal to the row index

      mysr <- gtkTreeSelectionGetSelectedRows(.local$ss)$retval 
      if(length(mysr)){
        first.row.idx <- gtkTreePathGetIndices(mysr[[1]])+1
        if(identical(first.row.idx, row.idx) && 
         !(col.idx%in%gsc) ) {

          if(!is.null(.local$entry)) {
            .local$entry$unrealize()
            .local$entry <- NULL
          }     
          .local$do.paint <- TRUE          
          .local$rectangles <- list()          
          .local$view$setCursorOnCell(info$path, info$column, renderer, FALSE)                                            
          UpdateColumnSelection(.local, selectedColumns)                     
          UpdateSelectionRectangle(.local)           
          return(TRUE)               
        }
       }
      }
      UpdateColumnSelection(.local, selectedColumns) 
    }  else if (typ == as.integer(5)) { # ignore dclick
      return(TRUE)
    }
   }
   
   ddf <- dim(.local$theFrame)
    # reset the cursor so it's visible  
  #.local$flash.cursor <- TRUE
    #  Flash the rectangles off when you click
  .local$do.paint <- TRUE        # 2-12-10
  if(length(.local$rectangles)){
#  print("resetting") 
    #.local$rectangles <- list() # 3-2-10
	widget$getBinWindow()$invalidateRect(NULL, FALSE) # 041012
  }
  return(FALSE)
}

# Update rectangle from .local$view selected rows.
# Don't draw if only 1 selected square
# Don't update if nothing changed.
UpdateSelectionRectangle <- function(.local, selection=NULL, widget=.local$view){
  allColumns <- .local$allColumns
  sw.va <- .local$sw.view.va
  ncolmax <- ncol(.local$theFrame)
  
  path2.index <- path1.index <- integer(0)
  
  row.idx <- col.idx <- NULL

  selections <- .local$selections

  t.cols <- integer(0)
  t.rows <- integer(0)

  if(length(selections) > 1){
    xors <- get.xor.selection(.local$selections)
    selections <- get.xor.rectangles(xors)
    t.rows <- as.integer(rownames(xors)[rowSums(xors)>0])
    t.cols <- as.integer(colnames(xors)[colSums(xors)>0])
  } else if (length(selections) == 1) {
    ls1 <- .local$selections[[1]]
    t.seq <- sort(c(ls1$start['row.idx'], ls1$end['row.idx']))
    t.from <- t.seq[1]
    t.to <- t.seq[2]
    if(!is.na(t.from) && !is.na(t.to)){
      t.rows <- seq(from=t.from, to=t.to)
      t.cols <- seq(from=.local$selections[[1]]$start['col.idx'],
        to =.local$selections[[1]]$end['col.idx'])
    }
  }

  UpdateColumnSelection(.local, t.cols)

    # row name rectangles
  if(0 && length(t.rows)){
    view.rn <- .local$view.rn
    va.value <- gtkAdjustmentGetValue(.local$sw.view.va)
    t.blocks <- MergeAdjacent(t.rows)-1
    rn.rectangles <- vector("list", nrow(t.blocks))
	  rn.col <- .local$rowname.column.obj$column
    for(ii in seq(length=nrow(t.blocks))){
	    path.start <-  gtkTreePathNewFromString(t.blocks[ii,1])
	    path.end <-  
        gtkTreePathNewFromString(t.blocks[ii,2])
  	  rn.rect.start <- 
        gtkTreeViewGetCellArea(view.rn, path.start, rn.col)$rect
  	  rn.rect.end <- 
        gtkTreeViewGetCellArea(view.rn, path.end, rn.col)$rect
	    rn.rect <- gdkRectangleUnion(rn.rect.start, rn.rect.end)$dest
	    rn.rect$y <- rn.rect$y + va.value # adjust for scroll
	    rn.rectangles[[ii]] <- rn.rect
    }
    .local$rn.rectangles <- rn.rectangles
  }

  if(length(selections)){
    va.value <- gtkAdjustmentGetValue(.local$sw.view.va)
    view <- .local$view 
    rectangles <- list()
	  for(ii in seq.int(length.out=length(selections))){
      sel.rect <- selections[[ii]]
      if(length(sel.rect$start) < 2 || length(sel.rect$end) < 2) break;
    row.block <- sort(c(sel.rect$start['row.idx'], sel.rect$end['row.idx']))
	  row.block.start <- row.block[1]
	  path.start <-  gtkTreePathNewFromString(row.block.start-1)

	  row.block.end <- row.block[2] 
    if(row.block.start == row.block.end)
      path.end <- path.start
    else
	    path.end <-  gtkTreePathNewFromString(row.block.end-1)

    col.block <- sort(c(sel.rect$start['col.idx'], sel.rect$end['col.idx']))
	  col.block.start <- col.block[1]
	  col.start <- allColumns[[col.block.start]]$column
	  rect <- gtkTreeViewGetCellArea(view, path.start, col.start)$rect

  	col.block.end <- col.block[2] 
    if(col.block.start != col.block.end || row.block.end != row.block.start){
  	  col.end <- allColumns[[col.block.end]]$column
  	  rect.end <- gtkTreeViewGetCellArea(view, path.end, col.end)$rect
  	  rect <- gdkRectangleUnion(rect, rect.end)$dest
    }
	  rect$y <- rect$y + va.value # adjust for scroll
	  rectangles[[length(rectangles) + 1]] <- rect
	}
	.local$rectangles <- rectangles        

  } else { # if length(...)
	  .local$rectangles <- list()    
  }

  DoUserCall("Selection", list(selections=selections), .local)
      
  gdkWindowInvalidateRect(.local$viewGetBinWindow, NULL, FALSE)      
  return(FALSE)   
}  

ViewButtonRelease <- function(widget, event, data){
  #if(.local$do.rendererEditingStarted){ # if we've started editing, don't do anything
  .local <- data
  .local$last.view.event <- event
  typ <- event[['type']]
  stat <- event[['state']]
  row.idx <- col.idx <- integer(0)
    # ignore this if there's an entry in progress
  if(is.null(.local$entry)){ # 041012 removed (if (1 || is.null ...
    allColumns <- .local$allColumns
    sw.va <- .local$sw.view.va
    info <- widget$getPathAtPos(event$x, event$y)
    selectedColumns <- GetSelectedColumns(.local)

      # This block deals with dragging outside the view rectangle
    if(is.null(info$column) || is.null(info$path)){
      view <- .local$view
      ptr <- view$getBinWindow()$getPointer()
      vr <- view$getVisibleRect()$visible.rect

      if(is.null(info$column))
        info$column <- widget$getPathAtPos(event$x, 0)$column
		  if(is.null(info$column)){
		    sw.ha <- .local$sw.ha
		    sw.ha.value <- sw.ha$value
		    direction <- ifelse (ptr$x - sw.ha.value <= vr$width/2, -1, 1)
		    if(direction == 1) col.idx <- length(allColumns)
		    if(direction == -1) col.idx <- 1
		  } else {
		    col.idx <- GetColIdx(info$column)
		  }

    if(is.null(info$path)){
      if(ptr$y < 0) ptr$y <- 0
      info$path <- widget$getPathAtPos(0, ptr$y)$path
    }

    if(is.null(info$path)) { # 
      sw.va <- .local$sw.view.va
      sw.va.value <- sw.va$value
      direction <- ifelse (sw.va$value+ptr$y <= vr$height/2, -1, 1)
      if(direction == 1) row.idx <- nrow(.local$theFrame)
      if(direction == -1) row.idx <- 1
    } else {
     row.idx <- info$path$getIndices()+1
    }
  } else {
     row.idx <- info$path$getIndices()+1
		 col.idx <- GetColIdx(info$column)
  }
   # restrict selection from last column
   # New code
    

   if (event[["button"]] == as.integer(1) && !is.null(.local$selections)){
     if(typ == GdkEventType['button-release']){ # single clicked
        if(as.flag(stat) & GdkModifierType['control-mask']){
          stopifnot(!is.null(.local$selections) 
            && length(.local$selections))
         .local$selections[[length(.local$selections)]]$end <- 
            c(row.idx=row.idx, col.idx=col.idx)
        } else if(as.flag(stat) & GdkModifierType['shift-mask']){
        } else {
          stopifnot(!is.null(.local$selections))
         .local$selections[[1]]$end <- c(row.idx=row.idx, col.idx=col.idx)
       }
     }
   }
   
    if(as.flag(stat) & GdkModifierType['control-mask']){
      #selectedColumns <- union(selectedColumns, col.idx)
    } else if (!is.null(.local$start.column.select) && length(col.idx)) {
      selectedColumns <- .local$start.column.select:col.idx
    } else if (length(col.idx)) {
      selectedColumns <- col.idx
    }# otherwise col.idx is NA


     # This may or may not be needed. 4-10-10
   .local$rectangles <- list() # reset our rectangle-drawing
  # 3-2-10
   UpdateSelectionRectangle(.local)
# split into rectangles


   # if we're still using the timeout to paint rectangles,
   # try to remove it
  if(!is.null(.local$do.paint.timeout)) {
    gSourceRemove(.local$do.paint.timeout) 
    .local$do.paint.timeout <- NULL
  }
  }
  return(FALSE)
}

  # Sync column selection with mouse-down
ViewMotionNotify <- function(widget, event){
#  while(gtkEventsPending())
#    gtkMainIteration()  
#  allColumns <- .local$allColumns
  #.local$last.coord <- 
  if(0){ # code for making pointer go to cross on drag handle
  pointer <- gdkWindowGetPointer(event[["window"]])
  pp <- gtkTreeViewGetPathAtPos(widget, event[["x"]], event[["y"]])
  dis <- gtkWidgetGetDisplay(widget)  
  win <- gtkWidgetGetWindow(widget)    
  if ((pp$cell.x-2)^2 + (pp$cell.y-2)^2 < 50){
    gdkWindowSetCursor(win, gdkCursorNewFromName(dis, "cross"))
  } else {
    gdkWindowSetCursor(win, gdkCursorNewFromName(dis, "arrow"))
  }
  }

#  if (as.flag(pointer$mask) & GdkModifierType["button1-mask"]){
#    info <- gtkTreeViewGetPathAtPos(widget, pointer[["x"]], pointer[["y"]])
#    if (info$retval){
#      col.idx <- GetColIdx(info$column)
#      if(is.null(.local$start.select.column)) .local$start.select.column <- col.idx
#		  new.sel <- sort(.local$start.select.column:col.idx)
#		  UpdateColumnSelection(.local, new.sel)     
#   }
#  } 
  return(FALSE)
}

Cell3rdButtonMenu <- function(.local, row.idx, col.idx){
  
  m <- gtkMenu()

  pasteItem <- gtkMenuItem("Paste (Ctrl-V)")
  m$append(pasteItem)
  gSignalConnect(pasteItem, "activate", function(...)         
    ReadFromClipboardWrapper(.local))      

  if(!is.null(col.idx) && !is.null(row.idx)){
	  cutItem <- gtkMenuItem("Cut")
	  copyItem <- gtkMenuItem("Copy (Ctrl-C)")
	  copyWithNamesItem <- gtkMenuItem("Copy With Names (Ctrl-Shift-C)")
	  #m$append(cutItem)
	  m$append(copyItem)
	  m$append(copyWithNamesItem)	  
	  #lapply(c(cutItem), gtkWidgetSetSensitive, FALSE)
	  gSignalConnect(copyItem, "activate", 
		function(...) {
		  w <- TransientWindow("Copying...", .local)
		  on.exit(w$destroy())    
      nf <- GetSelectedCells(.local)
      if(!is.null(nf))
        CopyToClipboard(nf, 
        do.rownames=FALSE, do.colnames=FALSE)
    })
	  gSignalConnect(copyWithNamesItem, "activate", 
		function(...) {
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())   
      nf <- GetSelectedCells(.local)
      if(!is.null(nf))
        CopyToClipboard(nf, 
        do.rownames=!RowNamesAreNull(.local$theFrame), do.colnames=!ColNamesAreNull(.local$theFrame))
    })
	  m$append(gtkSeparatorMenuItem())

	  view <- .local$view	
	  theFrame <- .local$theFrame


	  editFactorsItem <- gtkMenuItem("Edit Factors...")
	  randomizeItem <- gtkMenuItem("Randomize Selected")
	  fillItem <- gtkMenuItem("Fill Down")	  
	  fillCyclicItem <- gtkMenuItem("Fill In Cycles")
      m$append(editFactorsItem)
	  m$append(randomizeItem)
	  m$append(fillItem)
	  m$append(fillCyclicItem)
	  typ <- GetClasses(.local$theFrame)[col.idx+COLUMN_OFFSET]	  
	  if(typ != "factor")
		lapply(c(editFactorsItem, fillCyclicItem), gtkWidgetSetSensitive, FALSE)
	  gSignalConnect(editFactorsItem, "activate", function(...) {
	   # DoFactorEditor(theFrame, .local, col.idx + COLUMN_OFFSET))
		DoFactorEditor(theFrame, .local$toplevel, col.idx + COLUMN_OFFSET, 
		 FactorEditorHandler, data=.local)
	  })
	  gSignalConnect(randomizeItem, "activate", function(...){
		   sr <- GetSelectedRows(.local$view, .local)                            
		   if(length(sr)==0) sr <- 1:(dim(.local$theFrame)[1]-EXTRA_ROW)
		   if(length(sr)){
		     entry.frame <- theFrame[sr, GetSelectedColumns(.local)+COLUMN_OFFSET, drop=F]
		     entry.frame <- entry.frame[sample(1:(dim(entry.frame)[1])),,drop=F]
		     task <- list(list(func="ChangeCells", 
		        arg = list(nf=entry.frame, row.idx=sr, col.idx=GetSelectedColumns(.local)+COLUMN_OFFSET)))
		     DoTaskWrapper(.local, task)
		   }
		 })  
	  gSignalConnect(fillItem, "activate", function(...){
		   sr <- GetSelectedRows(.local$view, .local)
		   sc <- GetSelectedColumns(.local)+COLUMN_OFFSET
		   if(length(sr)==0) sr <- 1:(dim(.local$theFrame)[1]-EXTRA_ROW)  # None selected: fill the whole column
		   if(length(sr)==1) sr <- sr:(dim(.local$theFrame)[1]-EXTRA_ROW) # One selected: fill everything below
		   if(length(sr) && length(sc)){
		       # if the user's selected noncontiguous rows, we're still drawing just one selected rectangle
		       # so we want to fill the whole thing anyway           
		     sr <- range(sr)[1]:range(sr)[2]
		     entry.frame <- theFrame[sr, sc, drop=F]
		     for(jj in 1:length(sc))
		       entry.frame[,jj] <- entry.frame[1,jj]
		     task <- list(list(func="ChangeCells", 
		        arg = list(nf=entry.frame, row.idx=sr, col.idx=GetSelectedColumns(.local)+1)))
		     DoTaskWrapper(.local, task)
		   }
		 })  
	  gSignalConnect(fillCyclicItem, "activate", function(...) {
		   sr <- GetSelectedRows(.local$view, .local)
		   cc <- col.idx+COLUMN_OFFSET         
		   df1 <- dim(.local$theFrame)[1]
		   if(length(sr)==0) sr <- 1:(dim(.local$theFrame)[1]-EXTRA_ROW)  # None selected: fill the whole column
		   if(length(sr)==1) sr <- sr:(dim(.local$theFrame)[1]-EXTRA_ROW) # One selected: fill everything down              
		   if(length(sr) < 2 && df1 > 1) sr <- 1:(df1-EXTRA_ROW)
		            
		   if(length(sr)) {
		       # if the user's selected noncontiguous rows, we're still drawing just one selected rectangle
		       # so we want to fill the whole thing anyway           
		     sr <- range(sr)[1]:range(sr)[2]
		            
		     DoBlockSize(
		      theFrame[sr, cc],
		      .local$toplevel,
		      BlockSizeHandler,
		      data =  list(.local=.local, row.idx=sr, col.idx=cc), 
		      start.lvl=theFrame[sr[1], cc])
		     }
		 })
  } #row.idx, col.idx != null

  # if not editable, remove most options
  if(!.local$editable)	lapply(c(pasteItem, cutItem, editFactorsItem, randomizeItem, fillItem, fillCyclicItem), gtkWidgetSetSensitive, FALSE)

  lapply(c(randomizeItem, fillItem, fillCyclicItem), gtkWidgetSetSensitive, FALSE)

  return(m)
}

# Function to call when the data is sorted
SortHandler <- function(new.order, .local){
  dd <- dim(.local$theFrame)
  if(dd[1] > 1){
    theFrame <- .local$theFrame[new.order, ,drop=F]
    task <- list(
      list(func="ChangeRowNames", 
        arg = list(theNames = rownames(theFrame), row.idx=1:(dd[1]-EXTRA_ROW))),
      list(func="ChangeCells", 
        arg = list(nf=theFrame, row.idx=1:(dd[1]-EXTRA_ROW))))
    DoTaskWrapper(.local, task)
  }
}

OpenFileWrapper <- function(fileName, .local){
	stopifnot(!is.null(fileName) && !is.na(fileName) && nchar(fileName))
	ext <- strsplit(fileName, "[.]")[[1]]
	ext <- tolower(ext[length(ext)])
	if(identical(ext, "csv")){
	  sep <- ifelse(identical(options("OutDec")[[1]], "."), ",", ";")
	} else if (identical(ext, "txt")) {
	  sep <- "\t"
	} else {
	  sep <- ""
	}

	rfc <- ReadFromClipboard(fromFile=TRUE, fileName=fileName,
	  sep=sep)
	dat <- rfc$xx
	if(!is.null(dat) && length(dim(dat))==2){
	  dat <- MakeInternalDataFrame(rfc$xx)
	  ReplaceEditWindow(dat, .local, kill.history=TRUE)
    dataset.name <- basename(fileName)
    if(nchar(ext))
      dataset.name <- substr(dataset.name, 1, nchar(dataset.name)-nchar(ext)-1)
    .local$dataset.name <- dataset.name
	}
}

OpenFile <- function(menuitem, user.data){
  .local <- user.data
	fileName <- my.gfile(type="open", multi=F)
	if(!is.null(fileName) && !is.na(fileName) && nchar(fileName) && file.exists(fileName)){
	  OpenFileWrapper(fileName, .local) 
    DoUserCall("OnLoad", list(sourc=fileName, typ="file"), .local)
	}
}

OpenURL <- function(menuitem, user.data){
    .local <- user.data
    rv <- run.dialog(dlg.list=list(title="Enter URL", url.stringItem="", label = "Enter URL here"))
   if(!is.null(rv) && nchar(rv$args$url)){
	  OpenFileWrapper(rv$args$url, .local) 
    DoUserCall("OnLoad", list(sourc=rv$args$url, typ="url"), .local)
   }
}

Corner3rdButtonMenu <- function(.local){	
  theFrame <- .local$theFrame
  m <- gtkMenu()
  cutItem <- gtkMenuItem("Cut")
  copyItem <- gtkMenuItem("Copy")
  pasteItem <- gtkMenuItem("Paste...")
  m$append(cutItem)
  m$append(copyItem)	  
  m$append(pasteItem)

  m$append(gtkSeparatorMenuItem())
  gSignalConnect(copyItem, "activate", function(...) {
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())    
    CopyToClipboard(MakeExternalDataFrame(theFrame, .local), do.rownames=!RowNamesAreNull(theFrame), do.colnames=!ColNamesAreNull(theFrame))
    })
  gSignalConnect(pasteItem, "activate", function(...){
      ReadFromClipboardWrapper(.local)
  })      
  lapply(c(cutItem), 
    gtkWidgetSetSensitive, FALSE)

  renameItem <- gtkMenuItem("Rename Dataset")
  gSignalConnect(renameItem, "activate", function(...) {
     obj <- .local$rowname.column.obj$eventbox
     EditHeaderBox(obj, handler = function(obj, event, data){
        .local = data$.local
        box = data$box
        label = data$label
        getText <- obj$getText()
    	  if(nchar(getText) > 0) {
          #label$setText(getText)          
      	  #assign(getText, MakeExternalDataFrame(.local$theFrame, .local$dataset.class), envir=.GlobalEnv)
          oldname <- .local$dataset.name
      	  .local$dataset.name <- getText	  
             # 3-13-10
        to.external <- MakeExternalDataFrame(.local$theFrame, .local)
				nam <- .local$dataset.name
				# clean up spaces
				my.assign <- function(nam){
					 eval(parse(text=paste(paste(".GlobalEnv", nam, sep="$"), "<- to.external")))
           DoUserCall("OnRename", list(name=nam, oldname=oldname), .local)
           DoUserCall("OnLoad", list(sourc=nam, typ="rdata"), .local)
        }
				tryCatch({
					my.assign(nam)
				}, error = function(e) {  # User called it something strange
					my.assign(deparse(nam))
				})
  
          message(paste("RGtk2DfEdit: Creating dataset", .local$dataset.name, "in global environment."))
        }
        obj$destroy()
        box$packEnd(label, FALSE, FALSE, 0)       
        gtkWidgetSetState(.local$rowname.column.obj$button, as.integer(0))
        FALSE
      }, data=list(.local=.local, data=NULL))            
    }
  )   #..
  m$append(renameItem)  
  m$append(gtkSeparatorMenuItem())  
          
  openItem <- gtkMenuItem("Open File (Ctrl-O)")  
  m$append(openItem)	        
  gSignalConnect(openItem, "activate", OpenFile, data=.local)
  
  openURLItem <- gtkMenuItem("Open URL (Ctrl-U)")  
  m$append(openURLItem)	        
  gSignalConnect(openURLItem, "activate", OpenURL, data=.local)

  saveItem <- gtkMenuItem("Save As CSV... (Ctrl-S)")  
  m$append(saveItem)	        
  gSignalConnect(saveItem, "activate", function(...){
    tryCatch({	    
      theName <- paste(.local$dataset.name, "csv", sep=".")
      theFile <- my.gfile(initialfilename=theName, type="save", multi=F)
      if(length(theFile) > 0 && nchar(theFile) && !is.na(theFile)){
        write.csv(MakeExternalDataFrame(.local$theFrame, .local), theFile, quote=F, row.names=!RowNamesAreNull(.local$theFrame))
      }
    }, error = function(e){warning(e)})       
  })      
      
  m$append(gtkSeparatorMenuItem())  
    	  
  sortItem <- gtkMenuItem("Sort...")	  
  gSignalConnect(sortItem, "activate", function(...){
    dd <- dim(.local$theFrame)
    if (EXTRA_ROW) DoSortDialog(.local$theFrame[-dd[1], -dd[2],drop=F], SortHandler, .local)
    else 
      DoSortDialog(.local$theFrame[,-dd[2],drop=F], SortHandler, .local)
  })
  m$append(sortItem)
  m$append(gtkSeparatorMenuItem())
  # move dataset 1 column along
  ordinalItem <- gtkMenuItem("Default Rows")
  gSignalConnect(ordinalItem, "activate", function(...) {
    dd1 <- dim(.local$theFrame)[1]
    if(dd1 > 1){
      task <- list(
        list(func="ChangeRowNames", 
          arg = list(theNames = 1:(dd1-EXTRA_ROW), row.idx=1:(dd1-EXTRA_ROW))))
      DoTaskWrapper(.local, task)
    }
  })    

  m$append(ordinalItem)      
  m$append(gtkSeparatorMenuItem())
  ordinal2Item <- gtkMenuItem("Default Columns")
  gSignalConnect(ordinal2Item, "activate", function(...) {
    dd2 <- dim(.local$theFrame)[2]
    if(dd2 > 2){
      task <- list(
        list(func="ChangeColumnNames", 
          arg = list(theNames = DEFAULT_COLNAMES[1:(dd2-2)], col.idx=2:(dd2-1))))
      DoTaskWrapper(.local, task)
    }
  })

  m$append(ordinal2Item)      
  m$append(gtkSeparatorMenuItem())

  #lapply(renameItem, gtkWidgetSetSensitive, FALSE)
  
  m$append(gtkSeparatorMenuItem())
  aboutItem <- gtkMenuItem("About...")
  m$append(aboutItem)	  
  gSignalConnect(aboutItem, "activate", function(...){
	  dlg <- gtkAboutDialogNew(show=F)
    dlg["authors"] <- c("Tom Taverner <t.taverner@gmail.com>", 
"for the Department of Energy (PNNL, Richland, WA)",
"2010, Battelle Memorial Institute", 
"Contributions from John Verzani and Liviu Andronic", 
"RGtk2 by Michael Lawrence and Duncan Temple Lang",
"cairoDevice by Michael Lawrence")
    dlg["program-name"] <- "RGtk2DfEdit"
    dlg["comments"] <- "A spreadsheet data frame editor for the R environment"
    dlg["copyright"] <- "GPLv2"
    dlg["version"] <- VERSION_STRING
    gtkDialogRun(dlg)
    gtkWidgetDestroy(dlg)
  })

  if(!.local$editable)	lapply(c(pasteItem, cutItem, renameItem, sortItem, ordinalItem, ordinal2Item), gtkWidgetSetSensitive, FALSE)
    	  
  return(m)
}


Row3rdButtonMenu <- function(.local, row.idx){	  
  theFrame <- .local$theFrame
  m <- gtkMenu()
  cutItem <- gtkMenuItem("Cut")
  copyItem <- gtkMenuItem("Copy (Ctrl-C)")
  pasteItem <- gtkMenuItem("Paste (Ctrl-V)")

  m$append(cutItem)
  m$append(copyItem)
  m$append(pasteItem)
  m$append(gtkSeparatorMenuItem())
  gSignalConnect(copyItem, "activate", function(...) {
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())      
     gsr <- GetSelectedRows(.local$view.rn, .local)
     CopyToClipboard(
        theFrame[gsr, -c(1, dim(theFrame)[2]), drop=F], 
        do.rownames=!RowNamesAreNull(theFrame), do.colnames=F)
     })
  gSignalConnect(pasteItem, "activate", function(...) {
   dat <- ReadFromClipboard(.local, row.names=1, sep="\t", stringsAsFactors=F)$xx # character vector
  if(!is.null(dat)){
	if(!length(row.idx)) row.idx <- 1
	col.idx <- 1
	row.idx <- row.idx[1]
	task <- GetTaskPasteIn(.local$theFrame, dat, row.idx, col.idx+COLUMN_OFFSET, do.colnames=F, do.rownames = F)
	DoTaskWrapper(.local, task)
  }
  })                  

  renameItem <- gtkMenuItem("Rename Row")
  m$append(renameItem)  
  m$append(gtkSeparatorMenuItem())  
  gSignalConnect(renameItem, "activate", function(...) {
    info <- gtkTreeViewGetCursor(.local$view.rn)
    gtkTreeViewSetCursorOnCell(.local$view.rn, info$path, info$focus.column, info$focus.column$getCellRenderers()[[1]], TRUE)                                            
  })                  

  insertItem <- gtkMenuItem("Insert")
  gSignalConnect(insertItem, "activate", function(...) {
     task <- list(list(func="InsertNARows", 
                 arg = list(row.idx=row.idx)))
     DoTaskWrapper(.local, task)
  })
  deleteItem <- gtkMenuItem("Delete")
  gSignalConnect(deleteItem, "activate", function(...){	  
    sr <- GetSelectedRows(.local$view.rn, .local)
    if(length(sr)){
        # We have to block the selection-changed signal, otherwise we get a crash
        # This is because this fires the signal off for each deleted row!
	    task <- list(list(func="DeleteRows", arg = list(row.idx=sr)))
	    #gSignalHandlerBlock(.local$ss.rn, .local$ss.rn.changed.signal)
	    .local$rectangles <- list()	    
      DoTaskWrapper(.local, task)
      #gtkTreeSelectionUnselectAll(.local$ss.rn)      
	    #gSignalHandlerUnblock(.local$ss.rn, .local$ss.rn.changed.signal)      
    }
  })

  clearItem <- gtkMenuItem("Clear Contents")
  m$append(insertItem)
  m$append(deleteItem)
  m$append(clearItem)

  gSignalConnect(clearItem, "activate", function(...){	  	  
    row.idx <- GetSelectedRows(.local$view.rn, .local)
    RowNamesClearContents(row.idx, .local)
  })

  lapply(c(cutItem), gtkWidgetSetSensitive, FALSE)
  if(!length(GetSelectedRows(.local$view.rn, .local)))	  
	  lapply(c(deleteItem, pasteItem), gtkWidgetSetSensitive, FALSE)
  if(EXTRA_ROW && row.idx == dim(theFrame)[1])	  
	  lapply(c(cutItem, copyItem, deleteItem, pasteItem, clearItem, cutItem), 
      gtkWidgetSetSensitive, FALSE)  
      
  m$append(gtkSeparatorMenuItem())
  setAsNamesItem <- gtkMenuItem("To Column Names")
  m$append(setAsNamesItem)	  

  gSignalConnect(setAsNamesItem, "activate", function(...){
    theNames <- as.character(theFrame[row.idx,-1])
    theNames[is.na(theNames)] <- ""	    
    theNames <- make.unique(theNames)	    
    dd2 <- dim(.local$theFrame)[2]
    task <- list(
      list(func="ChangeColumnNames", 
        arg = list(theNames = theNames, col.idx=2:dd2)),
      list(func="DeleteRows", 
        arg = list(row.idx=row.idx))
     )
    DoTaskWrapper(.local, task)      
  })    	      	                	  

  if(!.local$editable)	lapply(c(pasteItem, cutItem, deleteItem, insertItem, clearItem, setAsNamesItem), gtkWidgetSetSensitive, FALSE)
  
  return(m)
}


ColumnClearContents <- function(col.idx, .local){
  if(length(col.idx)==0) return(FALSE)
  nf <- .local$model[,col.idx, drop=F]        
  stopifnot(ncol(nf) > 0)  
  nf[,] <- ""
  task <- list(list(func="ChangeCells", 
    arg = list(nf=nf, col.idx=col.idx)))
  DoTaskWrapper(.local, task)    
}

Column3rdButtonMenu <- function(.local, col.idx){	
  theFrame <- .local$theFrame
  typ <- GetClasses(theFrame)[col.idx+COLUMN_OFFSET]
  lastColumn <- col.idx == length(.local$allColumns)  # is this the last column
   
  theColumn <- theFrame[,col.idx+COLUMN_OFFSET,drop=F]
  m <- gtkMenu()
  cutItem <- gtkMenuItem("Cut")
  copyItem <- gtkMenuItem("Copy")
  pasteItem <- gtkMenuItem("Paste")
  m$append(cutItem)
  m$append(copyItem)
  m$append(pasteItem)
  m$append(gtkSeparatorMenuItem())
  gSignalConnect(copyItem, "activate", 
    function(...) {
      w <- TransientWindow("Copying...", .local)
      on.exit(w$destroy())        
      CopyToClipboard(.local$theFrame[-dim(.local$theFrame)[1],GetSelectedColumns(.local)+1,drop=F],
       do.colnames=T)
       })
  gSignalConnect(pasteItem, "activate", function(...) {
    #dat <- ReadFromClipboard(header=T, sep="\t", stringsAsFactors=F)$xx # character vector
    #task <- GetTaskPasteIn(.local$theFrame, dat, 
    #  1, col.idx+COLUMN_OFFSET, do.colnames= T)
    #DoTaskWrapper(.local, task)
	ReadFromClipboardWrapper(.local)
   })

  renameItem <- gtkMenuItem("Rename Column")
  gSignalConnect(renameItem, "activate", function(...) {
      obj <- .local$allColumns[[col.idx]]$eventbox
      EditHeaderBox(obj, handler = function(obj, event, data){
        col.idx = data$data
        box = data$box
        label = data$label
        .local <- data$.local  
        new.name <- obj$getText()
        task <- list(list(func="ChangeColumnNames", 
         arg = list(theNames = new.name, col.idx=col.idx+COLUMN_OFFSET)))
        obj$destroy()
        box$packEnd(label, FALSE, FALSE, 0)
        if (nchar(new.name) > 0 && new.name != colnames(.local$theFrame)[col.idx+COLUMN_OFFSET]){
          DoTaskWrapper(.local, task)
           #box$setTooltipText(new.name) 
box$setTooltipText(paste(DEFAULT_COLNAMES[col.idx], "\n", new.name, sep=""))
        }
        FALSE
      }, data=list(.local=.local, col.idx=col.idx))
    }
  )   #..

    
  insertItem <- gtkMenuItem("Insert")
  deleteItem <- gtkMenuItem("Delete")
  clearItem <- gtkMenuItem("Clear") 

  m$append(renameItem)
  m$append(gtkSeparatorMenuItem())    
  m$append(insertItem)
  m$append(deleteItem)
  m$append(clearItem)

  m$append(gtkSeparatorMenuItem())  
  gSignalConnect(clearItem, "activate", function(...) ColumnClearContents(GetSelectedColumns(.local)+COLUMN_OFFSET, .local))

  gSignalConnect(insertItem, "activate", function(...) {
     task <- list(list(func="InsertNAColumns", 
                 arg = list(col.idx=col.idx+COLUMN_OFFSET)))
     DoTaskWrapper(.local, task)                       
  })
  gSignalConnect(deleteItem, "activate", function(...) {
     sc <- GetSelectedColumns(.local)+COLUMN_OFFSET
     if(length(sc)){	    
       task <- list(list(func="DeleteColumns", 
                   arg = list(col.idx=GetSelectedColumns(.local)+COLUMN_OFFSET)))
       DoTaskWrapper(.local, task)	  
     }
 })  
  
  sortItem <- gtkMenuItem("Sort...")	  
  gSignalConnect(sortItem, "activate", function(...) {
    dd <- dim(.local$theFrame)
    DoSortDialog(.local$theFrame[, -dd[2],drop=F], SortHandler, .local)
  })
  m$append(sortItem)
  m$append(gtkSeparatorMenuItem())
    	  
  lapply(c(cutItem), gtkWidgetSetSensitive, FALSE)
  if(col.idx == length(.local$allColumns)) 
    lapply(c(deleteItem), gtkWidgetSetSensitive, FALSE)

    # If it's a factor, allow coercion to level type or numeric ordering
  if ("factor"%in%typ[[1]]) {
		item <- gtkCheckMenuItem("To Factor Levels")
		levelType <- class(levels(theFrame[,col.idx+COLUMN_OFFSET]))
		if(!is.atomic(levelType)) levelType <- "character" # default to character
		gSignalConnect(item, "button-release-event", function(item, evt, levelType){
	    task <- list(list(func="CoerceColumns", 
        arg = list(theClasses = levelType, col.idx=col.idx+COLUMN_OFFSET, to.levels=TRUE)))
      DoTaskWrapper(.local, task)  		  
			m$popdown()               
			return(TRUE)
		  }, levelType)
		m$append(item) 		
		
 		item <- gtkCheckMenuItem("To Factor Ordering")  		
		gSignalConnect(item, "button-release-event", function(item, evt, levelType){
	    task <- list(list(func="CoerceColumns", 
        arg = list(theClasses = "integer", col.idx=col.idx+COLUMN_OFFSET)))
      DoTaskWrapper(.local, task)  		  
			m$popdown()               
			return(TRUE)
		  }, levelType)
		m$append(item) 		 		
 		              
  }  else {  
    dataTypeNames <- list(Character="character", Integer="integer", Factor="factor", Logical="logical", Numeric="numeric")
    dataTypeItems <- list()
    for(theNewTypeName in names(dataTypeNames)){
  		item <- gtkCheckMenuItem(theNewTypeName)
  		item$setDrawAsRadio(TRUE)
  		dataTypeItems[[length(dataTypeItems)+1]] <- item
  		gSignalConnect(item, "button-release-event", function(item, evt, theNewTypeName){  		
  	    task <- list(list(func="CoerceColumns", 
          arg = list(theClasses = dataTypeNames[[theNewTypeName]], col.idx=GetSelectedColumns(.local)+COLUMN_OFFSET)))
        DoTaskWrapper(.local, task)  		  
  			m$popdown()               
  			return(TRUE)
  		  }, theNewTypeName)
  		if (dataTypeNames[[theNewTypeName]]%in%typ[[1]]) item$setActive(TRUE)  
  		m$append(item)
    }
  }
  m$append(gtkSeparatorMenuItem())
  editFactorsItem <- gtkMenuItem("Edit Factors...")  

  m$append(editFactorsItem)  
  m$append(gtkSeparatorMenuItem())
  abbreviateItem <- gtkMenuItem("Shorten")

  m$append(abbreviateItem)
  setAsNamesItem <- gtkMenuItem("To Row Names")

  m$append(setAsNamesItem)

  gSignalConnect(editFactorsItem, "activate", function(...) 
    DoFactorEditor(theFrame, .local$toplevel, col.idx + COLUMN_OFFSET, 
      FactorEditorHandler, data=.local))  
  if(typ != "factor") editFactorsItem$setSensitive(FALSE)
  
  gSignalConnect(abbreviateItem, "activate", function(...){
     abcol <- data.frame(X=cbind(abbreviate(as.character(theColumn[[1]]), minlength=10)))	     
     task <- list(
        list(func="CoerceColumns", 
          arg = list(theClasses = "character", idx=col.idx+COLUMN_OFFSET)),
         list(func="ChangeCells", 
         arg = list(nf=abcol, col.idx=col.idx+COLUMN_OFFSET))
    )
     DoTaskWrapper(.local, task)               
  })	  
  gSignalConnect(setAsNamesItem, "activate", function(...){
    theNames <- as.character(theColumn[[1]])
    theNames[is.na(theNames)] <- ""	    
    theNames <- make.unique(theNames)
    dd1 <- dim(.local$theFrame)[1]
    task <- list(
      list(func="ChangeRowNames", 
        arg = list(theNames = theNames, row.idx=1:dd1)),
      list(func="DeleteColumns", 
        arg = list(col.idx=col.idx+COLUMN_OFFSET))
    )
    DoTaskWrapper(.local, task)      
  })
  
 if(lastColumn) { # disable everything except for insert

   lapply(m$getChildren(), gtkWidgetSetSensitive, FALSE)
   lapply(list(insertItem), gtkWidgetSetSensitive, TRUE)
  }

 if(!.local$editable) { # disable everything except for insert
   lapply(m$getChildren(), gtkWidgetSetSensitive, FALSE)
   lapply(list(copyItem), gtkWidgetSetSensitive, TRUE)
  }

  return(m)
}

# just update the dataset name when we change the text
CornerBoxButtonPress <- function (obj, event, data){
  .local <- data$.local
  #gSignalHandlerBlock(.local$ss.rn, .local$ss.rn.changed.signal)
  #gtkWidgetGrabFocus(.local$view)
  #gSignalHandlerUnblock(.local$ss.rn, .local$ss.rn.changed.signal)

    # get out of editing or selection
  button <- event[['button']]
  typ <- event[['type']]
  stat <- event[['state']]
  col.idx <- data$col.idx
  if (button == as.integer(1)){
   if(typ == as.integer(4)){ # single clicked
     #gtkTreeViewGetSelection(.local$view.rn)$unselectAll()
     SelectAll(.local)
   } else if(0 && typ == as.integer(5)){ # turned off double clicked 
    gtkWidgetSetState(.local$rowname.column.obj$button, as.integer(1))
               
    EditHeaderBox(obj, handler = function(obj, event, data){
        .local = data$.local
        box = data$box
        label = data$labeld
        getText <- obj$getText()
    	  if(nchar(getText) > 0) {
          label$setText(getText)                                                      
      	  #assign(getText, MakeExternalDataFrame(.local$theFrame, .local$dataset.class), envir=.GlobalEnv)
      	  .local$dataset.name <- getText	  
             # 3-13-10
        to.external <- MakeExternalDataFrame(.local$theFrame, .local)
				nam <- .local$dataset.name
				# clean up spaces
				my.assign <- function(nam){
					 eval(parse(text=paste(paste(".GlobalEnv", nam, sep="$"), "<- to.external")))
          
        }
				tryCatch({
					my.assign(nam)
				}, error = function(e) {  # User called it something strange
					my.assign(deparse(nam))
				})  
          message(paste("RGtk2DfEdit: Creating dataset", .local$dataset.name, "in global environment."))
        }
        obj$destroy()
        box$packEnd(label, FALSE, FALSE, 0)       
        gtkWidgetSetState(.local$rowname.column.obj$button, as.integer(0))
        FALSE
      }, data=list(.local=.local, data=NULL))
  	} # end double clicked
  	} # end clicked
    if (button == as.integer(3)){ # our popup menu  
  	m <- Corner3rdButtonMenu(.local)
    gtkMenuPopupHack(m, button = event$GetButton(),
        activate.time = gdkEventGetTime(event))
    return(FALSE)
  }
  return(TRUE)
}

ModelChangeDatasetName <- function(.local, theName){
  .local$dataset.name <- theName
  .local$rowname.column.obj$eventbox$getChildren()[[1]]$getChildren()[[1]]$setText(theName)
}
   
# Sets up column header box for editing
# handler is function to call on focus out event
# Passes to handler: list(data=data, box=box, label=label, .local=.local)
EditHeaderBox <- function(obj, handler, data){
  .local <- data$.local
  col.idx <- data$col.idx

  box = obj$getChildren()[[1]]
  label = box$getChildren()[[1]]
  height = box$allocation$height-HEADER_BOX_MARGIN
  width = box$allocation$width
  box$remove(label)
  entry <- gtkEntryNew()
  entry$setText(label$getText())
  entry$setHasFrame(FALSE)
  #entry$modifyBase(as.integer(1), selectedColumnColor)    
  entry$modifyBase(as.integer(1), as.GdkColor(c(255,255,255)*256))  
  makeBlack <- function(y) entry$modifyText(y, as.GdkColor("black"))
  sapply(as.integer(0:1), makeBlack)    
  entry$modifyBase(as.integer(1), as.GdkColor(c(255,255,255)*256))  
  #entry$setAlignment(1)
  if(is.numeric(height) && length(height)==1 && is.numeric(width) && length(width)==1){
    entry$setSizeRequest(width, height) # 1 pixel margins I guess?
  }
  box$packEnd(child=entry, expand=TRUE, fill=TRUE, padding=0)
  entry$grabFocus()
  gSignalConnect(entry, "key-press-event", function(obj, event, data=col.idx){
    #if (event[["keyval"]]%in%myValidNavigationKeys) .local$view$grabFocus()
  if (event[["keyval"]]%in%c(GDK_Return, GDK_Escape)) .local$view$grabFocus()
    return(FALSE)
  })    
  #gSignalConnect(entry, "button-press-event", function(obj, event, data=col.idx){
  #  return(TRUE)
  #})      
  
  gSignalConnect(entry, "focus-out-event", handler, data=list(data=col.idx, box=box, label=label, .local=.local))
  #gSignalConnect(entry, "unrealize", handle)  
  
  #.local$entry <- entry  
  return(TRUE)
}       

# column header clicked
# obj is eventbox
ColumnHeaderButtonPress <- function (obj, event, data){
  .local <- data$.local
  .local$draw.cursor <-FALSE
  view <- .local$view
  if(!is.null(.local$entry)) {
    if(!is.null(.local$entry.focus.out))
      gSignalHandlerDisconnect(.local$entry, .local$entry.focus.out)

    gtkWidgetUnrealize(.local$entry)
    .local$entry <- NULL
    .local$entry.focus.out <- NULL
    gtkWidgetGrabFocus(view)
  }

  if(0 && !is.null(.local$entry)) {
    #tryCatch({
    if(gtkWidgetGetRealized(.local$entry))
      gtkWidgetUnrealize(.local$entry)
    .local$entry <- NULL
    #}, silent=F)
  }

  if(0){
  # get out of editing or selection    
  # To avoid condition where column header was just pressed and 
  # rows were previously selected, we block the changed signal
      # kill any active entries

  if(0 && !view$isFocus()) {
	gSignalHandlerBlock(.local$ss.rn, .local$ss.rn.changed.signal)
    view$grabFocus()
    gtkTreeSelectionUnselectAll(.local$ss.rn)    
	gSignalHandlerUnblock(.local$ss.rn, .local$ss.rn.changed.signal)
  }
  } # if 0
  view.rn <- .local$view.rn      
  model <- .local$model
  allColumns <- .local$allColumns
  

  # unselect the main view
  if(0 && length(.local$rectangles)){
    .local$rectangles <- list()
     gdkWindowInvalidateRect(.local$viewGetBinWindow, NULL, FALSE)      
  }      
  if(length(.local$selections)){
    
  }
   #for(tv in c(view, view.rn)){
  ##  pfc <- gtkTreeViewGetCursor(tv)
  #  gtkTreeSelectionUnselectAll(gtkTreeViewGetSelection(tv))
  #}                         
  #gtkTreeSelectionUnselectAll(gtkTreeViewGetSelection(view))
  #gtkTreeSelectionUnselectAll(gtkTreeViewGetSelection(view.rn))

    # prevents cursor showing up
  #gtkWidgetGrabFocus(gtkWidgetGetParent(.local$view))  
  button <- event[['button']]
  typ <- event[['type']]
  stat <- event[['state']]
  col.idx <- data$col.idx
  
  #gtkTreeViewSetCursorOnCell(view, gtkTreePathNewFromString(0))

    # ignore d-click on last column
  # ignore if it's the end column  
  lastColumn <- col.idx == length(allColumns)

  selectedColumns.new <- integer(0) # tom 092610
    # update: don't let double-click edit.    
  if (0 && !lastColumn && button == as.integer(1) && typ == as.integer(5) && stat == as.integer(0)){ # double clicked
    EditHeaderBox(obj, handler = function(obj, event, data){
      col.idx = data$data
      box = data$box
      label = data$label
      .local <- data$.local

      new.name <- obj$getText()
      task <- list(list(func="ChangeColumnNames", 
       arg = list(theNames = new.name, col.idx=col.idx+COLUMN_OFFSET)))
      obj$destroy()
      box$packEnd(label, FALSE, FALSE, 0)
      if (new.name != colnames(.local$theFrame)[col.idx+COLUMN_OFFSET])
        DoTaskWrapper(.local, task)        
      FALSE
    }, data=list(.local=.local, col.idx=col.idx))
  }

  if (button == as.integer(1) && typ == as.integer(4)){ # column clicked
   selections <- .local$selections
   theFrameNRows <- nrow(.local$theFrame)
   colSel <- list(start=c(row.idx=1, col.idx=col.idx),
        end=c(row.idx=theFrameNRows, col.idx=col.idx))

   selectedColumns <- GetSelectedColumns(.local)  	
   if (as.flag(stat) & GdkModifierType['shift-mask']) { # range
    if(is.null(selections) || !length(selections))
      selections <- list(list(start=colSel$start))
    # For the tricky case where we've selected a cell range
    # and shift-click - unsure of best behavior 
    selections[[length(selections)]]$start['row.idx'] <- 1
    selections[[length(selections)]]$end <- colSel$end

  } else if (as.flag(stat) & GdkModifierType['control-mask']) { # range
    # Simpler: Just add the column to whatever you clicked.
    selections <- .local$selections

    if(length(selections)){
        # we've only selected the start - bad news
       if( length(selections[[length(selections)]]) == 1)
         return(FALSE);
      #if(length(selections)[[length(selections)]] == 1)
      #  selections[[length(selections)]]$end <- 
      #    selections[[length(selections)]]$start

      xors <- get.xor.selection(selections)
      aicx <- as.integer(colnames(xors))
      selectedColumns <- aicx[colSums(xors)>0]
        # if the column index is not in the range 
        # of selected columns,
        # we want to add it
      whichidx <- which(aicx==col.idx)
      # if the entire column is already selected and TRUE
      # then set it to FALSE 
      if(col.idx%in%selectedColumns){
        #dont.add.flag <- TRUE
        if(all(xors[,whichidx])
          && nrow(xors)==theFrameNRows){
          xors[,whichidx] <- FALSE
          selections <- get.xor.rectangles(xors)
        } else {
          xors[,whichidx] <- FALSE
          selections <- get.xor.rectangles(xors)
          selections[[length(selections)+1]] <- colSel
        }
      # otherwise set the whole column to TRUE
      } else {
         selections[[length(selections)+1]] <- colSel
      }
    # if length(selections)
    } else {
       selections[[length(selections)+1]] <- colSel
    }

  } else {
    gtkTreeSelectionUnselectAll(.local$ss)    
    selections <- list()
    selections[[length(selections)+1]] <- colSel
 }
 .local$selections <- selections
  .local$do.paint <- TRUE

    UpdateSelectionRectangle(.local)
    gdkWindowInvalidateRect(.local$viewGetBinWindow, NULL, FALSE)

  #return(FALSE)
#  UpdateColumnSelection(.local, selectedColumns.new)  

  DoUserCall("ColumnClicked", list(idx = data$col.idx), .local)
   
  } # clicked
        
  if (button == as.integer(3)){ # our popup menu

  selectedColumns <- GetSelectedColumns(.local)
  if(!length(selectedColumns)){
    selectedColumns.new <- col.idx
    UpdateColumnSelection(.local, selectedColumns.new)  
  }

	  m <- Column3rdButtonMenu(.local, col.idx)
    gtkMenuPopupHack(m, button = event$GetButton(),
        activate.time = gdkEventGetTime(event))
    return(FALSE)
   }

  return(TRUE)
}

GetVerticalPosition <- function(.local){
  .local$sw.view.va$getValue()
}

SetVerticalPosition <- function(.local, value){
  .local$sw.view.va$setValue(value)
  .local$sw.rn.va$setValue(value)
}
# replace the entire gtktab
# new.df is the internal DF representation
# want to keep cursor position if it exists
ReplaceEditWindow <- function(theFrame, .local, kill.history=FALSE){
  old.v.pos <- GetVerticalPosition(.local)
  #old.h.pos <- .local$sw.ha$getValue()   

  UpdateColumnSelection(.local, integer(0))       
  if(kill.history) .local$undoStack <- NULL
  .local$theFrame <- theFrame
  gtktab.new <- MakeDFEditWindow(.local, theFrame) 
    
  .local$gtktab$destroy()  
  .local$gtktab <- gtktab.new
  .local$group.main$packStart(.local$gtktab, TRUE, TRUE, 0)
  .local$start.column.select <- NULL

  gTimeoutAdd(50, function(){SetVerticalPosition(.local, old.v.pos); return(FALSE)})
}

NewColumnToggle <- function(model, j, width=64, is.row = F, is.editable=T, resizable=T){
	renderer <- gtkCellRendererToggleNew()
	column <- gtkTreeViewColumnNew()  
    #column$packStart(renderer)
   gtkTreeViewColumnPackStart(renderer)
    #column$setTitle(as.character(j-1))
    gtkTreeViewColumnSetTitle(as.character(j-1))
	gtkTreeViewColumnSetFixedWidth(column, width)
	gtkTreeViewColumnSetSizing(column, GtkTreeViewColumnSizing['fixed'])
	gtkTreeViewColumnSetResizable(column, resizable)
	return(list(column=column, renderer=renderer, col.idx=j))
}

# Old code for pretty print:
#          txt <- .Call("R_getGObjectProps", cell, "text", PACKAGE = "RGtk2") # ~5 microseconds
#          if(txt == "1.#QNAN0" || txt == "-1.#IND00") # ~7 microseconds
#            .Call("R_setGObjectProps", cell, list(text="NA"), PACKAGE = "RGtk2")
#          else # ~13 microseconds, 22 if you suppress warnings
#            .Call("R_setGObjectProps", cell, list(text=sprintf(SPRINTF_FORMAT, as.numeric(txt))), PACKAGE = "RGtk2")    
# create a column with the title as the index
NewColumn <- function(model, j, width=64, is.row = F, is.editable=T, resizable=T, min.width=20, max.width=150){
	renderer <- gtkCellRendererTextNew()
	#gtkCellRendererTextSetFixedHeightFromFont(renderer, 1)
	gObjectSetData(renderer, "column", j-1)
	renderer['xalign'] <- 1
    #if(is.editable){
	renderer['editable-set'] <- TRUE
	renderer['editable'] <- is.editable
    #}
	column <- gtkTreeViewColumnNewWithAttributes(as.character(j-1), renderer, text = j-1)  
	gtkTreeViewColumnSetFixedWidth(column, width)
    gtkTreeViewColumnSetMinWidth(column, min.width) # for min width
    gtkTreeViewColumnSetMaxWidth(column, max.width) # for min width
	gtkTreeViewColumnSetSizing(column, GtkTreeViewColumnSizing['fixed'])
	#gtkTreeViewColumnSetResizable(column, TRUE)
	return(list(column=column, renderer=renderer, col.idx=j))
}

# Must be called after widget is mapped
MakeButtonAndEventBox <- function(col.obj, label.str,  handler, .local){
	label <- gtkLabelNew(label.str)
	box <- gtkVBoxNew()
	#gtkBoxPackEnd(box, label, FALSE, FALSE, 0)
	eventbox <- gtkEventBoxNew()
	gtkContainerAdd(eventbox, box)
	view.col <- col.obj$column
	gtkTreeViewColumnSetWidget(view.col, eventbox)
	alignment <- gtkWidgetGetParent(eventbox)
	gtkAlignmentSet(alignment, 0, 1, 1, 1)
	col.idx <- GetColIdx(view.col)

	#label2 <- gtkLabelNew(DEFAULT_COLNAMES[col.idx])
	#gtkContainerAdd(box, label2)

  box$packEnd(child=label, expand=TRUE, fill=TRUE, padding=0)

	#gtkContainerAdd(box, label)

	col.obj$eventbox <- eventbox
  gtkWidgetModifyBg(eventbox, as.integer(1), selectedColumnColor)
	col.obj$button <- gtkWidgetGetParent(gtkWidgetGetParent(alignment))
	gSignalConnect(eventbox, "button-press-event", handler, data=list(col.idx=col.idx, .local=.local))
  if(col.idx==0) {
    txt <- "Right-click for options"
  } else {
    if (DEFAULT_COLNAMES[col.idx] == label.str)
      txt <- DEFAULT_COLNAMES[col.idx]
    else
      txt <- paste(DEFAULT_COLNAMES[col.idx], 
        "\n", label.str, sep="")
  }
  gtkWidgetSetTooltipText(eventbox, txt)
	return(col.obj)
}


# update.object: if you update the data frame at the same time                              
UpdateDfEditor <- function(.local, theFrame, rows.changed=NULL){ 
  if (is.null(rows.changed)){
    if(ncol(theFrame) != ncol(.local$model)){
      ReplaceEditWindow(theFrame, .local)
    } else {
      cn <- colnames(theFrame)
      for(jj in which(cn != colnames(.local$theFrame)))
        ModelChangeColumnName(.local, jj, cn[jj])        
      .local$model$setFrame(theFrame)
    }
    .local$LAST_PATH <- MakeLastPath(theFrame)
  } else {
     .RGtkCall("R_r_gtk_data_frame_set", .local$model, theFrame, as.list(as.integer(rows.changed - 1)), F)
  }
  .local$theFrame <- theFrame
  if(identical(.local$update.frame, TRUE)){
	  to.external <- MakeExternalDataFrame(.local$theFrame, .local)
	 	#assign(.local$dataset.name, to.external, envir=.GlobalEnv)   
	#3-13-10
	# Change this to work with embedded lists - different choices for data names  
	  # if it's "iris", save it as iris
	  # if it's "xx$iris" and "xx" exists and is a list, save it as xx$iris
	  # if it's "abc$iris" and "abc" isn't a list, save it as "abc$iris"
		#3-13-10
	  nam <- .local$dataset.name
	  # clean up spaces
	  my.assign <- function(nam){
					 eval(parse(text=paste(paste(".GlobalEnv", nam, sep="$"), "<- to.external")))
        }
		tryCatch({
		my.assign(nam)
		}, error = function(e) {  # User called it something strange
		my.assign(deparse(nam))
	  })
  }
}

ModelChangeColumnName <- function(.local, idx, new.name)
  .local$allColumns[[idx-1]]$eventbox$getChildren()[[1]]$getChildren()[[1]]$setText(new.name)

 	  # function that takes a NewColumn returned object and adds 
    # eventboxes, buttons, etc
    # MUST BE DONE AFTER VIEW IS MAPPED
InitColumn <- function(.local, col.obj, nam){
  col.obj <- MakeButtonAndEventBox(col.obj, label.str=nam, handler=ColumnHeaderButtonPress, .local=.local)
  if(0) gSignalConnect(col.obj$renderer, "editing-started", function(renderer, editable, path, data=list(.local=.local)){
    #.local <- data$.local
#print("started")
    FALSE
  })

  #if(1) gSignalConnect(col.obj$renderer, "editing-started", function(renderer, entry, path, data) {

#  print("Here")
#       gSignalConnect(entry, "map-event", function(entry, event, data=.local){ 
#print("A")
#Sys.sleep(1)  
#print("B")
#    insertion_point <- 0
#    pixel_size <- entry$getLayout()$getPixelSize()
#    click_info <- .local$last_click_info
#    col_width <- click_info$column$getWidth()
#    text_width <- pixel_size$width
#    text_pos <- (click_info$cell.x - col_width + text_width)/(text_width)
#    if(length(text_pos) != 1 || text_pos < 0 || text_pos > 1) text_pos <- 0
#    insertion_point <- round(text_pos*nchar(entry$getText()))   
#    gtkEditableSelectRegion(entry, insertion_point, insertion_point)
#    return(TRUE)  
# })      
#  }, data=list(codl.idx=col.obj$col.idx, .local=.local))



  if(1) gSignalConnect(col.obj$renderer, "editing-started", after=T, RendererEditingStarted, data=list(col.idx=col.obj$col.idx, .local=.local))
  col.obj$state.set <- NULL    
  return(col.obj)
}  

InitAllColumns <- function(.local, model, allColumns){
  ds.colnames <- colnames(model)
  for(j in 1:length(allColumns)){
    allColumns[[j]] <- InitColumn(.local, allColumns[[j]], ds.colnames[j+1])
    allColumns[[j]]$col.idx <- j+COLUMN_OFFSET
  }
  # The last column is not editable
  allColumns[[j]]$renderer['editable'] <- FALSE 
  gtkTreeViewColumnSetResizable(allColumns[[j]]$column, FALSE)
  #gtkWidgetModifyBg(allColumns[[j]]$eventbox, as.integer(1), bgColor)    		  
  return(allColumns)
}


MakeAVerticalScrollbar <- function(.local, gtktab, va){
  vbar <- gtkVScrollbarNew(va)
  gtktab$attach(vbar, 1, 2, 0, 1, 0, 5)
    # scroll on this doesn't repaint, kill it
  vbar$setEvents(GdkEventMask["button-motion-mask"])	  
  gSignalConnect(vbar, "scroll-event", function(...) {
    return(TRUE)
  })
  gSignalConnect(vbar, "button-press-event", function(...) {
     if(!gtkWidgetIsFocus(.local$view) && !gtkWidgetIsFocus(.local$view.rn)) .local$view$grabFocus()
     return(FALSE)
  })
  gSignalConnect(vbar, "button-release-event", function(...) {
    .local$viewGetBinWindow$invalidateRect(NULL, FALSE)
    return(FALSE)
  })
  return(vbar)
}


RowNamesExpose <- function(widget, event=NULL, data=NULL){
  .local <- data
  if(.local$do.paint && DO_CAIRO){
    currentVadj <- gtkAdjustmentGetValue(.local$sw.view.va) # this 
    for(r in .local$rn.rectangles){
     r$y <- r$y - currentVadj
      selectedColorRgb <- c(49, 106, 197)
      cr <- gdkCairoCreate(.local$view.rnGetBinWindow)
      cairoSetSourceRgba(cr, selectedColorRgb[1]/256, 
        selectedColorRgb[2]/256, selectedColorRgb[3]/256, 0.2)    
      cairoRectangle(cr, r$x, r$y+1, r$width, r$height)
      cairoClip(cr)
      cairoPaint(cr)
    }
  }
  return(FALSE)
}


ViewExpose <- function(widget, event=NULL, data=NULL){
  .local <- data
  isWindows <- PLATFORM_OS_TYPE == "windows"
  if(.local$do.paint){  
  #cat("*")
    currentVadj <- gtkAdjustmentGetValue(.local$sw.view.va) # this gets called *before* value-changed

    for(r in .local$rectangles){
      r$y <- r$y - currentVadj
      gdkDrawRectangle(.local$viewGetBinWindow, .local$gc.invert, filled = F,
         ifelse(isWindows, r$x, r$x-1),
         ifelse(isWindows, r$y, r$y),
         r$width+2,
         ifelse(isWindows, r$height+2, r$height))
    }

   if(DO_CAIRO) {
    for(r in .local$rectangles){
      r$y <- r$y - currentVadj
      selectedColorRgb <- c(49, 106, 197)
      cr <- gdkCairoCreate(.local$viewGetBinWindow)
      cairoSetSourceRgba(cr, selectedColorRgb[1]/256, 
        selectedColorRgb[2]/256, selectedColorRgb[3]/256, 0.2)    
#      cairoRectangle(cr, r$x, r$y+1, r$width, r$height)
      cairoRectangle(cr,
        r$x,
        r$y+1,
        ifelse(isWindows, r$width+2, r$width),
        ifelse(isWindows, r$height+2, r$height))
      cairoClip(cr)
      cairoPaint(cr)
    }
  }

  } # display mode

  if(.local$draw.cursor){
  # Draw the cursor
  cursor.info <- gtkTreeViewGetCursor(widget)
  path <- cursor.info$path
  column <- cursor.info$focus.column
  if(is.null(.local$entry) && !is.null(path) && !is.null(column)){
	  rect <- gtkTreeViewGetCellArea(widget, path, column)$rect
    gdkDrawRectangle(.local$viewGetBinWindow, .local$gc.cursor, filled = F, rect$x, rect$y +1, rect$width, 
    ifelse(isWindows, rect$height-1, rect$height-2))
    gdkDrawRectangle(.local$viewGetBinWindow, .local$gc.cursor, filled = T, rect$x+rect$width-4, rect$y+rect$height-1-3, 4, 3)
  }
 }

  return(FALSE)
}


DoScrollUpdate <- function(obj, evt, data){
  .local <- data
#  if(!gtkWidgetIsFocus(obj)) gtkWidgetGrabFocus(obj)    
  while(gtkEventsPending())
    gtkMainIteration()
  dir <- evt[["direction"]]
  va <- .local$vbar$getAdjustment()
  if (dir == GdkScrollDirection["down"]) {
    new.val <- min(va$getValue() + va[["step_increment"]], va[["upper"]] - va[["page_size"]])    
    va$setValue(new.val)
  } else if (dir == GdkScrollDirection["up"]) {
    va$setValue(va$getValue() - va[["step_increment"]])
  }
  FALSE
}  

MakeDFEditWindow <- function(.local, theFrame, size.request=c(600, 300), col.width=64){  

  .local$allColumns	<-	NULL
  .local$allow.key.flag	<-	TRUE
  .local$col.width	<-	col.width
  .local$do.paint	<- TRUE
  .local$doing.scroll	<-	FALSE
  .local$doingEntryIntoFactor	<-	FALSE
  .local$doingHScroll	<-	FALSE
  .local$FIRST_PATH	<-	gtkTreePathNewFromString(0)
  .local$flash.cursor	<-	TRUE
  .local$last.time	<-	0
  .local$last.place.clicked	<-	NULL
  .local$LAST_PATH	<-	NULL
  .local$model	<-	NULL
  .local$SCROLL_ROW_TIMEOUT	<-	150
  .local$scrollID	<-	0
  .local$draw.cursor <- TRUE
  .local$theFrame <- theFrame
  

  sw.view <- gtkScrolledWindowNew()
  .local$sw.view <- sw.view
  sw.view$setPolicy(GtkPolicyType['never'], GtkPolicyType['never'])
  sw.view$setSizeRequest(size.request[1], size.request[2])
  
  .local$LAST_PATH <- MakeLastPath(theFrame)  
  dataset.name <- .local$dataset.name
  .local$model <- rGtkDataFrame()
  view <- gtkTreeViewNewWithModel(.local$model)  
  .local$model$setFrame(.local$theFrame)  # This is causing a memory leak...

  gtktab <- gtkTableNew(2,2, FALSE)   
  .local$mapped <- FALSE
  .local$view <- view
  gtkTreeViewSetFixedHeightMode(view, TRUE)
  allColumns <- vector("list", (dim(theFrame)[2]-1))                    

  view.rn <- gtkTreeViewNewWithModel(.local$model)
  .local$view.rn <- view.rn  
  view.rn$SetFixedHeightMode(TRUE)
  rowname.column.obj <- NewColumn(.local$model, 1, width=10, is.row=T, resizable=F, is.editable=.local$editable)
  .local$rowname.column.obj <- rowname.column.obj
  gtkTreeViewAppendColumn(view.rn, rowname.column.obj$column)

  col.sizing <- ifelse(identical(.local$autosize, TRUE), 'autosize', 'fixed')
  for(j in 2:(dim(.local$model)[2])){
    tmp <- NewColumn(.local$model, j, width=col.width, is.editable=T)#.local$editable)
    gtkTreeViewAppendColumn(view, tmp$column)
    gtkTreeViewColumnSetSizing(tmp$column, GtkTreeViewColumnSizing[col.sizing])
	gtkTreeViewColumnSetResizable(tmp$column, TRUE)
    allColumns[[j-1]] <- tmp
  }       

  ss <- view$getSelection()
  .local$ss <- ss  
  ss$setMode(as.integer(3)) # multiple
  view$setRubberBanding(TRUE)

  sw.view$add(view)
  sw.view.va <- sw.view$getVadjustment()
  .local$sw.view.va <- sw.view.va   
  sw.view.ha <- sw.view$getHadjustment()
  .local$sw.ha <- sw.view.ha    
  .local$va <- sw.view.va
  
  ss.rn <- view.rn$getSelection()
  # 09-16-10 - this is screwing up memory caching
  # 10-21-10 - this is causing a crash on row deletion
  if(0) .local$ss.rn.changed.signal <- gSignalConnect(ss.rn, "changed",  function(treeselection, data){
##print(.local$allow.key.flag)
#      .local$allow.key.flag <- FALSE
     #w <- TransientWindow("Updating Selection...", .local)
#     # on.exit(w$destroy())  
      #PaintRowSelectionOnTimeout(.local)
      PaintSelectionOnTimeout(.local, selection=treeselection, widget=.local$view.rn)
      #.local$allow.key.flag <- TRUE
#     if(.local$do.paint){
#       .local$viewGetBinWindow$invalidateRect(NULL, FALSE)      
#      .local$do.paint <- FALSE
#     }
    }
  )             

  .local$ss.rn <- ss.rn    
  ss.rn$setMode(as.integer(3)) # multiple
  
  sw.rn <- gtkScrolledWindowNew()
  sw.rn$add(view.rn)
  view.rn$setSizeRequest(-1, 10)
  sw.rn$setPolicy(GtkPolicyType['never'], GtkPolicyType['never'])
  sw.rn.va <- sw.rn$getVadjustment()
  .local$sw.rn.va <- sw.rn.va

#if(0){
  for(ii in as.integer(0:4))
  gtkWidgetModifyBase(view, ii,  whiteColor)
  #gtkWidgetModifyBase(view, GtkStateType['active'], whiteColor)
  gtkWidgetModifyText(view, GtkStateType['selected'], as.GdkColor("black"))
  gtkWidgetModifyText(view, GtkStateType['active'], as.GdkColor("black"))
  sapply(list(view, view.rn), function(x){
    sapply(as.integer(0:4), 
      function(y) gtkWidgetModifyFg(x, y, as.GdkColor("black")))
  })

  # It would be nice to change the background raised appearance
  # for selected rows

  #img_pixbuf = gdkPixbufNewFromFile(image_filename)
  #img_pixbuf <- gdkPixbufNew("GDK_COLORSPACE_RGB", FALSE, 8, 1, 1)
#img_pixmap = img_pixbuf.render_pixmap_and_mask()[0]
#for state in (gtk.STATE_NORMAL, gtk.STATE_ACTIVE, #gtk.STATE_PRELIGHT,
#              gtk.STATE_SELECTED, gtk.STATE_INSENSITIVE):
#    style$bg_pixmap[state] = img_pixmap
#treeview.set_style(style)

  for(ii in as.integer(0:4))
  gtkWidgetModifyBase(view.rn, ii, bgColor)

  gtkWidgetModifyText(view.rn, GtkStateType['selected'], as.GdkColor("black"))
  gtkWidgetModifyText(view.rn, GtkStateType['active'], as.GdkColor("black"))  
#}
  view.rn$setEnableSearch(FALSE)
  view$setEnableSearch(FALSE)  

  #view$setGridLines(as.integer(3))  
  #view.rn$setGridLines(as.integer(3))    

  .local$style <- view$getStyle

  paned <- gtkHPanedNew()
  paned$add1(sw.rn)
  paned$add2(sw.view)
  paned$setPosition(100)  
  gtktab$attach(paned, 0, 1, 0, 1, 5, 5)
  hbar <- gtkHScrollbarNew(sw.view.ha)
  gtktab$attach(hbar, 0, 1, 1, 2, 5, 0)
  
  vbar <- MakeAVerticalScrollbar(.local, gtktab, sw.view.va)
  .local$vbar <- vbar  
    
  ## we're doing all this after we map
  ## Notebooks map too!

  gSignalConnect(view, "map", after=T, data=.local, function(view, data){   
    .local <- data
    if(.local$mapped) return(TRUE)
    .local$mapped <- TRUE    
    .local$do.paint <- FALSE        
    #.local$rowname.column.obj <- MakeButtonAndEventBox(.local$rowname.column.obj, 
    # label.str=.local$dataset.name, handler=CornerBoxButtonPress, .local=.local)      
    # Change the displayed name to "" to avoid confusion.
    .local$rowname.column.obj <- MakeButtonAndEventBox(.local$rowname.column.obj, 
      label.str="", handler=CornerBoxButtonPress, .local=.local)            
             
    allColumns <- InitAllColumns(.local, .local$model, allColumns)
    .local$allColumns <- allColumns    
  
   gSignalConnect(view,"key-press-event", ViewKeyPress, data=.local)  
   #gSignalConnect(view,"key-release-event", ViewKeyRelease, data=.local)  

    gSignalConnect(view,"button-press-event", ViewButtonPress, data=.local)
    gSignalConnect(view,"button-release-event", after=T, ViewButtonRelease, data=.local)    
    #gSignalConnect(view,"motion-notify-event", ViewMotionNotify)  
    gSignalConnect(view.rn,"key-press-event", RowNamesKeyPress, data=.local)          
    gSignalConnect(view.rn,"button-press-event", RowNamesButtonPress, data=.local)

#gSignalConnect(view.rn,"button-press-event", function(...) {print("called"); return(FALSE)}, data=.local)

      # tree_view != NULL failed
  	gSignalConnect(view.rn, "focus-in-event", function(...) {
      .local$draw.cursor <- FALSE
       return(FALSE)
      #gtkTreeSelectionUnselectAll(.local$ss)
    })

  	gSignalConnect(view, "focus-in-event", function(...) {
     .local$draw.cursor <- TRUE
     gtkTreeSelectionUnselectAll(.local$ss.rn)
     return(FALSE)
    })
  
    #gSignalConnect(view, "cursor-changed", ViewCursorChanged, data=.local)
      

    rows.renderer <- view.rn$getColumns()[[1]]$getCellRenderers()[[1]] 
    .local$rows.renderer <- rows.renderer      
    gSignalConnect(rows.renderer,"edited", function(renderer, path, new.text){ 
      idx <- as.integer(path)+1
      if(new.text != row.names(.local$theFrame)[idx]){      
        task <- list(list(func="ChangeRowNames", 
          arg = list(theNames = new.text, row.idx=idx)))
        DoTaskWrapper(.local, task)
      }
      return(FALSE)
    })

    gSignalConnect(view, "expose-event", after=T, ViewExpose, data=.local)    
    #gSignalConnect(view.rn, "expose-event", after=T, RowNamesExpose, data=.local)    

    gSignalConnect(view, "scroll-event", after=T, DoScrollUpdate, data=.local)
    gSignalConnect(view.rn, "scroll-event", after=T, DoScrollUpdate, data=.local)    
    
    gSignalConnect(view, "leave-notify-event", AddHScrollTimeout, data=.local)
    lapply(c("enter-notify-event", "button-release-event", "focus-out-event"), 
      function(evt) gSignalConnect(view, evt, RemoveHScrollTimeout, data=.local)) 
      
    .local$viewGetBinWindow <- view$getBinWindow()		      


    .local$view.rnGetBinWindow <- view.rn$getBinWindow()		      

    .local$gc.invert <- gdkGCNew(.local$viewGetBinWindow) 
    .local$gc.invert$setFunction('GDK_INVERT')

    .local$gc.cursor <- gdkGCNew(.local$viewGetBinWindow) 
    .local$gc.cursor$setLineAttributes(line.width=1, line.style=integer(3), cap.style=integer(0), join.style=integer(0))


    toplevel <- .local$group.main$getToplevel()
    if (toplevel$flags() & GtkWidgetFlags["toplevel"]){
      .local$toplevel <- toplevel
    }  

      # Kill drag and drop, it causes a crash    
      # http://www.mail-archive.com/gtk-app-devel-list@gnome.org/msg11623.html      
    screen <- .local$group.main$getScreen()
    settings <- gtkSettingsGetForScreen(screen)
    settings[["gtk-dnd-drag-threshold"]] <- 10000
    checkPtrType(view, "GtkTreeView")              
    .local$scrollRowNames <- TRUE

    sw.va.value.ch <- function(obj, data){
      sw2 <- data$sw2
      .local <- data$.local        
      if(!is.null(.local$scrollRowNames) 
         && .local$scrollRowNames){
          # Take over the event loop
          # See http://wiki.laptop.org/go/PyGTK/Smooth_Animation_with_PyGTK
        while(gtkEventsPending())
          gtkMainIteration()
        gtkAdjustmentSetValue(sw2, gtkAdjustmentGetValue(obj))
        .local$allow.key.flag <- TRUE  # For paging, however, this has a problem...
      }
    }

    gSignalConnect(sw.view.va, "value-changed", sw.va.value.ch, data=list(sw2=sw.rn.va, .local=.local))
    
    # for some reason the last column doesn't get alignment set properly
  	alignment <- gtkWidgetGetParent(allColumns[[length(allColumns)]]$eventbox)
  	gtkAlignmentSet(alignment, 0, 1, 1, 1)    
     	
    return(TRUE)  
  }) # end of map event callback

       
  return(gtktab)
} # end MakeDFEditWindow

ConvertToDataObj <- function(items){
  .e <- new.env()
  .e$items <- items
  .e$dataset.attributes <- attributes(items)  

  if(is.ts(items)) {
    items <- data.frame(x=as.numeric(items), check.names=FALSE)
  } else {
    tryCatch(.e$items <- data.frame(.e$items, check.names=FALSE, stringsAsFactors=F), error = function(e)
      tryCatch(.e$items <- as.matrix(.e$items), error = function(e)
         stop(paste("Passed a", paste(class(.e$items), collapse=", "), "but expected a valid data object"))))
  }

  rv <- list(items=.e$items, dataset.attributes = .e$dataset.attributes)
  .e <- NULL
  return(rv)
}

#' Data frame editor for RGtk2
#'
#' @param items A data frame (?) to display graphically
#' @param dataset.name Name for data set
#' @param container An RGtk2 container object to place widget within. (Can this be missing?)
#' @return An object of class GtkDfEdit for which a few RGtk2-style methods are defined
#' @export
gtkDfEdit <- function(items, dataset.name = deparse(substitute(items)), 
  size.request=c(600, 300), col.width = 64,
  dataset.class="data.frame", editable = TRUE, 
  autosize = length(dim(items)) < 2 || ncol(items)<25,
  update=TRUE,
  envir=.GlobalEnv, ...){
   # our local environment
  #print(dim(.local$items))
  .local <- new.env()
  .local$items <- items
  .local$dataset.name <- dataset.name
  .local$dataset.class <- class(items)

  stopifnot(editable %in% c(TRUE, FALSE))
  .local$editable <- editable

  stopifnot(autosize %in% c(TRUE, FALSE))
  .local$autosize <- autosize

  stopifnot(update %in% c(TRUE, FALSE))
  .local$update.frame <- update
    
  .local$dataset.attributes <- attributes(items)
    
    #3-13-10
      # have we been passed a string referring to an object?
      # Might be "a", or "a$b$c"      
  if(!any(class(items)%in%DATA_OBJECTS)){
    if (is.character(items) && nchar(items)){
			.local$dataset.name <- items
      tryCatch({ # try eval(), then get() just in case it's a name using "$"
				.local$items <- safe.eval(items, envir=envir)
				.local$dataset.class <- class(.local$items)
				.local$dataset.attributes <- attributes(.local$items)
        }, 
        error = function(e) tryCatch({
          .local$items <- safe.eval(items, envir=envir)
  				.local$dataset.class <- class(.local$items)            
  				.local$dataset.attributes <- attributes(.local$items)  				
          }, error = function(e)
            stop("Passed a name that isn't a valid object") ) )
       if(!any(class(.local$items)%in%DATA_OBJECTS)) {
         .local$rv <- ConvertToDataObj(.local$items)
         .local$items <- .local$rv$items; .local$dataset.attributes <- .local$rv$dataset.attributes
       }   
    } else {
      .local$rv <- ConvertToDataObj(items)
      .local$items <- .local$rv$items; .local$dataset.attributes <- .local$rv$dataset.attributes
    } 
   }       
    
  .local$theFrame <- MakeInternalDataFrame(.local$items, add.rows=EXTRA_ROW)  

  .local$items <- NULL
      
  .local$theStack	<-	list()        
  .local$selectedColumns	<- integer(0)                             
            
  .local$gtktab <- MakeDFEditWindow(.local, .local$theFrame, size.request, col.width) 
  
  group.main	<-	gtkVBoxNew()
  .local$group.main <- group.main
  group.main$packStart(.local$gtktab, TRUE, TRUE, 0)
  group.main$setData(".local", .local) # use getData
  class(group.main) <- c("GtkDfEdit", class(group.main))

  gSignalConnect(group.main, "map", function(...)
    DoUserCall("OnLoad", list(sourc=dataset.name, typ="rdata"), .local)
  )
#  gSignalConnect(group.main, "destroy", function(obj, ...){
#    print("Destroyed")
#    #browser()
#    .local <- obj$getData(".local")
#    .local$model$setFrame()
#    rm(list=ls(env=.local), envir=.local)
#    obj$setData(".local", NULL)
##    print(gc()[2,2])
#  })

  return(group.main)
}

#' get selected row and column indices
#'
#' @param object The RGtk2DfEdit object
#' @return the 1-indexed selected rows
#' @export
gtkDfEditGetSelection <- function(object){
  .local <- object$getData(".local")
  dmm <- object$getDimension()
  columns=GetSelectedColumns(.local)
  rows <- GetSelectedRows(.local$view, .local)
  if(!length(rows) && !length(columns)){
    rows <- GetSelectedRows(.local$view.rn, .local)
    if(length(rows))
      if(dmm[2]) columns <- 1:dmm[2]
  } else if (!length(rows) && length(columns)){
    if(dmm[1]) rows <- 1:dmm[1]
  }
	return(list(rows=rows, columns=columns))
}


## JV ADD METHODS
## Should be able to call these via RGtk2 dispatch: obj$getModel() instead of gtkDfEditGetModel(obj)
#' get Model from object
#'
#' @param object The RGtk2DfEdit object
#' @return the RGtk2DataFrame that is the backend model for the widget
#' @export
gtkDfEditGetModel <- function(object) {
  object$getData(".local")$model
}

#' Return the dimensions (nrow, ncol) of the RGtk2DfEdit object
#'
#' @param object The RGtk2DfEdit object
#' @return Returns the number of rows and columns -- not counting row names
#' @export
gtkDfEditGetDimension <- function(object) {
  dim(object$getModel()) - c(EXTRA_ROW,2)       # rownames
}

#' Return the columns of the RGtk2DfEdit object
#'
#' @param object The RGtk2DfEdit object
#' @return Returns the column names for the current object
#' @export
gtkDfEditGetColumnNames <- function(object) {
  model <- object$getModel()
  nms <- colnames(model)[-1]
  return(nms)
}

#' Return the row names of the RGtk2DfEdit object
#'
#' @param object The RGtk2DfEdit object
#' @return Returns the row names for the current object
#' @export
gtkDfEditGetRowNames <- function(object) {
  model <- object$getModel()
  nms <- model[,1, drop=TRUE]
 return(nms)
}

#' Return a data frame from the RGtk2DfEdit object
#'

#' @param object The RGtk2DfEdit object
#' @return Returns the data frame with row names and column names
#' @export
gtkDfEditGetDataFrame <- function(object) {
  .local <- object$getData(".local")
  #dimnames(model) <- list(object$getRowNames(), object$getColumnNames())
  return(MakeExternalDataFrame(object$getModel(), .local))
}

#' S3 data extraction method
#'
#' Grabs data frame then passes onto [.data.frame method
#' @method [ RGtkDfEdit
#' @param x The RGtk2DfEdit object
#' @param i Row index
#' @param j Column indext
#' @param drop passed to extraction for data frame
#' @return The extracted entries
#' @export
"[.GtkDfEdit" <- function(x, i, j, drop = TRUE) {
   if(missing(drop))
     if (missing(i)) {
       drop <- TRUE
     } else {
       drop <- length(x$getDimension()[2]) == 1
     }
   df <- x$getDataFrame()
   df[i,j, drop=drop]
}

gtkDfEditSetActionHandler <- function(object, func.name, handler=NULL, data=NULL) {
  l <- object$getData(".local")
  if(!length(l$changed.handler)) l$changed.handler <- list()
  l$changed.handler[[func.name]]$func <- handler
  l$changed.handler[[func.name]]$data <- data
  object$setData(".local", l)
  invisible()
}

# needs missing i, j methods
"[<-.GtkDfEdit" <- function(x, i, j, value) {   
   l <- x$getData(".local")
   task <- list(
     list(func="ChangeCells", 
      arg = list(nf=value, row.idx=i, col.idx=ToInternalColIdx(j), do.coercion=T)))
   DoTaskWrapper(l, task)
   x
}

gtkDfEditDoTask <- function(x, task){
  # Correct col.idx to internal
  if(length(task)){
		for(jj in 1:length(task)){
		  stopifnot(all(c("func", "arg")%in%names(task[[jj]])))
		  if("col.idx"%in%names(task[[jj]]$arg)) 
		    task[[jj]]$arg$col.idx <- ToInternalColIdx(task[[jj]]$arg$col.idx)
      if(task[[jj]]$func=="InsertRows") 
        task[[jj]]$arg$nf <- MakeInternalDataFrame(task[[jj]]$arg$nf, add.rows=F)
      else if(task[[jj]]$func=="InsertColumns") 
        task[[jj]]$arg$nf <- MakeInternalDataFrame(task[[jj]]$arg$nf, add.columns=F)
		}
		DoTaskWrapper(x$getData(".local"), task)
  }
}


gtkDfEditUndo <- function(x){
  DoUndo(x$getData(".local"))
}    


gtkDfEditSort <- function(x){
  l <- x$getData(".local")
  dd <- dim(l$theFrame)
  DoSortDialog(l$theFrame[, -dd[2],drop=F], SortHandler, l)
}

gtkDfEditEditFactors <- function(x){
  .local <- x$getData(".local")
  DoFactorEditor(.local$theFrame, .local$toplevel, integer(0), 
  FactorEditorHandler, data=.local)
}

gtkDfEditCommandData <- function(x){
  .local <- x$getData(".local")
  sc <- GetSelectedColumns(.local)
  if(length(sr <- GetSelectedRows(.local$view, .local))==0) 
    sr <- 1:(dim(.local$theFrame)[1]-1)
  CommandData(.local, sr, sc)
}

gtkDfEditGetDatasetName <- function(x){
  x$getData(".local")$dataset.name
}

gtkDfEditSetDatasetName <- function(x, new.name){
  ModelChangeDatasetName(x$getData(".local"), new.name)
}
