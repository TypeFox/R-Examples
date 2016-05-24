### R code from vignette source 'ex-tcltk-table.Rnw'

###################################################
### code chunk number 1: ex-tcltk-table
###################################################
library(tcltk)

## helpers
quoteIt <- function(string) {           
  doQuote <- function(x) {
    xx <- strsplit(x, '"', fixe = TRUE)[[1]]
    paste(paste('"', xx, '"', sep = ""), collapse = '\'"\'')
  }
  if(!length(string)) return("")
  has_double_quote <- grep('"',string)
  if(!length(has_double_quote))
    return(paste('"',string,'"',sep = ""))
  if (!length(grep("([$`])", string))) {
    paste("\"", gsub("([\"!\\])", "\\\\\\1", string), 
          "\"", sep = "")
  } else sapply(string, doQuote)
}


## covert a data frame into a character based on
.toCharacter <- function(x,width,...) UseMethod(".toCharacter")
.toCharacter.default <- function(x,width,...) as.character(x)
.toCharacter.integer <- function(x,width,...) {
 if(missing(width)) width <- max(nchar(as.character(x))) + 2  
  format(x, justify = "right", width = width)
}
.toCharacter.numeric <- function(x,width,...) {
  if(missing(width)) width <- max(nchar(as.character(x))) + 2
  format(x,trim = FALSE, width = width, justify = "right")
}
.toCharacter.factor <- function(x,width,...) {
  if(missing(width)) width <- max(nchar(as.character(x))) + 2
  .toCharacter(as.character(x),width,...)
}
.toCharacter.logical <- function(x,width,...) {
  if(missing(width)) width <- 7
  format(as.character(x), justify = "centre", width = width)
}
.toCharacter.data.frame <- function(x,width =  10, ...) {
  nms <- dimnames(x)
  DF <- as.data.frame(lapply(x,function(i) .toCharacter(i, width = width)),
                      stringsAsFactors = FALSE)
  dimnames(DF) <- nms
  return(DF)
}

addScrollbars <- function(parent, widget) {
  xscr <- ttkscrollbar(parent, orient = "horizontal",
                       command = function(...) tkxview(widget, ...))
  yscr <- ttkscrollbar(parent, orient = "vertical",
                       command = function(...) tkyview(widget, ...))

  tkconfigure(widget,
              xscrollcommand = function(...) tkset(xscr,...),
              yscrollcommand = function(...) tkset(yscr,...))

  tkgrid(widget, row = 0, column = 0, sticky = "news")
  tkgrid(yscr,row = 0,column = 1, sticky = "ns")
  tkgrid(xscr, row = 1, column = 0, sticky = "ew")
  tkgrid.columnconfigure(parent, 0, weight = 1)
  tkgrid.rowconfigure(parent, 0, weight = 1)
}


###################################################
### code chunk number 2: ex-tcltk-table.Rnw:67-68
###################################################
DF <- getCRANmirrors()[, c(1,2,5,4)]


###################################################
### code chunk number 3: notShown
###################################################
window <- tktoplevel()
tkwm.title(window, "Choose a CRAN mirror")
frame <- ttkframe(window, padding = c(3,3,3,12))
tkpack(frame, expand = TRUE, fill = "both")


###################################################
### code chunk number 4: ex-tcltk-table.Rnw:79-85
###################################################
frame_0 <- ttkframe(frame); tkpack(frame_0, fill = "x")
label <- ttklabel(frame_0, text = "filter:")
tkpack(label, side = "left")
filter_var <- tclVar("")
filter_entry <- ttkentry(frame_0, textvariable = filter_var)
tkpack(filter_entry, side = "left")


###################################################
### code chunk number 5: makeTreeview
###################################################
frame_1 <- ttkframe(frame)
tkpack(frame_1, expand = TRUE, fill = "both")
treeview <- ttktreeview(frame_1, columns = 1:ncol(DF), 
                  displaycolumns = 1:(ncol(DF) - 1), 
                  show = "headings",     # not "tree" 
                  selectmode = "browse") # single selection
addScrollbars(frame_1, treeview)


###################################################
### code chunk number 6: configureColumns
###################################################
widths <- c(100, 75, 400)            # hard coded
nms <- names(DF)
for(i in 1:3) {
  tcl(treeview, "heading", i, text = nms[i])
  tcl(treeview, "column", i, width = widths[i], 
      stretch = TRUE, anchor = "w")
}


###################################################
### code chunk number 7: ex-tcltk-table.Rnw:125-136
###################################################
fillTable <- function(treeview, DF) {
  children <- as.character(tcl(treeview, "children", ""))
  for(i in children) 
    tcl(treeview, "delete", i)                 # out with old
  shade <- c("none", "gray")
  for(i in seq_len(nrow(DF))) 
    tcl(treeview, "insert", "", "end", tag = shade[i %% 2], 
        text = "",  
        values = unlist(DF[i,]))               # in with new
  tktag.configure(treeview, "gray", background = "gray95")
}


###################################################
### code chunk number 8: ex-tcltk-table.Rnw:140-141
###################################################
fillTable(treeview, DF)


###################################################
### code chunk number 9: ex-tcltk-table.Rnw:146-155
###################################################
cur_ind <- 1:nrow(DF)
tkbind(filter_entry, "<KeyRelease>", function(W, K) {
  val <- tclvalue(tkget(W))
  poss_vals <- apply(DF, 1, function(...) 
                    paste(..., collapse = " "))
  ind<- grep(val, poss_vals)
  if(length(ind) == 0) ind <- 1:nrow(DF)
  fillTable(treeview, DF[ind,])
})


###################################################
### code chunk number 10: ex-tcltk-table.Rnw:160-169
###################################################
tkbind(treeview, "<Double-Button-1>", function(W, x, y) {
  sel <- as.character(tcl(W, "identify", "row", x, y))
  vals <- tcl(W, "item", sel, "-values")
  URL <- as.character(vals)[4]          # not tclvalue
  repos <- getOption("repos")
  repos["CRAN"] <- gsub("/$", "", URL[1L])
  options(repos = repos)
  tkwm.withdraw(tkwinfo("toplevel", W))
})


