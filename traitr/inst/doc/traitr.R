### R code from vignette source 'traitr.Rnw'

###################################################
### code chunk number 1: traitr.Rnw:70-165
###################################################
## function to use booktabs format for tables
## colTypes -- l r c p ... for latex

## eg. cat(booktabs(df, colTypes=c("l","l","p{.65\\textwidth}")))
##     cat(booktabs(df, caption="Plot command"))
booktabs <- function(x,caption="",subtop="", label="",
                     colTypes=NULL, ...) UseMethod("booktabs")

booktabs.matrix <- function(x,caption="",subtop="",label="",colTypes=NULL,...) {
  dnms <- dimnames(x)
  hasDimnames <- ifelse(!is.null(names(dnms)), TRUE, FALSE)

  hasRowNames <- TRUE
  numRowNames <- try(as.numeric(dnms[[1]]), silent=TRUE)
  if(!inherits(numRowNames, "try-error") &&
     all( numRowNames == 1:nrow(x))) {
    ## skip numbers
    hasRowNames <- FALSE
  }

  ## make character matrix using format
  m <- matrix(rep("", (nrow(x) + hasRowNames) * ncol(x)), nrow= nrow(x))

  ## colTypes -- vector of 
  if(is.null(colTypes)) {
    colTypes <- character(ncol(x))
    for(i in 1:ncol(x)) {
      colTypes[i] <- switch(class(x[,i])[1],
                            "numeric" = "r",
                            "logical" = "c",
                            "l")
    }
  }
  
  ## fill
  if(hasRowNames) {
    if(hasDimnames) {
      m[,1] <- paste("\\quad ", dnms[[1]], sep="")
    } else {
      m[,1] <- dnms[[1]]
    }
  }
  
  ## fill in column by column
  for(j in 1:ncol(x))
    m[, j + hasRowNames] <- switch(class(x[,j])[1],
                                   "numeric"=format(x[,j], digits=4),
                                   as.character(x[,j]))

  ## write out
  out <- paste("\\begin{table}",
               "\\centering",
               if(label != "") {
                 paste("\\label{", label, "}", sep="")
               },
               if(caption != "") {
                 paste("\\caption{", caption, "}", sep="")
               },
#               if(subtop != "") {
#                 paste("\\subtop{", subtop, "}", sep="")
#               },
               paste("\\begin{tabular}{",
                     "@{}",rep("l@{\\quad}",hasRowNames),
                     paste(colTypes,collapse=""),
                     "@{}}",
                     sep=""),
               "\\toprule",
               if(!is.null(names(dnms)[2])) 
                 paste("\\multicolumn{", ncol(x) + hasRowNames, "}{c}{", names(dnms)[2], "}\\\\", sep=""),
               ## header row
               if(hasRowNames)  {
                 paste(paste(c(names(dnms)[1], dnms[[2]]), collapse="&"), "\\\\", sep="")
               } else {
                 paste(paste(dnms[[2]], collapse="&"), "\\\\", sep="")
               },
               "\\midrule",
               paste(apply(m, 1, function(i)
                           paste(i, collapse="&")), collapse="\\\\"),
               "\\\\ \\bottomrule",
               "\\end{tabular}",
               "\\end{table}",
               "",
               sep="\n")

  out
}
booktabs.data.frame <- booktabs.matrix
  

booktabs.table <- function(x, caption="",subtop="",label="",
                           colTypes=NULL, ...) {
  y <- x[,]
  dimnames(y) <- dimnames(x)
  booktabs(y, caption, subtop, label, ...)
}


###################################################
### code chunk number 2: traitr.Rnw:172-175
###################################################
options(prompt=" ")
options(continue=" ")
options(width=80)


###################################################
### code chunk number 3: require (eval = FALSE)
###################################################
## require(traitr)
## require(gWidgets)
## options(guiToolkit="RGtk2")


###################################################
### code chunk number 4: dlg (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  number=numericItem(0),
##                  string=stringItem("")
##                  )
##                )
## dlg$make_gui()                          # method call using $ notation.


###################################################
### code chunk number 5: dlg2 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  number=numericItem(0),
##                  string=stringItem("")
##                  ),
##                OK_handler=function(.) { # . is reference to dlg object
##                  values <- .$to_R()
##                  f <- function(number,string) 
##                    cat("number is", number, "string is", string, "\n")
##                  do.call(f, values)
##                }
##                )
## dlg$make_gui()


###################################################
### code chunk number 6: traitr.Rnw:334-335 (eval = FALSE)
###################################################
## dlg$get_number()


###################################################
### code chunk number 7: traitr.Rnw:345-355 (eval = FALSE)
###################################################
## basic.t.test <- function(mean, mu, sd, alternative=c("two.sided","less","greater"),
##                          n, ...) {
##   alternative <- match.arg(alternative)
##   obs <- (mean - mu)/(sd/sqrt(n))
##   switch(alternative,
##          "greater" = 1 - pt(obs, df=n-1),
##          "less" = pt(obs, df=n-1),
##          2*(1 - pt(abs(obs), df=n-1))
##          )
## }


###################################################
### code chunk number 8: traitr.Rnw:359-377 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  mean=numericItem(0),
##                  sd=numericItem(1),
##                  n=numericItem(5),
##                  mu=numericItem(0),
##                  alternative=choiceItem("two.sided", values=c("two.sided","less","greater")),
##                  output=stringItem("", label="p.value")
##                  ),
##                title="A basic t-test interface",
##                help_string="Enter in summarized values then click OK to get the p-value",
##                OK_handler=function(.) {
##                  lst <- .$to_R()
##                  lst$output <- NULL     # not really needed, ignored by function 
##                  out <- do.call("basic.t.test", lst)
##                  .$set_output(out)
##                }
##                )
## dlg$make_gui()


###################################################
### code chunk number 9: traitr.Rnw:415-416 (eval = FALSE)
###################################################
## view <- aContainer("mean","sd","n","mu","alternative","output")


###################################################
### code chunk number 10: traitr.Rnw:427-430 (eval = FALSE)
###################################################
## view <- aContainer(aFrame("mean","sd","n", label="Statistics"),
##                    aFrame("mu", "alternative", label="Hypotheses"),
##                    "output")


###################################################
### code chunk number 11: traitr.Rnw:436-440 (eval = FALSE)
###################################################
## view <- aContainer(aFrame(aContainer("mean","sd","n"), label="Statistics"),
##                    aFrame(aContainer("mu","alternative"),label="Hypotheses"),
##                    separatorItem(),
##                    "output")


###################################################
### code chunk number 12: traitr.Rnw:446-448 (eval = FALSE)
###################################################
## dlg1 <- dlg$instance() ## instead of copying the definition above.
## dlg1$make_gui(gui_layout=view)


###################################################
### code chunk number 13: traitr.Rnw:471-477 (eval = FALSE)
###################################################
## positive_value <- function(., rawvalue) {
##   value <- as.numeric(rawvalue)
##   if(!value > 0)
##     stop("value is not positive")
##   value
## }


###################################################
### code chunk number 14: traitr.Rnw:491-493 (eval = FALSE)
###################################################
## sd <- numericItem(1, name="sd")
## sd$validate <- positive_value


###################################################
### code chunk number 15: traitr.Rnw:499-503 (eval = FALSE)
###################################################
## sd <- dlg$get_item_by_name("sd")         # lookup and return item by name
## sd$validate <- positive_value            # assigns method to item
## dlg1 <- dlg$instance()
## dlg1$make_gui(gui_layout=view)


###################################################
### code chunk number 16: traitr.Rnw:519-546 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  x = numericItem(NA, eval=TRUE),
##                  y = numericItem(NA, eval=TRUE),
##                  alternative=choiceItem("two.sided", 
##                    values=c("two.sided", "less", "greater")),
##                  mu = numericItem(0),
##                  paired=trueFalseItem(FALSE),
##                  var.equal=trueFalseItem(FALSE)
##                  ),
##                title="GUI with some parts disabled"
##                )
## view <- aContainer("x",
##                    aContext("y", context=dlg,
##                             enabled_when=function(.) { # y depends on x
##                               ## . here is the context value, not the container object
##                               val <- .$to_R()$x
##                               !is.null(val) && !is.na(val)  && (nchar(val) > 0)
##                             }),
##                    "alternative",
##                    "mu",
##                    aContainer("paired","var.equal", context=dlg,
##                               enabled_when=function(.) {
##                                 val <- .$to_R()$y
##                                 !is.null(val) && !is.na(val) &&  (nchar(.$get_y()) > 0)
##                               })
##                    )
## dlg$make_gui(gui_layout=view)


###################################################
### code chunk number 17: traitr.Rnw:572-588 (eval = FALSE)
###################################################
## drawGraph <- function(n,..) hist(rnorm(n))
## dlg <- aDialog(items=list(
##                  n=rangeItem(10, from=1, to=100, by=1),
##                  graph=graphicDeviceItem()
##                  ),
##                title="Draw a graph",
##                help_string="Adjust slider or click OK to produce a new graph",
##                model_value_changed=function(.) {
##                  l <- .$to_R()
##                  l$graph <- NULL        # not really necessary
##                  do.call("drawGraph", l)
##                },
##                OK_handler=function(.) do.call("drawGraph",.$to_R())
##                )
##                
## dlg$make_gui()


###################################################
### code chunk number 18: traitr.Rnw:635-659
###################################################
df <- rbind(
            c("stringItem","For holding strings"),
            c("numericItem","For numbers"),
            c("integerItem","For integers"),
            c("expressionItem","For R expressions"),
            c("trueFalseItem","For Boolean values"),
            c("choiceItem","For choosing one or more values from a list of possible values"),
            c("rangeItem","To select a value from a range of values"),
            c("buttonItem","For adding a button"),
            c("labelItem","For adding a label"),
            c("dateItem","For editing a calendar date."),
            c("separatorItem","To add a visual separator"),
            c("dataframeItem","To select a data frame"),
            c("variableSelectorItem","To select a variable from a data frame"),
            c("graphicDeviceItem","(RGtk2 only) To embed a graphic device"),
            c("formulaItem", "For formula specification (to be written)"),
            c("dfeditItem", "To edit a data set (to be written)"),
            c("itemList", "An item that stores a list of other items (or itemgroups)")
            )
colnames(df) <- c("Constructor","Description")
cat(booktabs(df,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Table of item constructors.",
             label="tab:item-constructors"))


###################################################
### code chunk number 19: traitr.Rnw:677-685 (eval = FALSE)
###################################################
## dfi <- dataframeItem(value=".GlobalEnv", name="dfi", 
##                      editor_style="compact") # alternative editor style
## dlg <- aDialog(items=list(
##                  dfi,
##                  vsi=variableSelectorItem("", multiple=FALSE, dataframeItem=dfi, 
##                    attr=list(size=c(200,200)))
##                  ))
## dlg$make_gui()


###################################################
### code chunk number 20: traitr.Rnw:717-744 (eval = FALSE)
###################################################
## hyps <- anItemGroup(items=list(
##                       mu=numericItem(0), 
##                       alternative=choiceItem("two.sided", c("two.sided","less","greater"))
##                       ),
##                     gui_layout=aFrame("mu","alternative", label="Hypotheses")
##                     )
## ttestDialog <- aDialog(items=list(
##                          x=numericItem(NA, eval=TRUE),
##                          y=numericItem(NA, eval=TRUE),
##                          hyps$instance()
##                          ),
##                        OK_handler=function(.) {
##                          do.call("t.test",.$to_R())
##                        }
##                        )
## wilcoxDialog <- aDialog(items=list(
##                           x=numericItem(NA, eval=TRUE),
##                           y=numericItem(NA, eval=TRUE),
##                           hyps$instance()
##                           ),
##                         OK_handler=function(.) {
##                           do.call("wilcox.test", .$to_R())
##                         }
##                         )
## 
## ttestDialog$make_gui()
## wilcoxDialog$make_gui()                 # shares alt info!


###################################################
### code chunk number 21: traitr.Rnw:784-790 (eval = FALSE)
###################################################
## i <- numericItem(0, name="item1")
## i$get_item1()
## i$set_item1(3)
## i$get_item1()
## try(i$set_item1("c(1,2,3)"))            # fails validation, still stored in model
## i$get_item1()


###################################################
### code chunk number 22: traitr.Rnw:794-798 (eval = FALSE)
###################################################
## i <- numericItem(0, name="item2", eval=TRUE)
## i$set_item2("c(1,2,3)")                 # now okay
## i$get_item2()
## i$to_R()                                # coerced


###################################################
### code chunk number 23: traitr.Rnw:803-810 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  x=numericItem(0),
##                  y=stringItem("a")
##                  ))
## dlg$get_x()
## dlg$set_y("some string")
## dlg$get_y()


###################################################
### code chunk number 24: traitr.Rnw:824-834 (eval = FALSE)
###################################################
## ig <- anItemGroup(items=list(
##                   x=numericItem(1),
##                   y=choiceItem("a", values=letters[1:5])
##                     )
##                   )
## ig$get_y()
## i <- ig$get_item_by_name("y")
## i$get_y()                               # same as above
## i$get_values()                          # get values
## i$set_values(letters)                   # to set values


###################################################
### code chunk number 25: traitr.Rnw:853-870 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  x=numericItem(0),
##                  y=numericItem(0)
##                  )
##                )
## dlg1 <- aDialog(items=list(
##                   a=numericItem(0)
##                   ),
##                 property_x_value_changed=function(., value, old_value) {
##                   .$set_a(.$get_a() + value) # add value to a (assumes numeric)
##                 }
##                 )
## 
## dlg$add_observer(dlg1)
## dlg1$get_a()
## dlg$set_x(10)
## dlg1$get_a()                            # updated by x


###################################################
### code chunk number 26: traitr.Rnw:877-886 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  f=fileItem(""),
##                  i=imageItem("",attr=list(size=c(480,480)))
##                  ),
##                property_f_value_changed=function(., value, old_value) {
##                  .$set_i(value)
##                },
##                buttons="Cancel")
## dlg$make_gui()


###################################################
### code chunk number 27: traitr.Rnw:901-908 (eval = FALSE)
###################################################
## i <- numericItem(0, name="x")
## j <- numericItem(1, name="x")
## j$get_x()
## j$set_model(i)
## j$get_x()
## i$set_x(10)
## j$get_x()


###################################################
### code chunk number 28: traitr.Rnw:957-967 (eval = FALSE)
###################################################
## mb_l <- list(File=list(
##                New=gaction("new", icon="new", handler=function(h,...) print("New")),
##                Quit=gaction("quit", icon="quit", handler=function(h,...) dlg$close_gui())
##                ))
## tb_l <- list(Quit=gaction("quit", icon="quit", handler=function(h,...) dlg$close_gui()))
## dlg <- aDialog(items=list(x=stringItem("some value")),
##                menu_list=mb_l,
##                toolbar_list=tb_l,
##                title="Dialog with menu and toolbar")
## dlg$make_gui()


###################################################
### code chunk number 29: traitr.Rnw:990-1000 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  x=numericItem(0)
##                  ),
##                title="Change x to be able to close",
##                buttons="Close",
##                Close_handler=function(.) {
##                  if(.$get_x() != 0)
##                    .$close_gui()
##                  })
## dlg$make_gui()


###################################################
### code chunk number 30: traitr.Rnw:1017-1027 (eval = FALSE)
###################################################
## fun <- function() {
##   a <- 1
##   dlg <- aDialog(items=list(),
##                  a=2,
##                  meth1=function(.) a)
##   dlg$meth2 <- function(.) a
##   dlg[['meth3']] <- function(.) a
##   c("1"=dlg$meth1(), "2"=dlg$meth2(), "3"=dlg$meth3())
## }
## fun()


###################################################
### code chunk number 31: traitr.Rnw:1073-1088
###################################################
df <- rbind(
            c("aContainer","Basic container, uses tabular layout"),
            c("aTableLayout","Tabular layout with more than 2 columns"),
            c("aGroup","Box container to pack in children left to right or top to bottom"),
            c("aFrame","Box container with decorative frame"),
            c("anExpandGroup","Box container with trigger to hide"),
            c("aPanedGroup","Two pane container"),
            c("aNotebook", "Notebook container"),
            c("aContext", "Provide context for an item or items")
            )
colnames(df) <- c("Constructor","Description")
cat(booktabs(df,
             colTypes=c("l","p{0.7\\textwidth}"),
             caption="Table of view constructors.",
             label="tab:view-constructors"))


###################################################
### code chunk number 32: traitr.Rnw:1171-1177 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(x=numericItem(1)))
## g <- aGroup()                           # define outside view to access later
## view <- aContainer("x",g)
## dlg$make_gui(gui_layout=view, visible=FALSE) # postpone showing, but create containers
## l <- glabel("Look ma, a gWidgets label", cont = g$container) # how to find container
## dlg$visible(TRUE)


###################################################
### code chunk number 33: traitr.Rnw:1195-1204 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(x=numericItem(0)))
## g <- aGroup(visible_when=function(.) FALSE) # suppress showing
## view <- aContainer("x", g)
## dlg$make_gui(gui_layout=view)
## ## now to add to the GUI at g:
## ig <- anItemGroup(list(y=stringItem("a string")))
## ig$make_gui(container=g)
## g$visible_when <- function(.) TRUE
## dlg$update_ui()                        


###################################################
### code chunk number 34: traitr.Rnw:1225-1233 (eval = FALSE)
###################################################
## m <- getCRANmirrors(all = FALSE, local.only = FALSE)[,c(1,2,4)]
## setCran <- function(.,...) {
##   URL <- .$get_cran()
##   repos <- getOption("repos")
##   repos["CRAN"] <- gsub("/$", "", URL[1L])
##   options(repos = repos)
##   .$close_gui()
## }


###################################################
### code chunk number 35: traitr.Rnw:1236-1248 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  cran=choiceItem(value=NA, values=m,
##                    show_label=FALSE,  # suppress label
##                      attr=list(chosencol=3, size=c(400,500)) #chosencol is URL one, not first
##                    )
##                  ),
##                title="Choose a CRAN Mirror",
##                help_string="Click a mirror, then OK, or double click on mirror",
##                OK_handler=setCran,                 # OK button click
##                property_cran_value_changed=setCran # double click
##                )
## dlg$make_gui()


###################################################
### code chunk number 36: replot (eval = FALSE)
###################################################
##   replot <- function(.) {
##     l <- .$to_R()
##     f <- function(dist, kernel, n, bw,...) {
##       y <- switch(dist, "Normal"=rnorm(n), rexp(n))
##       plot(density(y,  bw=bw, kernel=kernel), xlim=range(y)+c(-2,2), main="Density example")
##       points(y, rep(0,n))
##     }
##     do.call(f,l)
##   }
## 


###################################################
### code chunk number 37: traitr.Rnw:1283-1300 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  dist=choiceItem("Normal", values=c("Normal","Exponential"),
##                    show_label=FALSE),
##                  kernel=choiceItem("gaussian",
##                    values=c("gaussian", "epanechnikov", "rectangular",
##                      "triangular", "cosine"),
##                    show_label=FALSE),
##                  n=choiceItem(50L, as.integer(c(50,100,200,300)),
##                    show_label=FALSE),
##                  bw=rangeItem(value=1, from=0.05, to=2.00, by=0.05,
##                    show_label=FALSE),
##                  out=graphicDeviceItem()
##                  ),
##                help_string="Adjust a parameter to update graphic",
##                title="tkdensity through traitr",
##                buttons="Cancel",
##                model_value_changed=replot)


###################################################
### code chunk number 38: traitr.Rnw:1305-1314 (eval = FALSE)
###################################################
## view <- aGroup(aContainer(aFrame("dist", label="Distribution"),
##                           aFrame("kernel", label="Kernel"),
##                           aFrame("n", label="Sample size"),
##                           aFrame("bw", label="Bandwidth")
##                           ),
##                "out",
##                horizontal=TRUE)
## dlg$make_gui(gui_layout=view)
## replot(dlg)                           # initial plot


###################################################
### code chunk number 39: traitr.Rnw:1332-1337 (eval = FALSE)
###################################################
## i1 <- itemList(items=list(),
##                items_names="x",
##                item_factory=function(.) numericItem(0)
##                )
## i1$make_ui(container=gwindow("Basic Use"))


###################################################
### code chunk number 40: traitr.Rnw:1341-1343 (eval = FALSE)
###################################################
## ## add some items offline
## for(i in 1:2) i1$append_item(i1$item_factory())


###################################################
### code chunk number 41: traitr.Rnw:1345-1346 (eval = FALSE)
###################################################
## i1$to_R()


###################################################
### code chunk number 42: traitr.Rnw:1354-1361 (eval = FALSE)
###################################################
## i2 <- itemList(items=list(),
##                items_names="Data frame rows",
##                item_factory=function(.) {
##                  ig <- anItemGroup(items=list(a=numericItem(0), b=stringItem("")))
##                  ig$to_string <- function(.) .$get_b()
##                  ig
##                })


###################################################
### code chunk number 43: traitr.Rnw:1367-1380 (eval = FALSE)
###################################################
## i2$to_R <- function(.) {
##   items <- .$get_value()
##   if(length(items) == 0) {
##     out <- as.data.frame(.$item_factory()$to_R(), stringsAsFactors=FALSE)[0,]
##   } else {
##     out <- as.data.frame(items[[1]]$to_R(), stringsAsFactors=FALSE)
##     if(length(items) > 1) {
##       for(i in 2:length(items))
##         out[i,] <- items[[i]]$to_R()
##     }
##   }
##   out
## }


###################################################
### code chunk number 44: traitr.Rnw:1384-1385 (eval = FALSE)
###################################################
## i2$make_ui(container=gwindow("Basic Use to make data frame"))


###################################################
### code chunk number 45: traitr.Rnw:1415-1449 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  file=fileItem("", attr=list(
##                                      filter=list("CSV or TXT"=list(
##                                                    patterns=c("*.csv","*.txt")
##                                                    ),
##                                        "All files" = list(patterns=c("*"))
##                                        ))),
##                  header=trueFalseItem(TRUE, tooltip=paste("Variable onfirst line?")),
##                  sep=stringItem("", tooltip="Field separator character"),
##                  quote=stringItem("", tooltip="Set of quoting characters"),
##                  dec=stringItem(".", tooltip="Character used for decimal points"),
##                  as.is=trueFalseItem(!default.stringsAsFactors(), 
##                    tooltip="Do not convert character to factor"),
##                  na.strings=stringItem("NA", tooltip="Strings to be NA", eval=TRUE),
##                  nrows=integerItem(-1, tooltip="Max number of rows to read"),
##                  skip=integerItem(0, tooltip="Number of lines to skip at header"),
##                  check.names=trueFalseItem(TRUE, tooltip="If TRUE ensure names are valid"),
##                  fill=trueFalseItem(TRUE, tooltip="Fill unequal length rows if TRUE"),
##                  strip.white=trueFalseItem(TRUE),
##                  blank.lines.skip=trueFalseItem(TRUE, tooltip="If TRUE, skip blank lines"),
##                  comment.char=stringItem("#", tooltip="Comment character"),
##                  allowEscapes=trueFalseItem(TRUE, tooltip="C-style escapes read verbatim?"),
##                  stringsAsFactors=trueFalseItem(default.stringsAsFactors(), 
##                    tooltip="Characters converted to factors"),
##                  fileEncoding=stringItem(""),
##                  encoding=stringItem("unknown"),
##                  ## our things
##                  assign.to=stringItem("df", label="Assign to:"),
##                  output=tableItem(attr=list(size=c(400,400)), show_label=FALSE),
##                  file.type=stringItem("")
##                  ),
##                title="Read in txt or csv file",
##                help_string="Select a file, then adjust parameters."
##                )


###################################################
### code chunk number 46: traitr.Rnw:1461-1478 (eval = FALSE)
###################################################
## view <- aGroup(aNotebook(
##                          aNotebookPage("file",
##                                        separatorItem(),
##                                        "header", "sep","quote",
##                                        "dec", "fill", "comment.char",
##                                        label="Main args"),
##                          aNotebookPage("as.is","na.strings","nrows","skip",
##                                        "check.names","fill","strip.white","blank.lines.skip",
##                                        "allowEscapes","stringsAsFactors",
##                                        separatorItem(),
##                                        "fileEncoding", "encoding",
##                                        label="Extra args")
##                          ),
##                aContainer("assign.to",
##                           aFrame("output", label="Preview")
##                           ), 
##                horizontal=TRUE)


###################################################
### code chunk number 47: traitr.Rnw:1485-1497 (eval = FALSE)
###################################################
## dlg$read_file <- function(., file.type, output, assign.to, ...) {
##   if(file.type != "") {
##     out <- try(do.call(sprintf("read.%s",file.type), list(...)), silent=TRUE)
##     if(inherits(out, "try-error")) {
##       cat("Error reading file of type,", file.type, "\n")
##       out <- data.frame(V1="")
##     }
##   } else {
##     out <- data.frame(V1="")
##   }
##   return(out)
## }


###################################################
### code chunk number 48: traitr.Rnw:1503-1518 (eval = FALSE)
###################################################
## dlg$model_value_changed <- function(.) {
##   fname <- .$get_file()
##   if(file.exists(fname)) {
##     for(i in c("txt","csv")) {
##       if(grepl(paste("\\.",i,sep=""), fname))
##         .$set_file.type(c(txt="table",csv="csv")[i])
##     }
##   }
##   switch(.$get_file.type(),
##          "csv"={.$set_sep(","); .$set_quote('\"')},
##          "table"={},
##          {}
##          )
##   .$set_output(.$do_call("read_file",.$to_R()))
## }


###################################################
### code chunk number 49: traitr.Rnw:1534-1543 (eval = FALSE)
###################################################
## dlg$OK_handler <- function(.) {
##   out <- .$do_call("read_file",.$to_R())
##   assign.to <- .$get_assign.to()
##   if(exists(assign.to, envir=.GlobalEnv)) {
##     if(!gconfirm(sprintf("Overwrite variable %s?", assign.to)))
##       return()
##   }
##   assign(assign.to, out, envir=.GlobalEnv)
## }


###################################################
### code chunk number 50: traitr.Rnw:1547-1548 (eval = FALSE)
###################################################
## dlg$make_gui(gui_layout=view)


###################################################
### code chunk number 51: traitr.Rnw:1563-1577 (eval = FALSE)
###################################################
## dlg$property_file_value_changed <- function(., value, old_value) {
##   if(file.exists(value)) {
##     for(i in c("txt","csv")) {
##       if(grepl(paste("\\.",i,sep=""), value))
##         .$set_file.type(c(txt="table",csv="csv")[i])
##     }
##   }
##   switch(.$get_file.type(),
##          "csv"={.$set_sep(","); .$set_quote('\"')},
##          "table"={.$set_output(.$do_call("read_file",.$to_R()))},
##          {}
##          )
## 
## }


###################################################
### code chunk number 52: traitr.Rnw:1585-1594 (eval = FALSE)
###################################################
## nms <- names(dlg$get_items())
## nms <- setdiff(nms, c("file","file.type","assign.to","output"))
## QT <- sapply(nms, function(i) {
##   assign(sprintf("property_%s_value_changed",i),
##          function(., ...) {
##            .$set_output(.$do_call("read_file",.$to_R()))
##          },
##          envir=dlg)
## })


###################################################
### code chunk number 53: traitr.Rnw:1619-1628 (eval = FALSE)
###################################################
## m <- mtcars
## nms <- names(m)
## make_model <- function(i) {
##   l <- list(i = integerItem(i))
##   for(j in 1:ncol(m)) {
##     l[[nms[j]]] <- numericItem(m[i,j])
##   }
##   l
## }


###################################################
### code chunk number 54: traitr.Rnw:1634-1654 (eval = FALSE)
###################################################
## dlg <- aDialog(items=make_model(1),
##                title="Data frame scroller",
##                help_string="Press buttons to scroll through data set",
##                buttons=c("<< previous","next >>","SPACE","Cancel"),
##                set_row=function(.,i) {
##                  if(i < 1)
##                    i <- nrow(m)
##                  if(i > nrow(m))
##                    i <- 1
##                  ig <- make_model(as.numeric(i))
##                  .$set_model(anItemGroup(items=ig))
##                },
##                previous_handler=function(.) {
##                  i <- as.integer(.$get_i())
##                  .$set_row(i-1)
##                },
##                next_handler=function(.) {
##                  i <- .$to_R()$i        # same but different
##                  .$set_row(i + 1)
##                })


###################################################
### code chunk number 55: traitr.Rnw:1657-1658 (eval = FALSE)
###################################################
## dlg$make_gui()


###################################################
### code chunk number 56: traitr.Rnw:1671-1692 (eval = FALSE)
###################################################
## model <- aDialog(items=list(
##                    a=stringItem(""),
##                    b=stringItem("")
##                  )
##                  )
## 
## dlg1 <- aDialog(buttons="Next",
##                 Next_handler=function(.) {
##                   dlg2$make_gui(gui_layout=view2)
##                   .$close_gui()
##                 })
## view1 <- aContainer("a", context=model)
## 
## dlg2 <- aDialog(buttons = c("Finished"),
##                 Finished_handler = function(.) {
##                   print(model$to_R())
##                   .$close_gui()
##                 })
## view2 <- aContainer("b", context=model)
## 
## dlg1$make_gui(gui_layout=view1)


###################################################
### code chunk number 57: SpotfireExample (eval = FALSE)
###################################################
## theDesc <- paste("<b>Spotfire</b>",
##                  "The Spotfire web player (http://spotfire.tibco.com)",
##                  "has several demos built around a somewhat similar set-up:",
##                  "a description page, a data set, a graphic display of the data, and a set",
##                  "of controls to filter out the data that is being displayed in the graphic.",
##                  "",
##                  "This example shows how <b>traitr</b> can be used to make a simple version of such.",
##                  "",
##                  "Click the <i>Explore</i> tab to begin.",
##                  sep="\n")


###################################################
### code chunk number 58: dataDisplay (eval = FALSE)
###################################################
## theData <- mtcars
## makeLabel <- function(nr) sprintf("<b>%s</b> cases",nr)
## dataDisplay <- anItemGroup(items=list(
##                              data = tableItem(theData, name="data", attr=c(expand=TRUE)),
##                              label = labelItem(makeLabel(nrow(theData)), attr=c(markup=TRUE, expand=FALSE))
##                              ),
##                            attr=c(expand=TRUE)
##                            )


###################################################
### code chunk number 59: traitr.Rnw:1750-1753 (eval = FALSE)
###################################################
## ## synchronize labe with data dimension
## dataDisplay$property_data_value_changed <- function(., value, old_value)
##   .$set_label(makeLabel(nrow(.$get_data())))


###################################################
### code chunk number 60: traitr.Rnw:1760-1765 (eval = FALSE)
###################################################
## dataDisplay$make_default_gui_layout <- function(.) {
##   aGroup("data",
##          aGroup(labelItem("", attr=c(expand=TRUE)),"label", horizontal=TRUE),
##          horizontal=FALSE)
## }


###################################################
### code chunk number 61: filterBy (eval = FALSE)
###################################################
## var <- "cyl"
## varLevels <- sort(unique(theData[, var]))
## cylFilter <- anItemGroup(name=var,
##                          items=list(
##                            choice=choiceItem(varLevels, varLevels,
##                              by_index=FALSE, multiple=TRUE, show_label=FALSE)
##                            ),
##                          data=theData[, var],
##                          get_selected = function(.) {
##                            choice <- .$get_item_by_name("choice")
##                            value <- choice$get_choice()
##                            values <- choice$get_values()
##                            vals <- values[values %in% value]
##                            .$data %in% vals
##                          },
##                          make_default_gui_layout=function(.) {
##                            aFrame("choice", label="Cylinders")
##                         })


###################################################
### code chunk number 62: traitr.Rnw:1799-1814 (eval = FALSE)
###################################################
## var1 <- "wt"
## rng <- range(theData[, var1])
## wtFilter <- anItemGroup(name=var1,
##                         items=list(
##                           weight=rangeItem(value=rng[1] - .2, from=rng[1], to=rng[2], by=.2,
##                             show_label=FALSE, label="Weight >")
##                           ),
##                         data=theData[,var1],
##                         get_selected=function(.) {
##                           .$data >= .$to_R()$weight
##                         },
##                         make_default_gui_layout=function(.) {
##                           aFrame("weight", label="Weight > ")
##                         }
##                         )


###################################################
### code chunk number 63: traitr.Rnw:1819-1826 (eval = FALSE)
###################################################
## filters <- anItemGroup(items=list(cylFilter, wtFilter))
## filters$make_default_gui_layout <- function(.) {
##   aFrame(var,
##          var1,
##          label="Filter by",
##          attr=c(size=c(300,-1)))
## }


###################################################
### code chunk number 64: traitr.Rnw:1834-1839 (eval = FALSE)
###################################################
## filters$model_value_changed <- function(.) {
##   items <- .$get_items()
##   ind <- as.logical(items[[1]]$get_selected() * items[[2]]$get_selected())
##   dlg$update_data(ind)
## }


###################################################
### code chunk number 65: traitr.Rnw:1847-1856 (eval = FALSE)
###################################################
## dlg <- aDialog(items=list(
##                  Description=labelItem(theDesc, attr=c(markup=TRUE)),
##                  gd = graphicDeviceItem(),
##                  filters,
##                  dataDisplay
##                  ),
##                title="Spotfire example",
##                buttons="Cancel"
##                )


###################################################
### code chunk number 66: specifyDialog (eval = FALSE)
###################################################
## dlg$property_data_value_changed <- function(., value, old_value) 
##   .$draw_graphic()


###################################################
### code chunk number 67: traitr.Rnw:1866-1873 (eval = FALSE)
###################################################
## dlg$update_data <- function(., index) {
##   if(missing(index))
##     data <- theData
##   else
##     data <- theData[index,]
##   dataDisplay$set_data(data)
## }


###################################################
### code chunk number 68: traitr.Rnw:1879-1888 (eval = FALSE)
###################################################
## dlg$draw_graphic <- function(.) {
##   data <- dataDisplay$get_data()
##   if(nrow(data)) {
##     hist(data[,"mpg"])
##     rug(data[, "mpg"])
##   } else {
##     plot.new()
##   }
## }


###################################################
### code chunk number 69: traitr.Rnw:1894-1906 (eval = FALSE)
###################################################
## ## Layout for main GUI
## lyt <- aNotebook(aNotebookPage("Description",
##                                label="About"
##                                ),
##                  aNotebookPage(aPanedGroup(
##                                            aPanedGroup("gd",
##                                                        dataDisplay,
##                                                        horizontal=FALSE),
##                                            filters),
##                                label="Explore"
##                                )
##                  )


###################################################
### code chunk number 70: makeGUI (eval = FALSE)
###################################################
## w <- loadingAnimation()
## dlg$make_gui(gui_layout=lyt)
## w$close()


