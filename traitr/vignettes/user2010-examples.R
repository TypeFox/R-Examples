##' These are the examples for the useR2010! talk on traitr

require(gWidgets)
options(guiToolkit="RGtk2")
require(traitr)


##################################################
### basic t-test summary dialog example
##################################################

ttest <- function(xbar, s, n, mu, alternative = c("two.sided", "less", "greater")) {
  tstat <- (xbar - mu)/(s/sqrt(n))
  switch(match.arg(alternative),
         "less"=pt(tstat, df=n-1),
         "greater"=pt(tstat, df=n-1, lower.tail=FALSE),
         2*pt(abs(tstat), df=n-1, lower.tail=FALSE))
}


example1 <- function(whichOne) {

  dlg <- aDialog(items=list(
                   xbar=numericItem(0),
                   s=numericItem(1),
                   n=integerItem(2),      # 2 or more
                   mu=numericItem(0),
                   alternative=choiceItem(value="two.sided", values=c("two.sided", "less", "greater"))
                   ),
                 title="T-test p-value from summary",
                 help_string="Adjust values then click 'OK'",
                 
                 OK_handler=function(.) print(do.call(ttest, .$to_R()))
                 )
  

  if(whichOne == 1) {
    dlg$make_gui()
  } else if(whichOne == 2) {
    ## validation
    dlg1 <- dlg$instance()
    tmp <- dlg1$get_item_by_name("s")
    tmp$validate <- function(., rawvalue) if(as.numeric(rawvalue) > 0) return(rawvalue) else stop("s must be positive")
    
    tmp <- dlg1$get_item_by_name("n")
    tmp$validate <- function(., rawvalue) if(as.integer(rawvalue) > 1) return(rawvalue) else stop("n must be 2 or more")
    
    dlg1$make_gui()
  } else {
    ## layout
    dlg2 <- dlg$instance()
    dlg2$set_xbar(NA)
    view <- aContainer(aFrame(label="Data:",
                              aContainer("xbar","s","n")),
                       aFrame(label="Hypotheses:",
                              enabled_when=function(.) !is.na(.$get_xbar()), ## dlg2$get_xbar() == "NA"
                              aContainer("mu","alternative"))
                       )
    dlg2$make_gui(gui_layout=view)
    
  }

}

example1RGtk2 <- function() {
  ## do the same in RGtk2
  library(RGtk2)
  w <- gtkWindow()
  w$setTitle("T-test summary")
  
  ## make widgets
  widgets <- list()
  widgets$xbar <- gtkEntryNew(); widgets$xbar$setText(0)
  widgets$s <- gtkEntryNew(); widgets$s$setText(1)
  widgets$n <- gtkEntryNew(); widgets$n$setText(2)
  widgets$mu <- gtkEntryNew(); widgets$mu$setText(0)
  widgets$alternative <- gtkRadioButton(group=NULL, label="two.sided")
  QT <- sapply(c("less","greater"), function(i) widgets$alternative$NewWithLabelFromWidget(i))
  alternativeGp <- gtkVBox()
  QT <- sapply(rev(widgets$alternative$getGroup()),
               function(i) alternativeGp$PackStart(i))
  
  ## add widget
  tbl <- gtkTableNew()
  w$Add(tbl)
  
  nms <- names(widgets)
  QT <- sapply(1:4, function(i) {
    tbl$attach(gtkLabel(nms[i]), 0, 1, i-1, i)
    tbl$attach(widgets[[i]], 1, 2, i-1, i)
  })
  tbl$attach((l <- gtkLabel("alternative")), 0, 1, 4, 5)
  l['yalign'] <- 0                        # align to top
  tbl$attach(alternativeGp, 1, 2, 4, 5)
  
  bg <- gtkHBox()
  tbl$attach(bg, 0, 2, 5,6)
  OK <- gtkButton("OK")
  QT <- gSignalConnect(OK, "clicked", function(w,...) {
    altVals <- sapply(widgets$alternative$getGroup(), function(i) i['label'])
    ind <- sapply(widgets$alternative$getGroup(), function(i) i$getActive())
    val <- ttest(as.numeric(widgets$xbar$getText()),
                 as.numeric(widgets$s$getText()),
                 as.numeric(widgets$n$getText()),
                 as.numeric(widgets$mu$getText()),
                 alternative=altVals[ind])
    print(val)
  })
  bg$PackStart(OK)
  
  Cancel <- gtkButton("Cancel")
  QT <- gSignalConnect(Cancel, "clicked", function(...) w$destroy())
  bg$PackStart(Cancel)
}

###################################################
###  tcltk-demo
###################################################

makePlot <-  function(dist, kernel, n, bw,...) {
  y <- switch(dist, "Normal"=rnorm(n), rexp(n))
  plot(density(y,  bw=bw, kernel=kernel), xlim=range(y)+c(-2,2), main="Density example")
  points(y, rep(0,n))
}

exampleTcltkDemo <- function() {
  
  modelItems <- list(
                     dist=choiceItem("Normal", values=c("Normal","Exponential"),
                       show_label=FALSE),
                     kernel=choiceItem("gaussian",
                       values=c("gaussian", "epanechnikov", "rectangular",
                         "triangular", "cosine"),
                   show_label=FALSE),
                 n=choiceItem(50L, as.integer(c(50,100,200,300)),
                   show_label=FALSE),
                 bw=rangeItem(value=1, from=0.05, to=2.00, by=0.05,
                   show_label=FALSE)
                   )
## modelItems a list of items already defined
modelItems$out <- graphicDeviceItem() ## New item type
dlg <- aDialog(
  items= modelItems, ## also dist, kernel, n, bw
  help_string="Adjust a parameter to update graphic",
  title="tkdensity through traitr",
  buttons="Cancel",
  make_plot=function(.) {
    do.call(makePlot, .$to_R())
  },
  model_value_changed=function(.) .$make_plot()
               )
## view
view <- aGroup(
  aContainer(aFrame("dist", label="Distribution"),
   aFrame("kernel", label="Kernel"),
   aFrame("n", label="Sample size"),
   aFrame("bw", label="Bandwidth")
  ),
  "out",
  horizontal=TRUE)

dlg$make_gui(gui_layout=view)
dlg$model_value_changed()

}
###################################################
### Filter example
###################################################

exampleFilter <- function() {
wt <- mtcars$wt
cyls <- sort(unique(mtcars$cyl))

do_find_ind <- function(., value, old_value) {
#  wt <- mtcars$wt
#  cyls <- sort(unique(mtcars$cyl))
  ind <-  (wt <= .$get_wt()) &  mtcars$cyl %in% .$get_cyl()
  .$set_tbl(.$data[ind,])
}

dlg <- aDialog(items=list(
                 tbl=tableItem(mtcars, attr=list(size=c(300,200))),
                 wt=rangeItem(max(wt), from=min(wt), to=max(wt), by=.1, label="Weight <=",
                   tooltip="Slide to adjust maximum weight for data"),
                 cyl=choiceItem(cyls, values=cyls, multiple=TRUE,
                   label="Cyl==",
                   tooltip="Restrict number of cylinders in data set")
                 ),
               data=mtcars,
               status_text=sprintf("%s cases",nrow(mtcars)),
               property_wt_value_changed=do_find_ind,
               property_cyl_value_changed=do_find_ind,
               property_tbl_value_changed=function(., value, old_value) {
                 .$set_status_text(sprintf("%s cases", nrow(value)))
               },
               title="Filter example"
               )
view <- aContainer("tbl",
                   aFrame(label="Filter:",
                          "wt", "cyl"))
dlg$make_gui(gui_layout=view)
}

##################################################
### How to select an example
##################################################

dlg <- aDialog(items=list(
                 e0=labelItem("Run an example"),
                 e1=buttonItem("RGtk2 ttest example", label="", action=function(...) {
                   example1RGtk2()
                 }),
                 e2 = buttonItem("traitr ttest example", label="", action=function(...) {
                   example1(1)
                 }),
                 e3 = buttonItem("traitr ttest example, validation", label="", action=function(...) {
                   example1(2)
                 }),
                 e4 = buttonItem("traitr ttest example, layout", label="", action=function(...) {
                   example1(3)
                 }),
                 e5 = buttonItem("traitr tkdensity example", label="", action=function(...) {
                   exampleTcltkDemo()
                 }),
                 e6 = buttonItem("traitr filter example", label="", action=function(...) {
                   exampleFilter()
                 })
                 ),
               buttons="Cancel",
               title="Run traitr examples"
               )

lyt = aContainer("e0",
  aFrame(label="Examples:", "e1", "e2", "e3", "e4", "e5", "e6"))
  
dlg$make_gui(gui_layout=lyt)
