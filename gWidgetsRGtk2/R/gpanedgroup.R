setClass("gPanedgroupRGtk",
         contains="gContainerRGtk",
         prototype=prototype(new("gContainerRGtk"))
         )

## TODO: method obj[1 or 2 ] <- replacewidget
setMethod(".gpanedgroup",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   widget1, widget2, horizontal=TRUE, container=NULL, ...) {
            ## add a paned group
            
            force(toolkit)
            
            if(horizontal) {
              panedWindow = gtkHPanedNew()
            } else {
              panedWindow = gtkVPanedNew()
            }

            obj <- as.gWidgetsRGtk2(panedWindow)
           
            if(!missing(widget1) && !is.null(widget1))
              add(obj, widget1)
            if(!missing(widget2) && !is.null(widget2))
              add(obj, widget2)
            
            if (!is.null(container)) {
              if(is.logical(container) && container == TRUE)
                container = gwindow(visible=TRUE)
              add(container, obj, ...)
            }
            
            return(obj)
          })

as.gWidgetsRGtk2.GtkHPaned <- as.gWidgetsRGtk2.GtkVPaned <-
  function(widget,...) {
    
    obj = new("gPanedgroupRGtk", block=widget, widget=widget,
      toolkit=guiToolkit("RGtk2"))
    
  if(is.null(tag(obj,"leftgroup"))) {
    ## left or right *or* top or bottom
    leftgroup = ggroup()
    rightgroup = ggroup()

    ## already a child?
    if(!is.null(child <- widget$GetChild1()))
      add(leftgroup,child,expand=TRUE, fill="both")
    if(!is.null(child <- widget$GetChild2()))
      add(rightgroup,child,expand=TRUE, fill="both")
    
    widget$Pack1(leftgroup@widget@block)#, resize=FALSE, shrink=FALSE)
    widget$Pack2(rightgroup@widget@block)#, resize=FALSE, shrink=FALSE)
    
    tag(obj,"leftgroup") <- leftgroup
    tag(obj,"rightgroup") <- rightgroup
  }
   
  return(obj)
     
}


## add -- use this rather than at construction time
setMethod(".add",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gPanedgroupRGtk", value="gWidgetRGtk"),
          function(obj, toolkit, value, ...) {
            ctr = tag(obj,"ctr")
            if(is.null(ctr))
              ctr = 0

            if(ctr == 0) {
              add(tag(obj,"leftgroup"), value, expand=TRUE, fill="both")
              ctr = 1
            } else if(ctr ==1) {
              add(tag(obj,"rightgroup"), value, expand=TRUE, fill="both")
              ctr = 2
            } else {
              gwCat(gettext("Can only add two widgets to a gpanedgroup\n"))
            }
            tag(obj,"ctr") <- ctr
            
          })


### methods
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gPanedgroupRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            panedWindow <- obj@widget
            min <- panedWindow['min-position']
            max <- panedWindow['max-position']
            position <- panedWindow['position']

            return((position - min)/(max - min))
            
          })

## svalue sets position
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gPanedgroupRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(0 <= value && value <= 1) {
                     panedWindow <- obj@widget
                     min <- panedWindow['min-position']
                     max <- panedWindow['max-position']
                     placement <- min + value * (max - min)
                     panedWindow['position'] <- placement
                   }
                   return(obj)
                 })
