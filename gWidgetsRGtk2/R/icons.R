## add to stock icons

## function to look up stock icons
## ie. ok returns "gtk-ok"

##stockIcons <- list();
stockIcons <- new.env()
updateStockIcons <- new.env(); updateStockIcons[['value']] <- TRUE
#assignInNamespace("stockIcons",list(), ns = "gWidgetsRGtk2")
#assignInNamespace("updateStockIcons",TRUE, ns = "gWidgetsRGtk2")

loadGWidgetIcons = function() {
  ## add the icons
  ## we use xpm icons gimp can convert
  iconFullNames = list.files(system.file("images", package="gWidgetsRGtk2"))
  iconFullNames = iconFullNames[grep("\\.xpm$",iconFullNames)] ## just xpm
  iconNames = gsub("\\.xpm$","",iconFullNames)
  ## Loop over all to add here
  iconFullNames = paste(iconNames,".xpm", sep="")
  iconFiles = sapply(iconFullNames, function(name) {
    system.file("images",name, package="gWidgetsRGtk2")
  })
  
  addToGtkStockIcons(iconNames, iconFiles)
}

## add stock icons from files
setMethod(".addStockIcons",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit, iconNames, iconFiles, ...) {
            updateStockIcons[['value']] <- TRUE
#            assignInNamespace("updateStockIcons",TRUE, ns = "gWidgetsRGtk2")
            addToGtkStockIcons(iconNames, iconFiles)
          })

addToGtkStockIcons = function(iconNames, iconFiles) {

  iconfactory = gtkIconFactoryNew()
  for(i in seq_along(iconNames)) {
    iconsource = gtkIconSourceNew()
    iconsource$SetFilename(iconFiles[i])
    
    iconset = gtkIconSetNew()
    iconset$AddSource(iconsource)
    
    stockName = paste("gWidgetsRGtk2-",iconNames[i],sep="")
    
    iconfactory$Add(stockName, iconset)
    
    items = list(test=list(stockName, iconNames[i],"","",""))
    gtkStockAdd(items)
  }
  
  iconfactory$AddDefault()
  invisible(TRUE)
}

## find the stock icons. This includes those added bia loadGWidgetIcons()
setMethod(".getStockIcons",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit) {
##           if(getFromNamespace("updateStockIcons", ns = "gWidgetsRGtk2")) {
            if(updateStockIcons[['value']]) {
             ## create icon list
             .stockicons <- list()
             for(i in unlist(gtkStockListIds())) {
               name <- sub("[a-zA-Z0-9]*-","",i)
               .stockicons[[name]] = i
             }
             stockIcons[["value"]] <- .stockicons
             updateStockIcons[["value"]] <- FALSE
#             assignInNamespace("stockIcons", .stockicons, ns = "gWidgetsRGtk2")
#             assignInNamespace("updateStockIcons",FALSE, ns = "gWidgetsRGtk2")
           }
           stockIcons[["value"]]
#           return(getFromNamespace("stockIcons", ns = "gWidgetsRGtk2"))
         })
                

## name can be a vector
## return NA, if not there
getstockiconname <- function(name=NULL) {
  .stockicons = getStockIcons(toolkit=guiToolkit("RGtk2"))         # cache?

  if(is.null(name))
    return(unlist(.stockicons))

  if(length(name) == 0)
    return(character(0))

  tmpfun = function(names) {
    sapply(names, function(name) {
      ## already a stock name?
      if(name %in% .stockicons)
        return(name)

      if(name %in% names(.stockicons)) {
        return(.stockicons[[name]])
      } else {
        return(NA)
      }
    })
  }
  return(tmpfun(name))
} 


#################################################
## functions to deal with icons
## class to icon translation -- return stock name
## with prefix
## find the stock icons. This includes those added bia loadGWidgetIcons()


setMethod(".stockIconFromClass",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,theClass, ...) {
            default = "symbol_star"
            
            if(is.null(theClass) ||
               is.na(theClass) ||
               length(theClass) == 0
               )
              return(NA)

            theClass = theClass[1]
            
            if(theClass %in% .models)
              return(getstockiconname("lines"))
            if(theClass %in% .ts)
              return(getstockiconname("ts"))
            if(theClass %in% .functions)
              return(getstockiconname("function"))
            
            ret = switch(theClass,
              "numeric"= "numeric",
              "integer"= "numeric",
              "logical" = "logical",
              "character"="select-font",
              "matrix" = "matrix",
              "data.frame" = "dataframe",
              "list" = "dataframe",
              "complex"="numeric",
              "factor"="factor",
              "recordedplot" = "plot",
              NA)
  
            return(getstockiconname(ret))
          })

setMethod(".stockIconFromObject",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,obj, ...) {
            .stockIconFromClass(class(obj)[1])
          })


##
## 
##loadGWidgetIcons()

