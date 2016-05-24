## ## a notebook for holding help pages
## setClass("gHelpRGtk",
##          contains="gComponentRGtk",
##          prototype=prototype(new("gComponentRGtk"))
##          )

## setMethod(".ghelp",
##           signature(toolkit="guiWidgetsToolkitRGtk2"),
##           function(toolkit,
##                    topic=NULL, package=NULL,
##                    container = NULL,
##                    ...) {                                # passed to gnotebook
##             force(toolkit)

##             ## check if newversion of R, if so, we con't do a thing but return a label
##             if(getRversion() >= "2.10.0" && getRversion() < "2.11.0") {
##               l <- .glabel(toolkit, "ghelp needs to be updated for your version of R. Sorry.", cont=container)
##               return(l)
##             }


            
##             group = ggroup(horizontal=FALSE, container = container, ...)
##             notebook = gnotebook(...)
##             add(group, notebook, expand=TRUE)

##             obj = new("gHelpRGtk", block=group, widget=notebook,
##               toolkit=toolkit)
##             ##  obj = list(ref=group, gnotebook = notebook, notebook = notebook$notebook)
##             ## class(obj) = c("gHelp",class(notebook))
  

##             if(!is.null(topic))
##               .add(obj,toolkit, value = list(topic=topic, package=package))

##             invisible(obj)

##           })

## ##################################################
## ## gHelp methods
## ## workhorse is add -- value is either
## ## just a topic (not a list), or a list with components topic, package
## setMethod(".add",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gHelpRGtk"),
##           function(obj, toolkit, value, ...) {
##             if(is.list(value)) {
##               topic = value$topic
##               package = value$package
##             } else if(length(grep(":",value)) > 0) { # "stats:::t.test" works here
##               tmp = unlist(strsplit(value, ":+"))
##               package = tmp[1]
##               topic = tmp[2]
##             } else {
##               topic = value
##               package = NULL
##             }
            
##             ## error check
##             if(!is.character(topic) || length(topic) > 1 || length(topic) == 0) {
##               warning("Adios, adding to ghelp needs a valid topic\n")
##               return()
##             }

##             if(getRversion() < "2.10.0") {            
##               ## if package is NULL, we find them
##               if(is.null(package)) {
##                 possiblePackages = getPossiblePackages(topic)
##                 if(length(possiblePackages) > 0) {
##                   package = possiblePackages
##                 } else {
##                   warning(Paste("Can't find a package containing ", topic,"\n"))
##                   return()
##                 }
##               }
##               ## add a page for each package
##               for(pkg in package) {
##                 helpPage = makeHelpPage(topic, pkg)
##                 tag(helpPage,"topic") <- topic
##                 tag(helpPage,"package") <- pkg
##               add(obj@widget, helpPage, label = Paste("Help on ",pkg,"::",topic))
##               }
##               return()
##             } else if(getRversion() >= "2.11.0") {
##               ## add a page for each package
##               l <- list(topic=topic)
##               if(!is.null(package))
##                 l$package <- package
##               out <- do.call("help", l)
##               pkgname <-  basename(dirname(dirname(out)))
##               temp <- tools::Rd2txt(utils:::.getHelpFile(out), out = tempfile("Rtxt"), package=pkgname)
##               x <- readLines(temp)
##               unlink(temp)
##               helpPage <- gtext(cont=obj@widget, label=topic)
##               ## add text to gtext widget
##               for(i in x) {
##                 if(grepl("^_\b",i)) {
##                   insert(helpPage, gsub("_\b","",i), font.attr=c(weight="bold"))
##                 } else {
##                   insert(helpPage, i)
##                 }
##               }
##             }
##             return()
##             })

## ## value returns the topic of the current page or the one give by index
## setMethod(".svalue",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gHelpRGtk"),
##           function(obj, toolkit, index=NULL, drop=NULL, ...) {
##             pageno = svalue(obj@widget)
##             widget = obj@widget[pageno]
##             topic = tag(widget,"topic")
##             package = tag(widget,"package")
##             return(list(topic=topic, package=package))
##           })

## setMethod(".length",
##           signature(toolkit="guiWidgetsToolkitRGtk2",x="gHelpRGtk"),
##           function(x, toolkit) {
##             length(x@widget)
##           })

## setMethod(".dispose",
##           signature(toolkit="guiWidgetsToolkitRGtk2",obj="gHelpRGtk"),
##           function(obj, toolkit, ...) {
##             dispose(obj@widget)
##           })

## ##################################################
## ## helpers

## ## Return gtext widget with help page
## makeHelpPage = function(topic, pkg) {
##   helpFile = help(topic, package=force(pkg), verbose=TRUE)
## #  helpFile = system.file("help",topic,package=pkg)
##   if(helpFile[1] != "") {
##     ## deal with windows issue here
##     helpFile = zip.file.extract(as.character(helpFile),"Rhelp.zip")
##     text = readLines(helpFile)
##     text = sapply(text, function(i) gsub("\\_\\\b","",i))
##     helpPage = gtext(text[1],font.attr=c(weight="bold"))
##     add(helpPage, text[2])
##     add(helpPage, text[3], font.attr=c(weight="bold",size="big",color="blue"))
##     add(helpPage, text[-(1:3)])
##     ## This gave troubles when there were more than a few pages open!
## ##     sapply(text[-(1:3)], function(x) {
## ##       if( length(grep("^\\w+:", x)) > 0) {
## ##         tmp = unlist(strsplit(x,":"))
## ##         add(helpPage,Paste(tmp[1],":"),font.attr=c("blue"), do.newline=FALSE)
## ##         add(helpPage,paste(tmp[-1], sep="", collapse=":"))
## ##       } else {
## ##         add(helpPage,x)
## ##       }
## ##     })
##   } else {
##     helpPage = gtext(paste("Page for ",topic," in package ",pkg," was not found.",collapse=" "))
##   }
##   return(helpPage)
## }

## getPossiblePackages = function(topic) {
##   possiblePackages = c()
##   ## find all packages
##   lib.loc <- .libPaths()
##   packages <- .packages(all.available = TRUE, lib.loc = lib.loc)
##   for (lib in lib.loc) {
##     for (pkg in packages) {
##       dir <- system.file(package = pkg, lib.loc = lib)
##       ## XXX This needs to be rewritten for R version 2.11.0 or later
##       l <- list(topic, dir, "AnIndex", "help")
##       path <- do.call("index.search",l)
##       if(path != "")
##         possiblePackages = c(possiblePackages, pkg)
##     }
##   }
  
##   if(length(possiblePackages) == 0) {
##     warning("Adios, can't find a package to match ",topic,"\n")
##     return()
##   }
##   return(possiblePackages)
## }


## ##################################################

## ## is this of class gHelp?
## is.ghelp = function(x) {
##   is(x,"gHelpRGtk")
## }

## ## return name of any gHelp instances in environment
## findHelpObjectName = function(envir=.GlobalEnv) {
##   x = ls(envir=envir)
##   x[sapply(1:length(x),function(i)
##            is.ghelp(get(x[i],envir=envir)))]
##   x[!sapply(1:length(x), function(i)
##             is.invalid(get(x[i],envir=envir)))]
## }

## ##################################################
## ## This just pops up a window to show the argument from a help page


## ## Hack to open up help page to the argument
## showHelpAtArgument = function(argument, topic, package=NULL,
##   width=600, height=250) {
##   if(missing(argument) || missing(topic))
##     return()

##   if(is.null(package)) {
##     possiblePackages = getPossiblePackages(topic)
##     if(length(possiblePackages) > 0) {
##       package = possiblePackages
##     } else {
##       warning(Paste("Can't find a package containing", topic,"\n"))
##       return()
##     }
##   }

##   ## the widget
##   win=gwindow(Paste("Help on argument: ",topic), visible=FALSE) # set to visible if one is found
##   size(win) <- c(width,height)
##   group = ggroup(horizontal=FALSE, container=win)
##   textwindow = gtext()
##   add(group, textwindow, expand=TRUE)

##   for(pkg in package) {
## ##    helpFile = system.file("help",topic,package=pkg)
##     helpFile = help(topic, package=force(pkg), verbose=TRUE)[1]
##     if(helpFile != "") {
##       text = readLines(helpFile)
##       text = sapply(text, function(i) gsub("\\_\\\b","",i))
##       argPosition = grep(Paste(argument,": "), text)
##       if(length(argPosition) == 0) {
##         next
##       } else {
##         argPosition = argPosition[1] - 1
##         ##Found one
##         visible(win) <- TRUE            # show window
##       }

##       add(textwindow,Paste("From package:",pkg), font.attr=c(weight="bold"))
##       ## add first line (it has a :)
##       add(textwindow,text[argPosition+1],font.attr=c(weight="bold",color="blue"))
##       ## add until a :
##       i = 2; n = length(text)
##       while(length(grep(":",text[argPosition+i])) == 0 &&
##             (argPosition + i) <= n
##             ) {
##         add(textwindow,text[argPosition+i],font.attr=c(weight="bold",color="blue"))
##         i = i + 1
##       }
##       add(textwindow,"\n")
##     }
##   }
##   ## close button
##   buttonGroup = ggroup(container=group)
##   addSpring(buttonGroup)
##   gbutton("cancel", container=buttonGroup,
##           handler = function(h,...) dispose(h$obj))
  
## }




## ##################################################
## ## build on ghelp widget to make a browser with search,
## ## simpler than old pmg.helpBrowser. Break that into components

## ## a notebook for holding help pages
## setClass("gHelpbrowserRGtk",
##          contains="gComponentRGtk",
##          prototype=prototype(new("gComponentRGtk"))
##          )


## setMethod(".ghelpbrowser",
##           signature(toolkit="guiWidgetsToolkitRGtk2"),
##           function(toolkit,
##                    title = "Help browser", maxTerms=100,
##                    width=550, height=600) {

##             force(toolkit)
            
##             win = gwindow("Help browser", v=T)
##             size(win) <-  c(width,height)

##             obj=new("gHelpbrowserRGtk",block=win,widget=win,toolkit=toolkit)
##             ## obj = list(ref=win)
##             ##class(obj) = c("gHelpBrowser","gWidget")

  
##             gp = ggroup(horizontal = FALSE, container = win, expand=TRUE)
##             toolbarGroup = ggroup(container = gp)
## ##  toolbar = list()
## ##  toolbar$quit$handler=function(h,...) dispose(win)
## ##  toolbar$quit$icon = "quit"
## ##  toolbar$examples$handler=function(h,...) {
## ##    ## run example of current topic
## ##    lst = svalue(help.notebook)
## ##    if(!is.null(lst$topic))
## ##      do.call("example",lst)
## ##  }
## ##  toolbar$examples$icon = "evaluate"
## ##  gtoolbar(toolbar, style="both-horiz", container = toolbarGroup, expand=TRUE)

## #  add(toolbarGroup, gbutton("Quit", handler = function(h,...) dispose(win)))
##             quitHandler = function(h,...) dispose(win)
##             quitButton = ggroup(container=toolbarGroup)
##             add(quitButton,gimage("quit",dirname="stock",handler=quitHandler))
##             add(quitButton, glabel("Quit",handler = quitHandler))
##             runExamples = function(h,...) {
##               lst = svalue(help.notebook)
##               if(!is.null(lst$topic))
##                 do.call("example",lst)
##             }
##             examplesButton = ggroup(container=toolbarGroup)
##             add(examplesButton,gimage("evaluate",dirname="stock",handler=runExamples))
##             add(examplesButton, glabel("run examples",handler = runExamples))
##             ##  add(toolbarGroup, gbutton("close"))
##             addSpring(toolbarGroup)
            
##             ## others?
##             searchOptionsList = list(
##               "Help on function:" = function(...) NULL,
##               "help.search: apropos"=function(...) searchResultsApropos(...),
##               "help.search: pattern"=function(...) searchResultsHelpSearch(...)
##               )
##             searchOptions = gdroplist(names(searchOptionsList),
##               container = toolbarGroup)
##             searchBox = gedit("", container = toolbarGroup)
            
##             ## search through packages
##             expgp = gexpandgroup("Browse package help pages:",container = gp)
##             packageNotebook = gnotebook()
##             size(packageNotebook) <- c(400,300)
##                                         #  size(packageNotebook) <- c(leftnotebookwidth,notebookheight)
##             add(expgp,packageNotebook, expand=TRUE)
##             addhandlerchanged(packageNotebook,function(h,...) {
##               dispose(h$obj, to.right=TRUE)
##             })                  # delete to right, when changed
##             visible(expgp) <- FALSE
            
##             allPackages = .packages(all=TRUE)
##             packageList = gtable(data.frame("Package names"=I(allPackages)))
##             add(packageNotebook, packageList,label="All packages")
            
##             addhandlerdoubleclick(packageList, handler = function(h,...) {
##               ## get contents, show with filter
##               package = svalue(h$obj)
##               contents = getContentsOfPackage(package)
              
##               page = ggroup(horizontal=FALSE)
              
##               ## objectList
##               if(ncol(contents) >=2)
##                 objectList <- gtable(contents, filter.column=2)
##               else
##                 objectList <- gtable(contents)
##               add(page, objectList, expand=TRUE)
##               ## add to packageNotebook
##               add(packageNotebook,page, label=Paste("Objects in ",package))
##               addhandlerdoubleclick(objectList,action=package,
##                                     handler=function(h,...) {
##                                       topic = svalue(h$obj)
##                                       package = h$action
##                                       svalue(statusBar) <- Paste("Getting help page for ",topic)
##                                       add(help.notebook,list(topic=topic, package=package))
##                                       svalue(statusBar)
##                                       svalue(nb) <- 1 # help page
##                                       visible(expgp) <- FALSE
##                                       return(FALSE)
##                                     })
##               return(FALSE)             # doubleclick return for no more propogation
##             })
            
## ##################################################
##             nb = gnotebook(tab.pos=3)
##             add(gp,nb, expand=TRUE)
            
##             help.notebook  = ghelp(tab.pos=1,closebuttons=TRUE)     # bottom tab
##             emptyDataFrame = data.frame(Title=c(""), Package=c(""),Descr=c(""))
##             for(j in 1:3) emptyDataFrame[,j] <- as.character(emptyDataFrame[,j])
##             search.results = gtable(emptyDataFrame, filter.column=2)
            
##             add(nb, help.notebook, label="Help pages")
##             add(nb, search.results, label="Search results")
##             svalue(nb) <-1                # help page first
            
##             statusBar = gstatusbar(container=gp)
##             svalue(statusBar) <- "Enter search term in box, click ENTER to begin"
##             ## actions
##             ## double click on search results
##             addhandlerdoubleclick(search.results,
##                                   handler = function(h,...) {
##                                     vals = svalue(search.results, drop=FALSE) # a data frame
##                                     topic = vals[,1,drop=TRUE]
##                                     package = vals[,2,drop=TRUE]
##                                     svalue(statusBar) <-
##                                       Paste("Getting help page for ",topic)
##                                     add(help.notebook, list(topic=topic, package=package))
##                                     svalue(statusBar) # pops statusbar?
##                                     ## swap tabs
##                                     svalue(nb) <- 1
##                                     return(FALSE) # no mas
##                                   })
##             ## make search resuslts -- return dataframe with title, package, description
##             ## as character vectors
##             searchResultsApropos = function(query) {
##               out = help.search(apropos=query, ignore.case = TRUE)
##               out = out$matches
##               if(nrow(out) > 0) {
##                 out = out[1:min(nrow(out),maxTerms),c(1,3,2), drop=FALSE]
##               } else {
##                 out = c("no matches","","")
##               }
##               colnames(out) = c("topic","Package","title")
##               out = as.data.frame(out)
##               for(j in 1:3) out[,j] <- as.character(out[,j]) # avoid factors
##               return(out)
##             }
##             searchResultsHelpSearch = function(query) {
##               out = help.search(pattern=query, ignore.case = TRUE)
##               out = out$matches
##               if(nrow(out) > 0) {
##                 out = out[1:min(nrow(out),maxTerms),c(1,3,2), drop=FALSE]
##               } else {
##                 out = c("no matches","","")
##               }
##               colnames(out) = c("topic","Package","title")
##               out = as.data.frame(out, stringsAsFactors = FALSE)    
              
##               return(out)
##             }
            
##             addhandlerchanged(searchBox, handler = function(h,...) {
##               searchType = svalue(searchOptions, index=TRUE)
##               svalue(statusBar) <- "Getting to work"
##               if(searchType == 1) {
##                 ## first one is show help page
##                 topic = svalue(h$obj)
##                 add(help.notebook,topic)
##               } else {
##                 df = searchOptionsList[[searchType]](svalue(h$obj))
##                 ## set value in widget
##                 search.results[,] <- df
##                 ## raise search box
##                 svalue(nb) <-2
##                 svalue(statusBar) <-"Double click line to show help page"
##               }
##               svalue(statusBar)                   # pops
##             })
            
            
##             return(obj)
##           })


## ##################################################
## ## these are from old version
## ## contents a matrix with entry, keywords, description and URL
## getContentsOfPackage = function(package=NULL) {
##   if(is.null(package)) {
##     warning("Empty package name")
##     return(NA)
##   }

##   if(getRversion() < "2.10.0") {
##     contents = read.dcf(system.file("CONTENTS",package=package))
    
##     return(data.frame(Entry=I(contents[,1]),Keywords=I(contents[,3]),
##                       Description=I(contents[,4])))
##   } else if(getRversion() >= "2.11.0") {
##     l <- list(package=package)
##     tmp <- do.call("help", l)
##     contents <- tmp$info[[2]]
##     ## need to strip of white space
##     contents <- contents[grepl("^[a-zA-Z]", contents)]
##     contents <- sapply(strsplit(contents,"\\s+"), function(i) i[1])
##     contents <- contents[!grepl("^==",contents)]
##     contents <- contents[!contents == ""]
    
##     return(data.frame(Entry=contents,  stringsAsFactors=FALSE))
##   }
## }

