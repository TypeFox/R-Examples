
##' Widget to provide interface to help system
setClass("gHelp",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor for help system widget
##'
##' @export
ghelp <- function(
                  topic = NULL, package = NULL, container = NULL, ... ,
                  toolkit=guiToolkit()){
  widget <- .ghelp (toolkit,
                    topic=topic, package=package, container=container ,...
                    )
  obj <- new( 'gHelp',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' API
##' add, character -- add page character=c(topic), c(topic,package) "package:::topic"
##' add, list -- add page, list(topic, package=NULL)
##' svalue: return list(topic=topic, package=package)
##' length: number of pages
##' dispose: remove current page


##' generic for toolkit dispatch
##' @alias ghelp
setGeneric( '.ghelp' ,
           function(toolkit,
                    topic = NULL, package = NULL, container = NULL, ... )
           standardGeneric( '.ghelp' ))


##################################################
## ANY imlementation
setClass("gHelpANY",
         contains="gComponentANY",
         prototype=prototype(new("gComponentANY"))
         )

setMethod(".ghelp",
          signature(toolkit="ANY"),
          function(toolkit,
                   topic=NULL, package=NULL,
                   container = NULL,
                   ...) {                                # passed to gnotebook
            force(toolkit)
            
            ## check if newversion of R, if so, we con't do a thing but return a label
            if(!getRversion() >= "2.11.0" && getRversion() < "2.11.0") {
              gwCat("Needs a new version of R to work.\n")
              glabel("ghelp", container=container)
            }


            nb <- gnotebook(container=container, closebuttons=TRUE, ...)
            
            obj <- new("gHelpANY", block=nb, widget=nb,
              toolkit=toolkit)

            if(!is.null(topic))
              .add(obj, toolkit, value = list(topic=topic, package=package))

            invisible(obj)

          })

##################################################
## gHelp methods
## workhorse is add -- value is either
## just a topic (not a list), or a list with components topic, package
setMethod(".add",
          signature(toolkit="ANY",obj="gHelpANY", value="character"),
          function(obj, toolkit, value, ...) {

            if(length(grep(":",value)) > 0) { # "stats:::t.test" works here
              tmp = unlist(strsplit(value, ":+"))
              package = tmp[1]
              topic = tmp[2]
            } else {
              topic = value
              package = NULL
            }
            .add(obj, toolkit, list(topic=topic, package=package))
          })


##' add a page to the help notebook
##' @param obj ghelp object
##' @param toolkit toolkit
##' @param value a list with component topic and value
setMethod(".add",
          signature(toolkit="ANY",obj="gHelpANY", value="list"),
          function(obj, toolkit, value, ...) {
            topic <- value$topic
            package <- value$package

            nb <- obj@widget

            ## error check
            if(!is.character(topic) || length(topic) > 1 || length(topic) == 0) {
              warning(sprintf("Adding to ghelp needs a valid topic. You tried %s.\n",topic))
              return(NULL)
            }

            l <- .findHelpPage(topic, package)

            x <- l$x
            topic <- l$topic; package <- l$package

            
            if(!is.null(x)) {
              ## are we already present?
              if(n <- .length(obj, toolkit)) {
                for(i in 1:n) {
                  l <- list(topic=tag(nb[i],"topic"), package=tag(nb[i],"package"))
                  if(l$topic == topic && (
                       (is.null(package) | is.null(l$package)) ||
                       l$package == package)) {
                    svalue(nb) <- i
                    return(NULL)
                  }
                }
              }

              ## good to go
              nb <- obj@widget
              t <- gtext(container=nb, label=topic, expand=TRUE)
              tag(t, "topic") <- topic
              tag(t, "pacakge") <- package
              svalue(nb) <- length(nb)
              .insertHelpPage(t, x)
            }
            
            return(NULL)
          })

##' returns list of topic and package of current page
setMethod(".svalue",
          signature(toolkit="ANY",obj="gHelpANY"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            nb <- obj@widget

            if(n <- length(nb) == 0)
              return(NULL)
            
            if(is.null(index))
              index <- svalue(nb)
            else
              index <- min(1, max(n, as.integer(index)))
            
            page <- nb[index]
            l <- list(topic=tag(page, "topic"),
                      package=tag(page, "package"))

            return(l)
          })

##' number of pages in notebook
setMethod(".length",
          signature(toolkit="ANY",x="gHelpANY"),
          function(x, toolkit) {
            length(x@widget)
          })

##' dispose of current page
setMethod(".dispose",
          signature(toolkit="ANY",obj="gHelpANY"),
          function(obj, toolkit, ...) {
            dispose(obj@widget)
          })

### Helper functions for "add"
##' return help page a set of lines
.findHelpPage <- function(topic, package=NULL) {
  l <- list(topic=topic)
  if(!is.null(package))
    l$package <- package
  out <- do.call("help", l)
  if(length(out) == 0) return(NULL)
  
  pkgname <-  basename(dirname(dirname(out)))

  ## thanks to Josef L for this
  help.txt <- "" ## keep R CMD check happy  
  help.con <- textConnection("help.txt", "w", local = TRUE)
  tools::Rd2txt(utils:::.getHelpFile(out), out=help.con, package=pkgname,
                width=80L)
  close(help.con)
  
  return(list(x=help.txt,topic=topic, package=pkgname))
}

##' insert help page into text object
##' makes bold if speedy enough
.insertHelpPage <- function(obj, x) {
  isSlow <- obj@toolkit@toolkit == "tcltk" || obj@toolkit@toolkit == "RGtk2"
  dispose(obj)       # clear
  
  out <- c()
  for(i in x) {
    if(grepl("^_\b",i)) {
      if(isSlow)
        out <- c(out, gsub("_\b","",i))
      else
        insert(obj, gsub("_\b","",i), font.attr=c(weight="bold"))
    } else {
      if(isSlow)
        out <- c(out,i)
      else
        insert(obj, i,font.attr=c(weight="normal"))
    }
  }
  if(isSlow)
    svalue(obj) <- out
  else
    insert(obj, "", do.newline=FALSE, where="beginning")              
}
            

##################################################
## helpers

getPossiblePackages = function(topic) {
  possiblePackages = c()
  ## find all packages
  lib.loc <- .libPaths()
  packages <- .packages(all.available = TRUE, lib.loc = lib.loc)
  for (lib in lib.loc) {
    for (pkg in packages) {
      dir <- system.file(package = pkg, lib.loc = lib)
      path = utils:::index.search(topic, dir, "AnIndex", "help")
      if(path != "")
        possiblePackages = c(possiblePackages, pkg)
    }
  }
  
  if(length(possiblePackages) == 0) {
    warning("Adios, can't find a package to match ",topic,"\n")
    return()
  }
  return(possiblePackages)
}


##################################################
## This just pops up a window to show the argument from a help page


## Hack to open up help page to the argument
showHelpAtArgument = function(argument, topic, package=NULL,
  width=600, height=250) {
  if(missing(argument) || missing(topic))
    return()

  if(is.null(package)) {
    possiblePackages = getPossiblePackages(topic)
    if(length(possiblePackages) > 0) {
      package = possiblePackages
    } else {
      warning(Paste("Can't find a package containing", topic,"\n"))
      return()
    }
  }

  ## the widget
  win=gwindow(Paste("Help on argument: ",topic), visible=FALSE) # set to visible if one is found
  group = ggroup(horizontal=FALSE, container=win)
  textwindow = gtext("", container=group, expand=TRUE)
  size(textwindow) <- c(width,height)

  for(pkg in package) {
##    helpFile = system.file("help",topic,package=pkg)
    helpFile = help(topic, package=force(pkg), verbose=TRUE)[1]
    if(helpFile != "") {
      text = readLines(helpFile)
      text = sapply(text, function(i) gsub("\\_\\\b","",i))
      argPosition = grep(Paste(argument,": "), text)
      if(length(argPosition) == 0) {
        next
      } else {
        argPosition = argPosition[1] - 1
        ##Found one
        visible(win) <- TRUE            # show window
      }

      add(textwindow,Paste("From package:",pkg), font.attr=c(weight="bold"))
      ## add first line (it has a :)
      add(textwindow,text[argPosition+1],font.attr=c(weight="bold",color="blue"))
      ## add until a :
      i = 2; n = length(text)
      while(length(grep(":",text[argPosition+i])) == 0 &&
            (argPosition + i) <= n
            ) {
        add(textwindow,text[argPosition+i],font.attr=c(weight="bold",color="blue"))
        i = i + 1
      }
      add(textwindow,"\n")
    }
  }
  ## close button
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  gbutton("cancel", container=buttonGroup,
          handler = function(h,...) dispose(h$obj))
  
}




##################################################
## build on ghelp widget to make a browser with search,
## simpler than old pmg.helpBrowser. Break that into components

##################################################
## ghelpbrowser
##' Widget to provide interface to help system
setClass("gHelpBrowser",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' help browser widget, stand alone window
##'
##' @export
ghelpbrowser <- function(
                         title = "Help browser", maxTerms = 100, width = 1000, height = 600 ,
                         ...,
                         toolkit=guiToolkit()) {
  widget <- .ghelpbrowser(toolkit,
                          title=title, maxTerms=maxTerms, width=width
                          )
  obj <- new( 'gHelpBrowser',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' API:
##' visible<-, logical. Display or hide sidebar
##' size<-, set size of window
##' size, size of window

##' generic for toolkit dispatch
##' @alias ghelpbrowser
setGeneric( '.ghelpbrowser' ,
           function(toolkit,
                    title = "Help browser", maxTerms = 100,
                    width = 1000, height = 600 )
           standardGeneric( '.ghelpbrowser' ))



## a notebook for holding help pages
setClass("gHelpbrowserANY",
         contains="gComponentANY",
         prototype=prototype(new("gComponentANY"))
         )


##################################################
## build on ghelp widget to make a browser with search,
## simpler than old pmg.helpBrowser. Break that into components

## a notebook for holding help pages
setClass("gHelpbrowserANY",
         contains="gComponentANY",
         prototype=prototype(new("gComponentANY"))
         )


setMethod(".ghelpbrowser",
          signature(toolkit="ANY"),
          function(toolkit,
                   title = "Help browser", maxTerms=100,
                   width=1000, height=600) {

            force(toolkit)

            ## Main widget
            helpBrowser <- gwindow(gettext("Help browser"), visible=FALSE)

            ## we need to check what toolkit
            toolkitType <- helpBrowser@toolkit@toolkit # hackery


            ##' layout for help search (apropos, pattern)
            helpSearch <- function(container, ...) {
              g <- ggroup(horizontal=FALSE, expand=TRUE, container=container, ...)
              sg <- ggroup(container=g, horizontal=TRUE, fill="x", ...)
              cb <- gcombobox(c("Apropos", "Pattern"), container=sg)
              e <- gedit("", container=sg, expand=TRUE, fill="x")
              
              sr <- gtable(data.frame("Function"=character(0), Package=character(0),
                                      Title=character(0), stringsAsFactors=FALSE),
                           container=g, expand=TRUE, fill="both")
              addHandlerClicked(sr, handler=function(h,...) {
                sel <- svalue(h$obj, drop=FALSE)
                if(!is.null(sel)) {
                  l <- list(topic=sel[[1]], package=sel[[2]])
                  add(helpWidget, l)
                }
              })
              
              searchResultsApropos = function(query) {
                out = help.search(apropos=query, ignore.case = TRUE)
                out = out$matches
                if(nrow(out) > 0) {
                  out = out[1:min(nrow(out),maxTerms),c(1,3,2), drop=FALSE]
                } else {
                  out = c("no matches","","")
                }
                colnames(out) = c("Function","Package","Title")
                out = as.data.frame(out)
                for(j in 1:3) out[,j] <- as.character(out[,j]) # avoid factors
                return(out)
              }
              ##' results for help.search
              searchResultsHelpSearch = function(query) {
                out = help.search(pattern=query, ignore.case = TRUE)
                out = out$matches
                if(nrow(out) > 0) {
                  out = out[1:min(nrow(out),maxTerms),c(1,3,2), drop=FALSE]
                } else {
                  out = c("no matches","","")
                }
                colnames(out) = c("Function","Package","Title")
                out = as.data.frame(out)    
                for(j in 1:3) out[,j] <- as.character(out[,j]) # avoid factors
                
                return(out)
              }
              
              addHandlerChanged(e, handler=function(h,...) {
                query <- svalue(h$obj)
                toolkitType <- svalue(cb, index=1)
                out <- switch(toolkitType,
                              searchResultsApropos(query),
                              searchResultsHelpSearch(query))
                sr[] <- out
              })
              
            }
            
            browsePackages <- function(container, ...) {
              getContentsOfPackage <- function(package) {
                ## return a data frame with entry keywords description
                path <- system.file("help", package = package)
                contents <- .readRDS(sub("/help", "/Meta/Rd.rds", path, fixed = TRUE))
                return(data.frame(Entry=contents[,'Name'],
                                  Keywords=sapply(contents[,"Keywords"], paste, collapse=", "),
                                  Description=contents[,'Title'],
                                  stringsAsFactors = FALSE))
              }
              emptyDf <- data.frame(Entry=character(0),
                                    Keywords=character(0),
                                    Description=character(0), stringsAsFactors=FALSE)
              
              allPackages <- .packages(all.available=TRUE)
              curPackage <- NULL
              
              g <- ggroup(container=container, horizontal=FALSE, expand=TRUE, ...)
              g1 <- ggroup(container=g)
              glabel("Package:", container=g1, anchor=c(1,0))
              e <- gedit("", container=g1, anchor=c(-1,0),
                         handler=function(h,...) {
                           val <- svalue(h$obj)
                           if(val %in% allPackages) {
                             curPackage <<- val
                             contents <- getContentsOfPackage(val)
                             fnList[] <- contents
                           } else {
                             curPackage <<- NULL
                             fnList[] <- emptyDf
                           }
                         })
              e[] <- allPackages
              
              fnList <- gtable(emptyDf,
                               container=g,
                               expand=TRUE)
              addHandlerClicked(fnList, handler=function(h,...) {
                topic <- svalue(h$obj)
                if(nchar(topic))
                  add(helpWidget, list(topic=topic, package=curPackage))
              })
            }
            
            
            ##' layout the search pane area
            layoutSearch <- function(container) {
              layoutNb <- gnotebook(container=container, expand=TRUE)
              helpSearch(container=layoutNb, label="Help search")
              browsePackages(container=layoutNb, label="Browse packages")
              svalue(layoutNb) <- 1     # first tab
            }
            
            
            ##' layout the help pane area
            layoutHelp <- function(container) {
              tb <- ggroup(container=container, horizontal=TRUE, fill="x")

              glabel("Help for:", container=tb, anchor=c(1,0))
              gedit("", container=tb, anchor=c(-1, 0),  handler=function(h,...) {
                val <- svalue(h$obj)
                add(helpWidget, val)
              })

              
              if(toolkitType == "tcltk") {
                ## dispose if tcltk, otherweise close buttons work
                gseparator(horizontal=FALSE, container=tb)
                d <- gbutton("dispose", container=tb, handler=function(h,...) {
                  dispose(helpWidget)
                })
              }
              
              if(!toolkitType == "Qt") {
                ## had errors with running withi handler
                gbutton("Example", container=tb, handler=function(h,...) {
                  page <- helpWidget[svalue(helpWidget)]
                  do.call("example", list(topic=tag(page, "topic"), package=tag(page, "package")))
                })
              }
              addSpring(tb)

              searchCb <<- gcheckbox("Search box", container=tb, handler=function(h,...) {
                visible(obj) <- svalue(h$obj)
              })
              
              helpWidget <<- ghelp(container=container, expand=TRUE, fill="both")
              
            }
            
            
            ## widgets
            pg <- gpanedgroup(container=helpBrowser, horizontal=TRUE)

            if(toolkitType == "RGtk2") {
              searchPane <- ggroup(horizontal=FALSE, container=pg)
              helpPane <- ggroup(horizontal=FALSE, container=pg)
            } else {
              helpPane <- ggroup(horizontal=FALSE, container=pg)
              searchPane <- ggroup(horizontal=FALSE, container=pg)
            }
            
            helpWidget <- NULL # defined in layoutHelp
            searchCb <- NULL
            
            layoutSearch(searchPane)
            layoutHelp(helpPane)

            ## show gwindow

            obj <- new("gHelpbrowserANY", block= helpBrowser, widget=helpBrowser, toolkit=toolkit)
            tag(obj, "pg") <- pg
            tag(obj, "searchCb") <- searchCb
            tag(obj, "toolkitType") <- toolkitType
            
            visible(helpBrowser) <- TRUE
            size(obj) <- c(width, height)
            visible(obj) <- FALSE       # hide sidebar            

            return(obj)
          })

##' toggle sidebar
setReplaceMethod(".visible", 
          signature(toolkit="ANY",obj="gHelpbrowserANY", value="logical"),
          function(obj, toolkit, ..., value) {
            cb <- tag(obj, "searchCb"); svalue(cb) <- value
            toolkitType <- tag(obj, "toolkitType")
            pg <- tag(obj, "pg")
            if(value) {
              val <- tag(pg, "lastPosition")
              if(is.null(val) || is.nan(val))
                val <- 0.6
              val <- max(min(0.75, val), 0.6)
              if(toolkitType == "RGtk2")
                val <- 1 - val
              svalue(pg) <- val
            } else {
              tag(pg, "lastPosition") <- svalue(pg)
              svalue(pg) <- ifelse(toolkitType=="RGtk2",0,1)
            }
            obj
          })

##' report widget size
setMethod(".size", 
          signature(toolkit="ANY",obj="gHelpbrowserANY"),
          function(obj, toolkit, ...) {
            w <- obj@widget
            size(w)
          })

##' Set widget size
setReplaceMethod(".size", 
                 signature(toolkit="ANY",obj="gHelpbrowserANY", value="numeric"),
                 function(obj, toolkit, ..., value) {
                   w <- obj@widget
                   size(w) <- value
                   obj
                 })

