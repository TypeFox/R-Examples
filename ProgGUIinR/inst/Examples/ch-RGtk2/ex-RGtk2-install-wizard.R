### R code from vignette source 'ex-RGtk2-install-wizard.Rnw'

###################################################
### code chunk number 1: installPackagesWizard
###################################################
## gtk Assistant example
require(RGtk2)


###################################################
### code chunk number 2: defineAssistant
###################################################
assistant <- gtkAssistant(show=FALSE)
assistant$setSizeRequest(500, 500)
gSignalConnect(assistant, "cancel", 
               function(assistant) assistant$destroy())


###################################################
### code chunk number 3: makePages
###################################################
pages <- lapply(1:5, gtkVBox, spacing = 5, homogeneous = FALSE)
page_types <- c("intro", rep("confirm", 3), "summary")
sapply(pages, gtkAssistantAppendPage, object = assistant)
sapply(pages, gtkAssistantSetPageType, object = assistant, 
       type=page_types)


###################################################
### code chunk number 4: sideLogo1
###################################################
image <- gdkPixbuf(filename = imagefile("rgtk-logo.gif"))[[1]]
sapply(pages, gtkAssistantSetPageSideImage, object = assistant, 
       pixbuf = image)


###################################################
### code chunk number 5: ex-RGtk2-install-wizard.Rnw:52-59
###################################################
populate_page <- list()                
gSignalConnect(assistant, "prepare", 
       function(assistant, page, data) {
         page_no <- which(sapply(pages, identical, page))
         if(!length(page$getChildren()))
           populate_page[[page_no]]()
       })


###################################################
### code chunk number 6: ex-RGtk2-install-wizard.Rnw:68-74
###################################################
assistant$setForwardPageFunc(function(page_index, data) {
  if(page_index == 0 && have_CRAN()) 
    2L 
  else 
    as.integer(page_index + 1)
}, data=NULL)


###################################################
### code chunk number 7: ex-RGtk2-install-wizard.Rnw:78-80
###################################################
CRAN_package <- NA
install_options <- list() #type, dependencies, lib


###################################################
### code chunk number 8: HelperFunctions
###################################################
## Helper functions
##' return value or NA
##'
gtkTreeViewGetSelectedValue <- function(object, column) {
  cur <- object$getSelection()$getSelected()
  if(cur$retval)
    with(cur, object$getModel()$getValue(iter, column -1 )$value)
  else
    NA
}


have_CRAN <- function() getOption("repos")["CRAN"] != "@CRAN@"

##' from getCRANmirrors
set_CRAN <- function(url) {
  if(is.null(url)) return()
  repos <- getOption("repos")
  repos["CRAN"] <- gsub("/$", "", url)
  options(repos=repos)
}


###################################################
### code chunk number 9: page1
###################################################
populate_page[[1]] <- function() {
  assistant$setPageTitle(pages[[1]], "Install a CRAN package")
  pages[[1]]$packStart(label <- gtkLabel())
  pages[[1]]$packStart(gtkLabel(), expand=TRUE) # a spring
  
  label$setMarkup(paste(
       "<span font='x-large'>Install a CRAN package</span>",
       "This wizard will help install a package from",
       "<b>CRAN</b>. If you have not already specified a",
       "CRAN repository, you will be prompted to do so.",
       sep="\n"))
  assistant$setPageComplete(pages[[1]], TRUE)
}


###################################################
### code chunk number 10: CRANMirror
###################################################
## Not shown
populate_page[[2]] <- function() {
  assistant$setPageTitle(pages[[2]], "Select a CRAN mirror")

  CRAN_mirrors <- getCRANmirrors(all = FALSE, local.only = FALSE)[, c(1,2,4)]
  nms <- names(CRAN_mirrors)
  d <- rGtkDataFrame(CRAN_mirrors)
  #
  view <- gtkTreeView()
  mapply(view$insertColumnWithAttributes, -1, nms[1:2], 
         list(gtkCellRendererText()), text = 0:1)
  view$setModel(d)
  view$getSelection()$unselectAll()     # no selection
  gSignalConnect(view$getSelection(), "changed", function(view, ...) {
    CRAN_repos <- view$getSelectedValue(3)
    set_CRAN(CRAN_repos)
    assistant$setPageComplete(pages[[2]], TRUE)
  }, data=view, user.data.first=TRUE)
  
  
  sw <- gtkScrolledWindow(); sw$add(view)
  sw$setPolicy("automatic", "automatic")
  
  pages[[2]]$packStart(gtkLabel("Select a CRAN mirror"), expand=FALSE)
  pages[[2]]$packStart(sw, expand=TRUE, fill=TRUE)

}


###################################################
### code chunk number 11: SelectPacakge
###################################################
## Not shown
populate_page[[3]] <- function() {
  assistant$setPageTitle(pages[[3]], "Select a CRAN package")
  #
  avail_packages <- available.packages()[, c(1,2)]
  nms <- colnames(avail_packages)
  avail_packages_store <- rGtkDataFrame(avail_packages)
  #
  view <- gtkTreeView()
  mapply(view$insertColumnWithAttributes, -1, nms, 
         list(gtkCellRendererText()), text = 0:1)
  view$setModel(avail_packages_store)
  view$getSelection()$unselectAll()     # no selection
  gSignalConnect(view$getSelection(), "changed", function(view, ...) {
    CRAN_package <<- view$getSelectedValue(1)
    assistant$setPageComplete(pages[[3]], TRUE)
  }, data=view, user.data.first=TRUE) 
  #
  sw <- gtkScrolledWindow(); sw$add(view)
  sw$setPolicy("automatic", "automatic")
  #
  pages[[3]]$packStart(gtkLabel("Select a package to install"), expand=FALSE)
  pages[[3]]$packStart(sw, expand=TRUE, fill=TRUE)
}


###################################################
### code chunk number 12: ex-RGtk2-install-wizard.Rnw:193-262
###################################################
populate_page[[4]] <- function() {
  assistant$setPageTitle(pages[[4]], "Install a CRAN package")
  ##
  get_desc <- function(pkgname) {
    o <- "http://cran.r-project.org/web/packages/%s/%s"
    x <- readLines(sprintf(o, pkgname, "DESCRIPTION"))
    f <- tempfile(); cat(paste(x, collapse="\n"), file=f)
    read.dcf(f)
  }
  desc <- get_desc(CRAN_package)
  #
  label <- gtkLabel()
  label$setLineWrap(TRUE)
  label$setWidthChars(40)
  label$setMarkup(paste(
    sprintf("Install package: <b>%s</b>", desc[1,'Package']),
    "\n",
    sprintf("%s", gsub("\\n", " ", desc[1,'Description'])),
    sep="\n"))
  
  pages[[4]]$packStart(label)
  ##
  table <- gtkTable()
  pages[[4]]$packStart(table, expand=FALSE)
  pages[[4]]$packStart(gtkLabel(), expand=TRUE)
  
  ##
  combo <- gtkComboBoxNewText()
  pkg_types <- c("source", "mac.binary", "mac.binary.leopard",
                 "win.binary", "win64.binary")
  sapply(pkg_types, combo$appendText)
  combo$setActive(which(getOption("pkgType") == pkg_types)-1)
  gSignalConnect(combo, "changed", function(combo, ...) {
    cur <- 1L + combo$getActive()
    install_options[['type']] <<- pkg_types[cur]
  })
  table$attachDefaults(gtkLabel("Package type:"), 0, 1, 0, 1)
  table$attachDefaults(combo, 1, 2, 0, 1)

  ##
  checkButton <- gtkCheckButton()
  checkButton$setActive(TRUE)
  gSignalConnect(checkButton, "toggled", function(ck_btn) {
    install_options$dependencies <<- ck_btn$getActive()
  })
  table$attachDefaults(gtkLabel("Install dependencies"),
                       0, 1, 1, 2)
  table$attachDefaults(checkButton, 1, 2, 1, 2)

  ##
  file_chooser <- gtkFileChooserButton("Select directory...", 
                                      "select-folder")
  file_chooser$setFilename(.libPaths()[1])
  gSignalConnect(file_chooser, "selection-changed", 
                 function(file_chooser) {
                   dir <- file_chooser$getFilename()
                   install_options[['lib']] <<- dir
                 })
  table$attachDefaults(gtkLabel("Where"), 0, 1, 2, 3)
  table$attachDefaults(file_chooser, 1, 2, 2, 3)
  ## align labels to right and set spacing
  sapply(table$getChildren(), function(child) {
    widget <- child$getWidget()
    if(is(widget, "GtkLabel"))  widget['xalign'] <- 1
  })
  table$setColSpacing(0L, 5L)
  ##
  assistant$setPageComplete(pages[[4]], TRUE)
}


###################################################
### code chunk number 13: ex-RGtk2-install-wizard.Rnw:271-290
###################################################
populate_page[[5]] <- function() {
  assistant$setPageTitle(pages[[5]], "Done")
  install_options$pkgs <- CRAN_package
  out <- try(do.call("install.packages", install_options), 
             silent=TRUE)

  label <- gtkLabel(); pages[[5]]$packStart(label)
  if(!inherits(out, "try-error")) {
    label$setMarkup(sprintf("Package %s was installed.", 
                            CRAN_package))
  } else {
    label$setMarkup(paste(sprintf("Package %s, failed install", 
                                  CRAN_package),
                          paste(out, collapse="\n"),
                          sep="\n"))
  }

  assistant$setPageComplete(pages[[5]], FALSE)
}


###################################################
### code chunk number 14: showAssistant
###################################################
populate_page[[1]]()
assistant$show()


