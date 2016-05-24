## install cran

## helper functions
pmg.chooseCRANmirror = function(widget = NULL, doing.first=FALSE,...) {
  ## copied from tcltk widget

  ## if doing.first == TRUE, then call pmg.installCRANPackage after click
  ## This is from packages.R in chooseCRANmirror
  m <- try(read.csv(url("http://cran.r-project.org/CRAN_mirrors.csv"),
                    as.is=TRUE))
  if(inherits(m, "try-error"))
    m <- read.csv(file.path(R.home("doc"), "CRAN_mirrors.csv"), as.is=TRUE)  
    
  window=pmgWC$new(title="Select CRAN site", visible=FALSE)
  size(window) <- c(500,400)
  group = ggroup(horizontal=FALSE, container = window)

  tbl = gtable(utils::getCRANmirrors(), chosencol=4, 
    filter.column=2,
    handler = function(h,...) {
      URL = svalue(tbl)
      repos <- getOption("repos")
      repos["CRAN"] <- gsub("/$", "", URL[1])
      options(repos = repos)
      dispose(window)
      ## now install
      if(doing.first) pmg.installCRANPackage()
    })
  add(group, tbl, expand=TRUE)
  status = gstatusbar("Double click site to select", container=group)
  visible(window) <- TRUE               # now show the window
}

## return a data frame with the CRAN packages
.empty.CRANPackages.data.frame = function() {
  tmp = data.frame(Package="", CRAN.version="", Installed.version="",Depends="",Suggests="")
  for(j in 1:5) tmp[,j] = as.character(tmp[,j])
  return(tmp)
}
pmg.getCRANPackages = function() {
  
    ## what is installed?
    x <- installed.packages()
    i.pkgs <- as.character(x[, 1])
    i.vers <- as.character(x[, 3])
    
    
#    y = CRAN.packages()
    y = available.packages()
    if(nrow(y) == 0) {                    # if empty
      return(.empty.CRANPackages.data.frame())
    }

    
    c.pkgs <- as.character(y[, 1])
    c.vers <- as.character(y[, 2])
    c.depends <- as.character(y[,5])
    c.suggest <- as.character(y[,7])
    idx <- match(i.pkgs, c.pkgs)
    vers2 <- character(length(c.pkgs))
    
    xx <- idx[which(!is.na(idx))]
    vers2[xx] <- i.vers[which(!is.na(idx))]
    i.vers <- vers2
    
    cranPkgs = data.frame(
      Package=c.pkgs,
      CRAN.version=c.vers,
      Installed.version=i.vers,
      Depends = c.depends,
      Suggets = c.suggest
      )
    ## make character -- not factor, gtable barks otherwise
    for(j in 1:5) cranPkgs[,j] = as.character(cranPkgs[,j])
    ## filter out NA values
    cranPkgs = cranPkgs[!is.na(cranPkgs[,1]),]

    return(cranPkgs)
}

needToChooseCRANMirror = function() {
  if(is.null(getOption("repos")) ||
     is.na(getOption("repos")) ||
     getOption("repos") == "@CRAN@" ||
     getOption("repos") == ""
     )
    return(TRUE)
  else
    return(FALSE)
}


pmg.installCRANPackage = function() {

  if(needToChooseCRANMirror()) {
    return(pmg.chooseCRANmirror(doing.first=TRUE))
  }

  
  ## add list of packages to packageList
  addPackageList = function() {
    ## update repos from box
    ## we add to packageList provided various things are satisfied
    if(needToChooseCRANMirror()) {
      svalue(statusBar) <- "Set the CRAN repository before continuing"
      return()
    }
    if(is.null(svalue(libBox))) {
      svalue(statusBar)  <- "Set the 'Install to' directory before continuing."
      return()
    }

    svalue(statusBar)  <- "Loading available CRAN packages from internet"
    ## okay lets load it up

    m = pmg.getCRANPackages()
    packageList[,] = m
##    delete(packageGroup, packageList)
##    packageList <<- gtable(m , filter.labels = c("",letters),
##                               filter.FUN=filter.FUN)
##
##    add(packageGroup, packageList, expand=TRUE)
    enabled(installButton) <- TRUE                 # was grayed out
    svalue(statusBar) <-  ""
  }
  

  
  ## start with the GUI
  win = pmgWC$new("Install CRAN packages",v=T)

  mainGroup = ggroup(horizontal=FALSE, container=win)
  table = glayout(container=mainGroup)
  table[1,1] = glabel("CRAN repository:")
  reposBox = gedit(getOption("repos")[1],
    handler = function(h,...) {
      repository = as.character(svalue(h$obj))
      if(!is.empty(repository)) {
        options("repos",repository)
        addPackageList()
      }
  })
##   ## this is really changed
##   addhandlerkeystroke(reposBox,handler = function(h,...) {
##     repository = as.character(svalue(h$obj))
##     if(nchar(repository) > 0)
##       options("repos",repository)
##     addPackageList()
##   })
  table[1,2] = reposBox
  reposButton = gbutton("preferences",dirname="stock")
  addhandlerclicked(reposButton,
                    handler = function(h,...) {
                      pmg.chooseCRANmirror(widget=reposBox)
                      svalue(statusBar) <- ""
                      ## can't update here as repoxBox isn't set by now
                    })
  table[1,3] = reposButton
  
  table[2,1] = glabel("Install to:")
  libBox = gdroplist(.libPaths(), editable=TRUE,
    handler = function(h,...) {
      svalue(statusBar) <- ""
      addPackageList()
    }
    )
  table[2,2] = libBox
  libButton = gbutton("preferences",dirname="stock")
  addhandlerclicked(libButton, handler = function(h,...) {
    gfile("Pick a directory...", type = "selectdir", handler = function(h,...) {
      svalue(libBox) <- svalue(h$obj)
      addPackageList()
    })
  })
  table[2,3] = libButton

  table[3,1] = glabel("Install dependencies?")
  dependenciesBox = gdroplist(c("TRUE","FALSE"))
  table[3,2] = dependenciesBox
  
  table[4,1] = glabel("Package type:")
  typeBox = gdroplist(c("source","mac.binary","win.binary"))
  table[4,2] = typeBox
  
  statusBar = gstatusbar("")
  installButton = gbutton("Install selected package(s)",
    handler=function(h,...) {
      thePackages = svalue(packageList)
      svalue(statusBar) <- Paste("installing package: ",thePackages)
      install.packages(
                       pkgs = thePackages,
                       lib = svalue(libBox),
                       type = svalue(typeBox)
                       )
      svalue(statusBar) <- "Packages were installed. (Versions not updated)"
                                        #          dispose(win)
    })
  enabled(installButton) <- FALSE
  closeButton = gbutton("cancel",handler = function(h,...) {
    dispose(win)
  })
  

  packageGroup = ggroup()
  add(mainGroup, packageGroup, expand=TRUE)

  firstLetter = function(x) tolower(unlist(strsplit(x,""))[1])
  filter.FUN = function(d, val ) {
    if(val == "") 
      return(rep(TRUE, dim(d)[1]))
    else
      sapply(d[,1],firstLetter) == val
  }

  ## start with an empty data frame
  m = .empty.CRANPackages.data.frame()
  packageList <- gtable(m , filter.labels = c("",letters),
                             filter.FUN=filter.FUN)
  size(packageList) <- c(400,300)
  
  add(packageGroup, packageList, expand=TRUE)
  visible(table) <-TRUE
  gseparator(container=mainGroup)
  buttonGroup = ggroup(container=mainGroup)
  addSpring(buttonGroup)
  add(buttonGroup,installButton)
  add(buttonGroup,closeButton)
  
  add(mainGroup, statusBar)
  
  ## now add, hopefull the thing has been drawn alread
  addPackageList()


}
