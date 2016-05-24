"packageAdd" <-
  function(pkg, files, path = ".", document = TRUE) {
    ## utility to find the package's env on the search list
    findPkg <- function(pkg) {
      ev = .GlobalEnv
      pkgName <- paste("package", pkg, sep=":")
      while(!is.null(ev)) {
        ev = parent.env(ev)
        if(identical(pkgName, attr(ev, "name")))
          break
      }
      ev
    }
    evPkg <- findPkg(pkg)
    if(is.null(evPkg)) {
      require(pkg, character.only = TRUE)
      evPkg <- findPkg(pkg)
      if(is.null(evPkg))
        stop("Package \"", pkg,
             "\" should be available when packageAdd() is called")
    }
    topenvPrev <- options("topLevelEnvironment")
    on.exit(options(topenvPrev))
    for(file in files) {
      ev = new.env(parent = evPkg)
      ## simulate the package name as inserted by library()
      ## Needed for classs, method definitions but much too low-level
      ## to be exposed to poor programmers.  Oh well!
      assign(".packageName", pkg, envir = ev)
      ## arrange for default topenv() to go to this ev
      options(topLevelEnvironment = ev)
      exprs <- parse(file)
      eval(exprs, envir = ev)
      sourceFile = basename(file)
      manDir = file.path(path, pkg, "man")
      what <- objects(ev, all.names=TRUE)
      if(document) {
        docCommon <- character() # objects documented together
        for(name in what) {
          obj <- get(name, envir = ev)
          if(is(obj, "MethodsList")) {
            fName <- metaNameUndo(name)
            packageSlot(fName) <- pkg
            fileName <- file.path(manDir, paste(.topicName("methods",fName),"Rd",sep="."))
            promptMethods(fName,  fileName, obj)
          }
          else if(is(obj, "classRepresentation")) {
            clName <- as.character(metaNameUndo(name, "C"))
            fileName <- file.path(manDir, paste(.topicName("class",clName),"Rd",sep="."))
            promptClass(clName, fileName, where = ev)
          }
          else
            docCommon <- c(docCommon, name)
        }
        if(length(docCommon)) {
          fileName <- file.path(manDir, paste(docCommon[[1]], "Rd", sep="."))
          promptAll(docCommon, fileName)
        }
      }
      sourceCopyFile <- file.path(path, pkg, "R", sourceFile)
      file.copy(file, sourceCopyFile)
      message("Copied file ", file, " to ", sourceCopyFile)
    }
  }

.topicName <- function(type, what)
  paste(what, type, sep="-")
