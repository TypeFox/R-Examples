getPackages = function(...) {
  allPackages = .packages(all.available=TRUE)
  loaded = allPackages %in% .packages()
  data.frame(Package=allPackages, loaded=loaded, stringsAsFactors=FALSE)
}


pmg.loadPackages = function(width=300, height=400) {

  win = pmgWC$new("Load or detach packages", v=T)
  size(win) <- c( width, height)
  group = ggroup(horizontal=FALSE, container=win, expand=TRUE)


  packageHandler = function(obj) {
    packages = svalue(packageList, drop=FALSE)
    for(i in 1:nrow(packages)) {
      package = packages[i,1]              # character
      installed = packages[i,2]                # logical
      if(installed == TRUE) {
        cat("Detach",package,"\n")
        svalue(status) <- Paste("detach package ",package)
        pkg = Paste("package:",package)       # see help on detach
        detach(pos = match(pkg, search()))
        svalue(status)
      } else {
        svalue(status) <- Paste("Load package ",package)
	fn <- get("require")
        res = fn(package, character.only=TRUE)
        if(res == FALSE)
          cat(Paste("Couldn't load package",package,"\n"))
        else
          cat("Loaded package ",package,"\n")
        svalue(status)
      }
    }
    ## updata package list, 
    packageList[,] = getPackages()
  }

  ## store package into a separate group -- can update
  packageGroup = ggroup(container=group, expand=TRUE)
  packageList = gtable(getPackages(),
    multiple=TRUE, sort.columns = 1:2,
    handler = function(h,...) {
      packageHandler(h$obj)
    },
    container=packageGroup, expand=TRUE)

  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
#  gbutton("ok",handler = function(h,...) packageHandler(h$action),
#          action=packageList, container=buttonGroup)
  gbutton("cancel",container=buttonGroup, handler = function(h,...) dispose(win))

  status = gstatusbar("Double click on  package to load/detach",container=group)
  
    
}
