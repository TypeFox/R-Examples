## this is for interface to help stuff within pmg
## much of this moved into ghelp.R via ghelpbrowser



## make a global
pmg.helppagebrowser = NA


### RSiteSearch Dialog
RSiteSearch.Dialog = function() {
  ## RSiteSearch
  win = pmgWC$new("RSiteSearch()", visible=TRUE)
  
  table = glayout()
  spacinggroup = ggroup(horizontal=FALSE)
  ##  sitegroup = ggroup(container=spacinggroup)
  ## iaddlabel(siteentry,"Search terms:", pos=2, container=spacinggroup)
  table[1,1] =  glabel("Search terms:")
  siteentry = gedit("", width=50)#,container=sitegroup)
  table[1,2:4] = siteentry
  
  
  restrict = c("","Rhelp02a","Rhelp01","docs","functions")
  restrictPopup = gdroplist(restrict,multiple=TRUE)#,container=spacinggroup)
                                        #  size(restrictPopup) <- c(200,200)
#  iaddlabel(restrictPopup, "Restrict to:",pos=2, container=spacinggroup)
  table[2,1] = glabel("Restrict to:")
  table[2:3,2] = restrictPopup
  
  matchesPerPage = gspinbutton(min=20,max=100,value=20, step=10)#,container=spacinggroup)
  table[4,1] = glabel("Matches per page:")
  table[4,2] = matchesPerPage
#  iaddlabel(matchesPerPage,"Matches per page", pos=2, container=spacinggroup)
#  addSpring(spacinggroup)
#  add(spacinggroup,sitebutton,expand=FALSE)
  visible(table) <-  TRUE
  gp = ggroup(horizontal=FALSE, container = win)
  glabel("RSiteSearch goes to the internet to find\nmatches to your queries. The results are shown\nin a browser window.",container=gp)
  add(gp,gseparator())
  add(gp, table)
  buttonGroup = ggroup(container=gp)
  addSpring(buttonGroup)
  sitebutton = gbutton("find", container=buttonGroup)
  add(buttonGroup, gbutton("cancel",handler=function(h,...) dispose(win)))
  
  addhandlerclicked(sitebutton,action=siteentry,handler=function(h,...) {
    restrictValues = svalue(restrictPopup)
    if(is.null(restrictValues) || restrictValues=="")
      restrictValues = restrict[-1] # all of them but ""
    matches = svalue(matchesPerPage)
    RSiteSearch(svalue(h$action), restrict=restrictValues, matchesPerPage=matches)
  })

}


### View Vignettes dialog
viewVignettes.Dialog = function() {
  allVignettes = getAllVignettes()

  defaultHandler = function(h,...) {
    tmp =  svalue(vignetteList, drop=FALSE)
    topic = tmp[1,2, drop=TRUE]
    package = tmp[1,1, drop=TRUE]
    cat("Show vignette for ",topic,"\n")
    print(do.call("vignette",list(topic=topic, package=package)))
  }
  
  vignetteList = gtable(as.data.frame(allVignettes),
    filter.column = 1,
    handler = defaultHandler
  )

  win = pmgWC$new("View vignette",v=T)
  gp = ggroup(horizontal=FALSE, container=win)
  add(gp, vignetteList, expand=TRUE)
  buttonGroup = ggroup(container=gp)
  addSpring(buttonGroup)
  gbutton("ok",container = buttonGroup, handler = defaultHandler)
  gbutton("cancel",container=buttonGroup, handler = function(h,...) dispose(win))
  size(win) <- c(450,300)
}

### Demos
viewDemos.Dialog = function() {
  allDemos = demo(package = .packages(all.available = TRUE))$results
  demoList = gtable(allDemos[,-2], chosencol = 2, filter.column=1)
  addhandlerdoubleclick(demoList, handler=function(h,...) {
    cat("Demo runs in command line area")
    item = svalue(h$obj, drop=FALSE)
    do.call("demo",list(topic=item[1,2,drop=TRUE], package=item[1,1,drop=TRUE]))
  })
  ## create widget
  win = pmgWC$new("View demo in command line",v=T)
  gp = ggroup(horizontal=FALSE, container=win)
  add(gp, demoList, expand=TRUE)
  buttonGroup = ggroup(container=gp)
  addSpring(buttonGroup)
  gbutton("ok",container = buttonGroup, handler = function(h,...) {
    cat("Demo runs in command line area")
    item = svalue(demoList, drop=FALSE)
    do.call("demo",list(topic=item[1,2,drop=TRUE], package=item[1,1,drop=TRUE]))
    })
  gbutton("cancel",container=buttonGroup, handler = function(h,...) dispose(win))
  size(win) <- c(450,300)
}

##################################################
## helper furntions
## list all packages
getAllPackages = function() .packages(all.available=TRUE)

#
getAllVignettes = function() {
  allVignettes = hack.as.data.frame.matrix(vignette()[[4]][,c(1,3)])
  ## just the Packages name and Item
  return(allVignettes)
}

## in ghelp of iWidgets
## contents a matrix with entry, keywords, description and URL
#getContentsOfPackage = function(package=NULL) {
#  if(is.null(package)) {
#    warning("Empty package name")
#    return(NA)
#  }
#  contents = read.dcf(system.file("CONTENTS",package=package))
  
#  return(data.frame(Entry=I(contents[,1]),Keywords=I(contents[,3]),
#                    Description=I(contents[,4])))
#}

