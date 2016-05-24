
checkForUpdatesGUI = function() {
  win = pmgWC$new("Check for updates", visible=FALSE)
  g = ggroup(horizontal = FALSE, container=win)
  l = glabel(".-.-.-", container=g, expand=TRUE)
  sb = gstatusbar("Checking for updates", container=g)
  visible(win) <- TRUE
  
  val = checkForUpdates()
  if(length(val) == 0)
    svalue(l) <- "All the main pmg packages are up to date"
  else
    svalue(l) <- paste("You can upgrade", paste(val,sep=", "),".",sep=" ")
  svalue(sb) <- ""
}

checkForUpdates = function() {
  ## find any packages needing updates

  if(getCranSiteIfNeeded()) {
  
    thePackages = c("pmg","gWidgets","gWidgetsRGtk2","cairoDevice")
    
    oldPackages = old.packages()
    
    updateThese = thePackages[which(thePackages %in% rownames(oldPackages))]
    
    if(length(updateThese) > 0) {
      return(oldPackages[updateThese,"Package"])
    } else {
      return(c())
    }
  } else {
    cat("You need to set a CRAN repository to proceed\n")
  }
  
}

getCranSiteIfNeeded = function() {
  repos = getOption("repos")
  if ("@CRAN@" %in% repos) {
  
    setCRAN <- function(URL) {
      repos = getOption("repos")
      repos["CRAN"] <- gsub("/$", "", URL)
      options(repos=repos)
    }
    
    
    handler = function(h,...) {
      URL <- svalue(tbl) # get value  widget
      cat("Set CRAN site to",URL,"\n")
      setCRAN(URL)       # set URL
    }

    g = ggroup(horizontal = FALSE, container= NULL)
    glabel("Select a site\nthen click 'OK'", container=g)
    tbl <- gtable(
                  items=utils::getCRANmirrors(),
                  chosencol=4,     
                  filter.column=2,
                  container=g,
                  )
    size(tbl) <- c(200,300)
    gbasicdialog(title="Select a CRAN site", widget=g, handler=handler)
  } else {
    return(TRUE)
  }
}

#########################################

pmg.about = function(container=NULL) {

## image is group pmg via www.geom.uiuc.edu/~dpvc
  
  group = ggroup(horizontal=FALSE,container=container)
  size(group) <-  c(500,500)
  theFactsMam = read.dcf(system.file("DESCRIPTION",package="pmg"))
  glabel(Paste(
               "<b> P M G</b>\n",
               "<i>",
               theFactsMam[1,'Title'],
               "</i>\n",
               "Version ",
               theFactsMam[1,"Version"],
               "\n\n",
               theFactsMam[1,'URL'],
               "\n",
               "Comments to pmgRgui@gmail.com\n",
               "\n\n",
               theFactsMam[1,"Author"],
               "\n\n",
               theFactsMam[1,"Description"],
               "\n"
               ), markup=TRUE, container=group)
  addSpring(group, 10)
  gbutton("Check for updates", container = group,
          handler = function(...) checkForUpdatesGUI())
  
  return(group)
}

pmg.about.R = function(container=NULL) {

## image is group pmg via www.geom.uiuc.edu/~dpvc
  
  group = ggroup(horizontal=FALSE,container=container)
  gimage(system.file("images","Rlogo.jpg",package="pmg"),  container=group)
  glabel(paste(
               "<b> R </b>",
               "is a free software environment for statistical\n",
               "computing and graphics.\n\n",
               "<i>http://www.r-project.org</i>\n\n",
               R.version.string,
#               "Version ",
#               paste(R.version$major,R.version$minor,sep="."),
               "\n\n",
               sep=" ", collapse=" "
               ), markup=TRUE, container=group)

  return(group)
}
