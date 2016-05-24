
## Make a widget holding all the dynamic test

dtestsDialog = function() {

  tests = list(
    "1-sample test of proportion" = "gui.prop.test",
    "1-sample exact test of proportion" = "gui.binom.test",
    "2-sample test of proportion" = "gui.prop.test.2sample",
    "2-sample exact test of proportion" = "gui.binom.test.2sample",
    ##
    "1-sample t-test"="gui.t.test",
    "1-sample signed rank test"="gui.wilcox.signed.rank.test",
    "2-sample t-test" =   "gui.t.test.2.sample",
    "2-sample t-test var. equal" = "gui.t.test.2.sample.var.equal",
    "2-sample rank sum test" = "gui.wilcox.rank.sum.test",
    ##
    "Oneway ANOVA" = "gui.oneway.test",
    "Kruska-Wallis test" = "gui.kruskal.test",
    ##
    "Correlation test" = "gui.cor.test",
    ##
    "Test of variances" = "gui.var.test",
    "Ansari test" = "gui.ansari.test",
    "Bartlett test" = "gui.bartlett.test",
    "Fligner test" = "gui.fligner.test",
    "Mood test" = "gui.mood.test",
    ##
    "2-sample Kolmogorov-Smirnov test" = "gui.ks.test",
    "Shapiro.test" = "gui.shapiro.test"
    
    )
  
  dialogList = list()
#  for(i in names(tests)) {
#    tmp = tests[[i]]
#    if(length(tmp) == 1)
#      dialogList[[i]] <- do.call(tmp,list())
#    else
#      dialogList[[i]] <- do.call(tmp[1],tmp[2]) # tmp a list with args
#  }
  dialogList[["FirstOne"]] = ilabel("Select a test from popup")
  
  
  win = iwindow("Tests",v=T)
  gp = igroup(horizontal=FALSE, window=win)
  popupGroup = igroup(window=gp)
  add.spring(popupGroup)
  testPopup = idroplist(c("",names(tests)), win=popupGroup)
  
  testWindow = igroup(window=gp)
  add(testWindow,dialogList[["FirstOne"]], expand=TRUE)
  setdata(testWindow,dialogList,"dialogList")
  setdata(testWindow,dialogList[["FirstOne"]],"currentTest")
  
  addhandlerchanged(testPopup, handler = function(h,...) {
    popupValue = get.value(testPopup)
    if(!is.empty(popupValue)) {
      delete(testWindow,getdata(testWindow,"currentTest"))
      dialogList = getdata(testWindow, "dialogList")
      if(is.null(dialogList[[popupValue]])) {
        dialogList[[popupValue]] <- do.call(tests[[popupValue]],list())
        setdata(testWindow, dialogList, "dialogList")
      }
      add(testWindow,dialogList[[popupValue]], expand=TRUE)
      setdata(testWindow,dialogList[[popupValue]],"currentTest")
  }
  })

}


pmg.dynamic = function() {
  ## show dynamic dialogs
  pmg.dw = function(FUN.name,title=NULL) {
    if(is.null(title))
      title = gsub("^gui\\.","",FUN.name)
    widget = do.call(FUN.name,list(window = iwindow(title=title,visible=TRUE)))
  }
    
  ## menubar
  pmg.menu = list()
  pmg.menu$File$"Source file.."$handler =
    function(h,...) ifile("Source file",type="open", action="source")
  pmg.menu$File$"Source file.."$icon="file"
  pmg.menu$File$"Save Workspace..."$handler =
    function(h,...) ifile("Save workspace",type="save", action="save.image")
  pmg.menu$File$"Save Workspace..."$icon = "save"
  pmg.menu$File$"Restore Workspace"$handler =
    function(h, ...) ifile("Restore workspace",type="open", action="load")
  pmg.menu$File$"Restore Workspace"$icon = "revert-to-saved"
  pmg.menu$File$"Load package..."$handler =
    function(h,...) pmg.loadPackages()
  pmg.menu$File$"Install CRAN package..."$handler =
    function(h,...) pmg.installCRANPackage()
  pmg.menu$File$"Install CRAN package..."$icon = "network"
  pmg.menu$File$"Install local package..."$handler =
    function(h,...) {
      old = options("repos")$repos; options("repos"=NULL);
      ifile("Select a package file...","open",action="install.packages")
      options("repos"=old)
    }
  pmg.menu$File$"Install local package..."$icon = "file"
  pmg.menu$File$"Set working directory..."$handler =
    function(h,...) ifile("Select a directory","selectdir",action="setwd")
  pmg.menu$File$"Set working directory..."$icon = "directory"
  
  pmg.menu$File$"pmg options..."$handler =
    function(h,...) pmg.options()
  pmg.menu$File$"pmg options..."$icon = "preferences"
  pmg.menu$File$"Exit pmg"$handler =
    function(h,...)  {
      dispose(win)
    }
  pmg.menu$File$"Exit pmg"$icon ="quit"
  ##
  ##
  pmg.menu$Data$browseEnv$handler =
    function(h,...) browseEnv()
  pmg.menu$Data$"Load data set..."$handler =
    function(h,...) pmg.viewDataSets()

  pmg.menu$Plot$"Lattice Explorer"$handler =
    function(h,...) ilatticeexplorer(win = iwindow("Lattice Explorer", visible=TRUE))
##################################################
##  Dynamic tests



pmg.menu$Dtests$"All tests widget"$handler = function(h,...) dtestsDialog()
  
pmg.menu$DTests$centers$"t.test - one.sample"$handler =
  function(h,...) pmg.dw("gui.t.test")
pmg.menu$DTests$centers$"t.test - two.sample"$handler =
  function(h,...) pmg.dw("gui.t.test.2.sample")
pmg.menu$DTests$centers$"t.test - two.sample, var. equal"$handler =
  function(h,...) pmg.dw("gui.t.test.2.sample.var.equal", "t.test, var.equal=TRUE")
pmg.menu$DTests$centers$"wilcox signed rank test"$handler =
  function(h,...) pmg.dw("gui.wilcox.signed.rank.test","wilcox.test, one sample")
pmg.menu$DTests$centers$"wilcox rank sum test"$handler =
  function(h,...) pmg.dw("wilcox.rank.sum.test","wilcox test, two sample")
pmg.menu$DTests$centers$"oneway.test"$handler =
  function(h,...) pmg.dw("gui.oneway.test")
pmg.menu$DTests$centers$"kruskal.test"$handler =
  function(h,...) pmg.dw("gui.kruskal.test")
#
pmg.menu$DTests$scales$"var.test"$handler =
  function(h,...) pmg.dw("gui.var.test")
pmg.menu$DTests$scales$"ansari.test"$handler =
  function(h,...) pmg.dw("gui.ansari.test")
pmg.menu$DTests$scales$"bartlett.test"$handler =
  function(h,...) pmg.dw("gui.bartlett.test")
pmg.menu$DTests$scales$"fligner.test"$handler =
  function(h,...) pmg.dw("gui.fligner.test")
pmg.menu$DTests$scales$"mood.test"$handler =
  function(h,...) pmg.dw("gui.mood.test")
#
pmg.menu$DTests$shape$"ks.test"$handler =
  function(h,...) pmg.dw("gui.ks.test")
Shapiro.test = function(x,...) shapiro.test(x)
pmg.menu$DTests$shape$"shapiro.test"$handler =
  function(h,...) pmg.dw("gui.shapiro.test")
#
pmg.menu$DTests$proportion$"prop.test, one sample"$handler =
  function(h,...) pmg.dw("gui.prop.test")
pmg.menu$DTests$proportion$"prop.test, two sample"$handler =
  function(h,...) pmg.dw("gui.prop.test.2sample","prop.test")
pmg.menu$DTests$proportion$"binom.test, one sample"$handler =
  function(h,...) pmg.dw("gui.binom.test")
pmg.menu$DTests$proportion$"binom.test, two sample"$handler =
  function(h,...) pmg.dw("gui.binom.test.2sample","binom.test")
#
#pmg.menu$DTests$counts$"chisq.test"$handler =
#  function(h,...) pmg.dw(chisq.test.list)
#pmg.menu$DTests$counts$"mantelhaen.test"$handler =
#  function(h,...) pmg.dw(mantelhaen.test.list)
#pmg.menu$DTests$counts$"mcnemar.test"$handler =
#  function(h,...) pmg.dw(mcnemar.test.list)
#

  pmg.menu$DTests$correlation$"cor.test"$handler =
  function(h,...) pmg.dw("gui.cor.test")
#

  ## rename!
  pmg.menu$Tests = pmg.menu$DTests
  pmg.menu$DTests <- NULL
  
  ### HELP
pmg.menu$Help$"About PMG"$handler =
  function(h,...) pmg.about(win=iwindow("About P M G",v=T))
  ##pmg.about(window=iwindow(v=TRUE))
  pmg.menu$Help$"About PMG"$icon="about"
  pmg.menu$Help$"R helpbrowser"$handler =
    function(h,...) ihelpbrowser(win=iwindow("P M G Help",v=T))
  pmg.menu$Help$"R Site Search"$handler = function(h,...) RSiteSearch.Dialog()
  pmg.menu$Help$"View vignettes"$handler = function(h,...) viewVignettes.Dialog()
  pmg.menu$Help$"View demos"$handler = function(h,...) viewDemos.Dialog()
  pmg.menu$Help$"View vignette"$handler = function(h,...) vignette("manual",package="pmg")
  
  ### Make the main window
  win = iwindow("P M G", visible=TRUE)
  set.size(win, 500,300)                #wide, short
  gp = igroup(horizontal = FALSE, window = win)
  add(gp, imenu(pmg.menu))              # menubar
  add(gp,ipanedgroup(ivarbrowser(), idfnotebook()), expand=TRUE)
  add(gp,(statusBar <- istatusbar()))
  set.value(statusBar,"P M G Dynamic")
}

pmg.dynamic()
