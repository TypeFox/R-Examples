
###################################################
### code chunk number 145: HelpBrowser
###################################################
qsetClass("HelpBrowser", Qt$QTabWidget, function(parent=NULL){
  super(parent)
  #
  this$tabsClosable <- TRUE
  qconnect(this, "tabCloseRequested", function(index) {
    this$removeTab(index)
  })
  this$movable <- TRUE; this$usesScrollButtons <- TRUE
  #
  this$browser <- getOption("browser")
  options("browser" =  function(url) openPage(url))
})


###################################################
### code chunk number 146: openPage
###################################################
qsetMethod("openPage", HelpBrowser, function(url) {
  tokens <- strsplit(url, "/")[[1]]
  tab_title <- sprintf("%s: %s", tokens[length(tokens)-2], 
                       tokens[length(tokens)])
  webview <- Qt$QWebView()
  webview$setUrl(Qt$QUrl(url))
  this$currentIndex <- addTab(webview, tab_title)
})


###################################################
### code chunk number 147: illustrateHelpBrowser
###################################################
help_browser <- HelpBrowser()
help_browser$windowTitle <- "Help Browser example"
help_browser$show()
help_browser$raise()
##
options("help_type"="html")
help("mean")
help("boxplot")

