# require(shinyBS) # run second R to serve app
# bsDemo(port = 7654, FALSE)


require(rDVR)
require(RSelenium)
# start the video server
startVideoServer(savedir = getwd())
DVR <- rDVR(saveDir = getwd())
Sys.sleep(5)
DVR$start(fileName = "secondVid")
RSelenium::startServer() # start selenium server
remDr <- remoteDriver()
remDr$open()
remDr$maxWindowSize()

remDr$navigate("http://localhost:7654/")
tabElems <- remDr$findElements("css selector", "#tabset a")
length(tabElems)
sapply(tabElems, function(x){x$getElementText()[[1]]})
# tab ids change from session to session
tabID <- sapply(tabElems, function(x){x$getElementAttribute("href")[[1]]})
tabID <- gsub(".*(#.*)", "\\1", tabID)
# Cycle thru each tab using highlightElement then retrun to introduction
out <- lapply(tabElems[c(2:10, 1)], function(x){
  x$highlightElement()
  x$clickElement()
  Sys.sleep(2)
})
# INTERACT WITH THE NAVBARS
tabElems[[2]]$highlightElement()
tabElems[[2]]$clickElement()
# customise the Brand
brandElem <- remDr$findElement("id", "nbBrand")
brandElem$highlightElement()
brandElem$clearElement()
brandElem$sendKeysToElement(list("BossNB", key = "enter"))
Sys.sleep(2)
# click the Link 5 times
linkElem <- remDr$findElement("id", "nbLink1")
out <- try(remDr$executeScript("arguments[0].scrollIntoView(true);", list(linkElem)), silent = TRUE) # scroll down if necessary
linkElem$highlightElement()
out <- lapply(seq(5), function(x){linkElem$clickElement()
                                  Sys.sleep(1)
})

# navigate the dropdown
ddElem <- remDr$findElement("id", "nbdd")
ddElem$highlightElement()
# complicated dropdown easiest to interact by changing class
ddElem$setElementAttribute("class", "dropdown sbs-dropdown shiny-bound-input open")
ddElems <- remDr$findElements("css selector", "#nbdd li a")
sapply(ddElems, function(x){x$getElementText()[[1]]})
out <- lapply(ddElems, function(x){
  ddElem$setElementAttribute("class", "dropdown sbs-dropdown shiny-bound-input open")
  x$highlightElement()
  x$clickElement()
  Sys.sleep(1)
})

# toggle link 2
linkElem <- remDr$findElement("id", "nbLink2")
linkElem$highlightElement()
out <- lapply(seq(5), function(x){linkElem$clickElement()
                                  Sys.sleep(1)
})

# text input and button
textElem <- remDr$findElement("id", "nbText")
buttonElem <- remDr$findElement("id", "nbButton")
textElem$highlightElement()
buttonElem$highlightElement()

out <- lapply(seq(5), function(x){
  textElem$clearElement()
  textElem$sendKeysToElement(list(paste("press button", x), key = 'enter'))
  buttonElem$clickElement()
  Sys.sleep(1)
})

# INTERACT WITH THE ALERTS
tabElems[[3]]$clickElement()
tabElems[[3]]$highlightElement()
Sys.sleep(1)
alertElem <- remDr$findElement("id", "alert_anchor")
out <- try(remDr$executeScript("arguments[0].scrollIntoView(true);", list(alertElem)), silent = TRUE) # scroll down if necessary
alertElem$highlightElement()

# clear checkboxs that are selected
cbElems <- remDr$findElements("css selector", paste(tabID[[3]], "label.checkbox input"))
cbLabelElems <- remDr$findElements("css selector", paste(tabID[[3]], "label.checkbox"))
out <- lapply(seq_along(cbElems), function(x){
  cbLabelElems[[x]]$highlightElement()
  enabled <- cbElems[[x]]$isElementSelected()[[1]]
  if(enabled){cbElems[[x]]$clickElement()}
})

ddElem <- remDr$findElement("css selector", paste(tabID[[3]], ".selectize-input"))
for(i in 1:4){
  # DOM gets rewritten each time so ddElems needs to be found each time
  ddElem$clickElement()
  ddElems <- remDr$findElements("css selector", paste(tabID[[3]], ".selectize-dropdown .option"))
  ddElems[[i]]$highlightElement()
  ddElems[[i]]$clickElement()
  Sys.sleep(1)
}

# INTERACT WITH THE PROGRESS BARS
tabElems[[4]]$clickElement()
tabElems[[4]]$highlightElement()
Sys.sleep(1)

pbElem <- remDr$findElement("id", "pb1")
out <- try(remDr$executeScript("arguments[0].scrollIntoView(true);", list(pbElem)), silent = TRUE) # scroll down if necessary
pbElem$highlightElement()



pbradio <- remDr$findElements("css selector", "#pbradio input")
pbradiolabel <- remDr$findElements("css selector", "#pbradio label")

lapply(c(2,3,1), function(x){
  pbradiolabel[[x+1]]$highlightElement()
  pbradio[[x]]$clickElement()
  Sys.sleep(3)
})

ddElem <- remDr$findElement("css selector", paste(tabID[[4]], ".selectize-input"))
for(i in 1:5){
  # DOM gets rewritten each time so ddElems needs to be found each time
  ddElem$clickElement()
  ddElems <- remDr$findElements("css selector", paste(tabID[[4]], ".selectize-dropdown .option"))
  ddElems[[i]]$highlightElement()
  ddElems[[i]]$clickElement()
  Sys.sleep(1)
}


remDr$close()
remDr$closeServer()
DVR$save()
DVR$closeServer()

