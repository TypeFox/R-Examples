## ---- eval=FALSE, loadlib------------------------------------------------
#  library(abbyyR)

## ---- eval=FALSE, install------------------------------------------------
#  "
#  Get the latest version from github:
#  
#  # install.packages('devtools')
#  devtools::install_github('soodoku/abbyyR')
#  "

## ---- eval=FALSE, setapp-------------------------------------------------
#  setapp(c("factbook", "7YVBc8E6xMricoTwp0mF0aH"))

## ---- eval=FALSE, getinfo------------------------------------------------
#  getAppInfo()

## ---- eval=FALSE, listasks-----------------------------------------------
#  tasklist <- listTasks()

## ---- eval=FALSE, checktasksMdelete--------------------------------------
#  listTasks(excludeDeleted="true")

## ---- eval=FALSE, listtasksdaterange-------------------------------------
#  listTasks(fromDate="2015-05-30T20:28:43Z", toDate="2015-05-31T20:28:43Z")

## ---- eval=FALSE, listfinished-------------------------------------------
#  listFinishedTasks()

## ---- eval=FALSE, gettaskstatus------------------------------------------
#  getTaskStatus(taskId="47f9b0d4-79a2-4aed-b656-2683a85ac203")

## ---- eval=FALSE, deletetask---------------------------------------------
#  #deleteTask(taskId="47f9b0d4-79a2-4aed-b656-2683a85ac203")

## ---- eval=FALSE, submittask---------------------------------------------
#  submitImage(file_path="t1.tif", pdfPassword="")

## ---- eval=FALSE, processimage-------------------------------------------
#  processImage(file_path="t1.tif")

## ---- eval=FALSE, processRemote------------------------------------------
#  processRemoteImage(img_url="https://raw.githubusercontent.com/soodoku/abbyyR/master/inst/extdata/t1.TIF")

## ---- eval=FALSE, processdoc---------------------------------------------
#  res <- listTasks()

## ---- eval=FALSE, processdoc2--------------------------------------------
#  processDocument(taskId=res$id[res$status=="Submitted"][1])

## ---- eval=FALSE, getresults---------------------------------------------
#  # getResults()

