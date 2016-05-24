## ---- eval=FALSE, install------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("soodoku/abbyyR", build_vignettes = TRUE)

## ---- eval=FALSE, setapp-------------------------------------------------
#  setapp(c("app_id", "app_password"))

## ---- eval=FALSE, get_appinfo--------------------------------------------
#  getAppInfo()

## ---- eval=FALSE, listTasks----------------------------------------------
#  listTasks(fromDate="yyyy-mm-ddThh:mm:ssZ",toDate="yyyy-mm-ddThh:mm:ssZ")

## ---- eval=FALSE, listFinishedTasks--------------------------------------
#  listFinishedTasks()

## ---- getTaskStatus, eval=FALSE------------------------------------------
#  getTaskStatus(taskId="task_id")

## ---- deleteTask, eval=FALSE---------------------------------------------
#      deleteTask(taskId="task_id")

## ---- submitImage, eval=FALSE--------------------------------------------
#  submitImage(file_path="file_path", taskId="task_id", pdfPassword="")

## ---- processImage, eval=FALSE-------------------------------------------
#  processImage(file_path="file_path", language="English", profile="documentConversion")

## ---- processRemoteImage, eval=FALSE-------------------------------------
#  processRemoteImage(img_url="img_url", language="English", profile="documentConversion")

## ---- processDocument, eval=FALSE----------------------------------------
#  processDocument(task_id="task_id")

## ---- processBusinessCard, eval=FALSE------------------------------------
#  processBusinessCard(file_path="file_path")

## ---- processTextField, eval=FALSE---------------------------------------
#  processTextField(file_path="file_path")

## ---- processBarcodeField, eval=FALSE------------------------------------
#  processBarcodeField(file_path="file_path")

## ---- processCheckmarkField, eval=FALSE----------------------------------
#  processCheckmarkField(file_path="file_path")

## ---- processFields, eval=FALSE------------------------------------------
#  processFields(file_path="file_path")

## ---- processMRZ, eval=FALSE---------------------------------------------
#  processMRZ(file_path="file_path")

## ---- processPhotoId, eval=FALSE-----------------------------------------
#  processPhotoId(file_path="file_path")

