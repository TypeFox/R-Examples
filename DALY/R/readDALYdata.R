## Open DALY RData file
## Load data into DALY Calculator

readDALYdata <-
function(file = NULL, example = NULL){
  ## Ask 'fileName' to user
  fileName <- file
  if (is.null(file) & is.null(example))
    fileName <- tclvalue(tkgetOpenFile(filetypes = "{{R Images} {.RData}}"))

  if (!is.null(fileName) | !is.null(example)){
  
    if (is.null(example)){
      ## Evaluate 'fileName': should end with '.RData'
      if (!grepl(".RData$", fileName, ignore.case = TRUE)){
        tkmessageBox(message = paste("The file you selected",
                                     "is not an '.RData' file"),
                     title = "DALY calculator",
                     icon = "error")
        stop("The file you selected is not an '.RData' file", call. = FALSE)
      }

      ## Try loading data into current frame
      tryCatch(load(fileName, envir = DALYenv()),
               silent = TRUE,
               error =
                 function(e){
                   tkmessageBox(message = "Error while loading file",
                                title = "DALY calculator", icon = "error")
                   stop(paste("Error while loading file\n      ",
                        "The file may be corrupt or of unsupported format"),
                        call. = FALSE)
                  })
      DALY_name <- "DALY_data"

    } else {
      ## Load example dataset
      examples <- c("Neurocysticercosis", "Toxoplasmosis")
      DALY_name <- paste("DALY_", examples[example], sep = "")
      data(list = DALY_name, envir = DALYenv())
    }

    ## Reset DALY Calculator
    reset()

    ## Retrieve data
    model <- DALYget(DALY_name)$model
    settings <- DALYget(DALY_name)$settings
    data <- DALYget(DALY_name)$data

    ## Update 'disease model'
    DALYupdate("diseaseName", model$diseaseName)
    for (i in seq(8))
      DALYupdate(paste("outcome", i, "Name", sep = ""),
                 model$outcomeNames[[i]])

    ## Update 'settings'
    DALYassign("pop", settings$pop); DALYupdate(".pop")
    DALYassign("LE", settings$LE); DALYupdate(".LE")
    DALYupdate(".aw", settings$aw)
    DALYupdate(".dr", settings$dr)

    ## Update 'data'
    for (i in seq(8)){  # 8 outcomes
      for (j in seq(8)){  # 8 parameters per outcome
        DALYassign(paste("dist", DALYget("txtLbl")[j], i, sep = ""),
                   data[[i]][[j]]$dist)
        DALYassign(paste("strat", DALYget("txtLbl")[j], i, sep = ""),
                   data[[i]][[j]]$strat)
        DALYassign(paste(DALYget("txtlbl")[j], i, sep = ""),
                   data[[i]][[j]]$param)
        DALYupdate(paste(".dist", DALYget("txtLbl")[j], i, sep = ""))
        DALYupdate(paste(".strat", DALYget("txtLbl")[j], i, sep = ""))
        DALYupdate(paste(".", DALYget("txtlbl")[j], i, sep = ""))
      }
    }
	
    ## Clean-up 'DALY' database
    rm(list = DALY_name,
       envir = DALYenv())

    ## Status message
    tkmessageBox(message = "Data successfully loaded",
                 title = "DALY calculator")
  }
}