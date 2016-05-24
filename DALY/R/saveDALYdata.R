## Save 'DALY' database data to RData file
## NOTE: the current structure will be transformed to S4 class

saveDALYdata <-
function(){
  ## Define 'fileName' - default is "DALY_<diseaseName>.RData"
  fileName <-
    tclvalue(tkgetSaveFile(initialfile =
               paste("DALY_",
                     DALYtclvalue("diseaseName"),
                     ".RData", sep = ""),
               filetypes = "{{R Images} {.RData}}"))

  if(fileName != ""){
    ## Add ".RData" suffix if missing
    if (substr(fileName, nchar(fileName) - 4, nchar(fileName)) != "RData")
      fileName = paste(fileName, ".RData", sep = "")

    ## Save 'disease model'
    outcomeNames <- list()
    for (i in seq(8))
      outcomeNames[[i]] <- DALYtclvalue(paste("outcome", i, "Name", sep = ""))
    model <- list(diseaseName = DALYtclvalue("diseaseName"),
                  outcomeNames = outcomeNames)

    ## Save 'settings'
    settings <- list(pop = DALYget("pop"),
                     LE = DALYget("LE"),
                     aw = DALYtclvalue(".aw"),
                     dr = DALYtclvalue(".dr"))

    ## Save 'data'
    data <- list()
    for (i in seq(8)){  # 8 outcomes
      data[[i]] <- list()
      for (j in seq(8)){  # 8 parameters per outcome
        data[[i]][[j]] <-
          list(DALYget(paste("dist", DALYget("txtLbl")[j], i, sep = "")),
               DALYget(paste("strat", DALYget("txtLbl")[j], i, sep = "")),
               DALYget(paste(DALYget("txtlbl")[j], i, sep = "")))
        names(data[[i]][[j]]) <- c("dist", "strat", "param")
      }
      names(data[[i]]) <- DALYget("txtlbl")
    }

    ## Save all to file
    DALYassign("DALY_data",
               list(model = model, settings = settings, data = data))
    save("DALY_data", file = fileName,
         envir = DALYenv())
	
	## Cleanup 'DALY' database
    rm(list = "DALY_data", envir = DALYenv())

    ## Exit message
    tkmessageBox(title = "DALY calculator",
                 message = paste("Data successfully saved to\n\"",
                                 fileName, "\"", sep = ""))
  }
}
