#############################################################################################
## Project: PKgraph
## File: handler.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
##        - interface
## Notes:
#############################################################################################

################################################################################
## Project handler
################################################################################
openDataHandler <- function(h,...)
{

    gtmp.win = gwindow("Open", horizontal=FALSE)
    #size(gtmp.win) = c(5,5)
    gtgroup1 = ggroup(cont=gtmp.win, horizontal=FALSE)
    #gtgroup2 = ggroup(cont=gtmp.win, horizontal=FALSE)
    gf1 <- gframe(text = "Configure", markup = FALSE, pos = 0, horizontal=TRUE, container = gtgroup1)
    tbl <- glayout(cont=gf1)
    tbl[1,1, anchor=c(-1,-1)] = glabel("Choose file types:")
    #datatype = gradio(selected = 1, horizontal = FALSE, items = c("Data","Bootstrap result", "Outlier and influential result"))
    datatype = gdroplist(items = c("PK data","General data"))
    tbl[1,2] = datatype
    tbl[2,1, anchor=c(-1,-1)] = glabel("File format: ")

    #filetype = gradio(selected = 1, horizontal = FALSE, items = c("txt","csv"))
    filetype = gdroplist(items = c("txt","csv"))
    tbl[2,2] = filetype

    tbl[3,1, anchor=c(-1,-1)] = glabel("Data start from line ")

    #filetype = gradio(selected = 1, horizontal = FALSE, items = c("txt","csv"))
    startline = gdroplist(items = c(1:10))
    tbl[3,2] = startline
    
    tbl[4,1, anchor=c(-1,-1)] = glabel("Data separated by ")

    #filetype = gradio(selected = 1, horizontal = FALSE, items = c("txt","csv"))
    sepline = gdroplist(items = c("SPACE", "\\t", ","))
    tbl[4,2] = sepline

    tbl[5,1, anchor=c(-1,-1)] = glabel("Data has column names ")

    #filetype = gradio(selected = 1, horizontal = FALSE, items = c("txt","csv"))
    fileheader = gdroplist(items = c("TRUE","FALSE"))
    tbl[5,2] = fileheader
    
    tbl[6,1, anchor=c(-1,-1)] = glabel("_________________________")
    tbl[7,1, anchor=c(-1,-1)] = glabel("Data from R environment ")
    dataRen = gdroplist(items = c("", ls(envir=.GlobalEnv)))
    tbl[7,2] = dataRen

    gb1 = gbutton(text="Click to open", horizontal=FALSE )
    addhandlerclicked(gb1, function(h,...)
      {
        if (svalue(dataRen) != "")
        {    
             tmp.data <- eval(parse(text= svalue(dataRen)), envir=.GlobalEnv)
             if ((!is.data.frame(tmp.data)) && (!is.matrix(tmp.data)))
             {
                 ErrorMessage("Data is neither DATA.FRAME nor MATRIX")
                 return(NULL)             
             }
                  myType <- svalue(datatype)  # 0603
                  if ( myType == "PK data") data.config <- 1     # 0603
                  else data.config <- 0
                  
                  #.pk$setDatasets(tmp.data, thisDatano)
                  #.pk$setCurrentDatano(thisDatano)
                   #browser()
                  ## use number as data name
                  filename <- svalue(dataRen)
                  filename <- filename[length(filename)]
                  
                  #thisDataName <- paste(getTotalDataLen() + 1, "_", svalue(datatype), sep="")
                  thisDataName <- paste(getTotalDataLen() + 1, "_", filename, sep="")

                  setDatasets(tmp.data, thisDataName) # use no as data name
                  setCurrentDataType(svalue(datatype), thisDataName)

                  dispose(gtmp.win)

                  ptable=gtable(tmp.data, multiple=TRUE, expand=TRUE)
                  pkmain.add(ptable, as.character(thisDataName), override.closebutton = TRUE)

                  ## for command area
                  #attach(tmp.data)

                  ## configure data


                  ## setup status
                  svalue(pmg.statusBar) <- "Data is loaded successfully."

                  if (myType == "PK data")    # 0603
                  {
                      ggobiPlotType()
                                      
                  }             
        }
        else
        {
            gfile("Select a file",type="open",
                #action="read.csv",
                handler = function(h,...)
                {
                  if (svalue(sepline) == "SPACE") svalue(sepline) <- ""
                  if (svalue(filetype) == "txt")
                  {
                     
                     #tmp.data <<- do.call("read.table",list(h$file, header=as.logical(svalue(fileheader)), sep=svalue(sepline),
                                        #skip=as.numeric(svalue(startline))-1))
                     test <<- try(do.call("read.table",list(h$file, header=as.logical(svalue(fileheader)), sep=svalue(sepline),
                                        skip=as.numeric(svalue(startline))-1)))

                     if (inherits(test, "try-error")) 
                     {
                        ErrorMessage("Data can NOT be read in! Please check sample data using data(pkdata)")
                        return(NULL)
                     }
                     else
                     {
                        tmp.data <<- test
                     }                                        
                  }
                  else
                  {
                      #tmp.data <<- do.call("read.csv",list(h$file, header=as.logical(svalue(fileheader)),
                      #                  skip=as.numeric(svalue(startline))-1))
                      test <<- try(do.call("read.csv",list(h$file, header=as.logical(svalue(fileheader)),
                                        skip=as.numeric(svalue(startline))-1)))
                     if (inherits(test, "try-error")) 
                     {
                        ErrorMessage("Data can NOT be read in! Please check sample data using data(pkdata)")
                        return(NULL)
                     }
                     else
                     {
                        tmp.data <<- test
                     }                      
                  }
                  
                  myType <- svalue(datatype)  
                  if ( myType == "PK data") data.config <- 1     
                  else data.config <- 0
                  

                  ## use number as data name
                  filename <- unlist(strsplit(h$file, "\\\\"))
                  filename <- filename[length(filename)]
                  
                  #thisDataName <- paste(getTotalDataLen() + 1, "_", svalue(datatype), sep="")
                  thisDataName <- paste(getTotalDataLen() + 1, "_", filename, sep="")

                  setDatasets(tmp.data, thisDataName) # use no as data name
                  setCurrentDataType(svalue(datatype), thisDataName)

                  dispose(gtmp.win)

                  ptable=gtable(tmp.data, multiple=TRUE, expand=TRUE)
                  pkmain.add(ptable, as.character(thisDataName), override.closebutton = TRUE)


                  ## setup status
                  svalue(pmg.statusBar) <- "Data is loaded successfully."

                  if (myType == "PK data")    # 0603
                  {
                      ggobiPlotType()
                  }

                })
            }
         })

    tbl[8,2] = gb1
   
}

cleanDataHandler <- function(h,...)
{
    cleanAll()
    if (length(pmg.dialog.notebook) > 0)
        for (i in 1: length(pmg.dialog.notebook)) dispose(pmg.dialog.notebook)
    else  ErrorMessage("No data to clean!")

}

setHandler<-function(h,...)
{    
    gfile("Select a directory", "selectdir", action=("setwd"),
        handler = function(h,...)
        {  #browser()
           do.call(h$action,list(h$file))
          svalue(pk.dirname) <- paste("Current directory: ", h$file)
          pk.dir[,] <- data.frame(dir())

        })
}


saveHandler <- function(h,...)
{
   # check data exist
    if(!checkDataExist())
    {
        ErrorMessage("No data is available!")
        return(invisible(NULL))
    }
    
    gfile("Save currrent data", type="save",
              action="dput", handler = function(h,...)
              {
                   currentMain <- svalue(pmg.dialog.notebook)
                   currentData <- getCurrentData(currentMain)
                   if (length(grep(".txt", h$file)) == 0) 
                   {
                      filename <- paste(h$file, ".txt", sep="")
                   }
                   else
                   {
                      filename <- h$file
                   }
                   write.table(currentData, file=filename, sep="\t")
                   
                   svalue(pmg.statusBar) <- "File has been saved successfully."
              })
}
saveProjectHandler <- function(h,...)
{
   # check data exist
    if(!checkDataExist())
    {
        ErrorMessage("No data is available!")
        return(invisible(NULL))
    }
    
    gfile("Save project", type="save",
              handler = function(h,...)
              {
                dir.create(h$file)
                old.dir <- getwd()
                setwd(h$file)
                on.exit(setwd(old.dir))
                
                ## save data set
                pkdata <- getDatasets()
                dput(pkdata, "pkgraphData.txt")

                ## get saving format
                saveFormat <- getSaveFormat()
                dput(saveFormat, "pkgraphSaveFormat.txt")

                ## get figure config
                figConfig <- getFigConfig()
                dput(figConfig, "pkgraphFigConfig.txt")
                
                svalue(pmg.statusBar) <- "Project has been saved successfully."
                
              }
        )
}

restoreHandler <- function(h,...)
{
    gfile("Save project", type="selectdir",
              handler = function(h,...)
         {
          file.list <- dir(h$file)
          require.file <- c("pkgraphData.txt","pkgraphSaveFormat.txt","pkgraphFigConfig.txt")
          if (all(require.file %in% file.list))
          {
               old.dir <- getwd()
               setwd(h$file)
               on.exit(setwd(old.dir))

               ## get data set
               pkdata <- dget("pkgraphData.txt")
               pkname <- names(pkdata)
               sapply(1:length(pkname), function(i)
                      {
                          setDatasets(pkdata[[pkname[i]]], pkname[i])

                          ptable=gtable(pkdata[[pkname[i]]], multiple=TRUE, expand=TRUE)
                          pkmain.add(ptable, as.character(pkname[i]), override.closebutton = TRUE)
                      })


               ## get saving format
               setSaveFormat(dget("pkgraphSaveFormat.txt"))

               ## get figure config
               setFigConfig(dget("pkgraphFigConfig.txt"))


                      svalue(pmg.statusBar) <- "Project has been saved successfully."
          }
          else
          {
              ErrorMessage("You need all project files: pkgraphData.txt, pkgraphSaveFormat.txt, pkgraphFigConfig.txt")
              return(invisible(NULL))
          }
      })
    
}

helpHandler<-function(h,...)
{
    print(vignette("PKgraph", package="PKgraph", all=TRUE))
}

exitHandler<-function(h,...)
{
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
    cleanPKCode()
    cleanPKGGobi()
    dispose(PKW)
}


################################################################################
## Config handler
################################################################################
configDirHandler <- function(h,...)
{   
    gfile("Select a directory", "selectdir", action=("setwd"),
        handler = function(h,...)
        {  #browser()
           do.call(h$action,list(h$file))
          svalue(pk.dirname) <- paste("Current directory: ", h$file)
          
          pk.dir[] <- dir()

          ## setup status bar
          svalue(pmg.statusBar) <- "Working directory is setup successfully."
        })
}

configDataHandler <- function(h,...)
{
    if(checkDataExist()) pk.dataConfig()
    else gmessage("No data is available for configuration!", icon=c("warning"), title="Warning")
    return(invisible(NULL))
}

configFormatHandler <- function(h,...)
{
    format1 <- c("bmp", "jpeg")
    format2 <- c("png", "tiff")
    format3 <- c("pdf", "win.metafile")

    gconfigsave.win = gwindow("Set saving formats for figures", horizontal=FALSE)

    gtgroup1 = ggroup(cont=gconfigsave.win, horizontal=FALSE)

    gf1 <- gframe(text = "Configure", markup = FALSE, pos = 0, horizontal=TRUE, container = gtgroup1)
    tbl <- glayout(cont=gf1)
    tbl[1,1, anchor=c(-1,-1)] = glabel("Saving format for figures:")

    f1type = gcheckboxgroup(format1, horizontal=TRUE)
    tbl[1,2] = f1type

    f2type = gcheckboxgroup(format2, horizontal=TRUE)
    tbl[2,2] = f2type
    
    f3type = gcheckboxgroup(format3, horizontal=TRUE)
    svalue(f3type) = "pdf"
    tbl[3,2] = f3type
    
    tbl[4,1, anchor=c(-1,-1)] = glabel("Figure width:")
    fig.width = gedit(text = "default")
    tbl[4,2] = fig.width
    
    tbl[5,1, anchor=c(-1,-1)] = glabel("Figure height:")
    fig.height = gedit(text = "default")
    tbl[5,2] = fig.height
    
    gb1 = gbutton(text="Click to configure", horizontal=FALSE )
    tbl[6,2] = gb1
    
    addhandlerclicked(gb1, function(h,...)
    {
        save.command <- c(svalue(f1type), svalue(f2type), svalue(f3type))
        if(length(save.command)==0)
        {
            ErrorMessage("Please choose saving format for figures!")
            return(invisible(NULL))
        }
        if (svalue(fig.width)=="default") mywidth <- NA
        else mywidth <- as.numeric(svalue(fig.width))
        
        if (svalue(fig.height)=="default") myheight <- NA
        else myheight <- as.numeric(svalue(fig.height))
        
        newformat <- list(command=save.command, width=mywidth,height=myheight)
        
        setSaveFormat(newformat)
        dispose(gconfigsave.win)
        svalue(pmg.statusBar) <- "Set saving format for figures SUCCESSFULLY"
    })

}

configFigureHandler <- function(h,...)
{
    #format1 <- c("red", "blue", "green", "pink", ")
    format2 <- c("identity line", "loess")


    gfigconfig.win = gwindow("Set figure configuration", horizontal=FALSE)

    gtgroup1 = ggroup(cont=gfigconfig.win, horizontal=FALSE)

    gf1 <- gframe(text = "General configure", markup = FALSE, pos = 0, horizontal=TRUE, container = gtgroup1)
    tbl <- glayout(cont=gf1)
    tbl[1,1, anchor=c(-1,-1)] = glabel("color for figures:")

    f1type = gedit(text = "royalblue")
    tbl[1,2] = f1type

    tbl[2,1, anchor=c(-1,-1)] = glabel("model diagnostics:")

    f2type = gcheckboxgroup(format2, horizontal=TRUE)
    tbl[2,2] = f2type

    gb1 = gbutton(text="Click to configure", horizontal=FALSE )
    tbl[3,2] = gb1

    addhandlerclicked(gb1, function(h,...)
    {
        figConfig <- list()
        figConfig$col <- svalue(f1type)
        if (length(svalue(f2type)) > 0 )
        {
            if ("identity line" %in% svalue(f2type)) figConfig$identity <- 1
            if ("loess" %in% svalue(f2type)) figConfig$loess <- 1
        }
        
        setFigConfig(figConfig)
        dispose(gfigconfig.win)
        svalue(pmg.statusBar) <- "Set saving format for figures SUCCESSFULLY"
    })

}

################################################################################
## PKmodel handler
################################################################################






