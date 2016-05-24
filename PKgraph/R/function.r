#############################################################################################
## Project: PKgraph
## File: function.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
##        - interface
## Notes:
##     For ggobiImageHandler
#         - name can not be "psn.outlier.ggobiimageHandler", only works with "outlier.imagehandler"
#         -     reason: only one dot is allowed?
#############################################################################################

pkmain.add = function(widget, label,...) {
  add(pmg.dialog.notebook, widget, label=label,...) # add near beginning
}


################################################################################
## Menu "Project"
################################################################################
project.restore <- function()
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
                          pkmain.add(pkdata[[pkname[i]]], as.character(pkname[i]), override.closebutton = TRUE)
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

pk.dataConfig <- function()
{
   myterm <- c("No match", "ID", "DV", "TIME", "PRED", "RES", "WRES",
                 "IPRE", "IDV", "COV", "ETA",
                  "PARAMETERS")
   currentMain <- svalue(pmg.dialog.notebook) 
   myNames <- colnames(getCurrentData(currentMain))
   mylist <- list()
   knownterm <- getTerm()
   if (nrow(knownterm) > 0)
   {
       sapply(1:length(myNames), function(i)
              {
                  myvar.ind <- which(knownterm$VarName == myNames[i])
                  defaultTerm <- knownterm[myvar.ind,]$TermName
                  myterm.ind <- which(myterm==defaultTerm)
                  mylist[[myNames[i]]] <<- list(type = "gdroplist", 
                          items =c(defaultTerm,myterm[-c(myterm.ind)]), expand=TRUE)
              })
   }
   else
   {
       sapply(1:length(myNames), function(i)
                  mylist[[myNames[i]]] <<- list(type = "gdroplist", 
                          items =myterm, expand=TRUE))   
   }
   

   
   configData = SimpleGUI$new(message="PK Data Configuration",
        widgetList = mylist
    )
    
    configData$okButtonHandler = simple.okButtonHandler
    configData$show()
}
################################################################################
## Menu "Data mangement"
################################################################################
pk.subset <- function()
{
  # check data first
  if(!checkDataExist())
  {
      ErrorMessage("No data is available!")
      return(invisible(NULL))
  }
  
  TableTest = TableGUI$new(message="Subset", frameMessage="Subset variables", buttonMessage="Subset")
  currentPage <- svalue(pmg.dialog.notebook)
  TableTest$table.data = colnames(getCurrentData(currentPage))
  TableTest$okButtonHandler = subset.okButtonHandler
  TableTest$show()


}


pk.factor <- function()
{
  # check data first
  if(!checkDataExist())
  {
      gmessage("No data is available!",
              icon=c("warning"), title="Warning")
      return(invisible(NULL))
  }

  TableTest = TableGUI$new(message="Factor/Unfactor", frameMessage="Factor/UnFactor variables", buttonMessage="Factor/Unfactor")
  currentPage <- svalue(pmg.dialog.notebook)
  TableTest$table.data <- colnames(getCurrentData(currentPage))
  TableTest$checkbox.data <- c("factor", "numeric", "character")

  TableTest$okButtonHandler = factor.okButtonHandler
  TableTest$show()
}

################################################################################
## Menu "Summary"
################################################################################

pk.uni.menu <- function()
{
    # check data exist
    if(!checkDataExist())
    {
        ErrorMessage("No data is available for configuration!")
        return(invisible(NULL))
    }
    
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    cleanDataSpecialPlot()

    currentMain <- svalue(pmg.dialog.notebook)  

   BGTest = BasicGUI$new(message="Univarate plot",
        widgetList = list(
     x = list(type = "gdroplist", items = colnames(getCurrentData(currentMain))),
     "number of bins" = list(type="gedit",text=""),
     main = list(type="gedit",text=" "),
     xlab = list(type="gedit",text=" "),
     ylab = list(type="gedit",text=" "),
     type = list(type = "gdroplist", items = c("percent", "count", "density"))

    ),
       saveList = list(
     cond = list(type = "gdroplist", items = c("", colnames(getCurrentData(currentMain)))),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="1") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.uni.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler      
    BGTest$ggobiImageHandler = summary.uni.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()
}
#######################
## bivariates
#######################
pk.bi.menu <- function()
{
    if(!checkDataExist())
    {
        ErrorMessage("No data is available for configuration!")
        return(invisible(NULL))
    }
    
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    cleanDataSpecialPlot()

    currentMain <- svalue(pmg.dialog.notebook)  
    
   BGTest = BasicGUI$new(message="Bivarate plots",
        widgetList = list(
     x = list(type = "gdroplist", items = colnames(getCurrentData(currentMain))),
     y = list(type = "gdroplist", items = colnames(getCurrentData(currentMain))),
     main = list(type="gedit",text=" "),
     xlab = list(type="gedit",text=" "),
     ylab = list(type="gedit",text=" "),
     type = list(type = "gdroplist", items = c("p", "l", "b", "ts")) #,
     
    ),
       saveList = list(
     cond = list(type = "gdroplist", items = c("", colnames(getCurrentData(currentMain)))),
     layout_x = list(type="gedit",text="5"),
     layout_y = list(type="gedit", text="5") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.bi.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = summary.bi.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()
}
#######################
## trivariates
#######################
pk.tri.menu <- function()
{
    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

    tmp.name <- colnames(getCurrentData())
   BGTest = BasicGUI$new(message="Trivariate plots",
        widgetList = list(
     x = list(type = "gdroplist", items = tmp.name ),
     y = list(type = "gdroplist", items = tmp.name),
     z = list(type = "gdroplist", items = tmp.name),
     main = list(type="gedit",text=" "),
     xlab = list(type="gedit",text=" "),
     ylab = list(type="gedit",text=" "),
     zlab = list(type="gedit",text=" "),
     type = list(type = "gdroplist", items = c("h", "p", "l", "b"))

    ),
       saveList = list(
     "conditional var" = list(type = "gdroplist", items = c("", colnames(.pk$getCurrentData()))),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="1"),
     graphics = list(type="gradio", items = c("lattice", "ggplot2(No support for 3D)"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.tri.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = summary.tri.ggobiImageHandler
    BGTest$show()
}

pk.time.menu <- function()
{
    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    
   BGTest = BasicGUI$new(message="Time plot",
        widgetList = list(
     x = list(type = "gdroplist", items = colnames(getCurrentData())),
     y = list(type = "gdroplist", items = colnames(getCurrentData())),
     main = list(type="gedit",text=" "),
     xlab = list(type="gedit",text=" "),
     ylab = list(type="gedit",text=" "),
     type = list(type = "gdroplist", items = c("p", "l", "b")),

     xlim_start = list(type="gedit",text=""),
     xlim_end = list(type="gedit",text=""),
     ylim_start = list(type="gedit",text=""),
     ylim_end = list(type="gedit",text="")

    ),
       saveList = list(
     cond = list(type = "gdroplist", items = c("", colnames(.pk$getCurrentData()))),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="1") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.bi.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = summary.bi.ggobiImageHandler
    BGTest$show()
}
pk.para.menu <- function()
{
    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()

    currentMain <- svalue(pmg.dialog.notebook) 
    
   BGTest = BasicGUI$new(message="Parallel coordinates plot",
        widgetList = list(

     x=list(type="gcheckboxgroup", items=colnames(getCurrentData(currentMain))),
     main = list(type="gedit",text=" "),
     horizontal = list(type="gradio",items=c("TRUE", "FALSE"), horizontal=TRUE)

    ),
       saveList = list(
     "conditional var" = list(type = "gdroplist", items = c("", colnames(getCurrentData(currentMain)))),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="1") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2(Not support)"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.para.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = summary.para.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()
}

pk.heat.menu <- function()
{

    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    
   BGTest = BasicGUI$new(message="Heatmap",
        widgetList = list(

     x=list(type="gcheckboxgroup", items=colnames(getCurrentData())),
     main = list(type="gedit",text=" "),
     "scale by"=list(type="gradio", items = c("column", "row"), horizontal=TRUE),
     "dendrogram for row"=list(type="gradio", items = c( "no", "yes"), horizontal=TRUE),
     "dendrogram for column"=list(type="gradio", items = c( "no", "yes"), horizontal=TRUE)
    ),
       saveList = list(

    )
    )
    BGTest$okButtonHandler = summary.heat.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$saveImageHandler = saveImageHandler
    #BGTest$ggobiImageHandler = summary.heat.ggobiImageHandler
    BGTest$show()
}

pk.matrix.menu <- function()
{

    if(!checkDataExist())
    {
        gmessage("No data is available for configuration!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }
    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()

    currentMain <- svalue(pmg.dialog.notebook) 

   BGTest = BasicGUI$new(message="Scatterplot matrix",
        widgetList = list(

     x=list(type="gcheckboxgroup", items=colnames(getCurrentData(currentMain))),
     main = list(type="gedit",text=" ")
    ),
       saveList = list(

     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    )
    )
    BGTest$okButtonHandler = summary.matrix.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler      
        
    BGTest$datagroup = 0
    BGTest$saveImageHandler = saveImageHandler.matrix
    BGTest$ggobiImageHandler = summary.matrix.ggobiImageHandler
    BGTest$show()
    
    cleanDataSpecialPlot()
}
################################################################################
## Menu "PK model"
################################################################################
pk.model.ind <- function()
{
    if (!checkforPKModel()) return(invisible(NULL))

    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()

    currentMain <- svalue(pmg.dialog.notebook) 
    
## specific check
   ## check var type for this function, particularly
    match.term <- getTerm()
    require.term <- c("ID")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData(currentMain))
   id.ind <- which(require.term==all.items)
   ind.items <- all.items[-id.ind]

   if (!config.check)
   {
       BGTest = BasicGUI$new(message="PK individual plot",
            widgetList = list(
         x = list(type = "gdroplist", items = c("", ind.items)),
         y = list(type = "gdroplist", items = c("", ind.items)),
         main = list(type="gedit",text=""),
         xlab = list(type="gedit",text=""),
         ylab = list(type="gedit",text=""),
         type = list(type = "gdroplist", items = c("p", "l", "b"))
         ),

            saveList = list(
         layout_x = list(type="gedit",text="5"),
         layout_y = list(type="gedit", text="5") ,
         graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                            )

        )
        BGTest$okButtonHandler = model.ind.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$figureGroup = 1
        BGTest$saveImageHandler = saveImageHandler
       BGTest$show()
    }
}

pk.model.gof <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))

    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    
   ## specific check
    match.term <- getTerm()
    require.term <- c("PRED", "IPRE", "DV", "IDV", "WRES")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData())
   pred.items <- match.term[match.term$TermName=="PRED",]$VarName
   ipre.items <- match.term[match.term$TermName=="IPRE",]$VarName
   dv.items <- match.term[match.term$TermName=="DV",]$VarName
   idv.items <- match.term[match.term$TermName=="IDV",]$VarName
   wres.items <- match.term[match.term$TermName=="WRES",]$VarName
   
   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Goodness of fit plots",
            widgetList = list(
         # |IWRES¦ vs IPRED, DV vs PRED/IPRE, Predictions vs IDV
         "DV vs PRED" = list(type="glabel", text=""),
         "DV_1" = list(type="gdroplist", items = c(dv.items)),
         "PRED_1" = list(type="gdroplist", items = c(pred.items)),
         "_________________________" = list(type="glabel", text=""),
         "DV vs IPRE"  = list(type="glabel", text=""),
         "DV_2" = list(type="gdroplist", items = c(dv.items)),
         "IPRE" = list(type="gdroplist", items = c(ipre.items)),
         "_________________________" = list(type="glabel", text=""),
         "WRES vs IDV"  = list(type="glabel", text=""),
         "WRES_3" = list(type="gdroplist", items = c(wres.items)),
         "IDV_3" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "PRED vs IDV" = list(type="glabel", text=""),
         "PRED_4" = list(type="gdroplist", items = c(pred.items)),
         "IDV_4" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "IPRE vs IDV" = list(type="glabel", text=""),
         "IPRE_5" = list(type="gdroplist", items = c(ipre.items)),
         "IDV_5" = list(type="gdroplist", items = c(idv.items))

        ),
       saveList = list(
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="1") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    ))
        BGTest$okButtonHandler = model.gof.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler
       BGTest$show()
   }
}

pk.model.struct <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))

    # clean code from other window
    cleanPKCode()
    cleanPKGGobi()
    
   ## specific check
    match.term <- getTerm()
    require.term <- c("PRED", "IPRE", "DV", "IDV", "WRES", "COV")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData())
   pred.items <- match.term[match.term$TermName=="PRED",]$VarName
   ipre.items <- match.term[match.term$TermName=="IPRE",]$VarName
   dv.items <- match.term[match.term$TermName=="DV",]$VarName
   idv.items <- match.term[match.term$TermName=="IDV",]$VarName
   wres.items <- match.term[match.term$TermName=="WRES",]$VarName
   cov.items <- match.term[match.term$TermName=="COV",]$VarName

   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Structural model diagnostics",
            widgetList = list(
         "PRED vs DV|IDV" = list(type="glabel", text=""),
         "PRED_1" = list(type="gdroplist", items = c(pred.items)),
         "DV_1" = list(type="gdroplist", items = c(dv.items)),
         "IDV_1" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "IPRE vs DV|IDV"  = list(type="glabel", text=""),
         "IPRE" = list(type="gdroplist", items = c(ipre.items)),
         "DV_2" = list(type="gdroplist", items = c(dv.items)),
         "IDV_2" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "WRES vs IDV" = list(type="glabel", text=""),
         "WRES_3" = list(type="gdroplist", items = c(wres.items)),
         "IDV_3" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "WRES vs IDV(bw)" = list(type="glabel", text=""),
         "WRES_4" = list(type="gdroplist", items = c(wres.items)),
         "IDV_4" = list(type="gdroplist", items = c(idv.items)),
         "_________________________" = list(type="glabel", text=""),
         "WRES vs PRED" = list(type="glabel", text=""),
         "WRES_5" = list(type="gdroplist", items = c(wres.items)),
         "PRED_5" = list(type="gdroplist", items = c(pred.items)),
         "_________________________" = list(type="glabel", text=""),
         "WRES vs PRED(bw)" = list(type="glabel", text=""),
         "WRES_6" = list(type="gdroplist", items = c(wres.items)),
         "PRED_6" = list(type="gdroplist", items = c(pred.items)),
         "_________________________" = list(type="glabel", text=""),
         "PRED vs DV|Covariates" = list(type="glabel", text=""),
         "PRED_7" = list(type="gdroplist", items = c(pred.items)),
         "DV_7" = list(type="gdroplist", items = c(dv.items)),
         "COV_7" = list(type="gdroplist", items = c(cov.items)),
         "_________________________" = list(type="glabel", text=""),
         "IPRED vs DV|Covariates" = list(type="glabel", text=""),
         "IPRE_8" = list(type="gdroplist", items = c(ipre.items)),
         "DV_8" = list(type="gdroplist", items = c(dv.items)),
         "COV_8" = list(type="gdroplist", items = c(cov.items))

        ),
       saveList = list(
     layout_x = list(type="gedit",text="5"),
     layout_y = list(type="gedit", text="5") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    ))
        BGTest$okButtonHandler = model.struct.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler
       BGTest$show()
   }
}

pk.model.resid <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))
    # clean code from other window
    cleanPKCode()
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
    cleanPKGGobi()
    
   ## specific check
    match.term <- getTerm()
    require.term <- c("WRES", "PRED", "COV", "IPRE")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData())
   wres.items <- match.term[match.term$TermName=="WRES",]$VarName
   pred.items <- match.term[match.term$TermName=="PRED",]$VarName
   cov.items <- match.term[match.term$TermName=="COV",]$VarName
   
   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Residual model diagnostics",
            widgetList = list(
         "Distribution of WRES:" = list(type="gdroplist", items = c(wres.items)),
         "_________________________" = list(type="glabel", text=""),
         "Distribution of WRES(QQ):"  = list(type="gdroplist", items = c(wres.items)),
         "_________________________" = list(type="glabel", text=""),
         "Individual distribution" = list(type="glabel", text=""),
         " of WRES:" = list(type="gdroplist", items = c(wres.items)),
         "_________________________" = list(type="glabel", text=""),
         "Individual distribution:" = list(type="glabel", text=""),
         " of WRES(QQ)" = list(type="gdroplist", items = c(wres.items)),
         "_________________________" = list(type="glabel", text=""),
         "|WRES|_1 vs PRED:" = list(type="glabel", text=""),
         "|WRES|_1" = list(type="gdroplist", items = c(wres.items)),
         "PRED_1" = list(type="gdroplist", items = c(pred.items)),
         "_________________________" = list(type="glabel", text=""),
         "Covariates vs |WRES|:" = list(type="glabel", text=""),
         "Covariates_2" = list(type="gdroplist", items = c(cov.items)),
         "|WRES|_2" = list(type="gdroplist", items = c(wres.items)),
         
         
         "_________________________" = list(type="glabel", text=""),
         "|WRES|_3 vs PRED|Covariates:" = list(type="glabel", text=""),
         "|WRES|_3" = list(type="gdroplist", items = c(wres.items)),
         "PRED_3" = list(type="gdroplist", items = c(pred.items)),
         "Covariates_3" = list(type="gdroplist", items = c(cov.items)),
         "_________________________" = list(type="glabel", text=""),
         "|IWRES| vs IPRE|Covariates:" = list(type="glabel", text=""),
         "|IWRES|" = list(type="gdroplist", items = c(wres.items)),
         "IPRE" = list(type="gdroplist", items = c(pred.items)),
         "Covariates_4" = list(type="gdroplist", items = c(cov.items)),
         "_________________________" = list(type="glabel", text=""),
         "Autocorrelation of WRES:" = list(type="gdroplist", items = c(wres.items))

        ),
       saveList = list(
     layout_x = list(type="gedit",text="5"),
     layout_y = list(type="gedit", text="5") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    ))
        BGTest$okButtonHandler = model.resid.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler        
       BGTest$show()
   }
   
   ## for abs value,
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
}
pk.model.para <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))
  
  cleanPKCode()
  cleanPKGGobi()
  
   ## specific check
    match.term <- getTerm()
    require.term <- c("PARAMETERS")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData())
   para.items <- match.term[match.term$TermName=="PARAMETERS",]$VarName

   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Parameter",
            widgetList = list(
         "Distribution of parameters:" = list(type="gdroplist", items = c(para.items)),
         "________________________________" = list(type="glabel", text=""),
         "Distribution of parameters (QQ):"  = list(type="gdroplist", items = c(para.items)),
         "________________________________" = list(type="glabel", text=""),
         "Scatterplot matrix of parameters:" = list(type="glabel", text=""),
         "________________________________" = list(type="glabel", text=""),
         "Parameter vs parameter:" = list(type="glabel", text=""),
         "Parameters_x:" = list(type="gdroplist", items = c(para.items)),
         "Parameters_y:" = list(type="gdroplist", items = c(para.items))

         #slider = list(type = "gslider", value = 10),
         # radio = list(type="gradio", items = 1:3, horizontal=FALSE)
        ),
       saveList = list(
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    ))
        BGTest$okButtonHandler = model.para.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler.pkmodel        
       BGTest$show()
   }
}
pk.model.cov <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))

  cleanPKCode()
  cleanPKGGobi()
  
   ## specific check
    match.term <- getTerm()
    require.term <- c("PARAMETERS", "ETA", "WRES", "COV")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term

    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")

   all.items <- colnames(getCurrentData())
   para.items <- match.term[match.term$TermName=="PARAMETERS",]$VarName
   cov.items <- match.term[match.term$TermName=="COV",]$VarName
   eta.items <- match.term[match.term$TermName=="ETA",]$VarName
   wres.items <- match.term[match.term$TermName=="WRES",]$VarName

   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Covariate model diagnostics",
            widgetList = list(
         "Scatterplot matrix of covariates" = list(type="glabel", text=""),
         "________________________________" = list(type="glabel", text=""),
         "Parameters vs covariates" = list(type="glabel", text=""),
         "Parameters:" = list(type="gdroplist", items = c(para.items)),
         "Cov_P:" = list(type="gdroplist", items = c(cov.items)),
         "________________________________" = list(type="glabel", text=""),
         "ETAs vs covariates" = list(type="glabel", text=""),
         "ETAS:" = list(type="gdroplist", items = c(eta.items)),
         "Cov_E:" = list(type="gdroplist", items = c(cov.items)),
         "________________________________" = list(type="glabel", text=""),
         "WRES vs covariates" = list(type="glabel", text=""),
         "WRES:" = list(type="gdroplist", items = c(wres.items)),
         "Cov_W:" = list(type="gdroplist", items = c(cov.items))

        ),
       saveList = list(
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)

    ))
        BGTest$okButtonHandler = model.cov.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler.pkmodel                
       BGTest$show()
   }
}

pk.model.random <- function()
{
## GENEARL CHECK
  if (!checkforPKModel()) return(invisible(NULL))

  cleanPKCode()
  cleanPKGGobi()
  
## specific check
   ## check var type for this function, particularly
    match.term <- getTerm()
    require.term <- c("ETA")
    PK.match <- match(require.term, match.term[,2]) ## 2 column is term
    
    ## config
    config.check <- FALSE
    if(length(PK.match[is.na(PK.match)]) > 0)
        config.check <- gmessage(paste(dQuote(require.term[is.na(PK.match)]), " is NOT defined! Please Define it in Configure Menu.", sep=""),
                 icon = c("error"),title="Error")
   
   if (!config.check)
   {
       BGTest = BasicGUI$new(message="Random effects diagnostics",
            widgetList = list(
         "Distribution of ETAS" = list(type="gdroplist",
                    items=c(match.term[match.term$TermName == require.term,]$VarName)),
         "Distribution of ETAs (QQ)" = list(type="gdroplist", items=c(match.term[match.term$TermName == require.term,]$VarName)),
         "Scatterplot matrix of ETAs" = list(type="glabel", text="")

        ),
       saveList = list(
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
    ))
        BGTest$okButtonHandler = model.random.okButtonHandler
        BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler              
        BGTest$ggobiImageHandler = model.ggobiImageHandler
        BGTest$saveImageHandler = saveImageHandler.pkmodel                        
       BGTest$show()
   }
}
################################################################################
## Menu "Model validation"
################################################################################
################################################################################
pk.outlier.psn <- function()
{
  # TODO: CHECK DATA?
    if(!checkDataExist())
    {
        gmessage("No data is available for outlier detection!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }
    
  cleanPKCode()
  cleanPKGGobi()
  
  currentPage <- svalue(pmg.dialog.notebook)
    BGTest = BasicGUI$new(message="Influence analysis summary (PsN)",  #parent.win = PKW,
                      widgetList = list(
                   "PsN Case deletion diagnostics:" = list(type="glabel", text=""),
                   "________________________________" = list(type="glabel", text=""),
                   "Result file:" = list(type = "gdroplist", items = c("", dir(pattern="csv"))),
                   "Note 1:" = list(type="glabel", text="raw_results1.csv (default file)"),
                   "________________________________" = list(type="glabel", text=""),
                   "Deleted ID file:" = list(type = "gdroplist", items = c("", dir(pattern="csv"))),
                   "Note 2:" = list(type="glabel", text="skipped_individuals1.csv (default file)")                   

                   ) ,
                       saveList = list(
                     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                    )
                  )


     BGTest$okButtonHandler = psn.outlier.okButtonHandler
     BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler           
     BGTest$ggobiImageHandler = psn.outlier.ggobiImageHandler

     BGTest$datagroup = 0
     BGTest$saveImageHandler = saveImageHandler
     BGTest$show()
     
}

pk.outlier.vis <- function()
{
  # TODO: CHECK DATA?
    if(!checkDataExist())
    {
        gmessage("No data is available for outlier detection!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

  cleanPKCode()
  cleanPKGGobi()

  currentPage <- svalue(pmg.dialog.notebook)
    BGTest = BasicGUI$new(message="Visualization for influence analysis",  #parent.win = PKW,
                      widgetList = list(
                   "Visualization for NONMEM runs:" = list(type="glabel", text=""),
                   "________________________________" = list(type="glabel", text=""),
                   "Target directory path:" = list(type="gedit", text=""),
                   "Simulation folder pattern:" =list(type="gedit", text="NM_run"),
                   "NONMEM result file name:" =list(type="gedit", text="CS1_IV1ESTFPDF-1.fit"),
                   "Patient ID:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage))),
                   "Plot variable:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage)))

                   ),
                       saveList = list(
                     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                    )
                  )


     BGTest$okButtonHandler = vis.outlier.okButtonHandler
     BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler           
     BGTest$ggobiImageHandler = vis.outlier.ggobiImageHandler

     BGTest$datagroup = 0
     BGTest$saveImageHandler = saveImageHandler
     BGTest$show()

}

pk.boot.vis <- function()
{
  # TODO: CHECK DATA?
    if(!checkDataExist())
    {
        gmessage("No original data is available for bootstrap!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

  cleanPKCode()
  cleanPKGGobi()

  currentPage <- svalue(pmg.dialog.notebook)
    BGTest = BasicGUI$new(message="Visualization for bootstrap",  #parent.win = PKW,
                      widgetList = list(
                   #"Option 1: PsN Case-delete-one results" = list(type="glabel", text=""),
                   #"________________________________" = list(type="glabel", text=""),
                   #"Result file:" = list(type = "gdroplist", items = c("", dir(pattern="csv"))),
                   #"Deleted ID file:" = list(type = "gdroplist", items = c("", dir(pattern="csv"))),
                   #"Option 2: NONMEM runs" = list(type="glabel", text=""),
                   "________________________________" = list(type="glabel", text=""),
                   "Target directory path:" = list(type="gedit", text=""),
                   "Bootstrap folder pattern:" =list(type="gedit", text="NM_run"),
                   "NONMEM result file name:" =list(type="gedit", text="CS1_IV1ESTFPDF-1.fit"),
                   "________________________________" = list(type="glabel", text=""),
                   "Bootstrap key table path:" =list(type="gedit", text=""),
                   "Bootstrap key table name:" =list(type="gedit", text="included_individuals1.csv"),
                   "Patient ID:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage))),
                   "Plot variable:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage)))
                   ),
                       saveList = list(
                     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                    )
                  )


     BGTest$okButtonHandler = psn.bootstrap.vis.okButtonHandler
     BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler           
     BGTest$ggobiImageHandler = boot.vis.ggobiImageHandler
     BGTest$datagroup = 0
     BGTest$saveImageHandler = saveImageHandler
     BGTest$show()

    cleanDataSpecialPlot()
    cleanDataLayoutPlot()

}

pk.boot.sum <- function()
{
  # TODO: CHECK DATA?
    if(!checkDataExist())
    {
        gmessage("No data is available for bootstrap detection!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

  cleanPKCode()

    BGTest = BasicGUI$new(message="Bootstrap summary (PsN)",  # parent.win = PKW,
                      widgetList = list(
                   "PsN result file:" = list(type = "gdroplist", items = c(dir(pattern="csv"))),
                   "Note 1:" = list(type="glabel", text="raw_results1.csv (default file)"),
                   "________________________________" = list(type="glabel", text=""),                   
                   "Bootstrap key file:" = list(type = "gdroplist", items = c(dir(pattern="csv"))),
                   "Note 2:" = list(type="glabel", text="included_individuals1.csv (default file)"),
                   "________________________________" = list(type="glabel", text=""),                                                         
                   "number of bins" = list(type="gedit",text="")

                   ),
                       saveList = list(
                     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                    )
                  )


     BGTest$okButtonHandler = psn.bootstrap.sum.okButtonHandler
     BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler           
     BGTest$ggobiImageHandler = psn.bootstrap.ggobiImageHandler
     BGTest$datagroup = 1
     BGTest$saveImageHandler = saveImageHandler
     BGTest$show()

}
################################################################################
pk.sim <- function()
{
  # TODO: CHECK DATA?
    if(!checkDataExist())
    {
        gmessage("No data is available for outlier detection!",
                icon=c("warning"), title="Warning")
        return(invisible(NULL))
    }

  cleanPKCode()
  cleanPKGGobi()

  currentPage <- svalue(pmg.dialog.notebook)
    BGTest = BasicGUI$new(message="Influential analysis",  #parent.win = PKW,
                      widgetList = list(
                   "Simulation runs:" = list(type="glabel", text=""),
                   "________________________________" = list(type="glabel", text=""),
                   "Directory path_Model 1:" = list(type="gedit", text=""),
                   "Simulation folder pattern_Model 1:" =list(type="gedit", text="NM_run"),
                   "NONMEM result file name_Model 1:" =list(type="gedit", text="CS1_IV1ESTFPDF-1.fit"),
                   "-------------" = list(type="glabel", text=""),
                   "Directory path_Model 2:" = list(type="gedit", text=""),
                   "Simulation folder pattern_Model 2:" =list(type="gedit", text="NM_run"),
                   "NONMEM result file name_Model 2:" =list(type="gedit", text="CS1_IV1ESTFPDF-1.fit"),
                   "________________________________" = list(type="glabel", text=""),
                   "Patient ID:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage))),
                   "Plot variable:" = list(type="gdroplist", items=colnames(getCurrentData(currentPage)))
                   ),
                       saveList = list(
                   "conditional var" = list(type = "gdroplist", items = c("", colnames(getCurrentData(currentPage)))),
                   layout_x = list(type="gedit",text="5"),
                   layout_y = list(type="gedit", text="5") ,
                     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
                    )
                  )


     BGTest$okButtonHandler = psn.sim.okButtonHandler
     BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler           
     BGTest$ggobiImageHandler = validation.ggobiImageHandler
     BGTest$datagroup = 1
     BGTest$saveImageHandler = saveImageHandler
     BGTest$show()

}


################################################################################
## Menu "model comparison"
################################################################################
pk.com.config <- function()
{
  labelMessage1 <- "Please move target dataset from left TABLE to right TABLE."
  labelMessage2 <- "After choosing data, click to map."
  winTitle <- "Configure datasets"
  statusMessage <- "Data is ready for model comparison."

  selectDataDialog(winTitle, labelMessage1, labelMessage2, statusMessage, menuOption=1)

}

pk.com.map <- function()
{

   data.name <- getComDataName()
   # check com data exists
   if (length(data.name) != 2)
   {
      ErrorMessage("You need to have TWO data sets for model comparison!")
      return(invisible(NULL))
   }

   data1 <- getCurrentData(data.name[1])
   data2 <- getCurrentData(data.name[2])
   
   name1 <- colnames(data1)
   name2 <- colnames(data2)
   
   if (length(name1)<1 || length(name2) <1)
   {
      ErrorMessage("Data does Not have column names. You need to have column names to map!")
      return(invisible(NULL))
   }

   widget.list <- list()
   sapply(1:length(name1), function(i)
          {
                widget.list[[name1[i]]] <<- list(type="gdroplist", items=c("", name2), expand=TRUE)
          })

   configData = SimpleGUI$new(message="Configure mapping for model comparison",
        widgetList = widget.list, 
        note="Note: The goal of this step is to merge two data. 
        Please choose data variables before fit.
        Yes - ID/TIME/DV/Covariates
        No - PRED/IPRED/WRES/Parameters" ## 1106
    )

    configData$okButtonHandler = com.map.okButtonHandler
    configData$show()
}

pk.com.hist <- function()
{
    # get current data type
    currentPage <- svalue(pmg.dialog.notebook)
    datatype <- getCurrentDataType(currentPage)
    
    if(datatype != modelComType)
    {
        ErrorMessage("You need to generate model comparison dataset from submenu: select datasets and configure mapping!")
        return(invisible(NULL))
    }

    ## prepare for save and ggobi. clean code
    # clean code from other window
    cleanPKCode()
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
    cleanPKGGobi()
    
    data.name <- getComDataName()
    mydata <- getCurrentData()
    myname <- colnames(mydata)
    name1 <- myname[grep(".x", myname)]
    name2 <- myname[grep(".y", myname)]
    
   BGTest = BasicGUI$new(message="Model comparison",
        widgetList = list(
     "Model 1: " = list(type="glabel", text=data.name[1]),
     x1 = list(type = "gdroplist", items = name1),
     "________________________________" = list(type="glabel", text=""),
     "Model 2: " = list(type="glabel", text=data.name[2]),
     x2 = list(type = "gdroplist", items = name2),
     "number of bins" = list(type="gedit",text=""),
     main = list(type="gedit",text=""),
     xlab = list(type="gedit",text=""),
     ylab = list(type="gedit",text=""),
     type = list(type = "gdroplist", items = c("percent", "count", "density"))
    ),
       saveList = list(
     "cond" = list(type = "gdroplist", items = c("", myname)),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="2") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
    )
    )
    BGTest$okButtonHandler = com.hist.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = model.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()

}
pk.com.scatter <- function()
{
    # get current data type
    currentPage <- svalue(pmg.dialog.notebook)
    datatype <- getCurrentDataType(currentPage)

    if(datatype != modelComType)
    {
        ErrorMessage("You need to generate model comparison dataset from submenu: select datasets and configure mapping!")
        return(invisible(NULL))
    }

    ## prepare for save and ggobi. clean code
    # clean code from other window
    cleanPKCode()
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
    cleanPKGGobi()

    data.name <- getComDataName()
    mydata <- getCurrentData()
    myname <- colnames(mydata)
    name1 <- myname[grep(".x", myname)]
    name2 <- myname[grep(".y", myname)]

   BGTest = BasicGUI$new(message="Model comparison",
        widgetList = list(
     "Model 1: " = list(type="glabel", text=data.name[1]),
     x = list(type = "gdroplist", items = name1),
     "________________________________" = list(type="glabel", text=""),
     "Model 2: " = list(type="glabel", text=data.name[2]),
     y = list(type = "gdroplist", items = name2),
     main = list(type="gedit",text=" "),
     xlab = list(type="gedit",text=" "),
     ylab = list(type="gedit",text=" "),
     type = list(type = "gdroplist", items = c("p", "l", "b"))

    ),
       saveList = list(
     "cond" = list(type = "gdroplist", items = c("", myname)),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="2") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
    )
    )
    BGTest$okButtonHandler = com.scatter.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = model.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()
}

pk.com.time <- function()
{
    # get current data type
    currentPage <- svalue(pmg.dialog.notebook)
    datatype <- getCurrentDataType(currentPage)
    
    if(datatype != modelComType)
    {
        ErrorMessage("You need to generate model comparison dataset from submenu: select datasets and configure mapping!")
        return(invisible(NULL))
    }

    ## prepare for save and ggobi. clean code
    # clean code from other window
    cleanPKCode()
    cleanDataSpecialPlot()
    cleanDataLayoutPlot()
    cleanPKGGobi()

    data.name <- getComDataName()
    mydata <- getCurrentData()
    myname <- colnames(mydata)
    name1 <- myname[grep(".x", myname)]
    name2 <- myname[grep(".y", myname)]
    name3 <- myname[!(myname %in% c(name1, name2))]

   BGTest = BasicGUI$new(message="Model comparison",
        widgetList = list(
     x = list(type = "gdroplist", items = name3),
     "________________________________" = list(type="glabel", text=""),
     "Model 1: " = list(type="glabel", text=data.name[1]),
     y1 = list(type = "gdroplist", items = name1),
     "Model 2: " = list(type="glabel", text=data.name[2]),
     y2 = list(type = "gdroplist", items = name2),
     "________________________________" = list(type="glabel", text=""),
     main = list(type="gedit",text=""),
     xlab = list(type="gedit",text=""),
     ylab = list(type="gedit",text=""),
     type = list(type = "gdroplist", items = c("p", "l", "b"))
    ),
       saveList = list(
     "Transform" = list(type = "gdroplist", items = c("y1/y2", "log(y1/y2)")),
     "cond" = list(type = "gdroplist", items = c("", myname)),
     layout_x = list(type="gedit",text="1"),
     layout_y = list(type="gedit", text="2") ,
     graphics = list(type="gradio", items = c("lattice","ggplot2"), horizontal=TRUE)
    )
    )
    BGTest$okButtonHandler = com.time.okButtonHandler
    BGTest$cleanFigureButtonHandler = cleanFigureButtonHandler          
    BGTest$ggobiImageHandler = model.ggobiImageHandler
    BGTest$saveImageHandler = saveImageHandler
    BGTest$show()
}
