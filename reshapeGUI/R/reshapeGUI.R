#' reshapeGUI: A GUI for the reshape2 and plyr packages
#' ; A tool for learning how to use the functions, \code{melt}, \code{acast/dcast},
#' and \code{ddply}.
#'
#' This graphical user interface (GUI) was built with the gWidgets package, 
#' under the RGtk2 toolkit.
#'
#' Some features of the GUI are as follows:
#' \itemize{
#'   \item{Use either data frames from your workspace, or read in
#'     .csv files from your hard drive.}
#'   \item{Utilize help buttons (look for ?'s) within GUI for quick hints on
#'     how to use specific parts of the GUI}
#'   \item{Learn the syntax of the \code{melt}, \code{acast/dcast},
#'     and \code{ddply} functions from the dynamically updating code previews}
#'   \item{Make changes to code within the GUI for added flexibility}
#'   \item{Preview both raw data frames, as well as ones created with \code{melt},
#'     \code{acast/dcast}, and \code{ddply} for quick feedback and results}
#'   \item{Execute code in GUI to export created data frames to your workspace}
#'   \item{Export outputted data frame to .csv files}
#' }
#'
#' NOTE: This GUI is meant as a tool for learning how to use the specified functions,
#' not as a replacement for their use. Much of the flexibility of the functions is lost
#' within the GUI. I designed features like the code preview with the intention that 
#' the user would use the GUI in a sort of "training wheels" approach, eventually 
#' transitioning to being able to write the code independently. Good luck!
#'
#' @author Jason Crowley \email{crowley.jason.s@@gmail.com}
#' @keywords reshape plyr melt cast GUI aggregate
#' @examples
#' # Bring two datasets included with reshape2 into the workspace
#' data(tips)
#' data(french_fries)
#'
#' # Run the GUI
#' \dontrun{reshapeGUI()}
reshapeGUI <- function() {
  #require(gWidgets)
  options("guiToolkit"="RGtk2")
  #require(reshape2)
  #require(plyr)

  ### CREATE GUI ###
  w <- gwindow("Data Reshaping & Aggregating", visible = FALSE)
  # Tab layout
  nb <- gnotebook(cont = w)
  tot.count <- 1

  # Change size of window
  size(nb) <- c(800,600)

  # Data Read-in Tab
  p1 <- ggroup(cont = nb,label = "Data Selection",horizontal=FALSE,expand=TRUE)

  ### Tab 1: Data Selection ###
  wkspaceObj <- get.ObjectClasses()
  wkspaceDfs <- names(wkspaceObj)[wkspaceObj == "data.frame"]

  #p1selectionLabelBox <- ggroup(cont=p1)
  #glabel("Current Data Frames in Workspace",cont=p1selectionLabelBox)
  #addSpring(p1selectionLabelBox)
  selectionBox <- ggroup(cont=p1)
  p1wkspaceBox <- ggroup(cont=selectionBox,expand=TRUE,horizontal=FALSE)
  wkspaceDfList <- gtable(wkspaceDfs, cont = p1wkspaceBox, expand=TRUE)
  p1wkspaceRefreshBox <- ggroup(cont=p1wkspaceBox)
  addSpring(p1wkspaceRefreshBox)
  p1wkspace.refresh <- gbutton("Refresh Workspace",cont=p1wkspaceRefreshBox)
  addSpring(p1wkspaceRefreshBox)
  names(wkspaceDfList) <- "Current_Data_Frames_in_Workspace"
  size(wkspaceDfList) <- c(155,150)

  entryButtonBox <- ggroup(cont = selectionBox, expand=TRUE, horizontal=FALSE)
  addSpring(selectionBox)
  glabel("",cont=entryButtonBox)

  p1csvLabelBox <- ggroup(cont = entryButtonBox)
  glabel("File Loading/Exporting",cont=p1csvLabelBox)
  loadcsvBox <- ggroup(cont = entryButtonBox)
  loadcsvButton <- gbutton("Load .csv File",cont=loadcsvBox)
  glabel("Object Name:",cont=loadcsvBox)
  p1csvObjTitle <- gedit("newData",cont=loadcsvBox,width=10)
  loadcsvHeaderCheck <- gcheckbox("Header?",checked=FALSE,cont=loadcsvBox)

  p1WriteBox <- ggroup(cont=entryButtonBox)
  p1WriteButton <- gbutton("Export to .csv File",cont=p1WriteBox)
  glabel("File Name:",cont=p1WriteBox)
  p1WriteName <- gedit("",cont=p1WriteBox,width=15)
  glabel("Include Row Names?",cont=p1WriteBox)
  p1WriteRowNames <- gcheckbox("",cont=p1WriteBox)

  glabel("",cont=entryButtonBox)

  p1ManipLabelBox <- ggroup(cont=entryButtonBox)
  glabel("Data Manipulation",cont=p1ManipLabelBox)
  addSpring(p1ManipLabelBox)
  p1ManipButtonBox <- ggroup(cont=entryButtonBox)
  p1MeltButton <- gbutton("melt",cont=p1ManipButtonBox)

  p1CastButton <- gbutton("cast",cont=p1ManipButtonBox)

  p1PlyButton <- gbutton("ddply",cont=p1ManipButtonBox)
  p1Manip.help <- gbutton("?",cont=p1ManipButtonBox)
  p1ManipLabel <- glabel("",cont=p1ManipButtonBox)

  p1DataLabelBox <- ggroup(cont=p1)
  p1DataLabel <- glabel("Data Preview (select data frame above)",cont=p1DataLabelBox)
  addSpring(p1DataLabelBox)
  p1DataPreview <- gdf(data.frame(items=TRUE),cont=p1,expand=TRUE)
  visible(p1DataPreview) <- FALSE

  # Handler for data manipulation help button
  addHandlerClicked(p1Manip.help, function(h,...) {
    gmessage("Melt: Convert a data frame into a molten data frame by identifying ID variables and measured variables. This is a flexible form of data that can then be cast into the desired form or aggregated, if necessary.\n\nCast: Transform a molten data frame into the desired form (either data frame or array) by choosing which variables to cast in each dimension.\n\nddply: Split a data frame by each unique level of one or a combination of ID variables, apply a function to each subset (mean, for instance), and combine back into a data frame.","Melt/Cast/ddply Help")
  })

  # Handler for selecting a data frame in the table (updates preview)
  addHandlerClicked(wkspaceDfList, function(h,...) {
    if(length(svalue(wkspaceDfList)) > 0) {
      svalue(p1ManipLabel) <- ""
      delete(p1,p1DataPreview)

      tmpData <- eval(parse(text=svalue(wkspaceDfList)))
      svalue(p1DataLabel) <- paste("Data Preview: ",svalue(wkspaceDfList)," (",nrow(tmpData)," rows, ",ncol(tmpData)," columns)",sep="")

      add(p1,(p1DataPreview <<- gdf(tmpData)),expand=TRUE)
    } else {
      svalue(p1ManipLabel) <- ""
    }
  })

  # Handler for clicking Load .csv File button (also adds to table)
  addHandlerClicked(loadcsvButton, function(h,...) {
    delete(p1,p1DataPreview)
    
    header <- svalue(loadcsvHeaderCheck)
    eval(parse(text=paste(svalue(p1csvObjTitle)," <<- read.csv(file.choose(),header=header)",sep="")))
    new.wkspaceObj <- get.ObjectClasses()
    wkspaceDfList[] <- names(new.wkspaceObj)[new.wkspaceObj == "data.frame"]
    
    add(p1,(p1DataPreview <<- gdf(data.frame(items=TRUE))),expand=TRUE)
    visible(p1DataPreview) <- FALSE
  })   

  # Handler to update the data frame list if the tab is changed
  addHandlerChanged(nb, function(h,...) {
    delete(p1,p1DataPreview)
    svalue(wkspaceDfList) <- character(0)
    svalue(p1DataLabel) <- "Data Preview (select data frame above)"
  
    new.wkspaceObj <- get.ObjectClasses()
    wkspaceDfList[] <- names(new.wkspaceObj)[new.wkspaceObj == "data.frame"]
  
    add(p1,(p1DataPreview <<- gdf(data.frame(items=TRUE))),expand=TRUE)
    visible(p1DataPreview) <- FALSE
  })

  # Handler for workspace refresh button
  addHandlerClicked(p1wkspace.refresh, function(h,...) {
    delete(p1,p1DataPreview)
    svalue(wkspaceDfList) <- character(0)
    svalue(p1DataLabel) <- "Data Preview (select data frame above)"
  
    new.wkspaceObj <- get.ObjectClasses()
    wkspaceDfList[] <- names(new.wkspaceObj)[new.wkspaceObj == "data.frame"]
  
    add(p1,(p1DataPreview <<- gdf(data.frame(items=TRUE))),expand=TRUE)
    visible(p1DataPreview) <- FALSE
  })

  # Handler for write.csv button
  addHandlerClicked(p1WriteButton, function(h,...) {
    fileName <- svalue(p1WriteName)
    if(!grepl(".csv",fileName,fixed=TRUE)) {
      fileName <- paste(fileName,".csv",sep="")
    }
    cat("File written to ",getwd(),"/",fileName,"\n",sep="")
    eval(parse(text=paste("write.csv(",svalue(wkspaceDfList),",'",fileName,"', row.names=",svalue(p1WriteRowNames),")",sep=""))) 
  })

  ### Melting Tab ###
# Melt To-Do List:
#  - Figure out why variable list doesn't show up at load sometimes.
#    - Found that it always shows up when that is the active window...why?
#  - Adjust size of code box so all code is shown
#    - will probably need to make window larger
#    - found how to make box taller, but how to force code onto new line?
#        and format it nicely?
#  - Look into text formatting for invalid variable combinations
#    - i.e. No ID variable
#    - if this doesn't work, show error messages (do this either way)
#      - gmessage? glabel?
#  - Guess at ID variable. Measure variables?
#    - First variable? First factor variable?
#  - Come up with some way to identify variable type
#    - Text formatting (color)? Type in parentheses?
#  - Warn user when they hit Execute (that they'll overwrite anything named
#      ____________.)?
#  - Change the outgoing arrows handlers so that the user can determine the
#      order of the variables in the melting statement (don't sort)
#  - Mess around with switching variables between boxes...I think there's
#      some bugs
#  - Melting doesn't work when you try to melt in a order different from what
#      the columns are in...test this.
#  - Come up with some way for the Row.Names in the tail preview to reflect
#      what the actual row numbers will be
#      - calculation 
  addHandlerClicked(p1MeltButton, handler = function(h,...) {
  if(length(svalue(wkspaceDfList)) == 0) {
    svalue(p1ManipLabel) <- "Error: No data frame has been selected"
  } else {
    svalue(p1ManipLabel) <- ""

    data_string <- svalue(wkspaceDfList)
    count <- tot.count + 1
    eval(parse(text=paste("p",count," <- ggroup(cont = nb,label = 'melt: ",data_string,"',horizontal=FALSE,expand=TRUE)",sep="")))

    eval(parse(text=paste("data",count," <- ",data_string,sep="")))
    eval(parse(text=paste("dfVars",count," <- data.frame(name=names(data",count,"),type=get.VarTypes(data",count,"))",sep="")))

    # id variables
    eval(parse(text=paste("p",count,"varboxes <- ggroup(cont = p",count,", expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"idvarbox <- ggroup(cont = p",count,"varboxes, horizontal=FALSE, expand=TRUE)",sep="")))
    #eval(parse(text=paste("glabel('ID Variables',cont=p",count,"idvarbox)",sep="")))
    eval(parse(text=paste("p",count,"idvars <- gtable('', cont = p",count,"idvarbox, expand=TRUE, multiple=TRUE)",sep="")))
    eval(parse(text=paste("names(p",count,"idvars) <- 'ID_Variables'",sep="")))

    # id arrows
    eval(parse(text=paste("p",count,"id.arr <- ggroup(cont = p",count,"varboxes, horizontal=FALSE)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"id.arr)",sep="")))
    eval(parse(text=paste("p",count,"id.help <- gbutton('?',cont=p",count,"id.arr)",sep="")))
    eval(parse(text=paste("p",count,"left.id.arr <- gbutton('<-',cont=p",count,"id.arr)",sep="")))
    eval(parse(text=paste("enabled(p",count,"left.id.arr) <- FALSE",sep="")))
    eval(parse(text=paste("p",count,"right.id.arr <- gbutton('->',cont=p",count,"id.arr)",sep="")))
    eval(parse(text=paste("enabled(p",count,"right.id.arr) <- FALSE",sep="")))
    eval(parse(text=paste("addSpring(p",count,"id.arr)",sep="")))

    # unused variables
    eval(parse(text=paste("p",count,"allvarbox <- ggroup(cont = p",count,"varboxes, horizontal=FALSE, expand=TRUE)",sep="")))
    #eval(parse(text=paste("glabel('Unused Variables',cont=p",count,"allvarbox)",sep="")))
    eval(parse(text=paste("p",count,"allvars <- gtable(names(data",count,"),cont=p",count,"allvarbox, expand=TRUE, multiple=TRUE)",sep="")))
    eval(parse(text=paste("names(p",count,"allvars) <- 'Unused_Variables'",sep="")))

    # measure arrows
    eval(parse(text=paste("p",count,"measure.arr <- ggroup(cont = p",count,"varboxes, horizontal=FALSE)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"measure.arr)",sep="")))
    eval(parse(text=paste("p",count,"measure.help <- gbutton('?',cont=p",count,"measure.arr)",sep="")))
    eval(parse(text=paste("p",count,"right.measure.arr <- gbutton('->',cont=p",count,"measure.arr)",sep="")))
    eval(parse(text=paste("enabled(p",count,"right.measure.arr) <- FALSE",sep="")))
    eval(parse(text=paste("p",count,"left.measure.arr <- gbutton('<-',cont=p",count,"measure.arr)",sep="")))
    eval(parse(text=paste("enabled(p",count,"left.measure.arr) <- FALSE",sep="")))
    eval(parse(text=paste("addSpring(p",count,"measure.arr)",sep="")))

    # measure variables
    eval(parse(text=paste("p",count,"measurevarbox <- ggroup(cont = p",count,"varboxes, horizontal=FALSE, expand=TRUE)",sep="")))
    #eval(parse(text=paste("glabel('Measure Variables',cont=p",count,"measurevarbox)",sep="")))
    eval(parse(text=paste("p",count,"measurevars <- gtable('', cont=p",count,"measurevarbox, expand=TRUE, multiple=TRUE)",sep="")))
    eval(parse(text=paste("names(p",count,"measurevars) <- 'Measure_Variables'",sep="")))

    # reset button
    eval(parse(text=paste("p",count,"resetbox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))
    eval(parse(text=paste("p",count,"resetButton <- gbutton('Reset',cont=p",count,"resetbox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))

    # code box / preview button / execute button
    eval(parse(text=paste("p",count,"codebox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("p",count,"code.label <- glabel('Code:',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.box <- gedit(paste(data_string,'.melt <- melt(data = ',data_string,')',sep=''),cont=p",count,"codebox,expand=TRUE)",sep="")))
    #size(p2code.box) <- c(1,50)
    eval(parse(text=paste("p",count,"code.preview <- gbutton('Preview',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.execute <- gbutton('Execute',cont=p",count,"codebox)",sep="")))

    # raw data preview
    eval(parse(text=paste("p",count,"previewbox <- ggroup(cont=p",count,",expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawbox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawlabel <- glabel('Raw Data',cont=p",count,"rawbox)",sep="")))
    eval(parse(text=paste("p",count,"rawd <- gdf(data",count,",cont=p",count,"rawbox,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"meltbox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))

    # melted data preview
    eval(parse(text=paste("p",count,"meltlabel <- glabel('Melted Data',cont=p",count,"meltbox)",sep="")))
    eval(parse(text=paste("numMeltedRows <- num.MeltedRows(data",count,",paste('melt(',data_string,')',sep=''))",sep="")))
    eval(parse(text=paste("add(p",count,"meltbox,(p",count,"meltd <- gdf(eval(parse(text = paste('head(melt(head(',data_string,')))',sep=''))))),expand=TRUE)",sep="")))
    eval(parse(text=paste("add(p",count,"meltbox,(p",count,"meltdlab <- glabel('...')))",sep="")))
    eval(parse(text=paste("add(p",count,"meltbox,(p",count,"meltd2 <- gdf(eval(parse(text = paste('tail(melt(tail(',data_string,')))',sep=''))))),expand=TRUE)",sep="")))
    eval(parse(text=paste("row.names(p",count,"meltd2) <- as.character((numMeltedRows-5):numMeltedRows)",sep="")))
    #eval(parse(text=paste("visible(p",count,"meltd) <- FALSE",sep="")))
    #eval(parse(text=paste("visible(p",count,"meltdlab) <- FALSE",sep="")))
    #eval(parse(text=paste("visible(p",count,"meltd2) <- FALSE",sep="")))

    # Create handlers for when users select values in the tables
    eval(parse(text=paste("addHandlerClicked(p",count,"idvars, function(h,...) {
      enabled(p",count,"right.id.arr) <- TRUE
      enabled(p",count,"left.id.arr) <- FALSE
      enabled(p",count,"right.measure.arr) <- FALSE
      enabled(p",count,"left.measure.arr) <- FALSE
    })",sep="")))

    eval(parse(text=paste("addHandlerClicked(p",count,"allvars, function(h,...) {
      enabled(p",count,"right.id.arr) <- FALSE
      enabled(p",count,"left.id.arr) <- TRUE
      enabled(p",count,"right.measure.arr) <- TRUE
      enabled(p",count,"left.measure.arr) <- FALSE
    })",sep="")))

    eval(parse(text=paste("addHandlerClicked(p",count,"measurevars, function(h,...) {
      enabled(p",count,"right.id.arr) <- FALSE
      enabled(p",count,"left.id.arr) <- FALSE
      enabled(p",count,"right.measure.arr) <- FALSE
      enabled(p",count,"left.measure.arr) <- TRUE
    })",sep="")))

    # Create handlers for when the user clicks the help buttons
      # id variable help
      # measure variable help
    eval(parse(text=paste("addHandlerClicked(p",count,"id.help, function(h,...) {
      gmessage('ID Variables identify a unit of measurement. Collectively, all of the ID variables should identify each individual measurement unit. Factor variables are typically ID variables.','Melting Help')
    })",sep="")))

    eval(parse(text=paste("addHandlerClicked(p",count,"measure.help, function(h,...) {
      gmessage('Measure Variables hold information recorded on individual units of measurement. Numerical variables are typically (but not always) measure variables.','Melting Help')
    })",sep="")))

    # Create handlers for when the user clicks an arrow button
      # Eliminate selected value
      # If vector we're adding it to is empty
      # If vector we're adding it to is not empty
      # Update code box
    eval(parse(text=paste("addHandlerClicked(p",count,"right.id.arr, function(h,...) {
      curVal <- svalue(p",count,"idvars)

      p",count,"idvars[] <- setdiff(p",count,"idvars[], curVal)

      if(length(p",count,"allvars[]) == 0 | (length(p",count,"allvars[]) == 1 & (any(p",count,"allvars[] == '') | any(is.na(p",count,"allvars[])))))
        p",count,"allvars[] <- curVal

      else
        p",count,"allvars[] <- sort(unique(c(p",count,"allvars[], curVal)))
      enabled(p",count,"right.id.arr) <- FALSE

      updateMeltCodeBox(p",count,"code.box, p",count,"idvars, p",count,"measurevars, '",data_string,"')
    })",sep="")))

    # Eliminate selected value
    # If vector we're adding it to is empty
    # If vector we're adding it to is not empty
    # Update code box
    eval(parse(text=paste("addHandlerClicked(p",count,"left.id.arr, function(h,...) {
      curVal <- svalue(p",count,"allvars)

      p",count,"allvars[] <- setdiff(p",count,"allvars[], curVal)

      if(length(p",count,"idvars[]) == 0 | (length(p",count,"idvars[]) == 1 & (any(p",count,"idvars[] == '') | any(is.na(p",count,"idvars[])))))
        p",count,"idvars[] <- curVal

      else
        p",count,"idvars[] <- unique(c(p",count,"idvars[], curVal))
      enabled(p",count,"left.id.arr) <- FALSE


      updateMeltCodeBox(p",count,"code.box, p",count,"idvars, p",count,"measurevars, '",data_string,"')
    })",sep="")))

    # Eliminate selected value
    # If vector we're adding it to is empty
    # If vector we're adding it to is not empty
    # Update code box
    eval(parse(text=paste("addHandlerClicked(p",count,"right.measure.arr, function(h,...) {
      curVal <- svalue(p",count,"allvars)

      p",count,"allvars[] <- setdiff(p",count,"allvars[], curVal)

      if(length(p",count,"measurevars[]) == 0 | (length(p",count,"measurevars[]) == 1 & (any(p",count,"measurevars[] == '') | any(is.na(p",count,"measurevars[])))))
        p",count,"measurevars[] <- curVal

      else
        p",count,"measurevars[] <- unique(c(p",count,"measurevars[], curVal))
      enabled(p",count,"right.measure.arr) <- FALSE

      updateMeltCodeBox(p",count,"code.box, p",count,"idvars, p",count,"measurevars, '",data_string,"')
    })",sep="")))

    # Eliminate selected value
    # If vector we're adding it to is empty
    # If vector we're adding it to is not empty
    # Update code box
    eval(parse(text=paste("addHandlerClicked(p",count,"left.measure.arr, function(h,...) {
      curVal <- svalue(p",count,"measurevars)

      p",count,"measurevars[] <- setdiff(p",count,"measurevars[], curVal)

      if(length(p",count,"allvars[]) == 0 | (length(p",count,"allvars[]) == 1 & (any(p",count,"allvars[] == '') | any(is.na(p",count,"allvars[])))))
        p",count,"allvars[] <- curVal

      else
        p",count,"allvars[] <- sort(unique(c(p",count,"allvars[], curVal)))
      enabled(p",count,"left.measure.arr) <- FALSE

      updateMeltCodeBox(p",count,"code.box, p",count,"idvars, p",count,"measurevars, '",data_string,"')
    })",sep="")))

    # Creates handler for Reset button
      # Update code box
      # Remove preview of melted data
    eval(parse(text=paste("addHandlerClicked(p",count,"resetButton, function(h,...) {
      p",count,"idvars[] <- ''
      p",count,"allvars[] <- names(data",count,")
      p",count,"measurevars[] <- ''
      enabled(p",count,"right.id.arr) <- FALSE
      enabled(p",count,"left.id.arr) <- FALSE
      enabled(p",count,"right.measure.arr) <- FALSE
      enabled(p",count,"left.measure.arr) <- FALSE

      svalue(p",count,"code.box) <- paste(data_string,'.melt <- melt(data = ',data_string,')',sep='')

      num <- num.MeltedRows(data",count,",paste('melt(',data_string,')',sep=''))

      delete(p",count,"meltbox,p",count,"meltd2)
      delete(p",count,"meltbox,p",count,"meltdlab)
      delete(p",count,"meltbox,p",count,"meltd)
      add(p",count,"meltbox,(p",count,"meltd <<- gdf(eval(parse(text = paste('head(melt(head(',data_string,')))',sep=''))))),expand=TRUE)
      add(p",count,"meltbox,(p",count,"meltdlab <<- glabel('...')))
      add(p",count,"meltbox,(p",count,"meltd2 <<- gdf(eval(parse(text = paste('tail(melt(tail(',data_string,')))',sep=''))))),expand=TRUE)
      row.names(p",count,"meltd2) <- as.character((num-5):num)

      #visible(p",count,"meltd) <- FALSE
      #visible(p",count,"meltdlab) <- FALSE
      #visible(p",count,"meltd2) <- FALSE
    })",sep="")))

    # Creates handler for Preview button
      # Adjust row.names to match the true length of the molten dataset
    eval(parse(text=paste("addHandlerClicked(p",count,"code.preview, function(h,...) {
      evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))
      s1 <- substr(evalExpr,1,regexpr('=',evalExpr)+1)
      s2 <- substr(evalExpr,regexpr('=',evalExpr)+2,regexpr(',',evalExpr)-1)
      s3 <- substr(evalExpr,regexpr(',',evalExpr),nchar(evalExpr))

      hd_evalExpr <- paste('head(',s1,'head(',s2,')',s3,')',sep='')
      tl_evalExpr <- paste('tail(',s1,'tail(',s2,')',s3,')',sep='')

      tl_preview <- eval(parse(text = tl_evalExpr))
      num <- num.MeltedRows(data",count,",evalExpr)

      delete(p",count,"meltbox,p",count,"meltd2)
      delete(p",count,"meltbox,p",count,"meltdlab)
      delete(p",count,"meltbox,p",count,"meltd)
      add(p",count,"meltbox,(p",count,"meltd <<- gdf(eval(parse(text = hd_evalExpr)))),expand=TRUE)
      add(p",count,"meltbox,(p",count,"meltdlab <<- glabel('...')))
      add(p",count,"meltbox,(p",count,"meltd2 <<- gdf(eval(parse(text = tl_evalExpr)))),expand=TRUE)
      row.names(p",count,"meltd2) <- as.character((num-5):num)
    })",sep="")))

    # Creates handler for Execute button
    eval(parse(text=paste("addHandlerClicked(p",count,"code.execute, function(h,...) {
      meltedDataName <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), 1, regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))-2)
      evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))
      cat(paste('Executed the code:\n> ',meltedDataName,' <- ',evalExpr,'\n',sep=''))
      eval(parse(text = paste(meltedDataName,' <<- ',evalExpr,sep='')))
    })",sep="")))

    tot.count <<- tot.count + 1

  }
  })


  ### Casting Tab ###
  # Include:
  #   dcast only? allow acast, but no preview?
  #    
  addHandlerClicked(p1CastButton, handler = function(h,...) {
  if(length(svalue(wkspaceDfList)) == 0) {
    svalue(p1ManipLabel) <- "Error: No data frame has been selected"
  } else {
    svalue(p1ManipLabel) <- ""

    data_string <- svalue(wkspaceDfList)
    count <- tot.count + 1
    eval(parse(text=paste("p",count," <- ggroup(cont = nb,label = 'cast: ",data_string,"',horizontal=FALSE,expand=TRUE)",sep="")))

    eval(parse(text=paste("data",count," <- ",data_string,sep="")))
    eval(parse(text=paste("dfVars",count," <- data.frame(name=names(data",count,"),type=get.VarTypes(data",count,"))",sep="")))

    # setup
    eval(parse(text=paste("p",count,"setupbox <- ggroup(cont = p",count,",horizontal=FALSE)",sep="")))
    eval(parse(text=paste("p",count,"outputbox <- ggroup(cont = p",count,"setupbox)",sep="")))
    eval(parse(text=paste("glabel('Output:',cont=p",count,"outputbox)",sep="")))
    eval(parse(text=paste("p",count,"output.select <- gradio(c('Data Frame','Array'),cont=p",count,"outputbox,horizontal=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"output.help <- gbutton('?',cont=p",count,"outputbox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"outputbox)",sep="")))

    eval(parse(text=paste("p",count,"aggLabelBox <- ggroup(cont=p",count,"setupbox)",sep="")))
    eval(parse(text=paste("glabel('Columns to aggregate over:',cont=p",count,"aggLabelBox)",sep="")))
    eval(parse(text=paste("p",count,"agg.help <- gbutton('?',cont=p",count,"aggLabelBox)",sep="")))

    eval(parse(text=paste("p",count,"aggBox <- glayout(cont=p",count,"setupbox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[1,2] <- glabel('First Variable:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[1,3] <- glabel('Second Variable:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[2,1] <- glabel('    First Dim.:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[2,2] <- (p",count,"agg11 <- gdroplist(c(names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[2,3] <- (p",count,"agg12 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[3,1] <- glabel('    Second Dim.:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[3,2] <- (p",count,"agg21 <- gdroplist(c(names(data",count,"),'(no variable)','(all ID vars)'),selected=2,cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[3,3] <- (p",count,"agg22 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[4,1] <- glabel('    Third Dim.:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[4,2] <- (p",count,"agg31 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[4,3] <- (p",count,"agg32 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[5,1] <- glabel('    Fourth Dim.:',cont=p",count,"aggBox)",sep="")))
    eval(parse(text=paste("p",count,"aggBox[5,2] <- (p",count,"agg41 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[5,3] <- (p",count,"agg42 <- gdroplist(c('n/a',names(data",count,"),'(no variable)','(all ID vars)'),cont=p",count,"aggBox))",sep="")))

    eval(parse(text=paste("p",count,"aggBox[1,4] <- glabel('        ')",sep="")))
    eval(parse(text=paste("p",count,"aggBox[1,5] <- glabel('Calculate Margins for:')",sep="")))
    eval(parse(text=paste("p",count,"aggBox[2:5,5] <- (p",count,"margins <- gtable(c('(none)','(all)',names(data",count,")[1:2]),multiple=TRUE))",sep="")))
    eval(parse(text=paste("svalue(p",count,"margins) <- '(none)'",sep="")))
    eval(parse(text=paste("p",count,"aggBox[1,6] <- glabel('        ')",sep="")))

    eval(parse(text=paste("p",count,"aggBox[1,7] <- glabel('Aggregating Function:')",sep="")))
    eval(parse(text=paste("p",count,"aggBox[2,7] <- (p",count,"AggFun <- gdroplist(c('(no aggreg.)','mean','length','median','sd','sum','<Enter Your Own>'),editable=TRUE))",sep="")))
    eval(parse(text=paste("p",count,"aggBox[3,7] <- glabel('Subset Expression:')",sep="")))
    eval(parse(text=paste("p",count,"aggBox[4,7] <- (p",count,"Subset <- gedit(initial.msg='ex: variable == \"total_bill\"',width=20))",sep="")))


    # reset button
    eval(parse(text=paste("p",count,"resetbox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))
    eval(parse(text=paste("p",count,"resetButton <- gbutton('Reset',cont=p",count,"resetbox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))

    # code box / preview button / execute button
    eval(parse(text=paste("p",count,"codebox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("p",count,"code.label <- glabel('Code:',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.box <- gedit('',cont=p",count,"codebox,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"code.preview <- gbutton('Preview',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.execute <- gbutton('Execute',cont=p",count,"codebox)",sep="")))

    # raw data preview
    eval(parse(text=paste("p",count,"previewbox <- ggroup(cont=p",count,",expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawbox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawlabel <- glabel('Raw Data',cont=p",count,"rawbox)",sep="")))
    eval(parse(text=paste("p",count,"rawd <- gdf(data",count,",cont=p",count,"rawbox,expand=TRUE)",sep="")))

    # melted data preview
    eval(parse(text=paste("p",count,"castbox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"castlabel <- glabel('Casted Data',cont=p",count,"castbox)",sep="")))
    eval(parse(text=paste("add(p",count,"castbox,(p",count,"castd <- gdf(NA,expand=TRUE)))",sep="")))
    eval(parse(text=paste("visible(p",count,"castd) <- FALSE",sep="")))



    # function to update the code box
    eval(parse(text=paste("codeUpdate <- function(h,...) {
      if(svalue(p",count,"output.select) == 'Data Frame') {
        if(svalue(p",count,"agg11) == '(no variable)') {
          codeHolder <- paste('",data_string,".cast <- dcast(data = ",data_string,", formula = .',sep='')
        } else if (svalue(p",count,"agg11) == '(all ID vars)') {
          codeHolder <- paste('",data_string,".cast <- dcast(data = ",data_string,", formula = ...',sep='')
        } else {
          codeHolder <- paste('",data_string,".cast <- dcast(data = ",data_string,", formula = ',svalue(p",count,"agg11),sep='')
        }

        if(svalue(p",count,"agg12) != 'n/a') {
          if(svalue(p",count,"agg12) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg12) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg12),sep='')
          }
        }

        if(svalue(p",count,"agg21) == '(no variable)') {
          codeHolder <- paste(codeHolder,' ~ .',sep='')
        } else if (svalue(p",count,"agg21) == '(all ID vars)') {
          codeHolder <- paste(codeHolder,' ~ ...',sep='')
        } else {
          codeHolder <- paste(codeHolder,' ~ ',svalue(p",count,"agg21),sep='')
        }

        if(svalue(p",count,"agg22) != 'n/a') {
          if(svalue(p",count,"agg22) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg22) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg22),sep='')
          }
        }

        if(svalue(p",count,"AggFun) != '(no aggreg.)') {
          codeHolder <- paste(codeHolder,', fun.aggregate = ',svalue(p",count,"AggFun),sep='')
        }

        if (length(svalue(p",count,"margins)) > 1) {
          margVars <- paste(svalue(p",count,"margins),collapse='\", \"')
          codeHolder <- paste(codeHolder,', margins = c(\"',margVars,'\")',sep='')
        } else if(length(svalue(p",count,"margins)) > 0) {
          if(svalue(p",count,"margins) != '(none)') {
            if(svalue(p",count,"margins) == '(all)') {
              codeHolder <- paste(codeHolder,', margins = TRUE',sep='')
            } else {
              codeHolder <- paste(codeHolder,', margins = \"',svalue(p",count,"margins),'\"',sep='')
            }
          }  
        }

        if(svalue(p",count,"Subset) != '') {
          codeHolder <- paste(codeHolder,', subset = .(',svalue(p",count,"Subset),')',sep='')
        }

        svalue(p",count,"code.box) <- paste(codeHolder,')',sep='')
      } else {
        if(svalue(p",count,"agg11) == '(no variable)') {
          codeHolder <- paste('",data_string,".cast <- acast(data = ",data_string,", formula = .',sep='')
        } else if (svalue(p",count,"agg11) == '(all ID vars)') {
          codeHolder <- paste('",data_string,".cast <- acast(data = ",data_string,", formula = ...',sep='')
        } else {
          codeHolder <- paste('",data_string,".cast <- acast(data = ",data_string,", formula = ',svalue(p",count,"agg11),sep='')
        }

        if(svalue(p",count,"agg12) != 'n/a') {
          if(svalue(p",count,"agg12) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg12) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg12),sep='')
          }
        }

        if(svalue(p",count,"agg21) == '(no variable)') {
          codeHolder <- paste(codeHolder,' ~ .',sep='')
        } else if (svalue(p",count,"agg21) == '(all ID vars)') {
          codeHolder <- paste(codeHolder,' ~ ...',sep='')
        } else {
          codeHolder <- paste(codeHolder,' ~ ',svalue(p",count,"agg21),sep='')
        }

        if(svalue(p",count,"agg22) != 'n/a') {
          if(svalue(p",count,"agg22) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg22) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg22),sep='')
          }
        }

        if(svalue(p",count,"agg31) != 'n/a') {
          if(svalue(p",count,"agg31) == '(no variable)') {
            codeHolder <- paste(codeHolder,' ~ .',sep='')
          } else if (svalue(p",count,"agg31) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' ~ ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' ~ ',svalue(p",count,"agg31),sep='')
          }
        }

        if(svalue(p",count,"agg32) != 'n/a') {
          if(svalue(p",count,"agg32) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg32) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg32),sep='')
          }
        }

        if(svalue(p",count,"agg41) != 'n/a') {
          if(svalue(p",count,"agg41) == '(no variable)') {
            codeHolder <- paste(codeHolder,' ~ .',sep='')
          } else if (svalue(p",count,"agg41) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' ~ ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' ~ ',svalue(p",count,"agg41),sep='')
          }
        }

        if(svalue(p",count,"agg42) != 'n/a') {
          if(svalue(p",count,"agg42) == '(no variable)') {
            codeHolder <- paste(codeHolder,' + .',sep='')
          } else if (svalue(p",count,"agg42) == '(all ID vars)') {
            codeHolder <- paste(codeHolder,' + ...',sep='')
          } else {
            codeHolder <- paste(codeHolder,' + ',svalue(p",count,"agg42),sep='')
          }
        }

        codeHolder <- paste(codeHolder,', fun.aggregate = ',svalue(p",count,"AggFun),sep='')

        if (length(svalue(p",count,"margins)) > 1) {
          margVars <- paste(svalue(p",count,"margins),collapse='\", \"')
          codeHolder <- paste(codeHolder,', margins = c(\"',margVars,'\")',sep='')
        } else if(length(svalue(p",count,"margins)) > 0) {
          if(svalue(p",count,"margins) != '(none)') {
            if(svalue(p",count,"margins) == '(all)') {
              codeHolder <- paste(codeHolder,', margins = TRUE',sep='')
            } else {
              codeHolder <- paste(codeHolder,', margins = \"',svalue(p",count,"margins),'\"',sep='')
            }
          }  
        }

        if(svalue(p",count,"Subset) != '') {
          codeHolder <- paste(codeHolder,', subset = .(',svalue(p",count,"Subset),')',sep='')
        }

        svalue(p",count,"code.box) <- paste(codeHolder,')',sep='')
      }
    }",sep="")))

    # function to update the margin table
    eval(parse(text=paste("marginUpdate <- function(h,...) {
      vars <- unique(c(svalue(p",count,"agg11),svalue(p",count,"agg12),svalue(p",count,"agg21),svalue(p",count,"agg22),svalue(p",count,"agg31),svalue(p",count,"agg32),svalue(p",count,"agg41),svalue(p",count,"agg42)))
      vars <- vars[!(vars %in% c('(no variable)','n/a'))]
      if('(all ID vars)' %in% vars) {
        p",count,"margins[] <- c('(none)','(all)',names(data",count,")[!(names(data",count,") == 'value')])
        svalue(p",count,"margins) <- '(none)'
      } else {
        p",count,"margins[] <- c('(none)','(all)',vars)
        svalue(p",count,"margins) <- '(none)'
      }
    }",sep="")))

    # add handlers to formula fields to update margin table
    eval(parse(text=paste("addHandlerChanged(p",count,"agg11, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg12, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg21, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg22, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg31, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg32, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg41, marginUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg42, marginUpdate)",sep="")))

    # add handlers to all change-able fields to update code box
    eval(parse(text=paste("addHandlerClicked(p",count,"output.select, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg11, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg12, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg21, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg22, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg31, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg32, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg41, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"agg42, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerClicked(p",count,"margins, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerChanged(p",count,"AggFun, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerKeystroke(p",count,"AggFun, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerKeystroke(p",count,"Subset, codeUpdate)",sep="")))

    # handler for help buttons
    eval(parse(text=paste("addHandlerClicked(p",count,"output.help, function(h,...) {
      gmessage('Casted datasets can be saved in either a data frame or an array format. Data frames can only have two output dimensions, while array output does not share this restriction. Note: Previews of array output in this GUI will be coerced into a data frame.','Cast Help')
    })",sep="")))

    eval(parse(text=paste("addHandlerClicked(p",count,"agg.help, function(h,...) {
      gmessage('Choose which variables you want cast over which dimensions. Note that you can cast multiple variables over the same dimension. Due to space concerns, there are dropdown menus for only two per dimension in this GUI, but you are not limited to this. Similarly, you are not limited to four dimensions, either.\n\nNOTE: The \"variable\" column should almost always be one of the columns you cast.','Cast Help')
    })",sep="")))

    # handler for Reset button
    eval(parse(text=paste("addHandlerClicked(p",count,"resetButton, function(h,...) {
      svalue(p",count,"output.select) <- 'Data Frame'
      svalue(p",count,"agg11) <- names(data",count,")[1]
      svalue(p",count,"agg12) <- 'n/a'
      svalue(p",count,"agg21) <- names(data",count,")[2]
      svalue(p",count,"agg22) <- 'n/a'
      svalue(p",count,"agg31) <- 'n/a'
      svalue(p",count,"agg32) <- 'n/a'
      svalue(p",count,"agg41) <- 'n/a'
      svalue(p",count,"agg42) <- 'n/a'
      p",count,"margins[] <- c('(none)','(all)',names(data",count,")[1:2])
      svalue(p",count,"margins) <- '(none)'
      svalue(p",count,"AggFun) <- 'length'
      svalue(p",count,"Subset) <- ''
      svalue(p",count,"code.box) <- ''
      visible(p",count,"castd) <- FALSE
    })",sep="")))

    # handler for Preview button
    eval(parse(text=paste("addHandlerClicked(p",count,"code.preview, function(h,...) {
      evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))
      evalExpr <- paste('as.data.frame(',evalExpr,')',sep='')

      if(svalue(p",count,"output.select) == 'Array') {
        gmessage('Warning: The preview shown is an array output coerced into a data frame using as.data.frame. Executing this code will still produce array output of appropriate dimensions.','Output Warning')
      }

      delete(p",count,"castbox,p",count,"castd)
      add(p",count,"castbox,(p",count,"castd <<- gdf(eval(parse(text = evalExpr)))),expand=TRUE)

    })",sep="")))

    # Creates handler for Execute button
    eval(parse(text=paste("addHandlerClicked(p",count,"code.execute, function(h,...) {
      meltedDataName <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), 1, regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))-2)
      evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))
      cat(paste('Executed the code:\n> ',meltedDataName,' <- ',evalExpr,'\n',sep=''))
      eval(parse(text = paste(meltedDataName,' <<- ',evalExpr,sep='')))
    })",sep="")))

    tot.count <<- tot.count + 1

  }
  })


  ### ddply Tab ###
  addHandlerClicked(p1PlyButton, handler = function(h,...) {
  if(length(svalue(wkspaceDfList)) == 0) {
    svalue(p1ManipLabel) <- "Error: No data frame has been selected"
  } else {
    svalue(p1ManipLabel) <- ""

    data_string <- svalue(wkspaceDfList)
    count <- tot.count + 1
    eval(parse(text=paste("p",count," <- ggroup(cont = nb,label = 'ddply: ",data_string,"',horizontal=FALSE,expand=TRUE)",sep="")))

    eval(parse(text=paste("data",count," <- ",data_string,sep="")))
    eval(parse(text=paste("dfVars",count," <- data.frame(name=names(data",count,"),type=get.VarTypes(data",count,"))",sep="")))

    # set up the widget groups
    eval(parse(text=paste("p",count,"varbox <- ggroup(cont=p",count,",expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"bybox <- ggroup(cont=p",count,"varbox,horizontal=FALSE)",sep="")))
    eval(parse(text=paste("p",count,"summbox <- ggroup(cont=p",count,"varbox,horizontal=FALSE)",sep="")))
    eval(parse(text=paste("size(p",count,"summbox) <- c(600,100)",sep="")))

    # summarize by variable group
    eval(parse(text=paste("p",count,"bylabelbox <- ggroup(cont=p",count,"bybox)",sep="")))
    eval(parse(text=paste("glabel('Summarize by:',cont=p",count,"bylabelbox)",sep="")))
    eval(parse(text=paste("p",count,"SummBy.help <- gbutton('?',cont=p",count,"bylabelbox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"bylabelbox)",sep="")))
    eval(parse(text=paste("p",count,"byvars <- gtable(names(data",count,"),cont=p",count,"bybox, expand=TRUE, multiple=TRUE)",sep="")))

    # All Variables group
    eval(parse(text=paste("p",count,"TitleBox <- ggroup(cont=p",count,"summbox)",sep="")))
    eval(parse(text=paste("glabel('Choose Method of Summary:',cont=p",count,"TitleBox)",sep="")))
    eval(parse(text=paste("p",count,"MethodOfSumm.help <- gbutton('?',cont=p",count,"TitleBox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"TitleBox)",sep="")))
    eval(parse(text=paste("p",count,"OverallSummBox <- ggroup(cont=p",count,"summbox)",sep="")))
    eval(parse(text=paste("p",count,"OverallSumm.select <- gcheckbox('On All Variables',checked=FALSE,cont=p",count,"OverallSummBox)",sep="")))
    eval(parse(text=paste("p",count,"OverallSumm.stat <- gdroplist(c('mean','median','sd','sum','nrow','<Enter Your Own>'),cont=p",count,"OverallSummBox,editable=TRUE)",sep="")))

    # Specific Summary group
    eval(parse(text=paste("p",count,"SpecSummBox <- ggroup(cont=p",count,"summbox,horizontal=FALSE)",sep="")))
    eval(parse(text=paste("p",count,"SpecSummSelectBox <- ggroup(cont=p",count,"SpecSummBox)",sep="")))
    eval(parse(text=paste("p",count,"SpecSumm.select <- gcheckbox('Specific Statistics',checked=TRUE,cont=p",count,"SpecSummSelectBox)",sep="")))
    eval(parse(text=paste("p",count,"SpecSumm.add <- gbutton('Add Another',cont=p",count,"SpecSummSelectBox)",sep="")))
    eval(parse(text=paste("p",count,"SpecSummArgBox <- ggroup(cont=p",count,"summbox,horizontal=FALSE,expand=TRUE,use.scrollwindow=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"SpecSummArg1 <- ggroup(cont=p",count,"SpecSummArgBox)",sep="")))
    eval(parse(text=paste("p",count,"SpecSummArg1.name <- gedit('',cont=p",count,"SpecSummArg1,initial.msg='ex. average.tip',width=15)",sep="")))
    eval(parse(text=paste("glabel('=',cont=p",count,"SpecSummArg1)",sep="")))
    eval(parse(text=paste("p",count,"SpecSummArg1.command <- gedit('',cont=p",count,"SpecSummArg1,initial.msg='ex. mean(tip)')",sep="")))

    # reset button
    eval(parse(text=paste("p",count,"resetbox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))
    eval(parse(text=paste("p",count,"resetButton <- gbutton('Reset',cont=p",count,"resetbox)",sep="")))
    eval(parse(text=paste("addSpring(p",count,"resetbox)",sep="")))

    # code preview group
    eval(parse(text=paste("p",count,"codebox <- ggroup(cont = p",count,")",sep="")))
    eval(parse(text=paste("p",count,"code.label <- glabel('Code:',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.box <- gedit('',cont=p",count,"codebox,expand=TRUE)",sep="")))
    #eval(parse(text=paste("p",count,"code.update <- gbutton('Update Code',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.preview <- gbutton('Preview',cont=p",count,"codebox)",sep="")))
    eval(parse(text=paste("p",count,"code.execute <- gbutton('Execute',cont=p",count,"codebox)",sep="")))

    # raw data preview
    eval(parse(text=paste("p",count,"previewbox <- ggroup(cont=p",count,",expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawbox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"rawlabel <- glabel('Raw Data',cont=p",count,"rawbox)",sep="")))
    eval(parse(text=paste("p",count,"rawd <- gdf(data",count,",cont=p",count,"rawbox,expand=TRUE)",sep="")))

    # ddply'd data preview
    eval(parse(text=paste("p",count,"plybox <- ggroup(cont=p",count,"previewbox,horizontal=FALSE,expand=TRUE)",sep="")))
    eval(parse(text=paste("p",count,"plylabel <- glabel('Summarized Data',cont=p",count,"plybox)",sep="")))
    eval(parse(text=paste("add(p",count,"plybox,(p",count,"plyd <- gdf(NA)),expand=TRUE)",sep="")))
    eval(parse(text=paste("visible(p",count,"plyd) <- FALSE",sep="")))

    # handlers for the help buttons
    eval(parse(text=paste("addHandlerClicked(p",count,"SummBy.help, function(h,...) {
      gmessage('Summary statistics will be calculated for each unique level of the variable(s) selected below','ddply Help')
    })",sep="")))

    eval(parse(text=paste("addHandlerClicked(p",count,"MethodOfSumm.help, function(h,...) {
      gmessage('On All Variables: The specified function will be applied to all variables in the dataset\n\nSpecific Statistics: The specified, named statistics will be calculated using the summarize function. The field before the equals sign corresponds to the variable name the column will have in the output, while the field after the equals sign corresponds to the actual calculation you want to be done (i.e. mean(tip) to calculate the mean of the tip variable). Click Add Another to add another variable to the output.','ddply help')
    })",sep="")))

    for(i in 2:30) {
      eval(parse(text=paste("add(p",count,"SpecSummArgBox, (p",count,"SpecSummArg",i," <- ggroup()))",sep="")))
      eval(parse(text=paste("add(p",count,"SpecSummArg",i,", (p",count,"SpecSummArg",i,".name <- gedit('',width=15)))",sep="")))
      eval(parse(text=paste("add(p",count,"SpecSummArg",i,", (p",count,"SpecSummArg",i,".command <- gedit('')))",sep="")))

      eval(parse(text=paste("delete(p",count,"SpecSummArg",i,", p",count,"SpecSummArg",i,".command)",sep="")))
      eval(parse(text=paste("delete(p",count,"SpecSummArg",i,", p",count,"SpecSummArg",i,".name)",sep="")))
      eval(parse(text=paste("delete(p",count,"SpecSummArgBox, p",count,"SpecSummArg",i,")",sep="")))
    }

    # function to update the code box
    eval(parse(text=paste("codeUpdate <- function(h,...) {
      if(svalue(p",count,"SpecSumm.select)) {
        specStats <- ''
        splitVars <- paste(svalue(p",count,"byvars),collapse=', ')
        for(i in 1:SpecSumm.count) {
          specStats <- paste(specStats,eval(parse(text=paste('svalue(p',count,'SpecSummArg',i,'.name)',sep=''))),' = ',eval(parse(text=paste('svalue(p',count,'SpecSummArg',i,'.command)',sep=''))),', ',sep='')
        }
        specStats <- substr(specStats,1,nchar(specStats)-2)
        codeHolder <- paste('",data_string,".summary <- ddply(.data = ",data_string,", .variables = .(',splitVars,'), .fun = summarize, ',specStats,')',sep='')
        codeHolder <- gsub(',  = ,', ',', codeHolder)
        codeHolder <- gsub(',  = )', ')', codeHolder)

        svalue(p",count,"code.box) <- codeHolder
      } else {
        splitVars <- paste(svalue(p",count,"byvars),collapse=', ')
        summStat <- svalue(p",count,"OverallSumm.stat)
        codeHolder <- paste('",data_string,".summary <- ddply(.data = ",data_string,", .variables = .(',splitVars,'), .fun = ',summStat,')',sep='')

        svalue(p",count,"code.box) <- codeHolder
      }
    }",sep="")))

    # update the code box when either of the two checkboxes are clicked
    eval(parse(text=paste("addHandlerClicked(p",count,"OverallSumm.select, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerClicked(p",count,"SpecSumm.select, codeUpdate)",sep="")))

    # update the code box when the variable selection box is changed
    eval(parse(text=paste("addHandlerClicked(p",count,"byvars, codeUpdate)",sep="")))

    # update the code box when the Overall Summary stat droplist is clicked or typed in
    eval(parse(text=paste("addHandlerClicked(p",count,"OverallSumm.stat, codeUpdate)",sep="")))
    eval(parse(text=paste("addHandlerKeystroke(p",count,"OverallSumm.stat, codeUpdate)",sep="")))

    # handler for 'Add Another' button
    SpecSumm.count <- 1
    eval(parse(text=paste("addHandlerClicked(p",count,"SpecSumm.add, function(h,...) {
      SpecSumm.count <<- SpecSumm.count + 1
      eval(parse(text=paste('add(p",count,"SpecSummArgBox, (p",count,"SpecSummArg',SpecSumm.count,' <<- ggroup()))',sep='')))
      eval(parse(text=paste('add(p",count,"SpecSummArg',SpecSumm.count,', (p",count,"SpecSummArg',SpecSumm.count,'.name <<- gedit(width=15,initial.msg=\"ex. avg.tip.pct\")))',sep='')))
      eval(parse(text=paste('add(p",count,"SpecSummArg',SpecSumm.count,', glabel(\"=\"))',sep='')))
      eval(parse(text=paste('add(p",count,"SpecSummArg',SpecSumm.count,', (p",count,"SpecSummArg',SpecSumm.count,'.command <<- gedit(initial.msg=\"ex. mean(tip/total_bill)\")))',sep='')))
      eval(parse(text=paste('add(p",count,"SpecSummArg',SpecSumm.count,', (p",count,"SpecSummArg',SpecSumm.count,'.remove <- gbutton(\"Remove\")))',sep='')))

      eval(parse(text=paste('addHandlerKeystroke(p",count,"SpecSummArg',SpecSumm.count,'.name, codeUpdate)',sep='')))

      eval(parse(text=paste('addHandlerKeystroke(p",count,"SpecSummArg',SpecSumm.count,'.command, codeUpdate)',sep='')))

      eval(parse(text=paste('addHandlerClicked(p",count,"SpecSummArg',SpecSumm.count,'.remove, function(h,...) {
        svalue(p",count,"SpecSummArg',SpecSumm.count,'.command) <- \"\"
        svalue(p",count,"SpecSummArg',SpecSumm.count,'.name) <- \"\"
        delete(p",count,"SpecSummArg',SpecSumm.count,', p",count,"SpecSummArg',SpecSumm.count,'.remove)
        delete(p",count,"SpecSummArg',SpecSumm.count,', p",count,"SpecSummArg',SpecSumm.count,'.command)
        delete(p",count,"SpecSummArg',SpecSumm.count,', p",count,"SpecSummArg',SpecSumm.count,'.name)
        delete(p",count,"SpecSummArgBox, p",count,"SpecSummArg',SpecSumm.count,')
      })',sep='')))

      eval(parse(text=paste('addHandlerClicked(p",count,"SpecSummArg',SpecSumm.count,'.remove, codeUpdate)',sep='')))

    })",sep="")))

    # handler for All Variables checkbox
    eval(parse(text=paste("addHandlerClicked(p",count,"OverallSumm.select, function(h,...) {
      if(svalue(p",count,"OverallSumm.select)) {
        svalue(p",count,"SpecSumm.select) <- FALSE
      } else {
        svalue(p",count,"SpecSumm.select) <- TRUE
      }
    })",sep="")))

    # handler for Specific Stat checkbox
    eval(parse(text=paste("addHandlerClicked(p",count,"SpecSumm.select, function(h,...) {
      if(svalue(p",count,"SpecSumm.select)) {
        svalue(p",count,"OverallSumm.select) <- FALSE
      } else {
        svalue(p",count,"OverallSumm.select) <- TRUE
      }
    })",sep="")))

    # handler for updating code box when arg1 name changed
    eval(parse(text=paste("addHandlerKeystroke(p",count,"SpecSummArg1.name, codeUpdate)",sep="")))

    # handler for updating code box when arg1 command changed
    eval(parse(text=paste("addHandlerKeystroke(p",count,"SpecSummArg1.command, codeUpdate)",sep="")))

    # Creates handler for Reset button
      # Update code box
      # Remove preview of melted data
    eval(parse(text=paste("addHandlerClicked(p",count,"resetButton, function(h,...) {
      svalue(p",count,"OverallSumm.select) <- FALSE
      svalue(p",count,"SpecSumm.select) <- TRUE
      svalue(p",count,"OverallSumm.stat) <- 'mean'
      svalue(p",count,"byvars) <- ''

      svalue(p",count,"SpecSummArg1.name) <- ''
      svalue(p",count,"SpecSummArg1.command) <- ''

      for(i in 2:30) {
        eval(parse(text=paste('svalue(p",count,"SpecSummArg',i,'.command) <- \"\"
        svalue(p",count,"SpecSummArg',i,'.name) <- \"\"
        ',sep='')))
      }

      visible(p",count,"plyd) <- FALSE

      svalue(p",count,"code.box) <- ''
      SpecSumm.count <<- 1

    })",sep="")))

    # handler for Preview button
    eval(parse(text=paste("addHandlerClicked(p",count,"code.preview, function(h,...) {
      evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))

      delete(p",count,"plybox,p",count,"plyd)
      add(p",count,"plybox,(p",count,"plyd <<- gdf(eval(parse(text = evalExpr)))),expand=TRUE)

    })",sep="")))

    # Creates handler for Execute button
    eval(parse(text=paste("addHandlerClicked(p",count,"code.execute, function(h,...) {
      if(svalue(p",count,"OverallSumm.select)) {
        meltedDataName <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), 1, regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))-2)
        evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))
        cat(paste('Executed the code:\n> ',meltedDataName,' <- ',evalExpr,'\n',sep=''))
        eval(parse(text = paste(meltedDataName,' <<- ',evalExpr,sep='')))
      } else {
        meltedDataName <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), 1, regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))-2)
        evalExpr <- substr(svalue(eval(parse(text=paste('p',count,'code.box',sep='')))), regexpr('<-',svalue(eval(parse(text=paste('p',count,'code.box',sep='')))))+3, nchar(svalue(eval(parse(text=paste('p',count,'code.box',sep=''))))))

        evalExpr1 <- substr(evalExpr,1,regexpr('summarize,',evalExpr)+10)
        evalExpr2 <- substr(evalExpr,regexpr('summarize,',evalExpr)+11,nchar(evalExpr))
        evalExpr2 <- paste('\n  ',evalExpr2,sep='')
        evalExpr2 <- gsub(', ',',\n  ',evalExpr2)
        evalExpr2 <- paste(substr(evalExpr2,1,nchar(evalExpr2)-1),'\n)',sep='')
        evalExpr <- paste(evalExpr1,evalExpr2,sep='')

        cat(paste('Executed the code:\n> ',meltedDataName,' <- ',evalExpr,'\n',sep=''))
        eval(parse(text = paste(meltedDataName,' <<- ',evalExpr,sep='')))
      }
    })",sep="")))

    tot.count <<- tot.count + 1

  }
  })

  ###

  visible(w) <- TRUE

}

# Helper functions for gui_reshape()

#' Updates the code preview box in the Melt tab
#'
#' @param code widget to update
#' @param id vector of ID variables
#' @param measure vector of measure variables
#' @param dataName character string of data frame
#' @keywords internal
updateMeltCodeBox <- function(code, id, measure, dataName) {
  if(length(measure[]) == 0 | (length(measure[]) == 1 & (any(measure[] == "") | any(is.na(measure[])))))
    codeHolder <- paste(dataName,".melt <- melt(data = ",
      dataName,", id.vars=c(",
      paste("'",paste(id[],collapse="','"),"'",sep=""),
      "))",sep="")
  else
    codeHolder <- paste(dataName,".melt <- melt(data = ",
      dataName,", id.vars=c(",
      paste("'",paste(id[],collapse="','"),"'",sep=""),"), measure.vars=c(",
      paste("'",paste(measure[],collapse="','"),"'",sep=""),
      "))",sep="")

  # Adjust for if the user moved a variable out of the ID Vars, leaving it as NA instead of ""
  if(length(id[]) == 1 & any(is.na(id[])))
    codeHolder <- sub('NA','',codeHolder)

  svalue(code) <- codeHolder
    
}

#' Returns a vector of variable types of the variables in the given df
#' ; (the vector's labels are the variable names)
#'
#' @param df data frame
#' @keywords internal
get.VarTypes <- function(df) {
  return(unlist(lapply(unclass(df),class)))
}

#' Returns a vector of classes for the objects in the user's workspace
#' ; used by the GUI to retrieve a list of data frames in the workspace
#'
#' @keywords internal
get.ObjectClasses <- function() {
  obj <- ls(name=1)
  classes <- rep(NA,length(obj))
  names(classes) <- obj
  for(i in 1:length(obj)) {
    classes[i] <- class(eval(parse(text=obj[i])))
  }
  return(classes)
}

#' Returns the number of rows in melted dataset 
#' given the data frame and the melting command (character string)
#'
#' @keywords internal
num.MeltedRows <- function(df,melt.comm) {
  # Returns the number of rows that will be in a molten data frame, given
  # the original data frame and the melting statement.
  n <- dim(df)[1]
  idvars.Stated <- grepl("id.vars",melt.comm)
  measurevars.Stated <- grepl("measure.vars",melt.comm)
  if(measurevars.Stated) {
    p.measure <- length(eval(parse(text = substr(melt.comm,regexpr("measure.vars",melt.comm)+13,nchar(melt.comm)-1))))
    p <- p.measure
  } else {
    if(idvars.Stated) {
      p.id <- length(eval(parse(text = substr(melt.comm,regexpr("id.vars",melt.comm)+8,regexpr(")",melt.comm)))))
      p <- dim(df)[2] - p.id
    } else {
      p <- dim(df)[2] - 1
    }
  }
  return(n*p)
}