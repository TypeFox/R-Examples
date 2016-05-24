#' @name gimmeInteractive
#' @aliases gimmeInteractive
#' @title Graphical User Interface designed for use with the gimme package
#' @description This function launches a graphical user interface which calls the main functions
#' for the gimme package. gimmeInteractive() requires the package gWidgetsRGtk2, which requires the GTK+ library. 
#' Please install this package before using gimmeInteractive().
#' @usage gimmeInteractive()
#' @author Hallie Pike
#' @export
gimmeInteractive <- function() {
  # added code to require the package RGtk2 and GTK+
  if (!requireNamespace("gWidgets2RGtk2", quietly = TRUE)) {
    stop("gWidgets2RGtk2 needed for this function to work. Please install it, along with the GTK+ library, before running gimmeInteractive().",
         call. = FALSE)
  }
  gui.env <- new.env()
  gui.env$Header <- TRUE
  gui.env$Plot <- TRUE
  gui.env$Ar <- FALSE
  gui.env$Subgroup <- FALSE
  gui.env$Paths <- NULL

  introwin <- gwindow(title="gimmeInteractive", visible=TRUE)
  gwin <- gwindow(title="gimme", visible=FALSE)
  iwin <- gwindow(title="indSEM", visible=FALSE)
  awin <- gwindow(title="aggregate", visible=FALSE)

  introgroup <- ggroup(horizontal=FALSE, container=introwin)
  agroup <- ggroup()
  bgroup <- ggroup(container=agroup)
  lyt <- glayout(container=agroup)

  gimmebutton <- gbutton("gimme",
                         container=introgroup)
  indSEMbutton <- gbutton("indSEM",
                          container=introgroup)
  aggbutton <- gbutton("aggSEM",
                       container=introgroup)

  size(gimmebutton) <- c(10,50)
  size(indSEMbutton) <- c(10,50)
  size(aggbutton) <- c(10,50)

  code1 <- lyt[18,1, anchor=c(1,0), expand=TRUE] <- glabel("data = ", container=lyt)
  code2 <- lyt[18,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code7 <- lyt[19,1, anchor=c(1,0), expand=TRUE] <- glabel("out = ", container=lyt)
  code8 <- lyt[19,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code3 <- lyt[20,1, anchor=c(1,0), expand=TRUE] <- glabel("sep = ", container=lyt)
  code4 <- lyt[20,2, anchor=c(-1,0), expand=TRUE] <- glabel(NULL, container=lyt)
  code5 <- lyt[21,1, anchor=c(1,0), expand=TRUE] <- glabel("header = ", container=lyt)
  code6 <- lyt[21,2, anchor=c(-1,0), expand=TRUE] <- glabel(TRUE, container=lyt)
  code9 <- lyt[22,1, anchor=c(1,0), expand=TRUE] <- glabel("ar = ", container=lyt)
  code10 <- lyt[22,2, anchor=c(-1,0), expand=TRUE] <- glabel(FALSE, container=lyt)
  code11 <- lyt[23,1, anchor=c(1,0), expand=TRUE] <- glabel("plot = ", container=lyt)
  code12 <- lyt[23,2, anchor=c(-1,0), expand=TRUE] <- glabel(TRUE, contaienr=lyt)
  code13 <- lyt[24,1, anchor=c(1,0), expand=TRUE] <- glabel("subgroup = ", container=lyt, visible=FALSE)
  code14 <- lyt[24,2, anchor=c(-1,0), expand=TRUE] <- glabel(FALSE, container=lyt, visible=FALSE)

  back <- gimage(stock.id="go-back",
                 dirname="stock",
                 container=bgroup,
                 size="large_toolbar",
                 anchor=c(0,1),
                 handler=function(h,...) {
                   #visible(introwin) <- TRUE
                   if (visible(gwin)==TRUE) delete(gwin, agroup)
                   if (visible(awin)==TRUE) delete(awin, agroup)
                   if (visible(iwin)==TRUE) delete(iwin, agroup)
                   visible(separ6) <- TRUE
                   visible(subgroupLabel) <- TRUE
                   visible(subgroup) <- TRUE
                   if (visible(iwin)==TRUE | visible(awin)==TRUE) {
                     visible(code13) <- FALSE
                     visible(code14) <- FALSE
                   }
                   visible(gwin) <- FALSE
                   visible(iwin) <- FALSE
                   visible(awin) <- FALSE
                 })

  datalabel <- lyt[2,1] <- glabel("Where are the data files located?",
                                  container=lyt)
  data <- lyt[2,2] <- gfilebrowse(text="Select a folder...",
                                  type="selectdir",
                                  quote=TRUE,
                                  container=lyt,
                                  action=code2,
                                  handler=function(h,...) {
                                    oldVal <- svalue(h$obj)
                                    str <- gsub("'", replacement="", oldVal, fixed=TRUE)
                                    svalue(h$action) <- dQuote(str)
                                    assign("Data", str, envir=gui.env)
                                  })
  separ1 <- lyt[3,1:2] <- gseparator(container=lyt)
  
  outLabel <- lyt[4,1] <- glabel("Where would you like results to be saved?", container=lyt)
  out <- lyt[4,2] <- gfilebrowse(text="Select a folder...",
                                 type="selectdir",
                                 quote=TRUE,
                                 container=lyt,
                                 action=code8,
                                 handler=function(h,...) {
                                   oldVal <- svalue(h$obj)
                                   str <- gsub("'", replacement="", oldVal, fixed=TRUE)
                                   svalue(h$action) <- dQuote(str)
                                   assign("Out", str, envir=gui.env)
                                 })
  
  separ3 <- lyt[5,1:2] <- gseparator(container=lyt)

  sepLabel <- lyt[6,1] <- glabel("What is the format of the data files?",
                                 container=lyt)
  sep <- lyt[6,2] <- gcombobox(c("", "Space Delimited", "Tab Delimited", "Comma Delimited"),
                               container=lyt,
                               action=code4,
                               handler=function(h,...) {
                                 oldVal <- svalue(h$obj)
                                 if (oldVal=="Space Delimited") str <- ""
                                 else if (oldVal=="Tab Delimited") str <- "/t"
                                 else if (oldVal=="Comma Delimited") str <- ","
                                 svalue(h$action) <- dQuote(str)
                                 assign("Sep", str, envir=gui.env)
                               })
  separ2 <- lyt[7,1:2] <- gseparator(container=lyt)

  headerLabel <- lyt[8,1] <- glabel("Do the data files have a header?",
                                    container=lyt)
  header <- lyt[8,2] <- gradio(c("Yes", "No"),
                               container=lyt,
                               action=code6,
                               horizontal=TRUE,
                               handler=function(h,...) {
                                 oldVal <- svalue(h$obj)
                                 if (oldVal=="Yes") bool <- TRUE
                                 else if (oldVal=="No") bool <- FALSE
                                 svalue(h$action) <- bool
                                 assign("Header", bool, envir=gui.env)
                               })

  separ5 <- lyt[9,1:2] <- gseparator(container=lyt)

  arLabel <- lyt[10,1] <- glabel("Should gimme start with autoregressive effects?",
                                 container=lyt)
  ar <- lyt[10,2] <- gradio(c("Yes", "No"),
                            container=lyt,
                            action=code10,
                            selected=2,
                            horizontal=TRUE,
                            handler=function(h,...) {
                              oldVal <- svalue(h$obj)
                              if (oldVal=="Yes") bool <- TRUE
                              else if (oldVal=="No") bool <- FALSE
                              svalue(h$action) <- bool
                              assign("Ar", bool, envir=gui.env)
                            })
  separ4 <- lyt[11,1:2] <- gseparator(container=lyt)

  plotLabel <- lyt[12,1] <- glabel("Would you like plots?", container=lyt)
  plot <- lyt[12,2] <- gradio(c("Yes", "No"),
                              container=lyt,
                              action=code12,
                              selected=1,
                              horizontal=TRUE,
                              handler=function(h,...) {
                                oldVal <- svalue(h$obj)
                                if (oldVal=="Yes") bool <- TRUE
                                else if (oldVal=="No") bool <- FALSE
                                svalue(h$action) <- bool
                                assign("Plot", bool, envir=gui.env)
                              })
  separ6 <- lyt[13,1:2] <- gseparator(container=lyt)

  subgroupLabel <- lyt[14,1] <- glabel("Should gimme subgroup individuals?",
                                       container=lyt)
  subgroup <- lyt[14,2] <- gradio(c("Yes", "No"),
                                  container=lyt,
                                  action=code14,
                                  selected=2,
                                  horizontal=TRUE,
                                  handler=function(h,...) {
                                    oldVal <- svalue(h$obj)
                                    if (oldVal=="Yes") bool <- TRUE
                                    else if (oldVal=="No") bool <- FALSE
                                    svalue(h$action) <- bool
                                    assign("Subgroup", bool, envir=gui.env)
                                  })
  separ7 <- lyt[15,1:2] <- gseparator(container=lyt)

  pathwayLabel <- lyt[16,1] <- glabel("Would you like to define pathways?",
                                      container=lyt)
  pathway <- lyt[16,2] <- gradio(c("Yes", "No"),
                                 container=lyt,
                                 horizontal=TRUE,
                                 handler=function(h,...){
                                   oldVal <- svalue(h$obj)
                                   if (oldVal=="Yes") bool <- TRUE
                                   else if (oldVal=="No") bool <- FALSE
                                   assign("Pathway", bool, envir=gui.env)
                                 })
  separ8 <- lyt[17,1:2] <- gseparator(container=lyt)

  srun <- lyt[25,1:2] <- gbutton("Define Structure",
                                 container=lyt,
                                 handler=function(h,...) {
                                   sswin <- gwindow("Lagged Effects", visible=TRUE)
                                   fgroup <- ggroup(horizontal=FALSE,
                                                    container=sswin)
                                   cgroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup)
                                   dgroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup)
                                   egroup <- ggroup(horizontal=FALSE,
                                                    container=fgroup,
                                                    visible=FALSE)
                                   backb <- gimage(stock.id="go-back",
                                                   dirname="stock",
                                                   container=cgroup,
                                                   size="large_toolbar",
                                                   anchor=c(-1,1),
                                                   handler=function(h,...) {
                                                     if (svalue(sswin)=="Lagged Effects") {
                                                       visible(sswin) <- FALSE
                                                     }
                                                     else if(svalue(sswin)=="Contemporaneous Effects") {
                                                       visible(egroup) <- FALSE
                                                       visible(dgroup) <- TRUE
                                                       svalue(sswin) <- "Lagged Effects"
                                                       delete(egroup, cg)
                                                       delete(egroup, ssrun)
                                                     }
                                                   })
                                   path <- paste(gui.env$Data, list.files(gui.env$Data)[1], sep="/")
                                   dir <- gsub("//", path, replacement="/", fixed=TRUE)
                                   dir <- gsub("\\", path, replacement="/", fixed=TRUE)
                                   assign("Dir", dir, envir=gui.env)
                                   values <- c(colnames(read.csv(gui.env$Dir, header=gui.env$Header)))
                                   ag <- glayout(container=dgroup, spacing=5)
                                   ag[1,1, anchor=c(1,0)] <- "Predicted:   "
                                   ag[1,2, anchor=c(-1,0)] <- "Predictor"
                                   for (i in 1:length(values)) {
                                     ag[i+1, 1, anchor=c(1,0)] <- paste(values[i], ":    ", sep="")
                                     ag[i+1, 2, anchor=c(-1,0)] <- gcheckboxgroup(values, horizontal=TRUE, container=ag)
                                   }
                                   cg <- glayout(spacing=5)
                                   cg[1,1, anchor=c(1,0)] <- "Predicted:   "
                                   cg[1,2, anchor=c(-1,0)] <- "Predictor"
                                   for (i in 1:length(values)) {
                                     cg[i+1, 1, anchor=c(1,0)] <- paste(values[i], ":    ", sep="")
                                     cg[i+1, 2, anchor=c(-1,0)] <- gcheckboxgroup(values[-i], horizontal=TRUE, container=cg)
                                   }
                                   crun <- gbutton("Next", container=dgroup)
                                   ssrun <- gbutton("Run gimme")

                                   pathhand <- addHandlerClicked(ssrun, handler=function(h,...){
                                     patha <- NULL
                                     pathc <- NULL
                                     for (i in 1:length(values)) {
                                       oldValA <- svalue(ag[i+1,2])
                                       oldValC <- svalue(cg[i+1,2])
                                       if (length(oldValA)!=0) {
                                         AVal <- paste(values[i], svalue(ag[i+1,2]), sep="~")
                                         patha <- c(patha, AVal)
                                         pathaa <- paste(patha, "lag", sep="")
                                         assign("patha", pathaa, envir=gui.env)
                                       }
                                       if (length(oldValC)!=0) {
                                         CVal <- paste(values[i], svalue(cg[i+1,2]), sep="~")
                                         pathc <- c(pathc, CVal)
                                         assign("pathc", pathc, envir=gui.env)
                                       }
                                     }
                                     oldVal <- c(gui.env$patha, gui.env$pathc)
                                     Val <- paste("\"",oldVal,"\"",collapse=",",sep="")
                                     ValA <- cat(Val)
                                     assign("Paths", c(ValA), envir=gui.env)
                                     enabled(ssrun) <- FALSE
                                     if (visible(gwin)==TRUE) {
                                       gimme(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
                                             out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot,subgroup=gui.env$Subgroup)
                                     }
                                     else if (visible(iwin)==TRUE) {
                                       indSEM(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
                                              out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot)
                                     }
                                     else if (visible(awin)==TRUE) {
                                       aggSEM(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
                                              out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot)
                                     }
                                   })

                                   addHandlerClicked(crun, action=ssrun, handler=function(h,...){
                                     svalue(sswin) <- "Contemporaneous Effects"
                                     visible(dgroup) <- FALSE
                                     visible(egroup) <- TRUE
                                     add(egroup, cg)
                                     add(egroup, ssrun)
                                     blockHandler(ssrun, pathhand)
                                     if (visible(iwin)==TRUE) svalue(h$action) <- "Run indSEM"
                                     else if (visible(awin)==TRUE) svalue(h$action) <- "Run aggSEM"
                                     unblockHandler(ssrun, pathhand)
                                   })

                                 })
  grun <- lyt[26,1:2] <- gbutton(NULL,
                                 container=lyt)

  enabled(srun) <- TRUE
  enabled(grun) <- FALSE

  grunhand <- addHandlerClicked(grun, handler=function(h,...){
    if (visible(gwin)==TRUE) {
      gimme(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
            out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot,subgroup=gui.env$Subgroup)
    }
    else if (visible(iwin)==TRUE) {
      indSEM(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
             out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot)
    }
    else if (visible(awin)==TRUE) {
      aggSEM(paths=gui.env$Paths,data=gui.env$Data,sep=gui.env$Sep,header=gui.env$Header,
             out=gui.env$Out,ar=gui.env$Ar,plot=gui.env$Plot)
    }
  })


  addHandlerSelect(pathway, handler=function(h,...) {
    if (svalue(pathway)=="Yes") {
      enabled(srun) <- TRUE
      enabled(grun) <- FALSE
    }
    else if (svalue(pathway)=="No") {
      enabled(srun) <- FALSE
      enabled(grun) <- TRUE
    }
  })

  addHandlerClicked(gimmebutton, handler=function(h,...){
    visible(gwin) <- TRUE
    add(gwin, agroup)
    visible(code13) <- TRUE
    visible(code14) <- TRUE
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run gimme"
    unblockHandler(grun, grunhand)
  }, action=grun)

  addHandlerClicked(indSEMbutton, handler=function(h,...){
    visible(iwin) <- TRUE
    visible(separ6) <- FALSE
    visible(subgroupLabel) <- FALSE
    visible(subgroup) <- FALSE
    visible(code13) <- FALSE
    visible(code14) <- FALSE
    add(iwin, agroup)
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run indSEM"
    unblockHandler(grun, grunhand)
  }, action=grun)

  addHandlerClicked(aggbutton, handler=function(h,...){
    visible(awin) <- TRUE
    visible(separ6) <- FALSE
    visible(subgroupLabel) <- FALSE
    visible(subgroup) <- FALSE
    visible(code13) <- FALSE
    visible(code14) <- FALSE
    add(awin, agroup)
    blockHandler(grun, grunhand)
    svalue(h$action) <- "Run aggSEM"
    unblockHandler(grun, grunhand)
  }, action=grun)
}


