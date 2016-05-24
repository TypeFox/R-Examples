library(SOMbrero) # Version 0.4

###############################################################################
## Global variables

# Max file input size :
options(shiny.maxRequestSize=30*1024^2)

# SOM training function
trainTheSom <- function(data, type, dimx, dimy, affectation, disttype, maxit,
                        varnames, rand.seed, scaling, eps0, init.proto, nb.save,
                        radiustype) {
  set.seed(rand.seed)
  
  if (type=="numeric")
    data <- data[, varnames]
  trainSOM(data, dimension=c(dimx,dimy), affectation=affectation, 
           dist.type=disttype, maxit=maxit, type=type, scaling=scaling, 
           init.proto=init.proto, nb.save=nb.save, radius.type=radiustype)
}

# List of somplot types options per SOM type and "what" :
all.somplot.types <- list("numeric"=
                            list("prototypes"=
                                   list("color", "3d", "lines", "barplot",
                                        "smooth distances"="smooth.dist",
                                        "polygon distances"="poly.dist",
                                        "grid distances"="grid.dist",
                                        "U matrix distances"="umatrix",
                                        "MDS"="mds", "radar"),
                                 "obs"=c("hitmap", "color", "lines", "barplot", 
                                          "names", "boxplot", "radar"),
                                 "energy"="Energy of backups"),
                          
                          "korresp"=
                            list("prototypes"=
                                   list("color", "3d", "lines", "barplot",
                                        "polygon distances"="poly.dist",
                                        "grid distances"="grid.dist",
                                        "U matrix distances"="umatrix",
                                        "MDS"="mds", "radar"),
                                 "obs"=c("hitmap", "names"),
                                 "energy"="Energy of backups"),
                          
                          "relational"=
                            list("prototypes"=
                                   list("lines", "barplot",
                                        "polygon distances"="poly.dist",
                                        "grid distances"="grid.dist",
                                        "U matrix distances"="umatrix",
                                        "MDS"="mds", "radar"),
                                 "obs"=c("hitmap", "names"),
                                 "energy"="Energy of backups"))
                                         
all.scplot.types <- list("numeric"=
                           list("prototypes"=
                                  list("dendrogram", "dendro3d", "color",
                                       "lines", "grid", "barplot", 
                                       "polygon distances"="poly.dist",
                                       "MDS"="mds", "radar"),
                                "obs"=c("hitmap", "color", "lines", "barplot", 
                                         "boxplot", "radar", "grid")),
                         "korresp"=
                           list("prototypes"=
                                  list("dendrogram", "color", "lines", "grid", 
                                       "barplot", 
                                       "polygon distances"="poly.dist",
                                       "MDS"="mds", "radar", "dendro3d"),
                                "obs"="hitmap"),
                         "relational"=
                           list("prototypes"=
                                  list("dendrogram", "lines", "barplot", "grid",
                                       "polygon distances"="poly.dist",
                                       "MDS"="mds", "radar", "dendro3d"),
                                "obs"="hitmap"))

###############################################################################

# Server
shinyServer(function(input, output, session) {

  #############################################################################
  ## Server variables
  server.env <- environment() # used to allocate in functions
  current.som <- NULL # this variable will contain the current SOM
  current.table <- NULL
  #############################################################################
  
  #### Panel 'About' (left hand side)
  ############################################################################## 
  # Output the sombrero logo :
  output$sombrero.logo <- renderImage(list(src="sombrero.png"), 
                                      deleteFile=FALSE)
  
  # Output the SAMM logo :
  output$samm.logo <- renderImage(list(src="samm.png"), deleteFile=FALSE)
  # Output the SAMM logo :
  output$miat.logo <- renderImage(list(src="miat.png"), deleteFile=FALSE)
  
  #### Panel 'Import data'
  ##############################################################################
  dInput <- reactive({
    in.file <- input$file1
    
    if (is.null(in.file))
      return(NULL)
    
    the.sep <- switch(input$sep, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec, "Period"=".", "Comma"=",")
    if (input$rownames) {
      the.table <- read.table(in.file$datapath, header=input$header, 
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else {
      the.table <- read.table(in.file$datapath, header=input$header, 
                              sep=the.sep, quote=the.quote, dec=the.dec)
    }
    
    # update the "input variables" checkbox (if somtype is numeric or integer)
    if (input$somtype =="numeric") {
      output$varchoice <- renderUI(
        checkboxGroupInput(inputId="varchoice", label="Input variables:",
                           choices=as.list(colnames(the.table)),
                           selected=as.list(colnames(the.table)[
                             sapply(the.table, class) %in%
                               c("integer", "numeric")])))
    } else output$varchoice <- renderText("")

    # update the map dimensions
    updateNumericInput(session, inputId="dimx", value={
      if (input$somtype =="korresp") {
        max(5,min(10,ceiling(sqrt((nrow(the.table)+ncol(the.table))/10))))
      } else max(5,min(10,ceiling(sqrt(nrow(the.table)/10))))
    })
    updateNumericInput(session, inputId="dimy", value={
      if (input$somtype =="korresp") {
        max(5,min(10,ceiling(sqrt((nrow(the.table)+ncol(the.table))/10))))
      } else max(5,min(10,ceiling(sqrt(nrow(the.table)/10))))
    })
    
    # update the max. iterations option
    updateNumericInput(session, "maxit", value=5 * nrow(the.table))

    # return the table
    server.env$current.table <- the.table
    the.table
  })
  
  # data preview table
  output$view <- renderTable({
    d.input <- dInput()
    if (is.null(d.input)) 
      return(NULL)
    if (ncol(d.input)>input$ncol.preview) 
      d.input <- d.input[,1:input$ncol.preview]
    head(d.input, n=input$nrow.preview) 
  })

  #### Panel 'Self-organize'
  #############################################################################

  # update the scaling option when input$somtype is changed
  output$scaling <- renderUI({
    selectInput(inputId="scaling", label="Data scaling:", 
                choices=switch(input$somtype,
                                "numeric"=c("unitvar", "none", "center"),
                                "korresp"=c("chi2"),
                                "relational"=c("none","cosine")),
                selected=switch(input$somtype, "numeric"="unitvar",
                                 "korresp"="chi2", "relational"="none"))
  })
  
  # update the scaling option when input$radiustype is changed
  output$disttype <- renderUI({
    selectInput(inputId="disttype", label="Distance scaling:", 
                choices=switch(input$radiustype,
                                "letremy"=c("letremy", "maximum", "euclidean",
                                             "manhattan", "canberra", "binary",
                                             "minkowski"),
                                "gaussian"=c("maximum", "euclidean", 
                                             "manhattan", "canberra", "binary",
                                             "minkowski")),
                selected=switch(input$radiustype, "letremy"="letremy",
                                 "gaussian"="euclidean"))
  })
  
  # update the initialization method when input$somtype is changed
  output$initproto <- renderUI({
    selectInput("initproto", label="Prototypes initialization method:", 
                choices=c("random","obs","pca"), 
                selected=switch(input$somtype, 
                                 "numeric"="random",
                                 "korresp"="random",
                                 "relational"="obs"))
  })

  # Train the SOM when the button is hit
  theSom<- function() {
    input$trainbutton
    server.env$current.som <- isolate(trainTheSom(current.table, input$somtype, 
                                                  input$dimx, input$dimy, 
                                                  input$affectation,
                                                  input$disttype, input$maxit, 
                                                  varnames=input$varchoice, 
                                                  rand.seed=input$randseed, 
                                                  scaling=input$scaling, 
                                                  eps0=input$eps0, 
                                                  init.proto=input$initproto, 
                                                  nb.save=input$nb.save,
                                                  radiustype=input$radiustype))
    
    updatePlotSomVar() # update variable choice for som plots
    updatePlotScVar() # update variable choice for sc plots

    # return the computed som
    server.env$current.som
  }

  # Render the summary of the SOM
  output$summary <- renderPrint({
    if (is.null(input$file1))
      return("First import a dataset.")
    if (input$trainbutton==0) {
      return("Hit the Train button to train the map.")
    }
    summary(theSom())
  })

  # Output the computed som object to be downloaded
  # TODO: output an error if map not trained
  output$som.download <- {
    downloadHandler(filename=function() {
        paste0("som ",format(Sys.time(),format="-%Y-%m-%d_%H:%M"),".rda",sep="")
      },
      content=function(file) {
        som.export <- server.env$current.som
        save(som.export, file=file)
      })
  }
  
  #### Panel 'Plot'
  ##############################################################################
  
  # Adapt plottype to the somtype and the "what" arguments
  observe({
    updateSelectInput(session, "somplottype", 
                      choices=all.somplot.types[[input$somtype]][[
                        input$somplotwhat]])
  })
  
  # update variables available for plotting
  updatePlotSomVar <- function() observe({
    tmp.names <- colnames(current.som$data)
    if (input$somtype =="korresp")
      tmp.names <- c(tmp.names, rownames(current.som$data))
    updateSelectInput(session, "somplotvar", choices=tmp.names)
    updateSelectInput(session, "somplotvar2", choices=tmp.names, 
                      selected=tmp.names[1:min(5,length(tmp.names))])
  })
  
  # Plot the SOM
  output$somplot <- renderPlot({
    if(is.null(current.table))
      return(NULL)
    if(input$trainbutton ==0)
      return(NULL)
    
    tmp.view <- NULL
    if (input$somtype =="korresp")
      tmp.view <- input$somplotrowcol
    
    if (input$somplotwhat =="energy")
      plot(current.som, what=input$somplotwhat)
    
    if (input$somplottype =="radar") {
      plot(x=current.som, what=input$somplotwhat, type=input$somplottype,
           variable=input$somplotvar, print.title=input$somplottitle,
           view=tmp.view, key.loc=c(-1,2), mar=c(0,10,2,0))
    } else {
      if (input$somplottype =="boxplot") {
        tmp.var <- (1:ncol(current.som$data))[colnames(current.som$data) %in% 
                                              input$somplotvar2]
      } else tmp.var <- input$somplotvar
    plot(x=current.som, what=input$somplotwhat, type=input$somplottype,
         variable=tmp.var, print.title=input$somplottitle,
         view=tmp.view)
    }
  })
  
  #### Panel 'Superclass'
  ##############################################################################
  # Input number of superclasses or cutting height
  output$scHorK <- renderUI(
    switch(input$sc.cut.choice, 
           "Number of superclasses"=
             numericInput("sc.k", "Number of superclasses:", 2, min=1,
                          max=input$dimx*input$dimy-1), 
           "Height in dendrogram"=
             numericInput("sc.h", "Height in dendrogram:", 10, min=0))
  )

  # Compute superclasses when the button is hit
  computeSuperclasses <- reactive({
    if (is.null(current.table))
      return(NULL)
    if (input$superclassbutton==0)
      return(NULL)
    
    isolate(switch(input$sc.cut.choice, 
                   "Number of superclasses"=
                     superClass(sommap=current.som, k=input$sc.k),
                   "Height in dendrogram"=
                     superClass(sommap=current.som, h=input$sc.h)))
  })

  output$sc.summary <- renderPrint({
    if (input$superclassbutton==0)
      return("Hit the 'Compute superclasses' button to see the results.")
    tmp.sc <- computeSuperclasses()
    summary(tmp.sc)
  })

  # Download the superclass classification
  # TODO: output an error if map not trained
  output$sc.download <- {
    
    
    downloadHandler(
      filename=function() {
        paste0("superclasses ", format(Sys.time(), format="%Y-%m-%d_%H:%M"),
              ".csv", sep="")
      },
      content=function(file) {
        the.sc <- computeSuperclasses()
        classes.export <- data.frame(obs=row.names(the.sc$som$data),
                                     cluster=the.sc$cluster[
                                       the.sc$som$clustering])
        write.csv(classes.export, file=file, row.names=FALSE)
      })
  }

  # Adapt scplottype to the somtype and the "what" arguments
  observe({
    updateSelectInput(session, "scplottype", 
                      choices=all.scplot.types[[input$somtype]][[
                        input$scplotwhat]])
  })
  
  # update variables available for plotting
  updatePlotScVar <- function() observe({
    tmp.names <- colnames(current.som$data)
    if (input$somtype =="korresp")
      tmp.names <- c(tmp.names, rownames(current.som$data))
    updateSelectInput(session, "scplotvar", choices=tmp.names)
    updateSelectInput(session, "scplotvar2", choices=tmp.names, 
                      selected=tmp.names[1:min(5,length(tmp.names))])
  })
  
  # Update SuperClass plot
  output$scplot <- renderPlot({
    if(is.null(current.table))
      return(NULL)
    
    the.sc <- computeSuperclasses()
    if(input$superclassbutton ==0)
      return(NULL)

    if (input$scplottype %in% c("grid", "dendrogram"))
      return(plot(the.sc, type=input$scplottype))

    tmp.view <- NULL
    if (input$somtype =="korresp")
      tmp.view <- input$scplotrowcol
    
    if (input$scplottype =="radar") {
      plot(x=the.sc, what=input$scplotwhat, type=input$scplottype,
                  variable=input$scplotvar, view=tmp.view, key.loc=c(-1,2),
                  mar=c(0,10,2,0))
    } else { 
      if (input$scplottype =="boxplot") {
        tmp.var <- (1:ncol(current.som$data))[colnames(current.som$data) %in% 
                                                input$scplotvar2]
      } else tmp.var <- input$scplotvar
      
      plot(x=the.sc, what=input$scplotwhat, type=input$scplottype,
           variable=tmp.var, view=tmp.view)
    }
  })

  #### Panel 'Combine with additional data'
  ##############################################################################
  
  # File input for additional variables
  dInputAdd <- reactive({
    in.file <- input$file2
    
    if (is.null(in.file))
      return(NULL)
    
    the.sep <- switch(input$sep2, "Comma"=",", "Semicolon"=";", "Tab"="\t",
                      "Space"="")
    the.quote <- switch(input$quote2, "None"="","Double Quote"='"',
                        "Single Quote"="'")
    the.dec <- switch(input$dec2, "Period"=".", "Comma"=",")
    
    if (input$rownames2) {
      the.table <- read.table(in.file$datapath, header=input$header2, 
                              sep=the.sep, quote=the.quote, row.names=1,
                              dec=the.dec)
    } else the.table <- read.table(in.file$datapath, header=input$header2, 
                                   sep=the.sep, quote=the.quote, dec=the.dec)
    
    updateAddPlotVar() # update variable selector    
    
    the.table
  })
  
  # additional data preview table
  output$addview <- renderTable({
    d.input <- dInputAdd()
    if (is.null(d.input)) 
      return(NULL)
    if (ncol(d.input)>input$ncol.preview.add) 
      d.input <- d.input[,1:input$ncol.preview.add]
    head(d.input, n=input$nrow.preview.add) 
  })
  
  # Adapt available variables from second file
  updateAddPlotVar <- function() observe({
    d.input <- dInputAdd()
    if(is.null(d.input))
      return(NULL)
    updateSelectInput(session, "addplotvar", choices=colnames(d.input), 
                      selected=colnames(d.input)[1])
    updateSelectInput(session, "addplotvar2", choices=colnames(d.input), 
                      selected=colnames(d.input)[1:min(5,ncol(d.input))])
  })
  
  # function to render Additional data Plot
  output$addplot <- renderPlot({
    d.input <- dInputAdd()
    if (is.null(d.input)) return(NULL)
    
    if (input$addplottype %in% c("pie","color","names")) {
      tmp.var <- input$addplotvar
    } else tmp.var <- input$addplotvar2
    
    if(input$addplottype =="radar") {
      plot(x=current.som, what="add", type=input$addplottype, 
           variable=d.input[,tmp.var], key.loc=c(-1,2), mar=c(0,10,2,0))
    } else if (input$addplottype !="graph") {
      plot(x=current.som, what="add", type=input$addplottype, 
           variable=d.input[,tmp.var])
    } else {
      adjBin <- as.matrix(d.input!=0)
      tmpGraph <- graph.adjacency(adjBin, mode="undirected")
      plot(current.som, what="add", type="graph", variable=tmpGraph)
    }
  })
})
