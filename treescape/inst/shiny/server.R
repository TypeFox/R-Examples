## DEFINE THE SERVER SIDE OF THE APPLICATION
shinyServer(function(input, output, session) {
  ## LOAD PACKAGES
  if(!require("ade4")) stop("package ade4 is required")
  if(!require("adegenet")) stop("package adegenet is required")
  if(!require("adegraphics")) stop("package ade4 is required")
  if(!require("ape")) stop("package ape is required")
  if(!require("distory")) stop("package distory is required")
  if(!require("fields")) stop("package fields is required")
  if(!require("htmlwidgets")) stop("package htmlwidgets is required")
  if(!require("MASS")) stop("package MASS is required")
  if(!require("phangorn")) stop("package phangorn is required")
  if(!require("treescape")) stop("package treescape is required")
 
  # suppress warning messages from creating temporary directories when 3d plotting
  suppressWarnings(warning("dir.create(dir)"))
  
  # the following resets the DensiTree plot every time the number of clusters changes - it was really slow without this
  rvs <- reactiveValues(showDensiTree=NULL)
  observeEvent(input$nclust, {
    rvs$showDensiTree <- NULL
  })
  observeEvent(input$selectedDensiTree, {
    rvs$showDensiTree <- 1  
  })
  
  
  
  ######################################
  ### Define main reactive functions
  ######################################
  
  getDataType <- reactive({
    input$datatype
  })
  
  getDataSet <- reactive({
    dataType <- getDataType()
    if(dataType=="exDengue"){
      return("Dengue")
    }
    if(dataType=="exWoodmice"){
      return("woodmiceTrees")
    }
    else {
      # extract file name
      strsplit(input$datafile$name, '[.]')[[1]][1]
    }
  })
  
  getSampleSize <- reactive({
    input$sampleSize
  })
  
  getRandSamp <- reactive({
    input$randSamp
  })
  
  ## GET DATA ##
  getData <- reactive({
    out <- NULL
    dataType <- getDataType()
    samp <- NULL
    
    ## data is a distributed dataset
    if(dataType=="exDengue"){
      if (!exists("DengueTrees")) { 
        data("DengueTrees", package="treescape", envir=environment()) }
      out <- get("DengueTrees")
    }
    if(dataType=="exWoodmice"){
      if (!exists("woodmiceTrees")) {
        data("woodmiceTrees", package="treescape", envir=environment()) }
      out <- get("woodmiceTrees")
    }
    
    ## data is an input file
    if(dataType=="file" && !is.null(input$datafile)){
      ## need to rename input file
      oldName <- input$datafile$datapath
      extension <- adegenet::.readExt(input$datafile$name)
      newName <- paste(input$datafile$datapath, extension, sep=".")
      file.rename(oldName, newName)
      
      if(tolower(extension) %in% c("rdata","rda")){
        out <- get(load(newName))
      }
      if(tolower(extension) %in% c("rds")){
        out <- readRDS(file=newName)
      }
      if(tolower(extension) %in% c("nex", "nexus")){
        if(!require(ape)) stop("ape is required to read in NEXUS (.nex, .nexus) files")
        out <- read.nexus(file=newName)
      }
      
      l <- length(out)
      
      ## fix potential bug with input of two trees
      validate(
        need(l>2, "treescape expects at least three trees. The function treeDist is suitable for comparing two trees.")
      )
      
      # get a manageable number of trees by sampling if necessary
      randSamp <- getRandSamp()
      if(randSamp == TRUE){
        sampleSize <- getSampleSize()
        if (l>sampleSize) {
          updateSliderInput(session, "sampleSize", "Size of random sample:", value=sampleSize, min=10, max=l, step=10)
          samp <- sample(1:l,sampleSize)
          out <- out[samp]
        }
        else{ # could only happen initially if <=10 trees supplied
          updateSliderInput(session, "sampleSize", "Size of random sample:", value=l, min=3, max=l, step=1)
        }
        
      }
      
      ## fix potential bug with tip labels - they need to match
      tipLabelProblem <- FALSE
      for (i in 1:length(out)) {
        if (!setequal(out[[i]]$tip.label,out[[1]]$tip.label)) {
          tipLabelProblem <- TRUE
          validate(
            need(!tipLabelProblem, "Trees must have identical tip labels for the current version of treescape")
          )
        } 
      }
      
    }
    
    validate(
      need(!is.null(out), "Waiting for data")
    )
    
    ## fix potential bug with names - they need to be defined and unique
    if(is.null(names(out))) {names(out) <- 1:length(out)}
    if(length(unique(names(out)))!=length(out)){
      warning("duplicates detected in tree labels - using generic names")
      names(out) <- 1:length(out)
    }
    
    ## return data
    # need to pass on the sample so that metaData can be sampled too
    if(is.null(samp)) samp <- 1:length(out)
    
    return(list(out=out,samp=samp))
  }) # end getData
  
  ## GET number of trees
  getLengthData <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(length(x))
  })
  
  ## GET tree names
  getTreeNames <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(names(x))
  })
  
  ## GET tip labels
  getTipLabels <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    return(x[[1]]$tip.label)
  })
  
  ## GET tree method
  getTreemethod <- reactive({
    input$treemethod
  }) # end getTreemethod
  
  ## GET number of axes retained
  getNaxes <- reactive({
    if(is.null(input$naxes)){
      naxes <- 3
    }
    else {
      naxes <- as.numeric(input$naxes)
      # when naxes changes we update the options available for the axes 
      # unfortunately I think they have to reset to their original 1,2,3 values
      # but at least they now only do this when naxes changes; they used to also do it for lambda etc.
      
      updateNumericInput(session,"xax", "Indicate the x axis", value=1, min=1, max=naxes)
      updateNumericInput(session,"yax", "Indicate the y axis", value=2, min=1, max=naxes)
      
      # (if relevant, update z axis selector too)
      dim <- getPlotDim()
      if (dim==3){
        updateNumericInput(session,"zax", "Indicate the z axis", value=3, min=1, max=naxes)
      }
      
    } 
    return(naxes)  
  }) # end getNaxes
  
  ## GET lambda
  getLambda <- reactive({
    l <- input$lambda
    ## the following removes the lambda error messages:
    validate(
      need(!is.null(l), "Loading data set")
    )	
    return(l)
  }) # end getLambda
  
  getTipsToEmphasise <- reactive({
    input$whichTips
  })
  
  getEmphWeight <- reactive({
    input$emphWeight
  })
  
  # GET the tree vectors as functions of lambda
  getKCtreeVecs <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    tips <- getTipsToEmphasise()
    weight <- getEmphWeight()
    df <- sapply(x, function(i) treeVec(i, return.lambda.function=TRUE, emphasise.tips=tips, emphasise.weight = weight)) 
  })
  
  
  # GET the tree vectors evaluated at lambda
  getKCtreeVecsAtLambda <- reactive({
    vectors <- getKCtreeVecs()
    l <- getLambda()
    validate(
      need(!is.null(vectors), "Analysing data")
    )
    t(sapply(vectors, function(i) i(l)))
  })
  
  
  ## GET KC matrix, evaluated at lambda
  getKCmatrix <- reactive({
    vls <- getKCtreeVecsAtLambda()
    as.dist(rdist(vls))
  }) # end getKCmatrix
  
  ## GET medTrees for all clusters
  getMedTreesList <- reactive({
    mat <- getKCtreeVecsAtLambda()
    groves <- getClusters()
    if(!is.null(groves$groups)){ # if clusters have been picked
      numGroups <- length(unique(groves$groups))
      med <- medTree(mat,groves$groups)
      lapply(1:numGroups, function(x) med[[x]]$treenumbers[[1]])
    }
    else{
      medTree(mat)$treenumbers[[1]]
    }
  })
  
  getMedTree <- reactive({
    data <- getData()
    x <- data$out
    whichClust <- input$selectedMedTree
    medList <- getMedTreesList()
    if(whichClust=="all"){
      x[[medList[[1]]]]
    }
    else{
      x[[medList[[as.numeric(whichClust)]]]]
    }
  })
  
  ## GET PCO analysis ##
  getPCO <- reactive({
    D <- getKCmatrix()
    naxes <- getNaxes()
    validate(
      need(!is.null(D), "Analysing data")
    )
    validate(
      need(!is.null(naxes), "Analysing data")
    )
    dudi.pco(D,scannf=FALSE,nf=naxes)
  }) # end getPCO
  
  ## GET ANALYSIS ##
  getAnalysis <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    
    naxes <- getNaxes()
    TM <- getTreemethod()
    
    ## select method used to summarise tree
    if(!is.null(TM)){
      if(TM %in% c("patristic","nNodes","Abouheif","sumDD")){
        ## run treescape
        res <- treescape(x, method=TM, nf=naxes)
      } 
      else if(TM=="metric"){
        ## don't actually need to call treescape here, to save on recomputation for varying lambda
        D <- getKCmatrix()
        pco <- getPCO()
        res <- list(D=D, pco=pco) 
      } 
      else if(TM=="RF"){
        # suppress unrooted warning
        D <- suppressWarnings(RF.dist(x))
        # suppress non-Euclidean distance warning
        pco <- suppressWarnings(dudi.pco(D, scannf=FALSE, nf=naxes))
        res <- list(D=D, pco=pco) 
      }
      else if(TM=="BHV"){
        D <- dist.multiPhylo(x)
        # suppress non-Euclidean distance warning
        pco <- suppressWarnings(dudi.pco(D, scannf=FALSE, nf=naxes))
        res <- list(D=D, pco=pco) 
      }
    }
    
    ## return results
    return(res)
  }) # end getAnalysis
  
  #################################################
  ### Little "get" functions to support getClusters
  #################################################
  
  getNclust <- reactive({
    if(!is.null(input$nclust)) {
      input$nclust
    } else {
      2
    }
  }) 
  
  getClustmethod <- reactive({
    input$clustmethod
  })
  
  
  ################
  ## GET CLUSTERS
  ################
  
  getClusters <- reactive({
    ## stop if clusters not required
    if(!input$findClusters) return(NULL)
    else if(input$clusterType=="meta") return(NULL)
    
    ## reset the densiTree plot to accommodate number of clusters available
    choices <- getClustChoices()
    updateSelectInput(session, "selectedDensiTree", "Choose collection of trees to view in densiTree plot", 
                      choices=choices, selected="")
    
    ## reset the median tree choices to accommodate number of clusters available
    updateSelectInput(session, "selectedMedTree", "Median tree from:", 
                      choices=choices, selected="all")
    
    ## get dataset
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    
    naxes <- getNaxes()
    TM <- getTreemethod()
    nclust <- getNclust()
    clustmethod <- getClustmethod()
    
    ## select method used to summarise tree
    if(!is.null(TM)){
      if(TM %in% c("patristic","nNodes","Abouheif","sumDD")){
        ## run findGroves
        res <- findGroves(x, method=TM, nf=naxes, nclust=nclust, clustering=clustmethod)
      } else if(TM %in% c("metric","BHV","RF")){
        res <- findGroves(getAnalysis(), nclust=nclust, clustering=clustmethod)
      } 
    }
    
    ## return results
    return(res)
    
  }) # end getClusters
  
  
  
  ## DYNAMIC UI COMPONENTS ##
  ## SELECTION OF MDS AXES
  output$naxes <- renderUI({
    if(!is.null(getLengthData())) {
      nmax <- getLengthData()
    } else {
      nmax <- 100
    }
    sliderInput("naxes", "Number of MDS axes retained:", min=2, max=nmax, value=3, step=1)
  })
  
  ## VALUE OF LAMBDA FOR METRIC
  output$lambda <- renderUI({
    ## if KC metric has been chosen
    TM <- getTreemethod()
    if(TM=="metric") {
      sliderInput("lambda", "Value of lambda", min=0, max=1, value=0, step=0.01)
    } else {
      NULL
    }
  })
  
  ## SELECTION OF NUMBER OF CLUSTERS
  output$nclust <- renderUI({
    if(!is.null(data <- getData())) {
      nmax <- length(data$out)
    } else {
      nmax <- 100
    }
    nmax <- min(20, nmax)
    sliderInput("nclust", "Number of clusters:", min=2, max=nmax, value=2, step=1)
  })
  
  ## SELECTION OF TIPS
  output$whichTips <- renderUI({
    # populate selection box with tip choices
    allTips <- getTipLabels()
    choices <- c("",allTips)
    names(choices) <- c("Type here to search tip names",allTips)
    selectInput("whichTips", "Select one or more tips to emphasise:", 
                choices=choices, selected=NULL, selectize=TRUE, multiple=TRUE)
  })
  
  ## GET METADATA ## for colouring trees by type
  getMetaData <- reactive({
    out <- NULL
    data <- getData()
    samp <- data$samp
    ## data is an input file
    if(input$clusterType=="meta" && !is.null(input$metadatafile)){
      ## need to rename input file
      oldName <- input$metadatafile$datapath
      extension <- adegenet::.readExt(input$metadatafile$name)
      newName <- paste(input$metadatafile$datapath, extension, sep=".")
      file.rename(oldName, newName)
      
      if(tolower(extension) %in% c("rdata","rda")){
        out <- get(load(newName))
        validate(
          need(class(out)%in%c("numeric","character","list","factor"), paste0("The class of the input is ", class(out), ". Please upload a single object of class list, numeric, factor or character, whose length is the same as the number of trees."))
        )
        
      }
      if(tolower(extension) %in% c("csv")){
        csvfile <- read.csv(file=newName, header=FALSE)
        out <- csvfile[,1]
        validate(
          need(class(out)%in%c("numeric","character","list","factor"), paste0("The first column of the csv file has been extracted. However, the class of the input is ", class(out), ". Please alter the entries so that it can be read by R as an object of class list, numeric, factor or character, whose length is the same as the number of trees."))
        )
      }
      
      if(class(out)=="list") {out <- unlist(out)}
      
      l <- getLengthData()
      out <- out[samp]
      validate(
        need(length(out)==l, paste0("The length of the metadata must be the same as the number of trees, which is ", l, ". However, the length of the input is ", length(out)))
      )
      
    }
    
    ## return metadata
    return(out)
  }) # end getMetaData
  
  
  ######################################################
  ### Little "get" functions to support getPlot
  ######################################################  
  
  getPalette  <- reactive({
    get(input$palette)
  })
  
  getLabcol <- reactive({
    ifelse(!is.null(input$labcol), input$labcol, "black")
  })
  
  getBgcol <- reactive({
    ifelse(!is.null(input$bgcol), input$bgcol, "white")
  })
  
  getXax <- reactive({
    input$xax
  })  
  
  getYax <- reactive({
    input$yax
  })  
  
  getZax <- reactive({
    input$zax
  })  
  
  getShowlabels <- reactive({
    input$showlabels
  })
  
  getLabelsize <- reactive({
    input$labelsize
  })
  
  getPointsize <- reactive({
    input$pointsize
  })
  
  getPlotFunction <- reactive({
    input$graphics
  })
  
  ##############  
  ## GET plot
  ##############
  
  ## GET whether plot is 2D (default) or 3D
  getPlotDim <- reactive({
    plotDim <- input$plot3D
    if(is.null(plotDim)) {2} # needed during startup
    else {return(plotDim)}
  })
  
  ## GET 2D plot
  getPlot <- reactive({
    
    res <- getAnalysis()
    pal <- getPalette()
    labcol <- getLabcol()
    groves <- getClusters()
    treeTypes <- getMetaData()
    showlabels <- getShowlabels()
    pointSize <- getPointsize()
    
    if(!is.null(treeTypes)) {
      groups <- treeTypes
      cols <- fac2col(1:length(unique(groups)),col.pal=pal)
    }
    else if (!is.null(groves)) {
      groups <- groves$groups
      cols <- fac2col(1:length(unique(groups)),col.pal=pal)
    }
    else {
      groups <- NULL
      n <- getLengthData()
      cols <- rep(labcol, n)
    }
    
    
    ## get aesthetics
    xax <- getXax()
    yax <- getYax()
    
    plotFunction <- getPlotFunction()
    
    if (plotFunction==1) {
      transitions <- input$transitions
      
      # labels and tree names
      treeNames <- getTreeNames()
      if (is.null(groups)) { tooltips <- paste0("Tree ", treeNames) }
      else { tooltips <- paste0("Tree ",treeNames,", cluster ",groups) }
      
      treeLabels <- NULL
      labelsize <- NULL
      
      if(showlabels==TRUE) {
        treeLabels <- getTreeNames()
        labelsize <- getLabelsize()
      }
      
      pointOpacity <- input$pointopacity
    
      plot <- plotGrovesD3(res$pco, xax=xax, yax=yax,
                 treeNames=treeLabels, labels_size=labelsize*5,
                 point_size = pointSize*40, point_opacity = pointOpacity,
                 groups=groups, colors=cols, col_lab="Cluster",
                 xlab=paste0("Axis ",xax), ylab=paste0("Axis ",yax),
                 tooltip_text = tooltips,
                 transitions=transitions, legend_width=50
      ) 
      # later could add:
      # other categories of variation e.g. metadata using symbols
    }
    
    else { # i.e. plotFunction==2
      bgcol <- getBgcol()
      scattertype <- input$scattertype
      screemds <- input$screemds
      optimlabels <- input$optimlabels
      labelsize <- getLabelsize()
      
      if(is.null(groves)){
        plot <- plotGroves(res$pco, groups=treeTypes, type=scattertype, xax=xax, yax=yax,
                   scree.posi=screemds, lab.optim=optimlabels,
                   lab.show=showlabels, lab.cex=labelsize,
                   lab.col=labcol,
                   point.cex=pointSize, bg=bgcol, col.pal=pal)
      } 
      else {
        ## plot with statistically identified groups
        plot <- plotGroves(groves, type=scattertype, xax=xax, yax=yax,
                   scree.posi=screemds, lab.optim=optimlabels,
                   lab.show=showlabels, lab.cex=labelsize,
                   lab.col=labcol,
                   point.cex=pointSize, bg=bgcol, col.pal=pal)
      }
    }
  return(plot)  
  })
  
  getDistPlot <- reactive({
    res <- getAnalysis()
    refTree <- input$selectedRefTree
    validate(
      need(refTree!="", "Select a reference tree")
    )
    groves <- getClusters()
    treeNames <- getTreeNames()
    pal <- getPalette()
    dists <- as.matrix(res$D)[refTree,] 
    g1 <- s1d.label(dists, labels=treeNames, poslabel="regular", p1d.horizontal=FALSE, p1d.reverse=TRUE, plot=FALSE)
    if(!is.null(groves$groups)){
      pal <- getPalette()
      nclusts <- getNclust()
      ordercols <- fac2col(1:nclusts, col.pal=pal)
      g2 <- s1d.boxplot(dists,fac=groves$groups, col=ordercols, p1d.horizontal=FALSE, plot=FALSE)
      ADEgS(c(g1, g2), layout = c(1, 2))
    }
    else{
      g1
    }
    
  })
  
  getPlotType <- reactive({
    input$plotType
  })
  
  ## TREESCAPE IMAGE ##
  output$treescapePlot <- renderUI({
    type <- getPlotType()
    if (type==1){ # i.e. full tree landscape
      plotFunction <- getPlotFunction()
      if (plotFunction==1) { # i.e. scatterD3
        scatterD3Output("scatterplotD3")
      }
      else { # i.e. adegraphics
        plotOutput("scatterplot", height = "800px") 
      }
    }
    else{ # i.e. distance from reference tree plot
      i <- input$stretch
      height <- as.character(paste0(i,"px"))
      plotOutput("DistPlot", height = height)  
    }
  })
  
  output$treescapePlot3D <- renderRglwidget({
    validate(
      need(packageVersion("rglwidget")>='0.1.1433',
           paste0("You are running version ",packageVersion("rglwidget")," of the package rglwidget, which contains a bug for 3D plotting. Please update to the latest version by running: install.packages('rglwidget', repos='http://R-Forge.R-project.org')")
      ))
    plot <- getPlot3d()
    plot
    rglwidget()
  })
  
  
  output$scatterplotD3 <- renderScatterD3({
    plotFunction <- getPlotFunction() # need to do this or you get an error when switching between plotGroves and plotGrovesD3
    if (plotFunction==1) { 
    withProgress(message = 'Loading plot',
                 value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                   }
                   myplot <- getPlot()
                   myplot
                 })
    }
  })
  
  output$scatterplot <- renderPlot({
    plotFunction <- getPlotFunction() # need to do this or you get an error when switching between plotGroves and plotGrovesD3
    if (plotFunction==2) {
    withProgress(message = 'Loading plot',
                 value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                   }
                   myplot <- getPlot()
                   myplot
                 })
    }
  }, res=120)
  
  output$DistPlot <- renderPlot({
      withProgress(message = 'Loading plot',
                   value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                     }
                     myplot <- getDistPlot()
                     plot(myplot)
                   })
  }, res=120)
  
  getPlot3d <- reactive({
    res <- getAnalysis()
    xax <- getXax()
    yax <- getYax()
    zax <- getZax()
    col <- getLabcol()
    
    # show clusters?
    clusts <- getClusters()
    treeTypes <- getMetaData()
    if (!is.null(clusts)){
      pal <- getPalette()
      cols3d <- fac2col(clusts$groups,col.pal=pal)
    }
    else if (!is.null(treeTypes)) {
      pal <- getPalette()
      cols3d <- fac2col(treeTypes,col.pal=pal)
    } 
    else{cols3d <- col}
    
    rgl::plot3d(res$pco$li[,xax],res$pco$li[,yax],res$pco$li[,zax], 
                type="s", size=getPointsize(),
                xlab="",ylab="",zlab="",
                col=cols3d, add=FALSE)
  })
  
  
  ## make Shepard plot
  getShep <- reactive({
    res <- getAnalysis()
    dim <- getPlotDim()
    xax <- getXax()
    yax <- getYax()
    
    if (dim==2) {  shep <- Shepard(res$D,as.matrix(res$pco$li[,xax],res$pco$li[,yax])) }
    
    else {
      zax <- getZax()
      shep <- Shepard(res$D,as.matrix(res$pco$li[,xax],res$pco$li[,yax],res$pco$li[,zax]))
    }
  })
  
  output$shepPlot <- renderPlot({
      withProgress(message = 'Loading Shepard plot',
                   value = 0, {
                     for (i in 1:15) {
                       incProgress(1/15)
                     }
                     shep <- getShep()
                     labcol <- getLabcol()
                     plot(shep, pch=19, cex=0.5, col=labcol, xlab="Distance in tree space", ylab="MDS distance")
                     
                   })
  }, res=120)
  
  output$shep <- renderUI({
    plotOutput("shepPlot", width="800px", height="800px")
  })
  
  
  ## make screeplot
  output$screePlot <- renderPlot({
    res <- getAnalysis()
    labcol <- getLabcol()
    barplot(res$pco$eig, col=labcol)
  }, res=120)
  
  output$scree <- renderUI({
    plotOutput("screePlot")
  })
  
  # get tree and aesthetics for plotting tree  
  getTreeChoice <- reactive({
    input$treeChoice
  })
  
  
  getTree <- reactive({
    data <- getData()
    x <- data$out
    validate(
      need(!is.null(x), "Loading data set")
    )
    treechoice <- getTreeChoice()
    if(treechoice=="med"){
      tre <- getMedTree()
    }
    else{
      g <- input$selectedGenTree
      validate(
        need(g!="", "Select tree to view")
      )
      treeNum <- as.numeric(g)
      tre <- x[[treeNum]]
    }
    
    # return tree
    if(!is.null(tre)){
      if(input$ladderize){
        tre <- ladderize(tre)
      }
      return(tre)   
    }
    else{
      NULL
    }
  })  
  
  ## PHYLOGENY ##
  output$tree <- renderPlot({
    tre <- getTree()
    if(!is.null(tre)){
      
      ## plot tree ##
      par(mar=rep(2,4), xpd=TRUE)
      plot(tre, type=input$treetype,
           use.edge.length=as.logical(input$edgelengths),
           show.tip.lab=input$showtiplabels, 
           font=as.numeric(input$tiplabelfont), 
           cex=input$tiplabelsize,
           direction=input$treedirection,
           edge.width=input$edgewidth,
           edge.color=input$edgecolor
      )
    }
  })
  
  ## DENSITREE
  
  # The slider bar is always at least 2 even when clusters haven't
  # been requested, so we can't just use getNclust.
  
  getNclustForDensiTree <- reactive({
    if(input$clusterType=="meta"){NULL}
    else{input$nclust}
  }) 
  
  getClustChoices <- reactive({
    nclust <- getNclustForDensiTree()
    if(is.null(nclust)){
      choices <- c("","all")
      names(choices) <- c("Choose one","All trees")
    }
    else{
      choices <- c("",1:nclust,"all")
      names(choices) <- c("Choose one",paste0("Cluster ",1:nclust),"All trees")
    }
    return(choices)
  })
  
  getDensiTree <- reactive({
    clusterNo <- input$selectedDensiTree
    if(clusterNo==""){
      NULL
    }
    else if(clusterNo=="all"){
      data <- getData()
      x <- data$out
      medList <- getMedTreesList()
      med <- x[[medList[[1]]]]
      return(list(trees=x,con=med))
    }
    else{
      data <- getData()
      x <- data$out
      clusts <- getClusters()
      clustTrees <- x[which(clusts$groups==as.numeric(clusterNo))]
      medList <- getMedTreesList()
      med <- x[[medList[[as.numeric(clusterNo)]]]]
      return(list(trees=clustTrees, con=med))
    }
  })  
  
  output$densiTree <- renderPlot({
    if(is.null(rvs$showDensiTree)) {NULL}
    else{
      withProgress(message = 'Loading densiTree plot',
                   detail = 'Note: the final stage of this process may take a while for large sets of trees',
                   value = 0, {
                     for (i in 1:30) {
                       incProgress(1/30)
                     }
                     clustTrees <- getDensiTree()
                     densiTree(clustTrees$trees, col=4, consensus=clustTrees$con, alpha=input$alpha, scaleX=input$scaleX)
                   })
    }
  })
  
  
  ## EXPORT TREES ##
  output$exporttrees <- downloadHandler(
    filename = function() { paste(getDataSet(), '.nex', sep='') },
    content = function(file) {
      if(!require(ape)) stop("ape is required to save trees into nexus file")
      data <- getData()
      x <- data$out
      if(!is.null(x) && inherits(x, "multiPhylo")) ape::write.nexus(x, file=file)
    })
  
  ## EXPORT ANALYSIS TO CSV ##
  output$exportrestocsv <- downloadHandler(
    filename = function() { paste(getDataSet(), "-analysis", '.csv', sep='') },
    content = function(file) {
      data <- getData()
      x <- data$out
      res <- getClusters()
      if(!is.null(res)){
        tab <- cbind.data.frame(res$groups, res$treescape$pco$li)
        names(tab) <- c("cluster", paste("PC", 1:ncol(res$treescape$pco$li), sep="."))
        row.names(tab) <- names(x)
      } else{
        res <- getAnalysis()
        tab <- res$pco$li
        names(tab) <- paste("PC", 1:ncol(tab), sep=".")
        row.names(tab) <- names(x)
      }
      if(!is.null(res)) write.csv(tab, file=file)
    })
  
  
  ## EXPORT ANALYSIS TO RDATA ##
  output$exportrestordata <- downloadHandler(
    filename = function() { paste(getDataSet(), "-analysis", '.RData', sep='') },
    content = function(file) {
      data <- getData()
      trees <- data$out
      analysis <- getClusters()
      if(is.null(analysis)) analysis <- getAnalysis()
      if(!is.null(analysis)) {
        save(trees, analysis, file=file)
      }
    })
  

  ## EXPORT 2D plotGroves MDS PLOT AS png ##
  output$downloadMDS <- downloadHandler(
    filename = function() { 
      paste0(getDataSet(),"scape2D.png") 
    },
    content = function(file) {
      myplot <- getPlot()
      png(file=file, width = 10, height = 10, units = 'in', res = 500)
      plot(myplot)
      dev.off()
      
      contentType = 'image/png'  
    }
  )
  
  ## EXPORT 2D plotGrovesD3 PLOT AS html ##
  output$downloadMDS2Dhtml <- downloadHandler(
    filename = function() { 
      paste0(getDataSet(),"scape2D.html") 
    },
    content = function(file) {
      htmlwidgets::saveWidget(
        getPlot(),
        file=file, 
        selfcontained = TRUE)
    },
    contentType = 'html'  
  )
  

  ## EXPORT 3D MDS PLOT AS html ##
  output$downloadMDS3Dhtml <- downloadHandler(
    filename = function() { paste0(getDataSet(),"scape3D.html") },
    content = function(file) {
      options(rgl.useNULL=FALSE)
      myplot <- getPlot3d()
      myplot
      rglwidget()
      rgl::writeWebGL(dir=getwd(), filename=file, snapshot=TRUE, width = 500, reuse=TRUE)
    },
    contentType = 'html'  
  )
  
  ## EXPORT SHEPARD PLOT AS PNG ##
  output$downloadShep <- downloadHandler(
    filename = function() { paste0(getDataSet(),"Shepard.png") },
    content = function(file) {
      shep <- getShep()
      labcol <- getLabcol()
      png(file=file, width = 10, height = 10, units = 'in', res = 500)
      plot(shep, pch=19, cex=0.5, col=labcol, xlab="Distance in tree space", ylab="Distance on MDS plot")
      dev.off()
    },
    contentType = 'image/png'
  )
  
  ## EXPORT SCREEPLOT AS PNG ##
  output$downloadScree <- downloadHandler(
    filename = function() { paste0(getDataSet(),"screeplot.png") },
    content = function(file) {
      res <- getAnalysis()
      labcol <- getLabcol()
      png(file=file, width = 5, height = 3, units = 'in', res = 500)
      barplot(res$pco$eig, col=labcol)
      dev.off()
    },
    contentType = 'image/png'
  )
  
  
  ## EXPORT TREE PLOT AS PNG ##
  output$downloadTree <- downloadHandler(
    filename = function() { paste0(getDataSet(),"SingleTree.png") },
    content = function(file) {
      tre <- getTree()
      png(file=file)
      plot(tre, type=input$treetype,
           show.tip.lab=input$showtiplabels, font=1, cex=input$tiplabelsize,
           direction=input$treedirection,
           edge.width=input$edgewidth)
      dev.off()
      contentType = 'image/png'
    }
  )
  
  ## EXPORT DENSITREE PLOT AS PNG ##
  output$downloadDensiTree <- downloadHandler(
    filename = function() { paste(getDataSet(), 'DensiTreeCluster',input$selectedDensiTree,'.png', sep='') },
    content = function(file) {
      clustTrees <- getDensiTree()
      png(file=file) 
      densiTree(clustTrees, col=4, alpha=input$alpha, scaleX=input$scaleX)
      dev.off()
      contentType = 'image/png'
    }
  )
  
  output$selectedGenTree <- renderUI({
    numTrees <- getLengthData()
    treeNames <- getTreeNames()
    choices <- c("",1:numTrees)
    names(choices) <- c("Choose one",treeNames)
    selectInput("selectedGenTree", "Choose individual tree", 
                choices=choices, selected="")
  })
  
  output$selectedRefTree <- renderUI({
    numTrees <- getLengthData()
    treeNames <- getTreeNames()
    choices <- c("",1:numTrees)
    names(choices) <- c("Choose one",treeNames)
    selectInput("selectedRefTree", "Select a reference tree", 
                choices=choices, selected="")
  })
  
 
  ## RENDER SYSTEM INFO ##
  output$systeminfo <- .render.server.info()
  
}) # end shinyServer