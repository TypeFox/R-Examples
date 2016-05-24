# library(EffectLiteR)
# library(lavaan)
# library(methods)
# library(shiny)
# library(foreign)
# library(ggplot2)
# library(nnet)
# library(lavaan.survey)

shinyServer(function(input, output, session) {
  
  ####### New analysis / reload button ########
  output$reload <- renderUI({
    if (input$newanalysis > 0) {
      tags$script("window.location.reload();")
    }
  })
  
  ######## Reactive Data Input ########
  dataInput <- reactive({
    inFile <- input$file1
    exdata <- input$exdata
    
    if(is.null(inFile)){      
      if(exdata==""){
        return(NULL)        
      }else if(exdata=="nonortho"){
        return(nonortho)  
      }else if(exdata=="example01"){
        return(example01)
      }else if(exdata=="example02lv"){
        return(example02lv)
      }else if(exdata=="example_multilevel"){
        return(example_multilevel)
      }                        
    }
    
    if(!is.null(inFile)){
      
      return(elrReadData(file=inFile$datapath,
                        name=inFile$name,
                        header=input$header,
                        sep=input$sep,
                        dec=input$dec,
                        na.strings=input$na.strings,
                        use.value.labels=input$vallabels))
      
    }
  })
  
  ###### Reactive Run Model #########
  model <- reactive({
        
    ## arguments for effectLite()
    d <- dataInput()
    mm <- mm()
    
    dv <- depv()
    x <- input$variablex
    
    k <- NULL; if(length(input$variablek) != 0){k <- input$variablek}
    fixed.cell <- FALSE; fixed.z <- FALSE
    if(input$fixed.cell == "fixed"){fixed.cell <- TRUE}
    if(input$fixed.cell == "fixed+e"){fixed.cell <- TRUE; fixed.z=TRUE}
    
    z <- NULL; if(length(input$variablez) != 0){z <- input$variablez}
    if(input$latentz & input$nlatentz > 0){z <- c(z,latentcov())}
    
    propscore <- NULL 
    if(length(input$propscore) != 0 & !input$propscoreformula){
      propscore <- input$propscore}
    if(input$prop.formula != "" & input$propscoreformula){
      propscore <- as.formula(input$prop.formula)}
    
    interactions <- input$interactions
    
    ids <- ~0
    if(input$ids != ""){ids <- as.formula(paste0(" ~ ", input$ids))}
    
    weights <- NULL
    if(input$weights != ""){weights <- as.formula(paste0(" ~ ", input$weights))}
    
    homoscedasticity <- input$homoscedasticity
    
    if(input$add.syntax == ""){
      add <- character()
    }else{
      add <- input$add.syntax
    }
    
    tryCatch(
      effectLite(y=dv, 
                 x=x,
                 k=k,
                 z=z,
                 data=d,
                 control=input$control,
                 measurement=mm,
                 missing=input$missing,
                 se=input$se,
                 bootstrap=input$bootstrap,
                 fixed.cell=fixed.cell,
                 fixed.z=fixed.z,
                 interactions=interactions,
                 propscore=propscore,                 
                 ids=ids,
                 weights=weights,
                 homoscedasticity=homoscedasticity,
                 add=add)
    )  
  })

  
  ######## Reactive zselect for Plot 2 ########
  zSelect <- reactive({
    zselect <- input$variablez
    
    if(!is.null(input$propscore)){
      d <- dataInput()
      x <- d[[input$variablex]]    
      ng <- length(unique(x))
      zselect <- c(zselect, paste0("logprop",1:(ng-1)))
    }
        
    return(zselect)
  })

  ######## Reactive gxselect for Plot 3 ########
  gxSelect <- reactive({
    
    d <- dataInput()
    x <- d[[input$variablex]]    
    ng <- length(unique(x))
    res <- paste0("g",2:ng-1)
    
    return(res)
  })
  
  
  ######## Reactive zselect2 for Plot 3 ########
  zSelect2 <- reactive({
    
    kstar <- NULL
    if(!is.null(input$variablek)){kstar <- "K"}
    zselect <- c(input$variablez, input$variablek, kstar, input$variablex)
    
    if(input$latentz & input$nlatentz > 0){
      nameslatentz <- c(input$name.etaz1, input$name.etaz2, input$name.etaz3,
                        input$name.etaz4, input$name.etaz5)
      nameslatentz <- nameslatentz[1:input$nlatentz]
      zselect <- c(nameslatentz, zselect)
    }
    
    d <- dataInput()
    x <- d[[input$variablex]]    
    ng <- length(unique(x))
    
    if(!is.null(input$propscore)){
      zselect <- c(zselect, input$propscore)
      zselect <- c(zselect, paste0("logprop",1:(ng-1)))
    }
    if(input$prop.formula != "" & input$propscoreformula){
      try({
        propscore <- as.formula(input$prop.formula)
        zselect <- c(zselect, all.vars(propscore[[3]]))
      }, silent=TRUE)
      zselect <- c(zselect, paste0("logprop",1:(ng-1)))
    }
    
    zselect <- c(zselect, paste0("g",1:(ng-1)))
    
    return(zselect)
  })

  
  ######## Reactive zselect3 (colour) for Plot 3 ########
  zSelect3 <- reactive({
    kstar <- NULL
    if(!is.null(input$variablek)){kstar <- "K"}
    zselect3 <- c("", input$variablek, kstar, input$variablex, input$variablez)
    
    d <- dataInput()
    x <- d[[input$variablex]]    
    ng <- length(unique(x))
      
    if(!is.null(input$propscore)){
      zselect3 <- c(zselect3, input$propscore)
      zselect3 <- c(zselect3, paste0("logprop",1:(ng-1)))
    }
    if(input$prop.formula != "" & input$propscoreformula){
      try({
        propscore <- as.formula(input$prop.formula)
        zselect3 <- c(zselect3, all.vars(propscore[[3]]))
      }, silent=TRUE)
      zselect3 <- c(zselect3, paste0("logprop",1:(ng-1)))
    }
    
    zselect3 <- c(zselect3, paste0("g",1:(ng-1)))
    
    return(zselect3)
  })
  
  
  ######## Reactive measurement model ########
  mm <- reactive({
    if(!input$latenty & !input$latentz){
      return(character())
      
    }else if(input$latenty | input$nlatentz > 0){
      
      ## determine number of cells
      d <- dataInput()
      ng <- length(unique(d[[input$variablex]]))
      nk <- 1
      if(length(input$variablek) != 0){
        for(i in 1:length(input$variablek)){
          tmpvar <- as.factor(d[[input$variablek[i]]])
          nk <- nk*length(levels(tmpvar))
        }        
      }
      ncells <- ng*nk
      
      names <- NULL
      indicators <- NULL
      mmodel <- NULL
      
      if(input$latenty){
        names$etay <- input$name.etay
        indicators$indicatorsy <- input$indicatorsy
        mmodel$mm.etay <- input$mm.etay
      }
      if(input$latentz & input$nlatentz > 0){
        names$etaz1 <- input$name.etaz1
        indicators$indicatorsz1 <- input$indicatorsz1
        mmodel$mm.etaz1 <- input$mm.etaz1
      }
      if(input$latentz & input$nlatentz > 1){
        names$etaz2 <- input$name.etaz2
        indicators$indicatorsz2 <- input$indicatorsz2
        mmodel$mm.etaz2 <- input$mm.etaz2
      }
      if(input$latentz & input$nlatentz > 2){
        names$etaz3 <- input$name.etaz3
        indicators$indicatorsz3 <- input$indicatorsz3
        mmodel$mm.etaz3 <- input$mm.etaz3
      }
      if(input$latentz & input$nlatentz > 3){
        names$etaz4 <- input$name.etaz4
        indicators$indicatorsz4 <- input$indicatorsz4
        mmodel$mm.etaz4 <- input$mm.etaz4
      }
      if(input$latentz & input$nlatentz > 4){
        names$etaz5 <- input$name.etaz5
        indicators$indicatorsz5 <- input$indicatorsz5
        mmodel$mm.etaz5 <- input$mm.etaz5
      }
      
      names <- unlist(names)
      mmodel <- unlist(mmodel)
      
      ## switch to default options if custom mmodel is not specified
      if(!input$custommeasmodels){mmodel <- NULL}
      
      mm <- generateMeasurementModel(
        names=names,
        indicators=indicators,
        ncells=ncells,
        model=mmodel
      )
      
      return(mm)
    }
  })
  
  ###### Reactive dependent variable ###########
  depv <- reactive({
    
    if(input$latenty){
      return(input$name.etay)
      
    }else{
      return(input$variabley)
      
    }
  })  

  ###### Reactive latent covariates ###########
  latentcov <- reactive({
    
    if(input$latentz == TRUE & input$nlatentz > 0){
        
      nameslatentcov <- NULL; nameslatentcov$etaz1 <- input$name.etaz1
      if(input$nlatentz > 1){nameslatentcov$etaz2 <- input$name.etaz2}
      if(input$nlatentz > 2){nameslatentcov$etaz3 <- input$name.etaz3}
      if(input$nlatentz > 3){nameslatentcov$etaz4 <- input$name.etaz4}
      if(input$nlatentz > 4){nameslatentcov$etaz5 <- input$name.etaz5}
      
      nameslatentcov <- unlist(nameslatentcov)
      return(nameslatentcov)
      
    }else{
      return(NULL)
      
    }
  })  
  
  
  ###### Update Variable Selectors UI ########
  observe({
    inFile <- input$file1
    exdata <- input$exdata
    
    if(is.null(inFile) & exdata=="")
      return(NULL)  
    
    d <- dataInput()
    
    updateSelectInput(session, "variabley", 
                      choices = c("", names(d)))
    updateSelectInput(session, "variablex", 
                      choices = c("", names(d)))
    updateSelectInput(session, "variablek", 
                      choices = c("", names(d)),
                      selected = "")
    updateSelectInput(session, "variablez", 
                      choices = c("", names(d)),
                      selected = "")
    updateSelectInput(session, "propscore", 
                      choices = c("", names(d)),
                      selected = "")
    updateSelectInput(session, "ids", 
                      choices = c("", names(d)),
                      selected = "")
    updateSelectInput(session, "weights", 
                      choices = c("", names(d)),
                      selected = "")
    updateSelectInput(session, "indicatorsy", choices = names(d))
    updateSelectInput(session, "indicatorsz1", choices = names(d))
    updateSelectInput(session, "indicatorsz2", choices = names(d))
    updateSelectInput(session, "indicatorsz3", choices = names(d))
    updateSelectInput(session, "indicatorsz4", choices = names(d))
    updateSelectInput(session, "indicatorsz5", choices = names(d))
  })


  ###### Update zselect for Plot 2 ########
  observe({
    zsel <- zSelect()
    updateSelectInput(session, "zselect", 
                    choices = zsel)  
  })

  ###### Update gxselect for Plot 3 ########
  observe({
    gxsel <- gxSelect()
    updateSelectInput(session, "gxselect", 
                      choices = gxsel)  
  })
  
  
  ###### Update zselect2 for Plot 3 ########
  observe({
    zsel <- zSelect2()
    updateSelectInput(session, "zselect2", 
                      choices = zsel)  
  })
  
  ###### Update zselect3 for Plot 3 ########
  observe({
    zsel3 <- zSelect3()
    updateSelectInput(session, "zselect3", 
                      choices = zsel3)  
  })
  
  ###### Update Control Group UI ########
  observe({
    inputx <- input$variablex
    
    if(inputx==""){
      return(NULL)        
    }else{      
      d <- dataInput()
      x <- as.factor(d[,inputx])
      
      updateSelectInput(session, "control", choices = levels(x))      
    }
  })  
  
  
  ###### Output Data Table #########  
  output$mytable1 = renderDataTable({ 
    d <- dataInput()
    dprint <- format(d, digits=3)
    dprint
  })

  ###### Output Conditional Effects Table #########
  output$helptextcondeffects <- renderPrint({
    if((input$variabley == "" & !input$latenty) || input$variablex == "" ){
           
      cat("Conditional effects are only shown if you have specified the dependent variable and the treatment variable.")
      
    }else{
      
      cat("This datatable shows the values and standard errors of the effect function for given values of the categorical and continuous covariates. Regression factor scores are used for latent covariates.")
    }
  })
    
  output$condeffs = renderDataTable({
    if((input$variabley == "" & !input$latenty) || input$variablex == "" ){
      return(NULL)
    }else{            
      m1 <- model()
      condprint <- format(m1@results@condeffects, digits=3)
      condprint
    }  
  })
  
  
  ###### Output Plot 1 #########
  output$helptextplot1 <- renderPrint({
    if(input$variabley == "" || input$variablex == "" || input$latenty == TRUE){
      
      cat("Plot 1 only works for a manifest dependent variable and a treatment variable.")
      
    }else{
      
      cat("Plot 1 shows a histogram of the dependent variable in each cell.")
    }
  })
    
  output$plot1 <- renderPlot({    
    
    if(input$variabley == "" || input$variablex == "" || input$latenty == TRUE){
      return(NULL)
    }else{
      
      m1 <- model()

      y <- m1@input@data[[input$variabley]]      
      cell <- m1@input@data[["cell"]]
      dp <- na.omit(data.frame(y,cell))
      binwidth <- (range(y, na.rm=TRUE)[2]-range(y, na.rm=TRUE)[1])/30

      p <- ggplot2::qplot(y, data=dp, geom="histogram",
                 binwidth=binwidth,
                 xlab=input$variabley,
                 main=paste0("Distribution of ", input$variabley, " in cells"))
      p <- p + ggplot2::facet_wrap( ~ cell)
      p <- p + ggplot2::theme_bw()
      print(p)
    }  
        
  })

  ###### Output Plot 2 #########  
  output$helptextplot2 <- renderPrint({
    if(input$variabley == "" || input$variablex == "" || 
         (is.null(input$variablez) & is.null(input$propscore)) || 
         input$latenty){
      
      cat("Plot 2 only works for a manifest dependent variable, a treatment variable, and at least one continuous covariate.")
      
    }else{
      
      cat("Plot 2 shows the regression of the dependent variable on the selected continuous covariate in each cell.")
    }
  })
  
  
  output$plot2 <- renderPlot({    
    
    if(input$variabley == "" || input$variablex == "" || 
         (is.null(input$variablez) & is.null(input$propscore)) || 
         input$latenty){
      
      return(NULL)
    }else{
      
      m1 <- model()
      
      y <- m1@input@data[[input$variabley]]
      zselected <- m1@input@data[[input$zselect]]
      cell <- m1@input@data[["cell"]]
      
      dp <- data.frame(y,cell,zselected)
      
      p <- ggplot2::qplot(y=y, x=zselected, data=dp, 
                 ylab=input$variabley,
                 xlab=input$zselect,                 
                 main=paste0("Regression of ", input$variabley, " on ", 
                             input$zselect, " in cells"))
      p <- p + ggplot2::facet_wrap( ~ cell)
      p <- p + ggplot2::geom_smooth(method = "lm")
      p <- p + ggplot2::theme_bw()
      print(p)                  
    }
        
  })

  
  ###### Output Plot 3 #########
  output$helptextplot3 <- renderPrint({
    if(input$variabley == "" || input$variablex == "" || 
         (is.null(input$variablez) & is.null(input$variablek) & !input$latentz & 
            is.null(input$propscore))  || 
         input$latenty || input$latentz){
      
      cat("Plot 3 is only shown if a dependent variable, a treatment variable, and at least one covariate is specified.")
      
    }else{
      
      cat("Plot 3 shows the regression of the selected effect function on the selected regressor.")
    }
  })
  
  output$plot3 <- renderPlot({    
    
    if((input$variabley == "" & !input$latenty) || input$variablex == "" || 
         (is.null(input$variablez) & is.null(input$variablek) & !input$latentz & 
            is.null(input$propscore))){
      return(NULL)
    }else{
      
      m1 <- model()
      condeffects <- m1@results@condeffects
      yselected <- round(condeffects[[input$gxselect]],4)    
      zselected <- condeffects[[input$zselect2]]
      colourselected <- condeffects[[input$zselect3]]
              
      g1label <- "(K,Z)"
      if(length(input$variablek) == 0){g1label <- "(Z)"}
      
      p <- ggplot2::qplot(y=yselected, x=zselected, 
                 data=condeffects,
                 ylab=paste0(input$gxselect,g1label),
                 xlab=input$zselect2,                 
                 main=paste0("Estimated regression of ",
                             paste0(input$gxselect,g1label), " on ", 
                             input$zselect2))
      p <- p + ggplot2::geom_smooth(method="loess")
      if(is.null(colourselected)){
        p <- p + ggplot2::geom_point(size=2.5)
      }else{
        p <- p + ggplot2::geom_point(ggplot2::aes(colour=colourselected),size=2.5)
        p <- p + ggplot2::guides(colour = ggplot2::guide_legend(input$zselect3))         }    
      p <- p + ggplot2::theme_bw()
      
      print(p)
      
    }
    
  })

  
  ###### Output EffectLiteR Summary #########
  output$summary <- renderPrint({
    
    if(input$variabley == "" & input$latenty == FALSE || input$variablex == ""){      
      
      cat("Please specify the outcome variable and the treatment variable")
      
    }else{
      
      m1 <- model()
      m1
    }
  })
  
  ###### Output ELR call #########
  output$elrcall <- renderPrint({
    if(input$variabley == "" & input$latenty == FALSE || input$variablex == ""){    
      cat("Please specify the outcome variable and the treatment variable")
    }else{
      dv <- depv()
      x <- input$variablex
      
      mm <- mm()
      printmm <- "character()"
      if(length(mm) != 0){printmm <- "mm"}
      
      printk <- "NULL"
      if(length(input$variablek) != 0){
        printk <- paste0("c(\"",
                         paste(input$variablek, collapse="\",\""), 
                         "\")")
      }
      
      fixed.cell <- FALSE; fixed.z <- FALSE
      if(input$fixed.cell == "fixed"){fixed.cell <- TRUE}
      if(input$fixed.cell == "fixed+e"){fixed.cell <- TRUE; fixed.z=TRUE}
      
      z <- NULL
      printz <- "NULL"
      if(length(input$variablez) != 0){z <- input$variablez}
      if(input$latentz & input$nlatentz > 0){z <- c(z,latentcov())}
      if(!is.null(z)){
        printz <- paste0("c(\"",
                         paste(z, collapse="\",\""), 
                         "\")")
      }
      
      propscore <- NULL 
      printpropscore <- "NULL"
      if(length(input$propscore) != 0 & !input$propscoreformula){
        propscore <- input$propscore
        printpropscore <- paste0("c(\"",
                                 paste(propscore, collapse="\",\""), 
                                 "\")")
      }
      if(input$prop.formula != "" & input$propscoreformula){
        printpropscore <- input$prop.formula
      }
      
      interactions <- input$interactions
      
      printids <- "~0"
      if(input$ids != ""){
        printids <- paste0(" ~ ", input$ids)
      }
      
      printweights="NULL"
      if(input$weights != ""){
        printweights <- paste0(" ~ ", input$weights)
        }
      
      homoscedasticity <- input$homoscedasticity
      
      printadd <- "character()"
      if(input$add.syntax != ""){printadd <- "add"}
      
      printbootstrap <- ""
      if(input$se == "boot"){
        printbootstrap <- paste0("bootstrap=",input$bootstrap, ", ")}
      
      tmp <- paste0("#### Call for effectLite #### \n\n",
                    "effectLite(",
                    "y=\"", dv, "\", ",
                    "x=\"", x, "\", ",
                    "k=", printk, ", ",
                    "z=", printz, ", ",
                    "data=data, ",
                    "control=\"", input$control, "\", ",
                    "measurement=", printmm, ", ",
                    "missing=\"", input$missing, "\", ",
                    "se=\"", input$se, "\", ",
                    printbootstrap,
                    "fixed.cell=", fixed.cell, ", ",
                    "fixed.z=", fixed.z, ", ",
                    "interactions=\"", interactions, "\", ",
                    "propscore=", printpropscore, ", ",
                    "ids=", printids, ", ",
                    "weights=", printweights, ", ",
                    "homoscedasticity=", homoscedasticity, ", ",
                    "add=", printadd,
                    ")")
      
      cat(tmp)  
    }  
  })
  

  ###### Output lavaan call #########
  output$lavcall <- renderPrint({
    if(input$variabley == "" & input$latenty == FALSE || input$variablex == ""){    
      cat("")
    }else{
      m1 <- model()
      cellabel <- paste0("c(\"",
                         paste(m1@input@vlevels$cell, collapse="\",\""), 
                         "\")")
      printbootstrap <- ""
      if(input$se == "boot"){
        printbootstrap <- paste0("bootstrap=",m1@input@bootstrap, ", ")}
      
      tmp <- paste0("#### Call for lavaan::sem #### \n\n",
                  "sem(", "model=model, ",
                  "group=\"", "cell", "\", ", 
                  "missing=\"", m1@input@missing, "\", ",
                  "se=\"", m1@input@se, "\", ", 
                  printbootstrap,
                  "group.label=", cellabel, ", ",
                  "data=data, ", 
                  "fixed.x=", m1@input@fixed.z, ", ",
                  "group.w.free=", !m1@input@fixed.cell, ", ",
                  "mimic=\"", "mplus", "\") ")
      cat(tmp, "\n")  
    }  
  })
  
  
      
  ###### Output Lavaan Syntax #########
  output$lavsyntax <- renderPrint({
    if(input$variabley == "" & input$latenty == FALSE || input$variablex == ""){    
      cat("")
    }else{
      m1 <- model()
      cat(m1@lavaansyntax@model)  
    }  
  })
  

  ###### Output Lavaan Results #########
  output$lavresults <- renderPrint({      
    
    if(input$variabley == "" & input$latenty == FALSE || input$variablex == ""){            
      
      cat("Please specify the outcome variable and the treatment variable")
      
    }else{
          
      m1 <- model()
      summary(m1@results@lavresults, fit.measures=TRUE)  
      ## maybe there was a reason I set fit.measures=FALSE in prior versions...
    }  
  })
  
  
  ###### Download Data (Conditional Effects Table) #######
  output$downloadConditionalEffects <- downloadHandler(
    filename = function() {
      paste('ELR-ConditionalEffects-', Sys.Date(), '.txt', sep='')
    },
    content = function(con) {
      m1 <- model()
      write.table(m1@results@condeffects, con, row.names=F, col.names=T, 
                  quote=F)
    }
  )
  
  ###### Download (transformed) Input Data #######
  output$downloadLavData <- downloadHandler(
    filename = function() {
      paste('ELR-data-', Sys.Date(), '.txt', sep='')
    },
    content = function(con) {
      m1 <- model()
      write.table(m1@input@data, con, row.names=F, col.names=T, 
                  quote=F)
    }
  )
  
  ##### Conditional effects II
  output$ui <- renderUI({
    
    m1 <- model()
    condeffects <- m1@results@condeffects
    vnamesk <- m1@input@vnames$k
    vlevelsk <- m1@input@vlevels$levels.k.original
    vnamesz <- m1@input@vnames$z
    numberks <- length(vnamesk)
    numberzs <- length(vnamesz)
    covnames <- c(vnamesk,vnamesz)
    uilist <- vector("list", length(covnames))
    
    if(numberks==0 & numberzs==0){return(NULL)}
    
    if(numberks>0){
    for(i in 1:numberks){
      uilist[[i]] <- selectInput(inputId = paste0('valk',i), 
                                 label = vnamesk[i], 
                                 choices = vlevelsk[i],
                                 width='90%')
    }}
    
    if(numberzs>0){
    for(i in 1:numberzs){
      uilist[[numberks+i]] <- numericInput(inputId = paste0('valz',i),
                label = vnamesz[i],
                value = round(mean(condeffects[[vnamesz[i]]], na.rm=T),3),
                width='90%')
    }}
    
    uilist
    
  })
  
  
  ## help text conditional effects 2
  output$helptextcondeffects2 <- renderPrint({
    if((input$variabley == "" & !input$latenty) || input$variablex == "" ){
      
      cat("Conditional effects are only available if you have specified the dependent variable and the treatment variable.")
      
    }else{
      
      cat("Conditional effects for user specified values of categorical and continuous covariates.")
    }
  })
  
  
  ## output conditional effects 2
  output$outputcondeffect2 <- renderPrint({
    
    m1 <- model()
    vnamesk <- m1@input@vnames$k
    vnamesz <- m1@input@vnames$z
    numberks <- length(vnamesk)
    numberzs <- length(vnamesz)
    covnames <- c(vnamesk,vnamesz)
    
    valuesk <- list()
    valuesz <- list()

    if(numberks>0){
      for(i in 1:numberks){
        valuesk <- c(valuesk, input[[paste0('valk',i)]])
      }}
    
    if(numberzs>0){
      for(i in 1:numberzs){
        valuesz <- c(valuesz, input[[paste0('valz',i)]])
    }}
    
    
    newdata <- data.frame(c(valuesk, valuesz))
    try({names(newdata) <- covnames}, silent=TRUE)
    
    try({indeff <- round(elrPredict(m1, newdata),3)}, silent=TRUE)
    try({print(indeff, row.names=F, print.gap=3)}, silent=TRUE)
    
  })
  
  
    
  
})