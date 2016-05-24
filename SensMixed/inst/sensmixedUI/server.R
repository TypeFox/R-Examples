# server.R
shinyServer(function(input, output, session) {
  tags$style(type="text/css", ".tab-content { overflow: visible; }")
  tags$head(
    tags$style(type="text/css", "html {overflow:hidden;}"))
  uploadData <- reactive({
    if(is.null(input$uploaddata) || input$uploaddata == 1){
      inFile <- input$file1
      
      if (is.null(inFile))
        {return()}
      
      
      return(read.delim(inFile$datapath, header=input$header, sep=input$sep, 
                        quote=input$quote,
                        dec = input$decimal, fileEncoding="UTF-8-BOM"))
    }
    else if(input$uploaddata == 2){
      TVbo$product <- interaction(TVbo$Picture, TVbo$TVset)
      return(TVbo)
    }
    else if(input$uploaddata == 3)
      return(ham)   
  })
  
  Data <- reactive({    
     input$goButton
  isolate({
     if (is.null(uploadData()))
       {return()}
    
    df.raw <- uploadData()     
    
    ## here the analysis of consumer/sensory data is sourced
    ## and saved in res variable
    #source(paste(system.file("sensmixedUI", package = "SensMixed"), "runAnalysis.R",  sep = "/"))
    source('runAnalysis.R', local=TRUE)
    
    return(res)
  })
  })
  
  ##### call utils functions ###################################################
  sensmixedPlot <- function(){
    if(input$analysis == "Consumer data") { return() }
    
    if(is.null(Data()))
      return("")
    
    if(class(Data()) == "consmixed") { return() }
    
    if(input$typeEffs == 1){
        return(plot(Data(), mult = input$representPlot, isFixed = FALSE,
                    isScaling = FALSE, cex = 2))
    }
    else if(input$typeEffs == 2){
      return(plot(Data(), mult = input$representPlot, isRand = FALSE, 
                  isScaling = FALSE, 
                  dprime = input$typePlot, cex = 2))
    }
    else{
      return(plot(Data(), mult = input$representPlot, isRand = FALSE, 
                  isFixed = FALSE, cex = 2))
    }
  } 
  
  ## here the step results are formatted using xtable
  #source(paste(system.file("sensmixedUI", package = "SensMixed"), "stepUtils.R",  sep = "/"))
  source('stepUtils.R', local=TRUE)
 # source(paste(system.file("sensmixedUI", package = "SensMixed"), "posthocUtils.R",  sep = "/"))
  source('posthocUtils.R', local=TRUE)
  source('MAMUtils.R', local=TRUE)
  ##############################################################################

  
  output$plotsSensMixed <- renderPlot({   
    sensmixedPlot()  
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste("plotSensmixed",input$typeEffs, 
                                  '.png', sep='') },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, 
                                                            height = height, 
                                                            res = 300, 
                                                            units = "in")
      ggsave(file, sensmixedPlot(), scale = input$scalePlot, device = device)
    }
  )

  output$downloadTable <- downloadHandler(
    filename = function() { paste("tableSensmixed", input$typeEffsTable, 
                                  '.doc', sep='') },
    content = function(file) {
      sink(file)
      saveToDoc(Data(), type = input$typetable2, typeEffs = input$typeEffsTable)
      sink()
    }, contentType = 'text/plain'
  )
  
  output$downloadStep <- downloadHandler(
    filename = function() { paste(getNameStep(), '.doc', sep='') },
    content = function(file) {
      sink(file)
      stepRandResult() 
      stepFixedResult()
      sink()
    }, contentType = 'text/plain'
  )
  
  output$downloadPosthocTable <- downloadHandler(
    filename = function() { paste(input$AttrPosthoc, input$whichPlot, 
                                  input$effsPlot, '.doc', sep='') },
    content = function(file) {
      sink(file)
      posthocResult()
      sink()
    }, contentType = 'text/plain'
  )
  
  output$downloadPosthocPlot <- downloadHandler(
    filename = function() { paste(input$AttrPosthoc, input$whichPlot, 
                                  input$effsPlot, '.png', sep='') },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, 
                                                            height = height, 
                                                            res = 300, 
                                                            units = "in")
      ggsave(file, posthocPlot(), device = device)
    }
  )
  
  output$downloadMAM <- downloadHandler(
    filename = function() { paste(input$AttrMAManalysis, " MAM results", 
                                  '.doc', sep='') },
    content = function(file) {
      sink(file)
      mamanova() 
      mamindiv()
      mamperf()
      mamposthoc()
      mamdiffmean()
      sink()
    }, contentType = 'text/plain'
  )

  
  output$tablesSensMixed <- renderPrint({
    if(is.null(uploadData())) { return() }
    if(input$analysis == "Consumer data") { return() }
    if(is.null(Data())){return()}
    if(class(Data()) == "consmixed") { return() }
    saveToDoc(Data(), type = input$typetable2, typeEffs = input$typeEffsTable)    
  })
  
  output$MAMtable <- renderPrint({
    mamanova()
  })
  
  output$MAMindiv <- renderPrint({
    mamindiv()
  })
  
  output$MAMperf <- renderPrint({
    mamperf()
  })
  
  
  output$MAMplotposthoc <- renderPlot({
    if(is.null(uploadData())) { return() }
    if(input$analysis == "Consumer data") { return() }
    if(is.null(Data()) || is.null(Data()$MAMan)){return()}
    if(class(Data()) == "consmixed") { return() }
    
    resposthoc <- Data()$MAMan[[5]][, , input$AttrMAManalysis]
    resci <- Data()$MAMan[[8]][, , input$AttrMAManalysis]
    
    tab <- cbind(resposthoc, resci)
    tab <- as.data.frame(tab)
    tab$levels <- rownames(tab)
    colnames(tab)[which(colnames(tab)=="Lower CI")] <- "lci"
    colnames(tab)[which(colnames(tab)=="Upper CI")] <- "uci"
    tab$col.bars <-  unlist(lapply(tab[,"Pval"], calc.cols2))
    
    ggplot(tab, aes(x=levels, y = Estimate, fill = col.bars)) + 
      geom_bar(position = "dodge", stat = "identity") +  
      geom_errorbar(aes(ymin = lci, ymax = uci ), 
                    colour="black", width=.1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
          axis.title.y = element_text(size = rel(1.4)), 
          axis.text = element_text(size = rel(1)), 
          legend.text = element_text(size = rel(1)), 
          legend.title = element_text(size = rel(1))) + 
      scale_fill_manual(values  = 
                          c(  "NS" = "grey", "p-value < 0.01" = "orange", 
                              "p-value < 0.05" = "yellow", 
                              "p-value < 0.001" = "red"), name="Significance")  
    
    })
  
  
  output$MAMposthoc <- renderPrint({
    mamposthoc()
  })
  
  output$MAMposthoc2 <- renderPrint({
    if(is.null(uploadData())) { return() }
    if(input$analysis == "Consumer data") { return() }
    if(is.null(Data()) || is.null(Data()$MAMan)){return()}
    if(class(Data()) == "consmixed") { return() }
    
    resposthoc <- Data()$MAMan[[5]][, , input$AttrMAManalysis]
    resci <- Data()$MAMan[[8]][, , input$AttrMAManalysis]
    resposthoc <- cbind(resposthoc, resci)
    resposthoc[ , "Pval"] <- 
      format.pval(resposthoc[, "Pval"], digits=3, eps=1e-3)
    resposthoc <- xtable(resposthoc, align = 
                           paste(c("l", rep("l", ncol(resposthoc))), 
                                 collapse = ""), 
                         display= c("s",rep("f", ncol(resposthoc))))
    caption(resposthoc) <- paste(" Pairwise product differences")
    print(resposthoc, caption.placement="top", table.placement="H", 
          html.table.attributes = getOption("xtable.html.table.attributes",
                                            "rules='groups' width='100%'"),
          type = "html")
  })
  
  output$MAMdiffmean <- renderPrint({
    mamdiffmean()
  })
  
  
  output$stepRand <- renderPrint({
    stepRandResult()  
  })
  
  output$stepFixed <- renderPrint({
    stepFixedResult()
  })
  
  output$posthocTable <- renderPrint({
    posthocResult()
  })

  output$posthocPlot <- renderPlot({
    posthocPlot()
  })
   
  output$contents <- renderDataTable({
    if(!is.null(uploadData()))
      return(uploadData())
    
  })

  output$helpprodstruct <- renderTable({
    helpprodstruct <- matrix(NA, nrow = 3, ncol = 1)
    rownames(helpprodstruct) <- c(1,2,3)
    colnames(helpprodstruct) <- "Explanations"
    helpprodstruct[1,1] <- "only main effects will enter the initial model"
    helpprodstruct[2,1] <- "main effects and 2-way interaction"
    helpprodstruct[3,1] <- "all main effects and all possible interaction"
    return(xtable(helpprodstruct))
  })

  output$helperrstruct <- renderTable({
    helperrstruct  <- matrix(NA, nrow = 3, ncol = 1)
    rownames(helperrstruct) <- c("No-Rep","2-WAY","3-WAY")
    colnames(helperrstruct) <- "Explanations"
    helperrstruct[1,1] <- "assessor effect and all possible interactions between assessor and product effects"
    helperrstruct[2,1] <- "No-Rep + replicate effect and replicate assessor interaction effect"
    helperrstruct[3,1] <- "assessor and replicate effect and interaction between them and interaction between them and Product_effects"
    return(xtable(helperrstruct))
  })
  
  output$helponeway <- renderTable({
    helponeway  <- matrix(NA, nrow = 2, ncol = 1)
    rownames(helponeway) <- c("No", "Yes")
    colnames(helponeway) <- "Explanations"
    helponeway[1,1] <- "considers multi-way product structure in the random part"
    helponeway[2,1] <- "considers just one product factor in the random part, where the product factor is chosen as the overall product factor combining each product-combination into a single factor with as many levels as there are different product combinations"    
    return(xtable(helponeway))
  })
   

  ## here the server part of the UI is sourced
  #source(paste(system.file("sensmixedUI", package = "SensMixed"), "serverUI.R",  sep = "/"))
  source('serverUI.R', local = TRUE)

  
  addTooltip(session, "plotsSensMixed", "title", placement = "bottom", 
             trigger = "click") 
  session$onSessionEnded(function() { stopApp() })

})