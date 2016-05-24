# load dataset soils
library(datacheck)

soil1 = read.csv(system.file("examples/soilsamples.csv", package="datacheck"))
recname = paste("Name",1:nrow(soil1),sep="_")
soil2 = cbind(recname,soil1)
soil2$ID[26:35]= - 123
soil2$Latitude[7:11] = 123


rule1=readLines(system.file("examples/soil_rules.R", package="datacheck"))
rule2=c("","sapply(recname, is.character)","!duplicated(recname)", rule1)

prof1 = datadict_profile(soil1, as_rules(rule1))
prof2 = datadict_profile(soil2, as_rules(rule2))

#############################################
shinyServer(function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "Soil_1" = soil1,
           "Soil_2" = soil2)
  })
  
  profileInput <- reactive({
    switch(input$dataset,
           "Soil_1" = prof1,
           "Soil_2" = prof2)
  })
  
  
  output$recLabels = renderUI({
    cn = names(datasetInput())
    un = rep(FALSE, length(cn))
    for(i in 1:length(cn)){
      un[i] = all(!duplicated(datasetInput()[,cn[i]]), !is.double(datasetInput()[,cn[i]]))
      #un[i] = all(!duplicated(soil2[,cn[i]]), !is.double(sol[,cn[i]]))
    }
    cn = cn[un]
    if(length(cn)>1) selectInput("labels", "Choose a variable for labeling records in heatmap:", cn)
    
  })
  
  output$recScores = renderUI({
    if(is_datadict_profile(profileInput())){
    scores = profileInput()$scores[1:(nrow(profileInput()$scores)-2),ncol(profileInput()$scores)]
    scores = sort(unique(scores))
    if(length(scores)>1){
    smin = min(scores)
    smax = max(scores)
    sval = max(smin, (smax-1))
    sliderInput("sscores","Restrict to low quality records till score of", min=smin, 
                max = smax, value = sval)
    }
    }
  })
  
  
  output$scoreSums = renderPlot({
    if(is_datadict_profile(profileInput())){
      profile = profileInput()
      score_sum(profile)
    }
  })
  
  # Show the first "n" observations
  output$view <- renderTable({
    
    sc = input$sscores
    dt = datasetInput()
    pi = profileInput()
    rs = dt
    if(!is.null(sc)){
      ds = pi$scores[1:nrow(dt),]
      ix = ds$Record.score <= sc
      rs = dt[ix,]
    }
    #head(rs, n = 300)
  })

  output$scores <- renderTable({
    sc = input$sscores
    dt = datasetInput()
    pi = profileInput()
    rs = as.data.frame("None", ncol=1, nrow=1)
    if(!is.null(sc)){
      ds = pi$scores[1:nrow(dt),]
      ix = ds$Record.score <= sc
      ds = ds[ix, ]
    }
    #head(ds, n = 300)
  }, digits = 0)
  
  
  output$heatmap <- renderPlot({
    if(is_datadict_profile(profileInput())){
      
      
      rm = 300
      title = paste("Heatmap of dataquality of up to first ",rm," records", sep="")
      
      heatmap_quality(profileInput(), input$labels, recMax = rm, scoreMax = input$sscores,
                      main = title
      )  
    }
    
  }, width = 700, height=1000)
  
  output$coverage <- renderPlot({
    if(is_datadict_profile(profileInput())){
      rule_coverage(profileInput())
    }
  })
  
  
  output$profile <- renderTable({
   as.data.frame(prep4rep(profileInput()$checks))  
  })
  
  output$descriptive <- renderTable({
    dt = datasetInput()
    if(names(dt)[1] != "V1"){
      short_summary(dt)
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { "scores.csv" },
    content = function(file) {
      write.csv(profileInput()$scores, file)
    }
  )
  
  
  
})