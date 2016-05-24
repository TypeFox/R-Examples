#### The file contains the server part of the UI (interactive UI)




## constructs tab panel for the input controls
tabPanel.input <- function(names.dd){
  if(input$analysis == "Consumer data")
    return(tabPanel("Input arguments",
                    selectInput("Response", "Select response", names.dd,
                                selected = "Liking"),           
                    selectInput("Consumers", "Select consumer", names.dd,
                                selected = "Consumer"),
                    selectizeInput("Products", "Select products", names.dd,  
                                   options = list(dropdownParent = 'body'),
                                   multiple = TRUE, selected = c("Product",
                                                                 "Information")),
                    selectizeInput("Consumerchar", 
                                   "Select consumer characteristics", 
                                   names.dd,  
                                   options = list(dropdownParent = 'body'),
                                   multiple = TRUE, selected = c("Gender",
                                                                 "Age")),
                    selectizeInput("Consumerfact", 
                                   "Consumer characteristics treated as factors", 
                                   names.dd,  
                                   options = list(dropdownParent = 'body'),
                                   multiple = TRUE, selected = "Gender")
    ))
  else {
    if(input$uploaddata == 2)
      return(tabPanel("Input arguments",
                    selectizeInput("Attributes", "Select attributes", 
                                   names.dd,
                                   options = list(dropdownParent = 'body'),
                                   multiple = TRUE, selected = 
                                     names.dd[5:(length(names.dd)-1)]),          
                    selectInput("Assessors", "Select assessor", 
                                names.dd, selected = 
                                  ifelse(input$uploaddata == 2, "Assessor",
                                         "")),
                    selectInput("Replications", "Select replications", 
                                names.dd, selected = 
                                  ifelse(input$uploaddata == 2, "Repeat",
                                         "")),
                    selectizeInput("Products", "Select products", names.dd,  
                                   options = list(dropdownParent = 'body'),
                                   multiple = TRUE, selected = c("product"))
      ))
    else
      return(tabPanel("Input arguments",
                      selectizeInput("Attributes", "Select attributes", 
                                     names.dd,
                                     options = list(dropdownParent = 'body'),
                                     multiple = TRUE),          
                      selectInput("Assessors", "Select assessor", 
                                  names.dd),
                      selectInput("Replications", "Select replications", 
                                  names.dd),
                      selectizeInput("Products", "Select products", names.dd,  
                                     options = list(dropdownParent = 'body'),
                                     multiple = TRUE))
      )
  }
}

## constructs tab panel for the modelling controls
tabPanel.model <- function(){
  if(input$analysis == "Consumer data")
    return(tabPanel("Modelling controls",
                    selectInput('struct', 'Select structure', 
                                c("1" = 1, "2" = 2, "3" = 3))             
    ))
  else
    return(tabPanel("Modelling controls",                
                     selectInput('struct', 
                                 'Select product structure', 
                                 c("1" = 1, "2" = 2, "3" = 3)),
                     bsPopover("struct", paste0("<p><b>1</b>: main effects</p>", 
                                                "<p><b>2</b>: main effects and ",
                                                "two way interactions</p>",
                                                "<p><b>3</b>: main effects and ",
                                                "all possible interactions</p>"), 
                               placement = "right", trigger = "hover"),
#                      bsCollapsePanel("Help product structure", 
#                                      tableOutput("helpprodstruct"), id="colll1", 
#                                      value="testlll1", style = "info"),
                     selectInput('errstruct', 
                                 'Select error structure', 
                                 c("No_Rep" = "No_Rep", 
                                   "2-WAY" = "2-WAY", 
                                   "3-WAY" = "3-WAY")),

                     bsPopover("errstruct", paste0("<p><b>No-Rep</b>: assessor " ,
                                                    "effect and all possible ",
                                                   " interactions between ",
                                                   "assessor and product ",
                                                   "effects</p>", 
                                                   "<p> <b>2-WAY</b>: <b>No-Rep</b> and",
                                                   " replicate effect and ",
                                                   "replicate assessor ",
                                                   "interaction effect</p>",
                                                   "<p> <b>3-WAY</b>: assessor ",
                                                   "and replicate effect and ",
                                                   "interaction between them ",
                                                   "and interaction between ",
                                                   "them and Product_effects</p>"), 
                               placement = "right", trigger = "hover"),

#                      bsCollapsePanel("Help error structure", 
#                                      tableOutput("helperrstruct"), 
#                                      id="col2", value="test2", style = "info"),
                     
#                     selectInput('oneway_rand', 'One-way product random part', 
#                                 c( "No" = FALSE, "Yes" = TRUE)),
#                     bsCollapsePanel("Help one-way product random part", 
#                                     tableOutput("helponeway"), 
#                                     id="col3", value="test3", style = "info"),
                     selectInput('MAM', 'Correct for scaling', c("Yes" = TRUE, 
                                                                 "No" = FALSE),
                                 selected = FALSE),
#                      bsPopover("MAM", paste0("<p> Mixed Assessor Model(MAM)",
#                                              " is applied as suggested by </p>", 
#                                              "<p>Brockhoff, Per Bruun, ",
#                                              "Pascal Schlich, and Ib Skovgaard.",
#                                              "Taking Individual Scaling ",
#                                              "Differences into Account by ",
#                                              "Analyzing Profile Data with the ",
#                                              " Mixed Assessor Model.” FOOD ",
#                                              " QUALITY AND PREFERENCE 39", 
#                                              "(2015): 156–166. Web.</p>"), 
#                                 placement = "right", trigger = "hover")
                     selectInput('multMAM', 'Mult-way scaling', c("No" = FALSE, 
                                                                  "Yes" = TRUE)),
                    bsPopover("multMAM", paste0("<p><b>No</b>: one scaling effect" ,
                                                  "</p>", 
                                                  "<p> <b>Yes</b>: multiple scaling",
                                                  " effects ",
                                                  "one for each product effect",
                                                  "</p>"), 
                              placement = "right", trigger = "hover")
    ))
}

## constructs tab panel for the analysis controls
tabPanel.an <- function(){
  if(input$analysis == "Consumer data")
    return(tabPanel("Analysis controls",
                    selectInput('alpharand', 
                                'Type 1 error for testing random effects', 
                                c("0.1" = 0.1, "0.2" = 0.2, "0.05" = 0.05)),
                    selectInput('alphafixed', 
                                'Type 1 error for testing fixed effects', 
                                c("0.05" = 0.05, "0.01" = 0.01, 
                                  "0.001" = 0.001)))
    )
  else
    return(tabPanel("Analysis controls",                     
                    selectInput('calc_post_hoc', 'Calculate post-hoc', 
                                c("Yes" = TRUE, "No" = FALSE)),
                    selectInput('simplerr', 'Simplification of error structure', 
                                c("Yes" = TRUE, "No" = FALSE)),
                    textInput("keep", label = "Effects to keep in a model", 
                              value = "Enter effects separated by space..."),
                    bsPopover("keep", paste0("<p>Assessor and interaction between Assessor and highest order product effects are always kept in the model</p>"), 
                              placement = "right", trigger = "hover"),
                    selectInput('alpharand', 
                                'Type 1 error for testing random effects', 
                                c("0.1" = 0.1, "0.2" = 0.2, "0.05" = 0.05)),
                    selectInput('alphafixed', 
                                'Type 1 error for testing fixed effects', 
                                c("0.05" = 0.05, "0.01" = 0.01, 
                                  "0.001" = 0.001))
    ))
}

output$antypeUI <- renderUI({ 
  #if(input$uploaddata == 1 || input$uploaddata == 2)
  #  antype <- "Sensory data"
  #else if(input$uploaddata == 3)
  #  antype <- "Consumer data"
  return(wellPanel(
              h4("Choose type of analysis"),
              radioButtons("analysis", "Analysis of", 
                           choices = c("Sensory data", "Consumer data"),
                           inline = TRUE)
             ))
})

output$AttrUI <- renderUI({ 
  if(is.null(uploadData()))
  {names.dd <- NULL}
  else{
    dd <- uploadData()
    names.dd <- colnames(dd)
  }
  tabsetPanel(
    tabPanel.input(names.dd),
    tabPanel.model(),
    tabPanel.an()
  )  
})


tabsCons <- function(){
  return(
    list(
      tabPanel("Data",
               h4("Choose data"),
               selectInput('uploaddata', '', 
                           c("Read CSV file from local drive" = 1, 
                             "TVbo data" = 2, "Ham data" = 3)),
               uiOutput("UploadUI")
      ),
      tabPanel("Step output",               
               sidebarLayout(
                 #                  conditionalPanel(condition =  "input.analysis == 'Sensory data'",
                 #                  sidebarPanel(
                 #                    uiOutput("AttrStepUI"))), 
                 sidebarPanel(
                   uiOutput("AttrStepUI")),
                 mainPanel(
                   htmlOutput("stepRand"), 
                   br(),
                   htmlOutput("stepFixed")
                 )
               )),
      tabPanel("Post-hoc",
               sidebarLayout(
                 sidebarPanel(
                   uiOutput("AttrPosthocUI"),
                   uiOutput("EffsPosthocUI")),                   
                 mainPanel(
                   plotOutput("posthocPlot"),
                   htmlOutput("posthocTable")
                 )
               ))
    )
  )
}


tabsSens <- function(){
  return(
    list(
      tabPanel("Data",
               h4("Choose data"),
               selectInput('uploaddata', '', 
                           c("Read CSV file from local drive" = 1, 
                             "TVbo data" = 2, "Ham data" = 3)),
               uiOutput("UploadUI")
               ),
      tabPanel("Plot output",
               helpText("Note: This output is only dedicated for analysis of sensory data"),
               conditionalPanel(
                 condition =  "input.analysis == 'Sensory data'",
                 sidebarLayout(
                   sidebarPanel(
                     selectInput('typeEffs', 'Plot effects', 
                                 c("random" = 1, "fixed" = 2, "scaling" = 3)),
                     selectInput('typePlot', 'Plot type', 
                                 c("F" = FALSE, "d-prime" = TRUE)),
                     selectInput('representPlot', 'Layout', 
                                 c("single" = FALSE, "multiple" = TRUE)),
                     numericInput('scalePlot', label = "Scale plot", value = 1),
                     downloadButton('downloadPlot', label = "Download Plot")   
                     ),               
                   mainPanel(plotOutput("plotsSensMixed"))
                   ), value = 1)),
      tabPanel("Table output",
           helpText("Note: This output is only dedicated for analysis of sensory data"),
           sidebarLayout(
             sidebarPanel(
               selectInput('typeEffsTable', 'Type of effects', 
                           c("random" = 1, "fixed" = 2, "scaling" = 3, 
                             "all" = 4)),
               selectInput("typetable2", "Type", c("html", "latex")),
               downloadButton('downloadTable', label = "Download Table")
             ),
             mainPanel(
               htmlOutput("tablesSensMixed")
             )
           ), 
           value = 2),
      tabPanel("Step output",               
           sidebarLayout(
             #                  conditionalPanel(condition =  "input.analysis == 'Sensory data'",
             #                  sidebarPanel(
             #                    uiOutput("AttrStepUI"))), 
             sidebarPanel(
               uiOutput("AttrStepUI")),
             mainPanel(
               htmlOutput("stepRand"), 
               br(),
               htmlOutput("stepFixed")
             )
           ),              
           value = 3),
      tabPanel("Post-hoc",
           sidebarLayout(
             sidebarPanel(
               uiOutput("AttrPosthocUI"),
               uiOutput("EffsPosthocUI")),                   
             mainPanel(
               plotOutput("posthocPlot"),
               htmlOutput("posthocTable")
             )
           ),              
           value = 4),
      tabPanel("MAM analysis",
               # Sidebar with a slider input
               
               sidebarPanel(
                 uiOutput("AttrMAManalysis"),
                 downloadButton('downloadMAM', label = "Download Table")
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 helpText("Note: Ouput only when the Correct for scaling = TRUE"),
                 
                 br(),
                 htmlOutput("MAMtable"),
                 br(),
                 htmlOutput("MAMindiv"),
                 br(),
                 htmlOutput("MAMperf"),
                 br(),
                 plotOutput("MAMplotposthoc"),
                 bsCollapsePanel("Table output", 
                                 htmlOutput("MAMposthoc"), id="colll1", 
                                 value="testlll1"),
                 br(),
                 htmlOutput("MAMdiffmean")
               ),              
               value = 5)
      ))
}

returnOutputs <- function(){
  #if(!is.null(input$uploaddata) && input$uploaddata == 2)
    
#  if( (is.null(input$analysis) && is.null(input$uploaddata)) 
#      || input$analysis == "Sensory data")
 #   return(tabsSens())
  if(is.null(input$analysis)){
   # if(!is.null(input$uploaddata) && input$uploaddata == 3)
      #return(tabsCons())
    return(tabsSens())
   # else
    #  return(tabsSens())
  }

  if(input$analysis == "Sensory data")
    return(tabsSens())
  else
    return(tabsCons())

}

output$theTabset <- renderUI({
  theOutputs <- returnOutputs()
  do.call(tabsetPanel, theOutputs)
})

output$AttrStepUI <- renderUI({
  if(is.null(Data())) {return()}    
  if(input$analysis == "Consumer data") {
    list(
      selectInput("typetable", "Type", c("html", "latex")),
      downloadButton('downloadStep', label = "Download Table")
    )
  } 
  else{
    list(
      selectInput("AttrStep", "Select attribute", names(Data()$step_res)),
      selectInput("typetable", "Type", c("html", "latex")),
      downloadButton('downloadStep', label = "Download Table"))
  }
})

output$AttrPosthocUI <- renderUI({
  if(is.null(Data()))    {return()}   
  if(input$analysis == "Consumer data") {
    selectInput("whichPlot", "Type of Plot", 
                c("LSMEANS" = "LSMEANS", 
                  "DIFF of LSMEANS" = "DIFF of LSMEANS"))
  }
  else{
    list(
      selectInput("AttrPosthoc", "Select attribute", names(Data()$step_res)),    
      selectInput("whichPlot", "Type of Plot", 
                  c("LSMEANS" = "LSMEANS", 
                    "DIFF of LSMEANS" = "DIFF of LSMEANS")))      
  }
  
  
})

output$AttrMAManalysis <- renderUI({
  if(is.null(Data()))    {return()}   
  if(class(Data()) == "consmixed") { return() }
    list(
      selectInput("AttrMAManalysis", "Select attribute", names(Data()$step_res))    
    )      
})

output$EffsPosthocUI <- renderUI({
  if(is.null(Data()))    {return()}   
  if(input$analysis == "Consumer data"){
    an.table <- Data()$anova.table
  }
  else{
    if(is.null(input$AttrPosthoc) || length(input$AttrPosthoc)>1)
    {return()}
    an.table <- Data()$step_res[[input$AttrPosthoc]]$anova.table
  }    
  
  if("elim.num" %in% colnames(an.table)){
    effs <- rownames(an.table[which(an.table$elim.num == "kept"), , 
                              drop = FALSE])
  }
  else
    effs <- rownames(an.table)
  list(
    selectInput("effsPlot", "Effects", effs),
    downloadButton('downloadPosthocTable', label = "Download Table"),
    downloadButton('downloadPosthocPlot', label = "Download Plot")
  )
})

output$UploadUI <- renderUI({    
  if(input$uploaddata == 1){ 
    verticalLayout(
      
      #tags$hr(),
      fileInput('file1', 
                'Choose CSV File from local drive, adjusting parameters if necessary',
                accept=c('text/csv', 'text/comma-separated-values,text/plain')),
      
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Semicolon=';',
                     Comma=',',
                     Tab='\t'),
                   'Semicolon'),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   'Double Quote'),
      radioButtons('decimal', 'Decimal',
                   c("Period" = ".", "Comma" = ",")),
      tags$head(tags$style(type="text/css",
                           "label.radio { display: inline-block; margin:0 10 0 0;  }",
                           ".radio input[type=\"radio\"] { float: none; }")),
      mainPanel(
        dataTableOutput('contents')
      )      
      
    )
  }
  else{
    verticalLayout(        
      mainPanel(
        dataTableOutput('contents')
      )         
    )
  }   
  
})