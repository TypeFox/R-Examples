# shiny microPlexer app user interface

# Kevin Keenan 2013

library("shiny")

# ui.R

# Define the user interface

shinyUI(pageWithSidebar(
  
  # App name
  headerPanel("microPlexer"),
  
  # Define sidebar contents
  
  sidebarPanel(
    
    fileInput('infile', 'Upload size range .csv file'),
      
    #textInput("chan", "How many fluorophore tags are you using?"),
    
#       sliderInput("obs", "Number of observations:", 
#                   min = 0, max = 1000, value = 500),
#       
#       # Show text when slider is set to 0
#       conditionalPanel(condition = "input.obs == 0",
#                        p("Conditional text here")),
    
    wellPanel(
      checkboxGroupInput("dyeCol", 
                         "Specify the colours of your fluorophores",
                         choices = list("yellow", "red", 
                                        "blue", "green",
                                        "orange")
                         )
    ),
        

    
    numericInput("proximity",  
                 "What is the minimum distance allowable between loci?",
                 value = 20, min = 0, max = 200, step = 10), 
                  
    radioButtons(inputId = "algorithm", label = "Grouping algorithm:",
                 choices = list("Maximum throughput" = "max",
                                "Balanced throughput" = "balanced")),

    
    wellPanel(
      sliderInput("maxLoci", "If you are using the balanced algorithm, \n
                  how many loci per multiplex would like, on average?",
                  min = 0, max = 100, step = 2, value = 10)
    ),
    
    submitButton(text = "Run microPlexer")
     
  ),
  
  mainPanel(
    downloadButton("dlplt"),
    plotOutput("plots", height = 3000)
  )
  
))