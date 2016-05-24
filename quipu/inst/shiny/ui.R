library(shiny)
library(quipu)

data(potato.quipu)
dat = potato.quipu

sort(unique(dat$accession_id))
ids = sort(unique(dat$accession_id))
# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title.
  headerPanel("Quipu demo"),
  
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view. The helpText function is also used to 
  # include clarifying text. Most notably, the inclusion of a 
  # submitButton defers the rendering of output until the user 
  # explicitly clicks the button (rather than doing it immediately
  # when inputs change). This is useful if the computations required
  # to render output are inordinately time-consuming.
  sidebarPanel(
   selectInput("acc_id", "Choose an accession:", 
                choices = ids ),
   radioButtons("layout", "Layout",list("full","no text"), "full"),  
   textInput("idLabel","Identifier label","ID label"),
   textInput("speciesName","Species name","Sample species"),
   textInput("setName","Set name","Sample set"),
  
   sliderInput("nodeG1", "Node color (group 1)",0.4, 2, 1.5, 0.1),
   sliderInput("nodeG2", "Node color (group 2)",0.4, 2, 1.2, 0.1),
   sliderInput("nodeG3", "Node color (group 3)",0.4, 2, 0.9, 0.1),
   sliderInput("nodeG4", "Node color (group 4)",0.4, 2, 0.6, 0.1),
   
   selectInput("colorG1", "Choose a color for group 1", choices = colors(), "red3"),
   selectInput("colorG2", "Choose a color for group 2", choices = colors(), "green"),
   selectInput("colorG3", "Choose a color for group 3", choices = colors(), "blue"),
   selectInput("colorG4", "Choose a color for group 4", choices = colors(), "grey50")
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations. Note the use of the h4 function to provide
  # an additional header above each output section.
  mainPanel(
    
    plotOutput("quipuPlot")#,
    
    #h4("Observations")
  )
))
