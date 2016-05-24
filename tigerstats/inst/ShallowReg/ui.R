library(shiny)

# Define UI for SlowGoodness application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Is the Regression Line too Shallow?"),
  
  # Sidebar
  sidebarPanel(
    
    sliderInput("n","Number of Points in the Cloud",min=100,max=2000,value=1200),
   
    sliderInput("rho","Target Correlation of Cloud",min=-1,max=1,step=.01,value=0.50),
    
    checkboxInput("showlines","Show SD Line and Regression Line",value=FALSE),
    helpText("For some values of the correlation, the regression line",
             "does not appear to pass through the cloud of points",
             "as well as the SD line does. Try various correlations",
             "and see for yourself!"),
    br(),
    helpText("But the regression line does a better job of predicting y-values from",
             "x-values.  To see this look at a vertical slice of the cloud:"),
    
    checkboxInput("showslice","Show a Slice of the Cloud",value=FALSE),
    
    helpText("The big point in the green slice indicates the mean y-value",
             "for all points in the slice.  Which line comes closer to it?",
             "Try some other slices, too."),
    
    sliderInput("slice","Choose a New Slice",min=1,max=10,value=4),
    
    helpText("To form an overall impression of how the regresion line",
             "is doing, look at the means for all of the slices at once:"),
  
    checkboxInput("showmeans","Show Means of All Slices",value=FALSE)

    
#     conditionalPanel(
#       condition="output.showslice==TRUE",
#       sliderInput("slice","Choose a New Slice",min=1,max=10,value=4)
#       ),
#     
#     conditionalPanel(
#       condition="output.showslice==TRUE",
#       checkboxInput("showmeans","Show Means of All Slices",value=FALSE)
#     )
  ),
  
  
  # Here comes the main panel
  
  mainPanel(
    
    plotOutput("cloud")
    
  )
  
))
