library(shiny)
#if (any(installed.packages()[,1]=="shinyIncubator")) library(shinyIncubator)
# runApp(".", launch.browser = TRUE)

# Define UI for STSRE application
shinyUI(pageWithSidebar(
  
  # STSRE title
  headerPanel(
    img(src="logo.png", height=126, width=626),
    tags$head(tags$style("body {background-color: #3078A6;}"))
  ),
  
  # Sidebar with tabs
  sidebarPanel(
    tabsetPanel(
      tabPanel(
        headerPanel(img(src="b0.png", height=35, width=35)),
        wellPanel(
          p(strong("The phenology package is a new tool to easily estimate density of animal present in a site partially during the year.")),
          tags$hr(),
          p("Phenology is developped by ",a("Marc Girondot", 
                                                    href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html",
                                                    target="_blank")),
          tags$hr(),         
          p(img(src="ACP_PDF.png", height=35, width=35), " The method has been published in ",a("Girondot, M. 2010. Estimating density of animals during migratory waves: application to marine turtles at nesting site. Endangered Species Research, 12, 85-105.", 
                                            href="http://www.int-res.com/articles/esr2010/12/n012p095.pdf",
                                            target="_blank"))
           ),
        
        value=0),
      
      
      
      tabPanel(
        headerPanel(img(src="b1.png", height=35, width=35)),
        wellPanel(
          p(strong("FIELD COUNTS")),
          fileInput("file1", "Please input a datasheet with counts", accept=c('text/csv', 
                    'text/comma-separated-values,text/plain', '.csv')),
          uiOutput("results"),
          helpText("Note: The input file must be as .csv or .txt. The first column must be the date and the second be the counts. An optional third column can be use to aggregate several datasets in a single analysis.")
        ),      
        value=1),
      
      id="conditionedPanels")),
  
  
  # Mainbar with tabs
  mainPanel(
#    progressInit(),
    conditionalPanel(
      condition="input.conditionedPanels==0",
      img(src="00.jpg", height=600, width=600)),
    
    
    
    conditionalPanel(
      condition="input.conditionedPanels==1",
      verbatimTextOutput(outputId="resultsInfo"),
      plotOutput(outputId="resultPlot", width="600px", height="400px")
    )
    
    )))