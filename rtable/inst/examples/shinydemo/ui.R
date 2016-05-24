
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("FlexTable show off"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("margin",
                  "Margin size:",
                  min = 0,
                  max = 20,
                  value = 4), 
      sliderInput("fontsize",
                  "Font size:",
                  min = 5,
                  max = 20,
                  value = 8), 
      checkboxInput("tablestyle", "Space columns:", TRUE )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      
      fluidRow(
        column(width = 8, align="center", 
          tags$h3("Font colors")
        )
      ),
      
      fluidRow(
        column(width = 2, align="center", 
               textInput("q1", label = "Q1", value = "#fee5d9") ), 
        column(width = 2, align="center", 
               textInput("q2", label = "Q2", value = "#fcae91") ),  
        column(width = 2, align="center", 
               textInput("q3", label = "Q3", value = "#fb6a4a") ),  
        column(width = 2, align="center", 
               textInput("q4", label = "Q4", value = "#cb181d") )
      ), 
      
      tableOutput("flextable"), 
      
      fluidRow(
        column(width = 8, align="center", 
               tags$h3("Download table")
        )
      ),
      
      fluidRow(
        column(width = 8, align="center", 
           radioButtons("filetype", "File type:", inline = TRUE, 
            choices = c("Word document"=".docx", "PowerPoint presentation" = ".pptx"))
        )
      ),
      
      fluidRow(
        column(width = 8, align="center", 
          downloadButton('downloadData', 'Download') 
        )
      )
      
    )
  )
))
