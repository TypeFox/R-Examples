library(bde)


# HEADER ------------------------------------------------------------------
header<-header<-list(h1('Bounded Density Estimation'),
                     img(src='./logos.png'))


# FOOTER ------------------------------------------------------------------

footer<-list()


# MAIN --------------------------------------------------------------------

# Common controls ---------------------------------------------------------
common<-list(
  uiOutput(outputId="main_opts")
)




# Density plot tab --------------------------------------------------------

density_tab<-list(plotOutput(outputId="plot_density"),
                  fluidRow(  
                    column(2,offset=1,
                           selectInput(inputId="density_plot_format",label="File format",
                                       choices=list("pdf" = "pdf",
                                                    "eps" = "eps",
                                                    "png" = "png",
                                                    "jpg" = "jpg"))),
                    column(2,
                           h5("Width"),
                           numericInput(inputId="density_plot_width",label="",
                                        value=1024,min=100,max=2500)),
                    column(2,
                           h5("Height"),
                           numericInput(inputId="density_plot_height",label="",
                                        value=1024,min=100,max=2500)),
                    column(2,
                           h5("Resolution (in dpi)"),
                           numericInput(inputId="density_plot_resolution", label="", 
                                                         value=200,min=50,max=600)),
                    column(2,br(),br(),
                           downloadButton(outputId="density_plot_save",label="Save plot"))))


# Distribution tab --------------------------------------------------------

distribution_tab<-list(plotOutput(outputId="plot_distribution"),
                      fluidRow(  
                        column(2,offset=1,
                               selectInput(inputId="distribution_plot_format",label="File format",
                                           choices=list("pdf" = "pdf",
                                                        "eps" = "eps",
                                                        "png" = "png",
                                                        "jpg" = "jpg"))),
                        column(2,
                               h5("Width"),
                               numericInput(inputId="distribution_plot_width",label="",
                                            value=1024,min=100,max=2500)),
                        column(2,
                               h5("Height"),
                               numericInput(inputId="distribution_plot_height",label="",
                                            value=1024,min=100,max=2500)),
                        column(2,
                               h5("Resolution (in dpi)"),
                               numericInput(inputId="distribution_plot_resolution", label="", 
                                            value=200,min=50,max=600)),
                        column(2,br(),br(),
                               downloadButton(outputId="distribution_plot_save",label="Save plot"))))

# Data tab ----------------------------------------------------------------

data_tab<-dataTableOutput("samplePoints")


# Extra -------------------------------------------------------------------

extra<-uiOutput(outputId="additional_opts")


# Main --------------------------------------------------------------------

main<-list(fluidRow(
  column(2,common),
  column(8,tabsetPanel(id="tabset",
                       tabPanel(title="Density plot",density_tab),
                       tabPanel(title="Distribution plot",distribution_tab),
                       tabPanel(title="Sample points",data_tab))),
  column(2,extra)))


# Load data ---------------------------------------------------------------

load<-fluidRow(column(2,fileInput(inputId="file",label="")))


# INTERFACE ---------------------------------------------------------------

shinyUI(fluidPage(title="bde GUI",theme="bde.css",
                  header,br(),
                  load,
                  main,
                  footer))