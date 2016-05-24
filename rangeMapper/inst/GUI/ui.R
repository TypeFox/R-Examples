library(shiny)
library(rangeMapper)

path = options('rangeMapper.path')$rangeMapper.path
con = rangeMap.open(path)
tabs = rangeMapProjInfo(con)
tabs[grep('^MAP_', tabs$tabname), 'rmap_type'] = 'map'
tabs[grep('^BIO_', tabs$tabname), 'rmap_type'] = 'bio'

M = tabs[which(tabs$rmap_type == 'map'), 'tabname']

cols =  c( list(rangeMapper = palette_rangemap('divergent')), rangeMapper:::brewer.pal.get () )


shinyUI(fluidPage(

  titlePanel( basename(path) ),


  sidebarLayout(
    sidebarPanel(
     selectInput('maps', 'Maps', M, selected = M[1] , multiple = TRUE),
      br(),
     selectInput('colorRamp', 'Color ramp', names(cols) )


    ),

    mainPanel(
      tabsetPanel(type = "tabs",
        tabPanel("Plot", plotOutput("plot")),
        tabPanel("Summary", tableOutput("summary"))
      )
    )
  )

))

