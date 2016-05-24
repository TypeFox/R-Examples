geochartPrereqs <- tagList(
  tags$head(
    tags$script(src="https://www.google.com/jsapi"),
    tags$script(src="geochart.js")
  )
)

geochart <- function(id, options=list()) {
  tags$div(id=id, class="shiny-geochart-output", `data-options`=RJSONIO::toJSON(options))
}

googleLineChart <- function(id, options=list()) {
  tags$div(id=id, class="google-linechart-output", `data-options`=RJSONIO::toJSON(options))
}

googleHistogram <- function(id, options=list()) {
  tags$div(id=id, class="google-histogram-output", `data-options`=RJSONIO::toJSON(options))
}

row <- function(...) {
	tags$div(class="row", ...)
}

col <- function(width, ...) {
	#tags$div(class=paste0("span", width), ...)
 tags$div(class=paste0("col-sm-", width), ...)
}

 
shinyUI(pageWithSidebar(
  headerPanel(paste("WPP", substr(wppExplorer:::wpp.data.env$package, 4, 8), "Explorer")),
  sidebarPanel(
    geochartPrereqs,
    tags$head(
         		#tags$style(type="text/css", ".jslider { max-width: 50px; }"),
         		#tags$style(type='text/css', ".well { padding: 0px; margin-bottom: 5px; max-width: 100px; }"),
		tags$style(type='text/css', ".span4 { max-width: 270px; }")
	),
    uiOutput('yearUI'),
    hr(),
    selectInput('indicator', h5('Indicator:'), wppExplorer:::wpp.data.env$indicators),
    conditionalPanel(condition=paste("input.indicator >", sum(attr(wppExplorer:::wpp.data.env$indicators, "settings")$by.age == FALSE)),
    	tags$head(tags$style(type="text/css", "#selagesmult { height: 150px; width: 70px}"),
    			  tags$style(type="text/css", "#selages { height:25px; width: 70px}"),
    			  tags$style(type="text/css", "#indsexmult { height: 50px; width: 90px}"),
    			  tags$style(type="text/css", "#indsex { height: 25px; width: 90px}")),
    	row(
			col(2, ''),
			col(4, uiOutput('sexselection')),
			col(1, ''),
    		col(3, uiOutput('ageselection'))
    	)
    ),
    textOutput('indicatorDesc'),
    hr(),
    selectInput('uncertainty', h5('Uncertainty:'), structure(as.character(1:3), names=c('80%', '95%', '+-1/2child')), 
    			multiple=TRUE, selected=1),
    textOutput('uncertaintyNote'),
    hr(),
    HTML("<p><small><b>Data Source:</b> United Nations, Department of Economic and Social Affairs, Population Division: <a href='http://esa.un.org/unpd/wpp' target='_blank'>World Population Prospects</a>. <a href='http://esa.un.org/unpd/ppp' target='_blank'>Probabilistic projections</a> based on <a href='http://www.pnas.org/content/early/2012/08/13/1211452109.abstract' target='_blank'>Raftery et al. (2012, PNAS)</a></small></p><p><small>Hana &#352;ev&#269;&#237;kov&#225;, <a href='https://www.csss.washington.edu' target='_blank'>CSSS</a>, University of Washington; <a href='http://bayespop.csss.washington.edu' target='_blank'>project website</a></small></p>"),
width=3
  ),
  mainPanel(
    tabsetPanel(
      tabPanel('Map',
      	 flowLayout(
			textOutput('mapyear'),
			checkboxInput('normalizeMapAndCountryPlot', 'Fixed scale over time', TRUE)
		),		
		hr(),
		geochart('map'),
		#htmlOutput('mapgvis'),
		hr(),
		conditionalPanel(condition='input.map_selection',
				plotOutput('countryPlot', height='300px'))
      ),
      #tabPanel('Data', 
		#textOutput('year1'),
		#checkboxInput('includeAggr1', 'Include Aggregations', FALSE),
      	#tableOutput('table')
      #),
      tabPanel('Sortable Data', 
      	flowLayout(
			textOutput('year2'),
			checkboxInput('includeAggr2', 'Include Aggregations', FALSE)
		),
		hr(),
      	DT::dataTableOutput('stable')
      ),
      tabPanel('Trends, Age profiles & Pyramids',
  		tags$head(
			tags$style(type="text/css", "#seltcountries { height: 450px; width: 150px}")
			),
			tags$div(
				class = "container",
				row(
					col(0.5, ''),
					col(2, uiOutput('cselection')),
				  	col(7, tabsetPanel(
				  				#tabPanel('Median',
				  				#	googleLineChart('trends', options=list(height=400, width=650)),
				  				#	checkboxInput('median.logscale', 'Log scale', FALSE)),
				  				tabPanel('Trends', 
				  					plotOutput('probtrends', height="400px", width="650px", 
				  							click = "probtrends_values", hover = "probtrends_values", 
				  							dblclick = "probtrends_zoom_reset", 
				  							brush = brushOpts(id = "probtrends_zoom", resetOnNew = TRUE)),
				  					flowLayout(
										checkboxInput('trend.logscale', 'Log scale', FALSE),
										textOutput("probtrends_selected")
										)
				  					),
				  				tabPanel('Age Profile', 
				  					googleLineChart('age.profileM', options=list(height=200, width=650)),
				  					googleLineChart('age.profileF', options=list(height=200, width=650)),
				  					checkboxInput('aprofile.logscale', 'Log scale', FALSE)),
				  				tabPanel('Pyramids', 
				  					plotOutput('pyramids', click = "pyramid_values", hover = "pyramid_values", 
				  						dblclick = "pyramid_zoom_reset",
				  						brush = brushOpts(id = "pyramid_zoom", resetOnNew = TRUE)
									),
				  					flowLayout(
				  						checkboxInput('proppyramids', 'Pyramid of proportions', FALSE),
				  						textOutput("pyramid_selected")
									)
				  				)	  			
				  			)
				  		)
 					),
 				row(
 					col(0.5, ''),
 					col(9, textOutput('trendstabletitle'),
 						   tableOutput('trendstable'))
 					)
 				)
 		),
 	tabPanel('Histogram',
 		#textOutput('year3'),
 		#hr(),
 		flowLayout(
      	checkboxInput('fiXscaleHist', 'Fixed x-axis over time', FALSE)
      	),
      	hr(),
      	#plotOutput('hist')
      	htmlOutput('ghist')
    ),
      tabPanel('Rosling Chart',
		htmlOutput('graphgvis'),
		row(col(1, "")),
		row(
			col(1,""),
			col(3, textOutput('AddIndicatorText')),
			col(1, actionButton("AddIndicator", "Add indicator"))
			)
      ),
      tabPanel("Help",
      	includeHTML("README.html")
      )
  ) #end tabsetPanel
  ) #end mainPanel
))