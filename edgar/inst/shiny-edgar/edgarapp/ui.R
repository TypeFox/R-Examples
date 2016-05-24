# create year list from 1994 to current year
tableyr <- as.character(c(1994:as.numeric(format(Sys.Date(), "%Y"))))
masteryr <- c("ALL",tableyr)

sidebar <- dashboardSidebar(
	sidebarMenu( 
	menuItem("Home", tabName = "homtab", icon = icon("home")),
	menuItem("Get Master Index", tabName = "getmaster", icon = icon("dashboard")),
	menuItem("Filings Info", icon = icon("th"), tabName = "getFilings"),
	menuItem("Download Filings", icon = icon("download"), tabName = "Filingsdownload"),
	menuItem("Get Daily EDGAR Info", icon = icon("th"), tabName = "getdailyinfo"),
	menuItem("Sentiment Analysis of 10K", icon = icon("cog"), tabName = "senti")
    ))

body <- dashboardBody(
	tabItems( 
		tabItem(tabName = "homtab",
            # apply css for shiny app
		        tags$head(
		          tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
		        ),
            actionButton("change.wd", "Change working directory"),
            br(),br(), tags$img(src = 'edgarimg.jpg', heigth = 700, width = 1100)
        ),
      
		tabItem(tabName = "getmaster",
			fluidRow(
				box(
					width = 6, status = "primary", solidHeader = TRUE,background = "yellow",
					title = "Download EDGAR Master ZIP Files",
					selectizeInput("getmasteryears", "Select Years:" ,choices = masteryr, 
									selected = NULL, multiple = TRUE, options = NULL),
					actionButton("down.master", "Download Master"),br(), br(),
					uiOutput("master.download.status")
				)
			)
		),
    
		tabItem(tabName = "getFilings",
		    fluidRow(
				box(
					width = 6, status = "primary", solidHeader = TRUE, background = "yellow",
					title = "Get Yearly EDGAR Filing Information",            
							selectInput("year.getinfo", "Select Year:",choices = c("",tableyr), selected = NULL),
				    actionButton("getinfo", "Get Info")
				),
				box(
				  width = 12, status = "primary", solidHeader = TRUE,
				  title = "Filing Information",            
				  uiOutput("tb"),
				  tags$head(tags$style("tfoot {display: table-header-group;}"))
				)
			)
				
		),

		tabItem(tabName = "Filingsdownload",
		    fluidRow(
				box(
					width = 6, status = "primary", solidHeader = TRUE, background = "yellow",
					title = "download EDGAR Filings", 
					selectInput("year3", "Select Years:", choices = c("",tableyr), selected = NULL),
					textInput("cik3", "Type CIK"),
					textInput("ftype3", "Type FORM TYPE"),
					actionButton("downFilings", "download Filings")
				),
				box(
				  width = 12, status = "primary", solidHeader = TRUE,
				  title = "Filing Download Status",            
				  uiOutput("stat.table.Filing")
				)
			)				
		),
    
		tabItem(tabName = "getdailyinfo", 
			fluidRow(
				box(
					width = 6, status = "primary", solidHeader = TRUE,background = "yellow",
					title = "Get Daily EDGAR Filing Information", 
					dateInput("dailydate", label = 'Select the date:', value = NULL),
					actionButton("get.dailyinfo", "Get Info")
				),
				box(
				  width = 12, status = "primary", solidHeader = TRUE,
				  title = "Daily Filing Information",            
				  uiOutput("daily.tb")
				)
			)
					
		),
		
		tabItem(tabName = "senti",
			fluidRow(
				box(
					width = 7, status = "primary", solidHeader = TRUE, background = "yellow",
					title = "Sentiment Analysis of 10K statement",
					actionButton("get.sentianalysis", "choose 10K Statement")
				),

				box(
					width = 6, status = "primary", solidHeader = TRUE,
					collapsible = TRUE, title = "Negative wordcloud",
					uiOutput("negwordcloud1")
				),
				box(
					width = 6, status = "primary", solidHeader = TRUE,  
					collapsible = TRUE, title = "Negative words histogram",
					uiOutput("neghist1")
				),
        
				box(
					width = 6, status = "primary", solidHeader = TRUE,
					collapsible = TRUE, title = "Positive wordcloud",
					uiOutput("poswordcloud1")
				),
        
				box(
					width = 6, status = "primary", solidHeader = TRUE, 
					collapsible = TRUE, title = "Positive words histogram",
					uiOutput("poshist1")
				),
				box(
					width = 6, status = "primary", solidHeader = TRUE, 	
					collapsible = TRUE,  title = "Polarity histogram",
					uiOutput("polarity.hist1")
				)
			)
		)
	)
)

# create dashboardPage
dashboardPage( skin="yellow",
		dashboardHeader(title = "EDGAR Filing mgmt"),
		sidebar,
		body
)
