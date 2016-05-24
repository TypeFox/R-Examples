### ui.R
#library(shiny)

message(paste0("\n***ATTENTION*** If web tool doesn't open automatically, ",
  "open your browser and point to the url above. ","Hit ESC or Ctrl+C ",
  "to stop the web tool."))


privateSales2 <- privateSales[order(privateSales$u5Sales, decreasing = TRUE),]
names2 <- privateSales2$Country
names2 <- paste(names2)
names2 <- names2[order(names2)]

for (i in 1:length(names2)){
	if (names2[i] == "Cote d'Ivoire"){
		names2[i] <- "Coté d'Ivoire"
	}
}

#print(names2)

# Define UI for slider demo application
shinyUI(pageWithSidebar(

	# Application title
	h2("Estimated Under-Five Malaria Deaths Due to Poor Quality Antimalarials"),

	# Sidebar with sliders that demonstrate various available options
	sidebarPanel(
    	h4("Input Selection"),

    	actionButton("goRefresh","Click to Refresh"),
      h4(""),

    	downloadButton('downloadInputs', 'Download Current Inputs as CSV'),
    	h4(""),

    	wellPanel(
    	# Simple integer interval
    	p("Latin Hypercube Sample Size"),
    	sliderInput("N", "",
        	       min=1000, max=10000, value=10000, step = 500),
       	helpText("Reduce sample size if you experience slow performance.")
       	),

    	wellPanel(
  		# Decimal interval with step value
  		p("Case Fatality Rate"),
    	sliderInput("CFR", "",
        	        min = 0, max = 0.05, value = c(0.002, 0.006), step= 0.001),
		helpText("Baseline Range = [0.002, 0.006]")
		),

		wellPanel(
  		# Specification of range within an interval
  		p("Prevalence of Poor Quality Antimalarials"),

    	conditionalPanel(
    		condition = "((input.countrySpecific != true) && (input.manuscript != true))",
    		sliderInput("Prev", "",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
   		),

  		conditionalPanel(
  		  condition = "input.manuscript != true",
  		  checkboxInput("countrySpecific","Check to specify different ranges for each country", value = FALSE)
  		),

    	conditionalPanel(
    		 condition = "input.countrySpecific == true",
    		 helpText("Click refresh button to update results after changing prevalence ranges.")
    	),

  		conditionalPanel(
  		  condition = "input.countrySpecific != true",
  		  checkboxInput("manuscript","Check to use the prevalence distributions from manuscript", value = TRUE)
  		),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		actionButton("goButton","Refresh")
    	),

   		conditionalPanel(
    		condition = "input.countrySpecific == true",
    		helpText("")
    	),

    	#conditionalPanel(
    	#	condition = "input.countrySpecific == true",
    	#	submitButton("Update View")
    	#),

    	########## START COUNTRY SPECIFIC CONDITIONAL PANELS ###########
    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev21", "Angola",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev23", "Benin",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev30", "Burkina Faso",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev12", "Burundi",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev29", "Cameroon",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev13", "Central African Republic",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev24", "Chad",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev9", "Congo",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev31", "Coté d'Ivoire",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev38", "Democratic Republic of the Congo",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev2", "Djibouti",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev3", "Equatorial Guinea",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev33", "Ethiopia",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev7", "Gabon",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev6", "Gambia",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev32", "Ghana",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),


    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev19", "Guinea",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev5", "Guinea-Bissau",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev35", "Kenya",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev11", "Liberia",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev25", "Madagascar",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev22", "Malawi",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev28", "Mali",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),


    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev8", "Mauritania",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev27", "Mozambique",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev4", "Namibia",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),


    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev26", "Niger",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev39", "Nigeria",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev20", "Rwanda",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev15", "Senegal",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	), # End Country 2 conditional Panel

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev18", "Sierra Leone",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev16", "Somalia",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev36", "Sudan",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev1", "Swaziland",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev34", "Tanzania",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev17", "Togo",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	), # End Country 2 conditional Panel

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev37", "Uganda",
    		            min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	), # End Country 2 conditional Panel

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev14", "Zambia",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	),

    	conditionalPanel(
    		condition = "input.countrySpecific == true",
    		sliderInput("Prev10", "Zimbabwe",
     		         min = 0, max = 1, value = c(0,0.4), step = 0.01)
    	)#,

  		########## END COUNTRY SPECIFIC CONDITIONAL PANELS ###########

    	) #End Prevalence well panel

  	),

 	mainPanel(
    	tabsetPanel(
    		tabPanel("Median Estimates",
    				h4("Median Death Estimates"),
    				downloadButton('downloadMedianEstimates', 'Download as PDF'),
    				plotOutput("raw"),
    				wellPanel(p("Median estimates of under-five malaria deaths caused by poor quality antimalarials",
    				"(error bars as interquartile range).")),

    				h4("Median Death Estimates as a Proportion of Total Malaria Related Deaths"),
    				downloadButton('downloadMalariaProp', 'Download as PDF'),
    				plotOutput("malariaProp"),
    				wellPanel(p("Median estimates of under-five malaria deaths caused by poor quality antimalarials",
    				"as a proportion of 2010 under-five malaria deaths (error bars as interquartile range)",
    				 a(href="http://apps.who.int/gho/data/node.main.GBDC-YEARS0-4?lang=en", "[Source: WHO Global Health Observatory Data Repository 2010]."))),

					h4("Median Death Estimates as a Proportion of Total All-Cause Deaths"),
					downloadButton('downloadDeathProp', 'Download as PDF'),
    				plotOutput("deathProp"),
    				wellPanel(p("Median estimates of under-five malaria deaths caused by poor quality antimalarials",
    				"as a proportion of 2012 under-five all-cause deaths (error bars as interquartile range)",
    				a(href="http://apps.who.int/gho/data/node.main.525", "[Source: WHO Global Health Observatory Data Repository 2012].")
    				 ))), #end tab panel

    		tabPanel("Summary Statistics",
    				h4("Summary Statistics"),
    				wellPanel(p("Estimated under-five malaria deaths caused by poor quality antimalarials",
    				"in 2013 (n = Latin Hypercupe Sample Size).")),
    				downloadButton('downloadSummary', 'Download as CSV'),
    				h4(""),
    				tableOutput("summaryTable")
    				),#end tab panel

    		tabPanel("Estimate Distributions",
    				h4("Estimate Distributions"),
       				selectInput("histogram1", "Choose a Country:", choices = names2),
    				h4(""),
    				downloadButton('downloadHist', 'Download as PDF'),
    		 		plotOutput("requestedHist"),
    		 		wellPanel(p("Distribution of estimated under-five malaria deaths caused by poor quality",
    				"antimalarials in 2013 (n = Latin Hypercube Sample Size). A",
    				"vertical red line is plotted at the median. A dotted vertical blue line is plotted at the mean.",
    				"Histogram bin sizes were calculated using Sturges' formula."))
    				),#end tab panel

# ## Sensitivity Analysis tab panel, currently commented out
#     		tabPanel("Sensitivity Analysis",
#     			h4("Sensitivity Analysis"),
#     			wellPanel(
#     			p("A sensitivity analysis was performed to see which input parameters were most",
#     			"responsible for the imprecision of our model output: the number of under-five malaria deaths",
#     			"(across the 39 sub-Saharan nations) caused by treatment with poor quality antimalarials.",
#     			"In our model there are three types of input parameters:"),
#     			p("1. The case fatality rate of under-five children who are treated with poor quality",
#     			"antimalarials (applies to all countries)"),
#     			p("2. Private sector antimalarial sales to malaria positive under-five children (country specific estimates)"),
#     			p("3. The proportion of private sector antimalarials that are poor quality (country specifc estimates)")),
#
#     			downloadButton('downloadPRCC', 'Download as CSV'),
#     			h4(""),
#     			tableOutput("PRCC")
#     			),#end tab panel

    		### input tab panel, currently commented out
    		# tabPanel(
    			# "Inputs",
    			# h4("Input Parameters"),
    			# tableOutput("InputParameters")
    			# ),#end tab panel

    		tabPanel("About",
    			h4("Methods"),
    			wellPanel(p("For each country we calculated the number of under-five deaths caused by malaria treatment failure due to consumption of poor quality antimalarials as the product of three inputs: the number of private sector antimalarials consumed by malaria-positive children in 2013, the proportion of private sector antimalarials consumed that are poor quality, and the case fatality rate (CFR) of under-five malaria-positive children who receive poor quality antimalarials. We selected the 39 sub-Saharan nations included in our analysis because antimalarial consumption estimates were available. Probability distributions were constructed for each input parameter, and an uncertainty analysis was conducted according to the Latin hypercube sampling method. Please see the publication for more details and references.")),
    			h4("Authors"),
    			wellPanel(p("J. Patrick Renschler, Kelsey Walters, Paul Newton, Ramanan Laxminarayan. \"Estimated under-five deaths associated with poor-quality antimalarials in sub-Saharan Africa\".  2014. Paper submitted.")),
				h4("Help"),
				wellPanel(p("For help refer to the github repository README:"),a(href="https://github.com/renschler/pqantimalarials", "https://github.com/renschler/pqantimalarials"),p("\n"),p("or email:"), a(href= "mailto:patrick.renschler@gmail.com", "patrick.renschler@gmail.com"))
    			)#end tab panel
		)
  	)
))
