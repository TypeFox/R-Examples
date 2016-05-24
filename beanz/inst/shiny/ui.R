##
##  date    : 06/15/2014
##  contact : CWANG68@JHMI.EDU
##
library(shinythemes)

fluidPage(theme = shinytheme("united"),
          ##render math symbol in dynamic ui
          includeScript('www/beanz.js'),
          tags$head(tags$title("BEANZ"),
                    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
                    tags$link(rel = "stylesheet", type = "text/css", href = "//fonts.googleapis.com/css?family=Oswald"),
                    tags$link(rel = "stylesheet", type = "text/css", href = "//fonts.googleapis.com/css?family=Lora")
                    ##,tags$script(src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML",
                    ##             type="text/javascript")
                    ),
          ##title box
          withTags({
              div(class="cheader", "Bayesian Analysis of Heterogeneity of Treatment Effect")
          }),

          ## wellPanel(
          ##     fluidRow(br(),
          ##              img(src="beans_logo.png", width="80%"),
          ##              align="center"),
          ##     hr(),
          ##     style="min-height=700px"
          ##     )

          ##main page
          uiOutput("mainpage"),

          ##foot
          ## withTags({
          ##     div(class="cfooter",
          ##         "A",
          ##         a("PCORI", href="http://www.pcori.org/"), "Project.")
          ## })
          includeHTML("www/beanz_footer.html")
          )
