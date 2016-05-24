library(shiny)
library(shinyAce)



shinyUI(navbarPage("IRTShiny", windowTitle="IRTShiny v1.1",

# ########## Adding loading message #########
#
# tags$head(tags$style(type="text/css", "
# #loadmessage {
# position: fixed;
# top: 0px;
# left: 0px;
# width: 100%;
# padding: 10px 0px 10px 0px;
# text-align: center;
# font-weight: bold;
# font-size: 100%;
# color: #000000;
# background-color: #CCFF66;
# z-index: 105;
# }
# ")),
#
# conditionalPanel(condition="$('html').hasClass('shiny-busy')",
# tags$div("Loading...",id="loadmessage")),
#
# ########## Added up until here ##########



#         navbarMenu("Dichotomous Models", icon = icon("chevron-right", lib = "font-awesome"),
        tabPanel("Dichotomous Models", icon = icon("chevron-right", lib = "font-awesome"),
                 sidebarPanel(
                   p(' '),
                   radioButtons("type", strong("Model selection:"),
                                list("One-parameter logistic model (1PL)" = "1PL",
                                     "Two-parameter logistic model (2PL)" = "2PL",
                                     "Three-parameter logistic model (3PL)" = "3PL"
                                ), selected = "2PL"
                   ),
                   submitButton("Update View"),
                   helpText("Click here to update your results, you need to do this after you change the data, model, or any of the settings"),
                   br(),
                   actionButton("quit", "Quit", icon("sign-out")),
                   helpText("Press Quit to exit the application")

                 ),
                 br(),


        br(),

        p('Note: Input values (either numeric or character) must be separated by tabs. Copy and paste from Excel.'),
        p('If you need to create 1-0 data first, use the Binary (1-0) Data Converter.'),
        p("Please make sure that your data includes the header (variable names) in the first row."),

        aceEditor("text1", value="Q01\tQ02\tQ03\tQ04\tQ05\tQ06\tQ07\tQ08\tQ09\tQ10\tQ11\tQ12\tQ13\tQ14\tQ15\tQ16\tQ17\tQ18\tQ19\tQ20\tQ21\tQ22\tQ23\tQ24\tQ25\tQ26\tQ27\tQ28\tQ29\tQ30\tQ31\tQ32\tQ33\tQ34\tQ35\tQ36\tQ37\tQ38\tQ39\tQ40\tQ41\tQ42\tQ43\tQ44\tQ45\tQ46\tQ47\tQ48\tQ49\tQ50\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t0\t1\t0\t0\t0\t1\t1\t1\t0\t0\t0\t0\t0\n1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t0\t0\t0\t0\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t0\n1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t0\t1\t1\t0\t0\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\n0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\n0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t0\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t0\t0\t1\t0\t1\t0\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t0\t0\t1\n0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t0\t0\t0\t0\t0\t1\t0\t0\t1\t1\t1\t0\t0\t0\t1\t1\t1\t0\t0\t0\t0\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t0\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\n0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t0\t1\t0\t1\t1\t1\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t0\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t0\t1\t0\t0\t1\t1\t0\t1\t0\n1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t1\t1\t1\t1\t1\n0\t0\t0\t1\t0\t0\t1\t0\t0\t1\t0\t0\t1\t1\t0\t1\t1\t0\t0\t0\t1\t1\t0\t1\t1\t1\t1\t0\t0\t0\t1\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t1\t1\t0\t1\t0\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t0\t1\t0\t1\t0\t0\t0\t0\t1\t1\t1\t0\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\n0\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t0\t0\t0\t0\t1\t0\t1\t1\t0\t0\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\n1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0\t0\t1\t1\t1\t0\t0\t1\t0\t0\t0\t1\t0\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t1\t1\t0\t0\t0\t1\t0\t1\t1\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t0\t0\t0\t0\t1\t0\t1\t1\t1\t0\t1\t0\t1\t1\t0\t0\t1\t1\t1\t0\t0\t0\n0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t0\t0\t1\t1\t0\t1\t1\n1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t1\t0\t1\t1\t1\t1\t1\t0\t0\t0\t1\t0\t0\t0\n0\t1\t0\t1\t0\t0\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t0\t0\t1\t0\t1\t1\t0\t0\t0\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t1\t1\t0\t0\t0\t0\t0\t0\t0\n0\t1\t0\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t1\t0\t0\t0\t0\t1\t0\t1\t1\t0\t0\t1\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t0\n1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\t0\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t0\t0\t0\t0\t1\t0\t1\t1\t0\t1\t1\t0\t0\t0\t1\t1\t0\t0\t0\n0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\n0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t1\t0\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t0\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t0\t0\t0\t0\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t0\t0\t0\t1\t0\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\n1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t0\t0\t1\t1\t0\t0\t1\t0\t1\t0\t0\t0\t1\t1\t1\t1\t0\t0\n1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t1\t0\t0\t1\t1\n0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t0\t0\t0\t1\t1\t0\t0\t1\t0\t0\t0\t1\t0\t1\t1\t0\t0\t1\t0\t1\t0\t1\t0\t1\t1\t0\t1\t1\t1\n0\t1\t0\t0\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t0\t0\t1\t0\t0\t0\t1\t1\t0\t0\t1\t0\t1\t0\t1\t0\t1\t0\t0\t0\t0\t0\t0\t0\n0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\t0\t1\t0\t0\t1\t1\t0\t0\t0\t1\t0\t1\t0\t0\t1\t0\t0\t0\t0\t1\t1\n0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t0\t0\t1\t1\t1\t0\t1\t1\t1\t0\t0\t0\t0\t0\t1\t1\t0\t1\n1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t0\t1\t0\t0\t0\t1\t0\t1\t0\t1\t0\t0\t1\t0\t0\t1\t1\t1\t0\t1\n1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t1\t1\t0\t1\t1\t0\t0\t0\t1\t0\t0\t1\t0\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t0\t0\n1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t1\t1\t0\t1\t0\t0\t0\t0\t0\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t1\t1\t0\t0\t0\t1\t1\t1\t0\t1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t0\t1\t1\t0\t1\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t0\t1\t1\t0\t1\t0\t0\t1\t0\t1\t1\t0\t0\t1\t0\t0\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t0\t1\t1\t0\t0\t1\t1\t1\t0\t1\t0\t0\t0\t0\t0\t0\t1\t1\t0\t0\t1\n0\t1\t1\t0\t1\t0\t1\t0\t0\t0\t0\t1\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0\t0\t1\t0\t1\t1\t0\t0\t0\t0\t0\t0\t1\t0\t1\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t1\t0\t1\t1\n1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t0\t0\t0\t0\t1\t0\t0\t0\t0\t1\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t0\t1\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t0\t0\t0\t1\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t0\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t1\t0\t0\t1\t1\t0\t0\t1\t0\t0\t0\t0\t1\t0\t0\t1\t0\t1\t0\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t0\t1\t1\t0\t0\t1\t0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\n1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t1\t0\t1\t0\t0\t1\t1\t0\t0\t1\t1\t1\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t1\t0\t0\t1\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t0\t0\n1\t1\t0\t0\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t0\t0\t0\t0\t1\t1\t0\t1\t0\t0\t0\t0\t1\t0\t1\t0\t1\t0\t0\t1\t1\t0\t0\t1\n0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t0\t1\t1\t0\t0\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\n1\t1\t0\t0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t1\t1\t0\t0\t0\t1\t0\t1\t0\t0\t0\t1\t0\t0\t1\t0\t0\n1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\t0\t1\t1\t0\t1\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t0\n0\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t0\t0\t0\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\t1\t1\n1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t0\t1\t0\t1\t0\t1\t1\t1\t0\t1\t0\t0\t0\t0\t0\t1\t0\t0\t0\t1\t1\t1\t0\t0\t1\t1\t0\t0\t0\t0\t0\n0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t0\t0\t0\t1\t1\t1\t1\t1\t1\t0\t0\t1\t0\t1\t1\t1\t1\t0\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t0\t1\t1\t0\t0\t0\n1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t1\t0\t0\t1\t1\t1\t1\t0\t1\t0\t0\t1\n0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t0\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t1\t0\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t1\t0\t1\t1\t1\t0\t0\t1\t0\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t1\t1\t1\t1\t1\t0\t1\t0\t0\t0\t1\t0\t1\t1\t1\t0\t1\t0\t0\t1\t1\t0\t1\t1\t0\t1\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\n1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t1\t1\t1\t0\t1\t0\t0\t0\t0\t1\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t0\t1\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t0\t0\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t0\t1\t0\t1\t0\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t0\t0\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t0\n1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t1\t0\t1\t0\t0\t1\t0\t0\t0\t0\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t0\t1\t0\t1\t0\t1\t1\t1\t1\t1\t1\t1\t0\t1\t0\t1\t1\t1\t0\t0\t1\n1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t0\t1\t0\t0\t1\t0\t1\t1\t1\t1\t1\t0\t1\t0\t1\t0\t1\t0\t0\t1\t1\t0\t1\t1\n1\t1\t1\t1\t1\t0\t0\t0\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t1\t0\t0\t1\t1\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t1\t1\t1\t1\t1\t0\t0\t1\t0\t0\t0\t1\t0\t0\t0\t0\n0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t0\t1\t1\t1\t0\t0\t1\t1\t1\t0\t1\t1\t1\t0\t1\t1\t1\t1\t0\t0\t1\t0\t0\t0\t1\t1\t1\t1\t0\t0\t1\t0\t1\t0\t0\t1\t1\t0\t0\t1",
            mode="r", theme="cobalt", height="400px"),

            br(),

            h3("Basic statistics (Using the raw score)"),
            verbatimTextOutput("textarea.out"),

            br(),

            h3("Histogram of the total score (Using the raw score)"),
           downloadButton('downloaddistPlot', 'Download the plot as pdf'),
            plotOutput("distPlot"),

            br(),

            h3("Boxplot (Using the raw score)"),
          downloadButton('downloadboxPlot', 'Download the plot as pdf'),
            plotOutput("boxPlot"),

            br(),

            h3("Cronbach's coefficient alpha (Using the raw score)"),
            verbatimTextOutput("alpha.result.out"),

            p('Drop if: Cronbach alpha when the item is removed', br(),
            'r dropped: item-total correlation without the item', br(),
            'r: item-total (point-biserial) correlation', br(),
            'IF: item facility or item mean (proportion correct answers)', br(),
            'ID: item discrimination (upper 1/3 - lower 1/3)', br()
            ),

            br(),

            h3("Item estimate"),
            verbatimTextOutput("item.est.out"),

            br(),

            h3("Person estimate (Theta)"),
            verbatimTextOutput("person.est.out"),

            br(),
            br(),

            h3("Plots"),

            plotOutput("ICC"),
            plotOutput("IIC"),
            plotOutput("TIC"),

            br(),
            br(),

            strong('R session info'),
            verbatimTextOutput("info1.out")

            ),

#           tabPanel("Model options", icon = icon("line-chart", lib = "font-awesome"),
#
#                    radioButtons("cormeasures", strong("Correlation model measures"),
#                                 c("raw correlation coefficient" = "COR",
#                                   "raw correlation coefficient corrected for its slight negative bias" = "UCOR",
#                                   "Fisherr-to-z transformed correlation coefficient" = "ZCOR"
#                                 ), selected = "ZCOR"),
#                    p(h6('Fishes r-to-z transformed correlation coefficient is the default estimator for the metafor package')),
#
#                    verbatimTextOutput('cormeasures.out')
#
#           )),



        tabPanel("Polytomous Models", icon = icon("chevron-right", lib = "font-awesome"),

        br(),

        p('Note: Input values (either numeric or character) must be separated by tabs. Copy and paste from Excel/Numbers.'),

        p(HTML("<b><div style='background-color:#FADDF2;border:1px solid black;'>Please make sure that your data includes the header (variable names) in the first row.</div></b>")),

        aceEditor("text2", value="i01\ti02\ti03\ti04\n2\t1\t3\t3\n2\t1\t3\t1\n2\t1\t2\t1\n5\t5\t5\t5\n4\t1\t4\t2\n3\t3\t3\t4\n3\t2\t4\t1\n4\t3\t4\t3\n4\t2\t4\t3\n3\t2\t3\t2\n3\t3\t4\t3\n3\t1\t4\t3\n3\t4\t4\t5\n4\t3\t3\t3\n3\t2\t3\t3\n2\t2\t3\t3\n3\t3\t5\t3\n3\t3\t4\t2\n1\t1\t1\t1\n1\t1\t5\t2\n3\t3\t4\t4\n3\t1\t4\t1\n6\t6\t6\t6\n3\t1\t4\t3\n2\t2\t4\t2\n3\t3\t3\t2\n4\t3\t4\t3\n3\t3\t5\t5\n1\t1\t3\t1\n1\t1\t1\t1\n3\t3\t3\t3\n4\t3\t3\t3\n2\t1\t4\t3\n2\t2\t4\t3\n4\t5\t6\t4\n4\t2\t4\t3\n3\t3\t3\t3\n4\t3\t4\t3\n1\t1\t1\t1\n2\t3\t2\t3\n4\t3\t5\t3\n2\t2\t2\t2\n3\t3\t4\t3\n4\t4\t4\t4\n5\t3\t5\t5\n4\t3\t4\t6\n2\t2\t4\t6\n5\t3\t5\t5\n4\t4\t5\t4\n2\t1\t3\t1\n1\t3\t4\t4\n5\t3\t5\t3\n2\t2\t4\t1\n5\t4\t3\t3\n3\t2\t4\t4\n1\t1\t2\t3\n3\t1\t2\t3\n2\t3\t5\t5\n3\t3\t4\t4\n5\t3\t4\t4\n2\t2\t5\t3\n2\t1\t3\t3\n2\t2\t4\t4\n4\t3\t3\t4\n1\t2\t2\t2\n4\t2\t4\t4\n3\t3\t3\t2\n2\t1\t2\t1\n5\t2\t5\t3\n3\t3\t3\t4\n5\t2\t4\t5\n2\t1\t1\t1\n3\t2\t5\t3\n4\t3\t4\t4\n4\t2\t4\t4\n1\t1\t1\t1\n4\t4\t4\t2\n3\t1\t3\t1\n1\t1\t3\t1\n3\t1\t4\t2\n3\t1\t3\t2\n1\t1\t3\t2\n4\t2\t3\t3\n1\t3\t2\t1\n3\t3\t4\t3\n3\t1\t5\t1\n4\t3\t5\t4\n2\t2\t2\t2\n1\t1\t1\t1\n4\t1\t2\t2\n3\t4\t5\t3\n4\t3\t4\t3\n3\t2\t2\t2\n3\t2\t3\t2\n3\t3\t5\t4\n3\t1\t3\t2\n3\t3\t4\t2\n3\t2\t3\t1\n3\t4\t5\t5\n2\t2\t3\t3\n4\t4\t4\t4\n3\t4\t6\t5\n1\t1\t3\t1\n1\t1\t1\t1\n3\t3\t4\t3\n2\t2\t4\t1\n1\t1\t2\t1\n2\t1\t2\t1\n3\t3\t4\t2\n2\t2\t4\t2\n2\t1\t3\t1\n1\t2\t3\t1\n4\t3\t4\t3\n1\t2\t3\t3\n2\t2\t3\t3\n3\t4\t4\t4\n1\t1\t1\t2\n1\t1\t2\t2\n2\t1\t3\t2\n3\t3\t5\t4\n4\t4\t4\t4\n3\t3\t3\t2\n4\t3\t4\t4\n1\t3\t4\t3\n3\t3\t3\t3\n2\t1\t2\t2\n2\t4\t5\t5\n1\t1\t4\t2\n2\t2\t5\t2\n4\t3\t4\t3\n2\t2\t3\t2\n4\t3\t4\t4\n1\t1\t1\t2\n3\t4\t4\t3\n3\t3\t4\t3\n4\t3\t4\t3\n3\t3\t5\t3\n2\t1\t2\t2\n3\t2\t4\t2\n5\t3\t5\t3\n3\t1\t1\t1\n4\t3\t5\t3\n3\t2\t3\t2\n2\t1\t2\t1\n2\t2\t2\t2\n3\t3\t4\t2\n1\t1\t1\t1\n1\t1\t1\t1\n4\t3\t3\t2\n3\t2\t2\t2\n4\t2\t3\t2\n4\t4\t4\t5\n3\t3\t5\t4\n3\t4\t4\t3\n3\t3\t4\t4\n3\t3\t3\t1\n3\t3\t4\t3\n5\t4\t5\t5\n3\t2\t3\t3\n4\t3\t5\t4\n1\t1\t3\t1\n2\t2\t4\t3\n3\t3\t5\t4\n3\t3\t4\t2\n2\t1\t2\t2\n3\t2\t4\t2\n1\t1\t4\t1\n4\t4\t4\t4\n4\t3\t3\t2\n2\t2\t4\t2\n3\t3\t4\t3\n4\t3\t5\t3\n5\t3\t4\t4\n3\t3\t5\t4\n3\t3\t4\t3\n2\t4\t4\t3\n3\t2\t4\t3\n4\t4\t4\t4\n3\t2\t5\t3\n3\t1\t3\t4\n4\t3\t4\t4\n3\t2\t4\t3\n4\t3\t4\t4\n2\t2\t3\t4\n2\t1\t3\t3\n3\t3\t5\t4\n3\t2\t4\t4\n5\t4\t5\t4\n4\t2\t5\t4\n4\t2\t4\t2\n2\t2\t5\t4\n3\t3\t4\t1\n2\t2\t5\t3",
            mode="r", theme="cobalt", height="400px"),

            br(),
            br(),

            h3("Basic statistics (Using the raw score)"),
            verbatimTextOutput("bs.poly.out"),

            br(),

            h3("Histogram"),
            radioButtons("meantotal1", "",
                list("Average" = "mean1",
                     "Total" = "total1"), selected = "mean1"),
            plotOutput("distPlot2"),
        downloadButton('downloaddistPlot2', 'Download the plot as pdf'),
            br(),

            h3("Box plot with individual data points"),

            radioButtons("meantotal2", "",
                list("Average" = "mean2",
                     "Total" = "total2"), selected = "mean2"),
            plotOutput("boxPlot2"),
        downloadButton('downloadboxPlot2', 'Download the plot as pdf'),
            br(),

            h3("Cronbach's coefficient alpha (Using the raw score)"),
            verbatimTextOutput("alpha.result2.out"),

            p('Drop if: Cronbach alpha when the item is removed', br(),
            'r dropped: item-total correlation without the item', br(),
            'r: item-total (point-biserial) correlation', br()
            ),

            br(),

            plotOutput("screePlot"),

            br(),
            br(),

            h3("Item estimate"),

            radioButtons("model", strong("Model selection:"),
                list("Generalized Parial Credit Model" = "gpcm.mdl",
                    "Graded Response Model" = "gr.mdl"
                ), selected = "gpcm.mdl"
            ),

            verbatimTextOutput("polyitem.est.out"),

            br(),

            h3("Person estimate"),

            radioButtons("person", p("Select:"),
                list("No estimation of theta (person's ability)" = "no.theta",
                    "Show estimation of theta" = "show.theta"
                ),
            ),

            verbatimTextOutput("polyperson.est.out"),

            br(),
            br(),

            h3("Plots"),
            numericInput("plot.item", "Indicate which item to plot:", 1),

            plotOutput("polyICC"),


            br(),

            plotOutput("polyIIC"),

            br(),

            plotOutput("polyTIC"),

            br(),
            br(),

            strong('R session info'),
            verbatimTextOutput("info2.out")

            ),



        tabPanel("Binary (1-0) Data Converter", icon = icon("table", lib = "font-awesome"),
                 strong('Option:'),


                 checkboxInput("colname", label = "The input data includes variable names (the header) in the first row.", TRUE),

                 br(),

                 p('Note: Input values (either numeric or character) must be separated by tabs. Copy and paste from Excel/Numbers.'),

                 aceEditor("text", value="i01\ti02\ti03\ti04\ti05\ti06\ti07\ti08\ti09\ti10\ti11\ti12\ti13\ti14\ti15\ti16\ti17\ti18\ti19\ti20\nA\tB\tB\tB\tB\tC\tB\tC\tB\tD\tD\tC\tA\tB\tA\tD\tB\tD\tA\tC\nC\tD\tA\tD\tC\tB\tD\tB\tD\tA\tD\tD\tA\tB\tC\tC\tC\tA\tD\tC\nB\tD\tC\tD\tA\tB\tA\tC\tB\tD\tB\tA\tA\tD\tD\tA\tB\tC\tB\tB\nC\tC\tD\tD\tD\tA\tA\tD\tD\tD\tA\tB\tC\tB\tD\tB\tC\tB\tC\tA\nA\tA\tA\tD\tA\tA\tD\tB\tA\tC\tA\tD\tC\tC\tC\tC\tA\tA\tA\tB\nA\tA\tB\tC\tC\tA\tA\tA\tA\tA\tB\tC\tC\tC\tC\tB\tD\tC\tD\tD\nA\tA\tB\tA\tA\tA\tD\tB\tC\tC\tB\tC\tD\tA\tB\tD\tB\tB\tB\tD\nA\tC\tA\tD\tC\tA\tD\tA\tA\tA\tD\tD\tC\tC\tB\tA\tD\tC\tA\tD\nD\tB\tA\tD\tD\tA\tD\tB\tB\tA\tB\tB\tB\tC\tA\tA\tD\tA\tC\tB\nC\tC\tA\tC\tB\tC\tD\tC\tA\tA\tD\tD\tA\tA\tB\tC\tB\tB\tC\tC\nD\tA\tC\tB\tD\tA\tD\tB\tD\tA\tA\tD\tD\tC\tA\tC\tD\tC\tA\tD\nA\tD\tA\tD\tC\tA\tA\tA\tC\tD\tB\tB\tB\tA\tA\tC\tC\tD\tC\tC\nD\tC\tB\tA\tD\tA\tD\tB\tB\tA\tB\tD\tC\tC\tC\tC\tD\tA\tB\tC\nD\tC\tC\tD\tA\tA\tD\tB\tD\tB\tD\tD\tC\tB\tB\tB\tD\tB\tC\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tC\tC\tC\tB\tC\tB\tA\tA\tB\nB\tC\tB\tC\tB\tB\tD\tD\tA\tA\tB\tA\tC\tD\tD\tB\tA\tA\tB\tC\nD\tA\tA\tA\tD\tD\tB\tC\tB\tA\tB\tA\tD\tB\tD\tC\tC\tD\tD\tC\nC\tD\tB\tC\tC\tA\tD\tC\tC\tB\tC\tA\tC\tA\tA\tC\tC\tB\tA\tC\nD\tB\tA\tC\tD\tB\tD\tB\tD\tC\tD\tD\tA\tC\tA\tD\tC\tD\tD\tD\nA\tC\tA\tD\tB\tC\tC\tD\tD\tC\tA\tB\tC\tC\tA\tD\tB\tC\tA\tB\nD\tC\tA\tD\tB\tA\tD\tA\tA\tA\tA\tD\tC\tC\tC\tC\tB\tB\tD\tD\nB\tA\tB\tB\tB\tA\tD\tD\tD\tD\tD\tD\tB\tB\tB\tD\tB\tD\tC\tC\nB\tC\tC\tA\tC\tC\tD\tB\tD\tA\tC\tD\tD\tA\tA\tA\tC\tB\tD\tC\nA\tC\tA\tA\tA\tB\tD\tB\tA\tA\tC\tC\tC\tC\tA\tB\tA\tB\tB\tB\nD\tD\tD\tB\tC\tA\tD\tB\tC\tA\tD\tD\tB\tD\tB\tC\tB\tA\tB\tA\nA\tD\tA\tB\tA\tA\tB\tC\tB\tB\tA\tA\tB\tA\tC\tA\tD\tB\tD\tB\nB\tD\tD\tA\tC\tB\tD\tD\tA\tA\tB\tD\tC\tA\tA\tD\tD\tA\tD\tD\nD\tC\tD\tB\tA\tA\tB\tC\tC\tB\tC\tD\tC\tC\tB\tD\tA\tB\tC\tB\nD\tB\tD\tB\tD\tA\tD\tC\tD\tD\tC\tD\tC\tD\tD\tB\tC\tC\tB\tD\nB\tB\tD\tC\tD\tA\tB\tB\tB\tB\tB\tC\tA\tC\tC\tA\tB\tA\tB\tA\nD\tC\tC\tC\tC\tA\tD\tA\tB\tA\tC\tD\tC\tC\tB\tA\tA\tC\tB\tB\nB\tA\tA\tB\tD\tA\tA\tD\tD\tA\tA\tC\tD\tA\tD\tA\tA\tC\tB\tA\nB\tA\tA\tB\tB\tC\tB\tA\tC\tA\tA\tC\tD\tD\tB\tD\tA\tB\tB\tB\nD\tA\tC\tB\tB\tA\tC\tB\tB\tD\tB\tD\tC\tA\tA\tA\tC\tD\tD\tD\nB\tC\tB\tA\tD\tA\tD\tB\tA\tA\tD\tB\tC\tA\tC\tB\tA\tB\tC\tC\nA\tA\tC\tB\tA\tA\tC\tB\tB\tC\tD\tA\tB\tC\tA\tD\tA\tB\tB\tD\nD\tC\tA\tC\tA\tA\tD\tB\tC\tB\tA\tD\tC\tB\tD\tC\tD\tD\tC\tA\nD\tC\tC\tB\tA\tA\tD\tB\tC\tB\tD\tD\tB\tA\tD\tA\tD\tC\tD\tD\nB\tC\tB\tC\tA\tA\tD\tB\tA\tA\tA\tC\tD\tB\tD\tC\tB\tD\tC\tC\nD\tC\tA\tB\tD\tA\tD\tA\tB\tB\tC\tC\tC\tA\tA\tA\tC\tC\tA\tC\nC\tC\tD\tC\tB\tD\tA\tC\tD\tA\tC\tB\tA\tD\tD\tA\tB\tD\tA\tC\nC\tD\tB\tD\tA\tA\tC\tB\tC\tA\tB\tD\tC\tD\tC\tA\tD\tC\tB\tA\nC\tD\tC\tB\tB\tA\tD\tB\tD\tA\tD\tD\tC\tB\tD\tB\tB\tA\tA\tB\nB\tB\tA\tB\tB\tA\tD\tB\tC\tC\tD\tB\tC\tC\tD\tA\tD\tC\tD\tC\nD\tB\tD\tC\tC\tA\tD\tD\tB\tA\tC\tA\tC\tC\tA\tB\tC\tC\tA\tA\nD\tA\tD\tA\tA\tA\tD\tB\tD\tB\tB\tB\tA\tB\tB\tC\tC\tB\tB\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tB\tD\tC\tC\tB\tA\tD\tA\tB\tB\nD\tC\tA\tD\tB\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tA\tC\tB\nD\tC\tD\tA\tB\tA\tD\tB\tD\tD\tA\tD\tC\tD\tB\tC\tD\tA\tA\tB\nA\tB\tA\tB\tA\tA\tD\tB\tA\tC\tC\tB\tC\tD\tB\tD\tC\tA\tD\tA\nB\tB\tA\tB\tC\tA\tD\tA\tB\tA\tC\tD\tC\tC\tC\tB\tB\tA\tD\tD\nD\tC\tC\tD\tD\tA\tD\tB\tD\tA\tC\tD\tC\tD\tD\tA\tB\tB\tA\tC\nD\tC\tD\tD\tD\tA\tD\tC\tB\tA\tA\tD\tC\tD\tB\tD\tD\tC\tB\tA\nA\tB\tA\tD\tC\tA\tB\tB\tC\tC\tC\tB\tD\tA\tA\tA\tB\tB\tC\tD\nB\tC\tA\tC\tA\tA\tD\tA\tB\tC\tD\tA\tC\tB\tC\tB\tC\tC\tB\tC\nC\tA\tB\tD\tC\tA\tD\tA\tA\tD\tC\tB\tA\tA\tB\tB\tB\tA\tC\tD\nA\tD\tA\tB\tC\tB\tD\tB\tB\tC\tB\tA\tC\tC\tA\tA\tA\tC\tC\tD\nB\tA\tA\tB\tD\tA\tD\tC\tD\tA\tA\tD\tC\tC\tD\tA\tD\tA\tA\tB\nD\tC\tA\tB\tC\tA\tD\tD\tA\tA\tD\tA\tC\tB\tC\tB\tA\tA\tA\tC\nA\tB\tD\tC\tC\tA\tD\tA\tD\tA\tD\tD\tC\tC\tA\tB\tC\tB\tB\tD\nD\tC\tA\tD\tD\tA\tD\tB\tC\tA\tA\tD\tA\tC\tB\tC\tD\tA\tD\tB\nB\tA\tA\tC\tC\tC\tD\tB\tB\tC\tA\tA\tA\tD\tB\tB\tD\tD\tD\tC\nD\tC\tA\tB\tB\tA\tD\tA\tB\tB\tA\tD\tC\tC\tA\tC\tD\tB\tA\tB\nD\tC\tA\tB\tC\tA\tD\tB\tA\tA\tA\tD\tC\tC\tB\tC\tC\tA\tD\tB\nC\tB\tA\tC\tC\tA\tB\tA\tA\tD\tB\tA\tC\tA\tA\tA\tA\tB\tC\tB\nC\tC\tD\tA\tD\tC\tB\tB\tA\tB\tC\tD\tC\tB\tC\tC\tD\tA\tA\tA\nB\tA\tB\tA\tB\tA\tA\tA\tC\tC\tC\tB\tD\tD\tB\tB\tB\tB\tA\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tB\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tC\tD\tA\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tB\tD\tA\tA\tB\nA\tC\tB\tA\tD\tA\tA\tB\tB\tB\tA\tD\tC\tB\tB\tA\tD\tD\tA\tD\nA\tC\tA\tB\tC\tB\tD\tB\tC\tD\tA\tD\tC\tC\tA\tD\tD\tD\tB\tA\nD\tA\tA\tA\tA\tA\tD\tB\tB\tB\tA\tA\tA\tC\tB\tC\tD\tD\tA\tB\nA\tC\tC\tD\tB\tA\tC\tD\tA\tA\tA\tA\tC\tC\tB\tA\tD\tC\tB\tB\nD\tC\tA\tD\tD\tA\tD\tA\tA\tD\tD\tD\tC\tC\tC\tD\tA\tA\tA\tB\nD\tC\tD\tB\tD\tA\tC\tA\tC\tA\tA\tD\tC\tC\tD\tD\tC\tC\tB\tD\nD\tC\tA\tD\tA\tA\tD\tB\tA\tB\tA\tD\tD\tC\tB\tB\tD\tD\tD\tC\nB\tC\tA\tB\tC\tA\tD\tB\tA\tD\tA\tD\tC\tC\tA\tC\tB\tB\tA\tB\nD\tC\tB\tD\tB\tA\tD\tC\tD\tA\tB\tD\tC\tC\tB\tC\tD\tD\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tC\tD\tC\tC\tD\tA\tD\tA\tB\tD\nD\tC\tA\tB\tD\tA\tD\tB\tD\tB\tA\tD\tC\tC\tB\tD\tA\tB\tA\tD\nC\tC\tA\tC\tD\tA\tD\tB\tD\tB\tA\tD\tA\tC\tD\tC\tC\tD\tA\tB\nB\tC\tA\tD\tC\tB\tB\tD\tB\tA\tB\tA\tD\tC\tD\tC\tA\tD\tD\tB\nC\tA\tD\tA\tA\tA\tB\tC\tA\tC\tB\tC\tC\tD\tB\tB\tC\tC\tD\tB\nB\tB\tA\tC\tA\tA\tD\tB\tC\tC\tD\tA\tB\tA\tC\tD\tD\tB\tC\tA\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tD\tD\tA\tA\tC\nC\tC\tA\tB\tC\tA\tD\tB\tD\tD\tA\tD\tA\tC\tD\tA\tC\tC\tB\tB\nA\tC\tB\tD\tD\tA\tD\tB\tA\tA\tA\tD\tC\tC\tD\tC\tA\tB\tA\tD\nD\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tA\tD\tA\tA\tB\nD\tC\tA\tD\tD\tA\tD\tB\tB\tD\tA\tD\tA\tC\tC\tC\tD\tA\tC\tB\nD\tB\tA\tB\tD\tA\tA\tC\tC\tD\tC\tC\tB\tD\tC\tC\tB\tB\tD\tB\nD\tA\tD\tB\tA\tA\tA\tC\tB\tA\tC\tA\tC\tD\tB\tC\tB\tC\tD\tD\nB\tD\tA\tB\tB\tA\tC\tD\tC\tB\tB\tD\tD\tA\tB\tA\tB\tD\tC\tC\nD\tA\tB\tB\tC\tA\tD\tB\tC\tA\tA\tD\tC\tC\tB\tC\tD\tC\tD\tB\nA\tA\tB\tB\tD\tA\tA\tD\tB\tA\tB\tD\tC\tC\tD\tB\tB\tA\tB\tD\nD\tC\tA\tD\tD\tA\tD\tB\tB\tC\tB\tD\tA\tA\tC\tB\tD\tA\tD\tB\nD\tC\tC\tD\tD\tA\tD\tB\tC\tA\tA\tD\tC\tC\tB\tB\tD\tD\tB\tB\nD\tC\tA\tD\tB\tA\tD\tB\tA\tA\tA\tD\tC\tC\tB\tB\tD\tA\tA\tB\nB\tB\tA\tC\tD\tA\tD\tB\tC\tA\tD\tB\tC\tC\tD\tC\tD\tC\tB\tD\nD\tA\tA\tD\tD\tA\tD\tB\tB\tA\tA\tD\tB\tC\tC\tC\tD\tA\tA\tB",
                           mode="r", theme="cobalt", height="400px"),

                 p("Input answer keys (Either numeric or character, separated by tabs.):"),

                 aceEditor("anskey", value="D\tC\tA\tD\tD\tA\tD\tB\tD\tA\tA\tD\tC\tC\tB\tC\tD\tA\tA\tB", mode="r", theme="chrome", height="50px"),

                 br(),

                 h3("Checking the 1-0 converted data"),
                 p('Only the first 10 observations are displayed.'),
                 tableOutput("check"),

                 br(),

                 downloadButton('downloadData', 'Download the 1-0 data file'),

                 br(),
                 strong('R session info'),
                 verbatimTextOutput("info3.out"),

                 p(br())

        ),
# tabPanel("Model Options", icon = icon("cog", lib = "font-awesome"),
#          strong('Option:'),
#
# #           radioButtons("optimmethod", strong("Optimization Method"),
# #                        c("Broyden-Fletcher-Goldfarb-Shanno (BFGS)" = "BFGS",
# #                          "Conjugate Gradients" = "CG",
# #                          "Limited-memory Broyden-Fletcher-Goldfarb-Shanno (BFGS)" = "L-BFGS-B"
# #                        ), selected = "BFGS"),
# #           p(h6('Broyden-Fletcher-Goldfarb-Shanno (BFGS) is the default method for the ltm package')),
# #
# #           verbatimTextOutput('optimmethod.out'),
#
#          br(),
#          p(br())
#
# ),

        tabPanel("About", icon = icon("dot-circle-o", lib = "font-awesome"),
                 strong('IRT with Shiny'),
                 p("IRTShiny Version 1.1"),
                 p("The goal of this project is to help students and researchers run IRT analysis as easily as possible."),
                 p('This application is developed with',
                   a("Shiny.", href="http://www.rstudio.com/shiny/", target="_blank"),
                   ''),
                 p('The code for this application is available at this',
                   a('GitHub.', href='https://github.com/kylehamilton/IRTShiny', target="_blank")),
                 p("If you're having problems with IRTShiny feel free to refer to our GitHub wiki or the documentation available on CRAN."),
                 a("CRAN page for IRTShiny", href="http://cran.r-project.org/web/packages/IRTShiny/index.html", target="_blank"),
                 br(),
                 a("GitHub Wiki page for IRTShiny", href="https://github.com/kylehamilton/IRTShiny/wiki", target="_blank"),
                 br(),
                 p("As always you are more than welcome to contact the project maintainer at kyle.hamilton@gmail.com"),



                 br(),

                 strong('List of Packages Used'), br(),
                 code('library(shiny)'),br(),
                 code('library(shinyAce)'),br(),
                 code('library(psych)'),br(),
                 code('library(CTT)'),br(),
                 code('library(ltm)'),br(),
                 code('library(beeswarm)'),br(),


                 br(),

                 h4('Acknowledgments and Authors'),

                 strong('Acknowledgments'),

                 p('William Kyle Hamilton would like to thank the ',
                   a("Health Communications and Interventions Lab at UC Merced", href="http://cameronhcilab.com/", target="_blank"),
                   'for their comments and beta testing efforts on this application'),

                 p('Atsushi Mizumoto would like to thank',
                   a("Dr. Luke Plonsky", href="http://oak.ucc.nau.edu/ldp3/", target="_blank"), 'and',
                   a("Dr. Yo In'nami", href="https://sites.google.com/site/yoinnami/", target="_blank"),
                   'for their support and feedback to create this web application.'),

                 br(),



                 h5('Authors'),

                 HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/kyle80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                 p(a("William Kyle Hamilton - University of California, Merced", href="http://www.kylehamilton.com", target="_blank")),
                 p("William Kyle Hamilton maintains this application and has authored new features."),

                 br(),
                 HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/atsushi80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                 p(a("Atsushi Mizumoto, PhD - Kansai University", href="http://mizumot.com", target="_blank"),br(),
                   p("Atsushi Mizumoto wrote the first version of this application; this application is a fork of the original which can be found", a("here", href="https://github.com/mizumot/irt", target="_blank"))


                 ),
                 br(),
                 strong('Bug Reports'),

                 p("If you discover a problem with IRTShiny please submit it to the project GitHub page",
                   a("https://github.com/kylehamilton/IRTShiny/issues", href="https://github.com/kylehamilton/IRTShiny/issues", target="_blank"),br()),

                 p("IRTShiny is an Open Source project, you are more than welcome to submit patches or features and help the project grow."),
                 strong('License'),

                 p("IRTShiny: Item Response Theory with Shiny"),
                 p(" Copyright 2015  William Kyle Hamilton and Atsushi Mizumoto"),

                 p(" This program is free software you can redistribute it and or modify
           it under the terms of the GNU General Public License as published by
           the Free Software Foundation either version 3 of the License or
           at your option any later version."),

                 p("This program is distributed in the hope that it will be useful,
           but WITHOUT ANY WARRANTY; without even the implied warranty of
           MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
           GNU General Public License for more details."),

                 p("You should have received a copy of the GNU General Public License
           along with this program.  If not, see", a("http://www.gnu.org/licenses/gpl.html", href="http://www.gnu.org/licenses/gpl.html", target="_blank"),br()),
                 img(src = "http://www.gnu.org/graphics/gplv3-127x51.png", seamless=NA),


                 br(),

                 strong('Futher Infomation'),
                 p("If you would like to learn more about the GNU General Public License and what it means tl'dr legal has a simple explaination which can be found here", a("https://www.tldrlegal.com/l/gpl-3.0", href="https://www.tldrlegal.com/l/gpl-3.0", target="_blank"),br()),


                 p(br())

        )

        # ))

)
)
