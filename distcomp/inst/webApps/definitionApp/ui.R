shinyUI(fluidPage(

  ## tags$head(tags$script("
  ##       window.onload = function() {
  ##           $('#mynavlist a:contains(\"Data Check\")').parent().addClass('disabled');
  ##           $('#mynavlist a:contains(\"Dry Run\")').parent().addClass('disabled');
  ##           $('#mynavlist a:contains(\"Output\")').parent().addClass('disabled');
  ##       };

  ##       Shiny.addCustomMessageHandler('activeNavs', function(nav_label) {
  ##           $('#mynavlist a:contains(\"' + nav_label + '\")').parent().removeClass('disabled');
  ##       });
  ##  ")),
  titlePanel(h3("Propose a new Distributed Computation")),
  navlistPanel(selected="Propose Computation", id='mynavlist',
               tabPanel("Propose Computation", icon=img(src="checkmark.png"),

                        textInput("nameIn", label = h5("Project name (50 characters max.)"),
                                  value = "Enter text..."),

                        textInput("descIn", label = h5("Description (250 characters max.)"),
                                  value = "Enter text..."),

                        br(),
                        actionButton("printProjectSummary", "Project Summary"),
                        h5(textOutput('projectSummary')),

                        textOutput('name'),

                        textOutput('desc'),

                        selectInput("compType",
                                    label = h5("Type of computation"),
                                    choices = lapply(availableComputations(),
                                      function(x) x$desc),
                                    selected = availableComputations()[1]$desc,
                                    width="250px"),
                        br(),
                        conditionalPanel(
                          condition = "input.gotoDataInputs != 0",
                          textOutput('compType')
                        ),
                        br(),
                        actionButton("gotoDataInputs", "Continue")

                        ),
               "-----",
               tabPanel("Help-FAQ"
                      , h5("here in help")
                        )

               )
))


