source("helper.R")

shinyUI(fluidPage(
    #theme="bootstrap.css",
    withMathJax(),
    titlePanel("Estimate Rayleigh \\(\\sigma\\) parameter from range statistics"),
    sidebarLayout(
        #####-------------------------------------------------------------------
        ## sidebar
        #####-------------------------------------------------------------------
        sidebarPanel(width=3,
           #####----------------------------------------------------------------
           ## Rayleigh sigma from range statistics
           conditionalPanel(condition="/Range statistics/.test(input.task)",
                radioButtons("rangeStatType1", label=h5("Range statistic"),
                             rangeStat, selected="1"),
                textInput("rangeStatStats", h5("Measured values"),
                          value=c("1 2 3")),
                numericInput("rangeStatN", h5("Shots per group"),
                             min=2, max=100, step=1, value=5),
                numericInput("rangeStatNGroups", h5("Number of groups"),
                             min=1, max=10, step=1, value=1),
                sliderInput("rangeStatCILevel", label=h5("CI level"),
                            min=0, max=1, value=0.9, step=0.05),
                numericInput("dstTrgt1", h5("Distance to target"),
                             min=0, step=1, value=100),
                selectInput("unitDstTrgt1", h5("Measurement unit distance"),
                            choices=unitsDst, selected=2),
                selectInput("unitXY1", h5("Measurement unit coordinates"),
                            choices=unitsXY, selected=3)
            ),

           #####----------------------------------------------------------------
           ## efficiency - required number of groups
           conditionalPanel(condition="/number/.test(input.task)",
                radioButtons("effStatType1", label=h5("Measured group statistic"),
                             rangeStatSig, selected="1"),
                textInput("effN1", h5("Shots per group"),
                          value=c("3 5 10 20")),
                sliderInput("effCILevel1", label=h5("CI level"),
                            min=0, max=1, value=0.9, step=0.05),
                sliderInput("effCIWidth1", label=h5("CI width (=2*E)"),
                            min=0, max=1, value=0.2, step=0.05)
            ),

           #####----------------------------------------------------------------
           ## efficiency - required number of groups
           conditionalPanel(condition="/CI width/.test(input.task)",
                radioButtons("effStatType2", label=h5("Range statistic"),
                             rangeStatSig, selected="1"),
                #numericInput("effN2", h5("Shots per group"),
                #             min=0, max=100, step=1, value=5),
                textInput("effN2", h5("Shots per group"),
                          value=c("3 5 10 20")),
                numericInput("effNGroups2", h5("Number of groups"),
                             min=1, step=1, value=1),
                sliderInput("effCILevel2", label=h5("CI level"),
                            min=0, max=1, value=0.9, step=0.05)
            ),

            #####---------------------------------------------------------------
            ## about
            conditionalPanel(condition="input.task == 'About'",
                h4("Background information")
            )
        ),

        #####-------------------------------------------------------------------
        ## main output area
        #####-------------------------------------------------------------------
        mainPanel(
            #tags$head(tags$style(type="text/css", ".container-fluid { max-width: 12600px; }")),
            tabsetPanel(
                #####-----------------------------------------------------------
                ## Range statistics -> Rayleigh sigma
                tabPanel(div("Range statistics", icon("arrow-right", lib="glyphicon"),"Rayleigh \\(\\sigma\\)"),
                    h6("Background information"),
                    p("Assuming a circular bivariate normal shot distribution, this
                      panel estimates the Rayleigh \\(\\sigma\\) parameter from
                      measured range statistics in a given number of shots per
                      group, averaged over a given number of groups.", br(),
                      "Distance to target, and information on the measurement
                      unit for distance and range statistic is only used for
                      the conversion to",
                      a("angular size", href="http://dwoll.shinyapps.io/shotGroupsAngular"),
                      "."),
                    verbatimTextOutput("range2sigma"),
                    p("Based on a Monte Carlo simulation with 2 million repetions of
                      2, ..., 100 shots per group in 1, ..., 10 groups.", br(),
                      "For more information, see the",
                      a("Ballistipedia entry on range statistics",
                        href="http://ballistipedia.com/index.php?title=Range_Statistics"),
                      ".")
                ),

                #####-----------------------------------------------------------
                ## Efficiency - required number of groups
                tabPanel("Efficiency: number of groups",
                    h6("Background information"),
                    p("Assuming a circular bivariate normal shot distribution, this
                      panel calculates the number of groups required to achieve
                      the desired coverage (level) for the confidence interval (CI)
                      for the Rayleigh \\(\\sigma\\) parameter, given the number
                      of shots per group, the desired CI width, and the type of
                      measured group statistic.", br(),
                      "The CI width is 2*E, where E is
                      the width as a fraction of the mean on either side."),
                    verbatimTextOutput("effNGroups"),
                    HTML("<ul>
                          <li>n - given number of shots per group</li>
                          <li>nGroupsReq - required number of groups with n shots each (including fractional groups)</li>
                          <li>nGroupsReqCeil - required number of groups with n shots each (full groups only)</li>
                          <li>nShotsReq - required total number of shots, assuming we're shooting groups of size n each (including fractional groups)</li>
                          <li>nShotsReqCeil - required total number of shots, assuming we're shooting groups of size n each (full groups only)</li>
                          <li>CIlevel - desired CI level (coverage probability)</li>
                          <li>CIwidth - desired CI width (as a fraction of the mean)</li>
                          </ul>"),

                    p("Based on a Monte Carlo simulation with 2 million repetions of
                      2, ..., 100 shots per group in 1, ..., 10 groups.", br(),
                      "For more information, see the",
                      a("Ballistipedia entry on range statistics",
                        href="http://ballistipedia.com/index.php?title=Range_Statistics"),
                      ".")
                ),

                #####-----------------------------------------------------------
                ## Efficiency - achievable CI width
                tabPanel("Efficiency: CI width",
                    h6("Background information"),
                    p("Assuming a circular bivariate normal shot distribution, this
                      panel calculates the width of the confidence interval (CI)
                      for the Rayleigh \\(\\sigma\\) parameter, given the number of
                      shots per group, the number of groups, the desired coverage
                      (level) for the CI, and the type of measured group statistic.", br(),
                      "The CI width is 2*E, where E is the width as a fraction of
                      the mean on either side."),
                    verbatimTextOutput("effCIWidth"),
                    HTML("<ul>
                          <li>n - given number of shots per group</li>
                          <li>nGroups - given number of groups with n shots each</li>
                          <li>CIlevel - desired CI level (coverage probability)</li>
                          <li>CIwidth - achievable CI width (as a fraction of the mean)</li>
                          </ul>"),
                    p("Based on a Monte Carlo simulation with 2 million repetions of
                      2, ..., 100 shots per group in 1, ..., 10 groups.", br(),
                      "For more information, see the",
                      a("Ballistipedia entry on range statistics",
                        href="http://ballistipedia.com/index.php?title=Range_Statistics"),
                      ".")
                ),

                #####-----------------------------------------------------------
                ## about
                tabPanel("About",
                    h6("About shotGroups"),
                    p("The", a("shotGroups", href="http://CRAN.R-project.org/package=shotGroups"),
                      "package for", a("R", href="http://www.r-project.org/"),
                      "provides functions to read in, plot,
                      statistically describe, analyze, and compare shooting data with respect
                      to group shape, precision, and accuracy. This includes graphical methods,
                      descriptive statistics, and inference tests using standard, but also
                      non-parametric and robust statistical methods. The data can be imported
                      from files produced by", a("OnTarget PC and OnTarget TDS",
                                                 href="http://ontargetshooting.com/tds/"), ", ",
                      a("TARAN", href="http://taran.ptosis.ch/"),
                      "or from custom data files in text format with a similar structure.
                      For further explanations and an example walkthrough, see the",
                      a("package vignette",
                        href="http://cran.rstudio.com/web/packages/shotGroups/vignettes/shotGroups.pdf"),
                      "."),
                    p("shotGroups and this web application are written by:", br(),
                      "Daniel", HTML("Wollschl&auml;ger"),
                      a("<dwoll@kuci.org>", href="mailto:dwoll@kuci.org"), br(),
                      "Source code shotGroups:",
                      a("http://github.com/dwoll/shotGroups/",
                        href="http://github.com/dwoll/shotGroups/"), br(),
                      "Source code for this Shiny app:",
                      a("http://github.com/dwoll/shotGroupsRangeStat/",
                        href="http://github.com/dwoll/shotGroupsRangeStat/")),

                    h6("More shotGroups web applications"),
                    p("Comprehensive shot group analysis:",
                      a("http://dwoll.shinyapps.io/shotGroupsApp/",
                        href="http://dwoll.shinyapps.io/shotGroupsApp/"), br(),
                      "Absolute", icon("resize-horizontal", lib="glyphicon"),
                      "angular size conversion:",
                      a("http://dwoll.shinyapps.io/shotGroupsAngular/",
                        href="http://dwoll.shinyapps.io/shotGroupsAngular/"), br(),
                      "Region", icon("resize-horizontal", lib="glyphicon"),
                      "hit probability calculations:",
                      a("http://dwoll.shinyapps.io/shotGroupsHitProb/",
                        href="http://dwoll.shinyapps.io/shotGroupsHitProb/")),

                    h6("Acknowledgements"),
                    p("Thanks to David Bookstaber for testing, feedback and data.")
                ),

                id="task"
            )
        )
    )
))
