source("helper.R")

shinyUI(fluidPage(
    #theme="bootstrap.css",
    withMathJax(),
    titlePanel(div("Angular", icon("resize-horizontal", lib="glyphicon"),
                   "absolute size conversion using shotGroups")),
    sidebarLayout(
        #####-------------------------------------------------------------------
        ## sidebar
        #####-------------------------------------------------------------------
        sidebarPanel(width=3,
            #####---------------------------------------------------------------
            ## angular size
            conditionalPanel(condition="input.task == 'Math'",
                h4("Background math"),
                img(src="AnglesCircle.png", width=150),
                p("Angle \\(\\varphi\\) (in degree) with corresponding arc length
                  \\(x\\) (in radian) in the unit circle.")
            ),

            conditionalPanel(condition="/angular/.test(input.task)",
                numericInput("dstTrgt1", h5("Distance to target"),
                             min=0, step=1, value=100),
                selectInput("unitDst1", h5("Measurement unit distance"),
                            choices=unitsDst, selected=2),
                textInput("angszeAbs1", h5("Absolute size"),
                          value=c("1 2 3")),
                selectInput("angszeUnitAbs1", h5("Unit absolute size"),
                            choices=unitsAbs, selected=6),
                h4("Output"),
                selectizeInput("angszeUnitAngOut1", h5("Unit angular size"),
                               choices=unitsAng, multiple=TRUE, selected=1:6)
            ),

            conditionalPanel(condition="/absolute/.test(input.task)",
                             #"input.task == 'Angular -> absolute size'",
                numericInput("dstTrgt2", h5("Distance to target"),
                             min=0, step=1, value=100),
                selectInput("unitDst2", h5("Measurement unit distance"),
                            choices=unitsDst, selected=2),
                textInput("angszeAng2", h5("Angular size"),
                          value=c("1 2 3")),
                selectInput("angszeUnitAng2", h5("Unit angular size"),
                            choices=unitsAng, selected=2),
                h4("Output"),
                selectizeInput("angszeUnitAbsOut2", h5("Unit absolute size"),
                               choices=unitsAbs, multiple=TRUE, selected=c(2, 6))
            ),

            conditionalPanel(condition="/\\+/.test(input.task)",
                             #"input.task == 'Abs+ang size -> distance'",
                textInput("angszeAbs3", h5("Absolute size"),
                      value=c("1 2 3")),
                selectInput("angszeUnitAbs3", h5("Unit absolute size"),
                            choices=unitsAbs, selected=1),
                textInput("angszeAng3", h5("Angular size"), value="1"),
                selectInput("angszeUnitAng3", h5("Unit angular size"),
                            choices=unitsAng, selected=2),
                h4("Output"),
                selectizeInput("angszeUnitDstOut3", h5("Unit distance"),
                               choices=unitsAbs, multiple=TRUE, selected=c(1, 4))
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
            #####---------------------------------------------------------------
            tabsetPanel(
                tabPanel("Math",
                    h6("Background math"),
                    p("In addition to absolute length units, group size is often reported in terms
                      of its angular diameter. Angles can be measured equivalently either in degree
                      or in radian. If \\(x\\) is the angular measurement in radian, and \\(\\varphi\\)
                      the angular measurement in degree for the same angle, then
                      \\(\\frac{x}{2 \\pi} = \\frac{\\varphi}{360}\\) such that conversion between
                      degree and radian is given by \\(x = \\frac{2 \\pi}{360} \\cdot \\varphi\\) and
                      \\(\\varphi = \\frac{360}{2 \\pi} \\cdot x\\)."),
                    p("The angular size of an object with absolute size \\(s\\)
                      is its angular diameter at a given distance \\(d\\). This
                      is the angle \\(\\alpha\\) subtended by the object with the
                      line of sight centered on it."),
                    img(src="AnglesCenterH.jpg", width=300),
                    p("Angular diameter of object with absolute size \\(s\\) at
                      distance to target \\(d\\). Right triangle formed by \\(d\\)
                      and object of size \\(s/2\\). \\(s\\) corresponds to angle
                      \\(\\alpha\\) (degree) and arc length \\(x\\) (radian)."),
                    p("The following measures for angular size are supported:"),
                    HTML("<ul>
                          <li>deg = degree. The circle is divided into 360 degrees.</li>
                          <li>MOA = minute of angle = arcmin. 1 MOA = \\(1/60\\) degree
                              such that the circle has \\(360 \\cdot 60 = 21600\\) MOA.</li>
                          <li>SMOA = Shooter's MOA = Inches Per Hundred Yards IPHY.
                              1 inch at 100 yards = 1 SMOA.</li>
                          <li>rad = radian. 1 radian is 1 unit of arc length on
                              the unit circle which has a circumference of \\(2 \\pi\\).
                              The circle circumference is divided into \\(2 \\pi\\) rad.</li>
                          <li>mrad = milliradian = 1/1000 radian. The circle circumference is
                              divided into \\(2 \\pi \\cdot 1000 \\approx 6283.19\\) mrad.</li>
                          <li>mil: 1 mil = \\(\\frac{2 \\pi}{6400}\\) radian -
                              the circle circumference is divided into 6400 mils.</li>
                          </ul>")
                ),

                #####-----------------------------------------------------------
                ## angular size
                tabPanel(div("Absolute", icon("arrow-right", lib="glyphicon"), "angular"),
                    h6("Convert absolute to angular size"),
                    verbatimTextOutput("abs2Ang"),
                    h6("Math"),
                    p("The angle \\(\\alpha\\) subtended by an object of size
                      \\(s\\) at distance \\(d\\) can be calculated from the
                      right triangle with hypotenuse \\(d\\) and cathetus \\(s/2\\):
                      \\(\\tan\\left(\\frac{\\alpha}{2}\\right) = \\frac{s}{2}
                      \\cdot \\frac{1}{d}\\), therefore \\(\\alpha = 2 \\cdot
                      \\arctan\\left(\\frac{s}{2 d}\\right)\\)."),
                   p("Assuming that the argument for \\(\\tan(\\cdot)\\) and the
                     result from \\(\\arctan(\\cdot)\\) are in radian, and that
                     distance to target \\(d\\) and object size \\(s\\) are measured
                     in the same unit, this leads to the following formulas for
                     calculating \\(\\alpha\\) in MOA, SMOA as well as \\(x\\) in
                     mrad and NATO mil based on \\(d\\) and \\(s\\)"),

HTML("<ul>
<li>Angle \\(\\alpha\\) in degree: \\(\\alpha = \\frac{360}{2 \\pi} \\cdot 2
    \\cdot \\arctan\\left(\\frac{s}{2 d}\\right) =
    \\frac{360}{\\pi} \\cdot \\arctan\\left(\\frac{s}{2 d}\\right)\\)</li>
<li>Angle \\(\\alpha\\) in MOA: \\(\\alpha = 60 \\cdot \\frac{360}{\\pi}
    \\cdot \\arctan\\left(\\frac{s}{2 d}\\right) = \\frac{21600}{\\pi} \\cdot
    \\arctan\\left(\\frac{s}{2 d}\\right)\\)</li>
<li>Angle \\(\\alpha\\) in SMOA: By definition, size \\(s=1\\) inch at
    \\(d = 100\\) yards (\\(= 3600\\) inch) is 1 SMOA.<br />
    Conversion factors to/from MOA are \\(\\frac{21600}{\\pi} \\cdot
    \\arctan\\left(\\frac{1}{7200}\\right) \\approx 0.95493\\),
    and \\(\\frac{\\pi}{21600} \\cdot
    \\frac{1}{\\arctan(1/7200)} \\approx 1.04720\\).<br />
    \\(\\alpha = \\frac{\\pi}{21600} \\cdot \\frac{1}{\\arctan(1/7200)} \\cdot
    \\frac{21600}{\\pi} \\cdot \\arctan\\left(\\frac{s}{2 d}\\right) =
    \\frac{1}{\\arctan(1/7200)} \\cdot \\arctan\\left(\\frac{s}{2 d}\\right)\\)</li>
<li>Arc length \\(x\\) in rad:  \\(x = 2 \\cdot \\arctan\\left(\\frac{s}{2 d}\\right)\\).</li>
<li>Arc length \\(x\\) in mrad: \\(x = 2000 \\cdot \\arctan\\left(\\frac{s}{2 d}\\right)\\).<br />
    Conversion factors to/from MOA are \\(\\frac{21600}{2000 \\pi} \\approx 3.43775\\)
    and \\(\\frac{2000 \\pi}{21600} \\approx 0.29089\\).</li>
<li>Arc length \\(x\\) in NATO mil: \\(x = \\frac{6400}{\\pi} \\cdot
    \\arctan\\left(\\frac{s}{2 d}\\right)\\).<br />
    Conversion factors to/from MOA are \\(\\frac{21600}{6400} = 3.375\\)
    and \\(\\frac{6400}{21600} \\approx 0.2962963\\).
</ul>")
                ),

                ## angular size
                tabPanel(div("Angular", icon("arrow-right", lib="glyphicon"), "absolute"),
                    h6("Convert angular diameter to absolute size"),
                    verbatimTextOutput("ang2Abs"),
                    h6("Math"),
                    p("Absolute object size \\(s\\) can be calculated from angular
                      diameter and distance to target \\(d\\), assuming that the
                      argument for \\(\\tan(\\cdot)\\) and the
                      result from \\(\\arctan(\\cdot)\\) are in radian, and that
                      distance to target \\(d\\) and object size \\(s\\) are measured
                      in the same unit:"),
HTML("<ul>
<li>From angle \\(\\alpha\\) in degree: \\(s = 2 \\cdot d \\cdot
    \\tan\\left(\\alpha \\cdot \\frac{\\pi}{360}\\right)\\)</li>
<li>From angle \\(\\alpha\\) in MOA: \\(s = 2 \\cdot d \\cdot
    \\tan\\left(\\alpha \\cdot \\frac{\\pi}{60 \\cdot 360} \\right) = 2 \\cdot d \\cdot
    \\tan\\left(\\alpha \\cdot \\frac{\\pi}{21600}\\right)\\)</li>
<li>From angle \\(\\alpha\\) in SMOA: \\(s = \\frac{21600}{\\pi} \\cdot
    \\arctan\\left(\\frac{1}{7200}\\right) \\cdot 2 \\cdot d \\cdot
    \\tan\\left(\\alpha \\cdot \\frac{\\pi}{21600}\\right)\\)</li>
<li>From arc length \\(x\\) in rad: \\(s = 2 \\cdot d \\cdot \\tan\\left(x
    \\cdot \\frac{1}{2}\\right)\\)</li>
<li>From arc length \\(x\\) in mrad: \\(s = 2 \\cdot d \\cdot \\tan\\left(x
    \\cdot \\frac{1}{2000}\\right)\\)</li>
<li>From arc length \\(x\\) in NATO mil: \\(s = 2 \\cdot d \\cdot \\tan\\left(x
    \\cdot \\frac{\\pi}{6400}\\right)\\)</li>
</ul>")
                ),

                #####-----------------------------------------------------------
                ## angular size
                tabPanel(div("Absolute+Angular", icon("arrow-right", lib="glyphicon"), "distance"),
                    h6("Calculate distance from absolute and angular size"),
                    verbatimTextOutput("angAbs2Dist"),
                    h6("Math"),
                    p("Distance to target \\(d\\) can be calculated from absolute
                      object size \\(s\\) and angular size, assuming that the
                      argument for \\(\\tan(\\cdot)\\) and the
                      result from \\(\\arctan(\\cdot)\\) are in radian, and that
                      distance to target \\(d\\) and object size \\(s\\) are measured
                      in the same unit:"),
HTML("<ul>
<li>From angle \\(\\alpha\\) in degree: \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(\\alpha \\cdot \\pi/360)}\\)</li>
<li>From angle \\(\\alpha\\) in MOA: \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(\\alpha \\cdot \\pi/21600)}\\)</li>
<li>From angle \\(\\alpha\\) in SMOA: \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(\\alpha \\cdot \\arctan(1/7200))}\\)</li>
<li>From arc length \\(x\\) in rad:  \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(x / 2)}\\)</li>
<li>From arc length \\(x\\) in mrad: \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(x / 2000)}\\)</li>
<li>From arc length \\(x\\) in NATO mil: \\(d = \\frac{s}{2} \\cdot
    \\frac{1}{\\tan(x \\cdot \\pi / 6400)}\\)
</ul>")
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
                      "Source code shiny app:",
                      a("http://github.com/dwoll/shotGroupsAngular/",
                        href="http://github.com/dwoll/shotGroupsAngular/")),

                    h6("More shotGroups web applications"),
                    p("Comprehensive shot group analysis:",
                      a("http://dwoll.shinyapps.io/shotGroupsApp/",
                        href="http://dwoll.shinyapps.io/shotGroupsApp/"), br(),
                      "Region", icon("resize-horizontal", lib="glyphicon"),
                      "hit probability calculations:",
                      a("http://dwoll.shinyapps.io/shotGroupsHitProb/",
                        href="http://dwoll.shinyapps.io/shotGroupsHitProb/"), br(),
                      "Estimate Rayleigh \\(\\sigma\\) from range statistics:",
                      a("http://dwoll.shinyapps.io/shotGroupsRangeStat/",
                        href="http://dwoll.shinyapps.io/shotGroupsRangeStat/")),

                    h6("Calculations"),
                    p("For details of the calculations used in this app, see the documentation for",
                      a("getMOA()",
                        href="http://www.rdocumentation.org/packages/shotGroups/functions/getMOA"),
                      "as well as for",
                      a("fromMOA()",
                        href="http://www.rdocumentation.org/packages/shotGroups/functions/fromMOA"),
                      "and the",
                      a("shotGroups vignette",
                        href="http://cran.rstudio.com/web/packages/shotGroups/vignettes/shotGroups.pdf"),
                      "section 3.5."),
                    h6("Acknowledgements"),
                    p("Thanks to David Bookstaber for testing, feedback and data.")
                ),

                id="task"
            )
        )
    )
))
