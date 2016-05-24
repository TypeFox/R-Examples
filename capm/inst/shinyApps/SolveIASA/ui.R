shinyUI(fluidPage(
  
  titlePanel('IASA model: Immigration, Abandonment, Sterilization and Adoption'),
  
  sidebarPanel(
    
    conditionalPanel(
      condition = 'input.conditionedPanels == 1', value = 1,
      HTML('<p><b>Initial conditions</b></p>'),
      br(),
      HTML('<p>Owned population</p>'),
      fluidRow(
        column(6,
               numericInput('f1', 'Intact females (f1)',
                            value = 33425, min = 0),
               numericInput('m1', 'Intact males (m1)',
                            value = 38039, min = 0)
               ),
        column(6,
               numericInput('fs1', 'Sterilized females (fs1)',
                            value = 10865, min = 0),
               numericInput('ms1', 'Sterilized males (ms1)',
                            value = 6008, min = 0)
               )
        ),
            
      br(),
      HTML('<p>Stray population</p>'),
      fluidRow(
        column(6,
               numericInput('f2', 'Intact females (f2)',
                            value = 3343, min = 0),
               numericInput('m2', 'Intact males (m2)',
                            value = 3804, min = 0)
               ),
        column(6,
               numericInput('fs2', 'Sterilized females (fs2)',
                            value = 109, min = 0),
               numericInput('ms2', 'Sterilized males (ms2)',
                            value = 68, min = 0))
        ),
      
      tags$hr(),
      HTML('<p><b>Parameters</b></p>'),
      br(),
      HTML('<p>Owned population</p>'),
      fluidRow(
        column(6,
               numericInput('b1', 'Births per year (b1)',
                            value = 21871, min = 0),
               numericInput('dm1', 'Male death rate (dm1)',
                            value = 0.098, min = 0, step = 0.01),
               numericInput('sm1', 'Male steril. rate (sm1)',
                            value = 0.028, min = 0, step = 0.01),
               numericInput('h1', 'Mean harem size (h1)',
                            value = 1, min = 0, step = 0.1)
               ),
        column(6,
               numericInput('df1', 'Female death rate (df1)',
                            value = 0.104, min = 0, step = 0.01),
               numericInput('sf1', 'Female steril. rate (sf1)',
                            value = 0.069, min = 0, step = 0.01),
               numericInput('k1', 'Carrying capacity (k1)',
                            value = 98050, min = 0),
               numericInput('ab', 'Abandonment rate (ab)',
                            value = 0.054, min = 0, step = 0.01)
               )
        ),
      fluidRow(
        column(5,
               numericInput('v', 'Immigration rate (v)',
                            value = 0.2, min = 0, step = 0.01)
               ),
        column(7,
               numericInput('z', 'Prop. of steril. immigrants (z)',
                            value = 0.1, min = 0, step = 0.01)
               )
        ),
      
      br(),
      HTML('<p>Stray population</p>'),
      fluidRow(
        column(6,
               numericInput('b2', 'Births per year (b2)',
                            value = 4374, min = 0),
               numericInput('dm2', 'Male death rate (dm2)',
                            value = 0.118, min = 0, step = 0.01),
               numericInput('sm2', 'Male steril. rate (sm2)',
                            value = 0.03, min = 0, step = 0.01),
               numericInput('h2', 'Mean harem size (h2)',
                            value = 0.5, min = 0, step = 0.1)
               ),
        column(6,
               numericInput('df2', 'Female death rate (df2)',
                            value = 0.125, min = 0, step = 0.01),
               numericInput('sf2', 'Female steril. rate (sf2)',
                            value = 0.05, min = 0, step = 0.01),
               numericInput('k2', 'Carrying capacity (k2)',
                            value = 8055, min = 0),
               numericInput('ad', 'Adoption rate (ad)',
                            value = 0.1, min = 0, step = 0.01)
               )
        ),
      
      tags$hr(),    
      numericInput('time', strong('Simulation time (initial time is 0)'), value = 15),
      numericInput('t_steps', 'Time steps', value = 1),
      
      tags$hr(),
      selectInput('output_var', strong('Choose a population'), 
                  choices = list(
                    'Owned intact animals (n1)',
                    'Owned sterilized animals (ns1)',
                    'Stray intact animals (n2)',
                    'Stray sterilized animals (ns2)',
                    'Owned animals (N1)',
                    'Stray animals (N2)',
                    'Total population (N)',
                    'Owned intact females (f1)',
                    'Owned sterilized females (fs1)',
                    'Owned intact males (m1)',
                    'Owned sterilized males (ms1)',
                    'Stray intact females (f2)',
                    'Stray sterilized females (fs2)',
                    'Stray intact males (m2)',
                    'Stray sterilized males (ms2)'),
                  'Owned sterilized animals (ns1)'
      )
    ),
    
    conditionalPanel(
      condition = 'input.conditionedPanels == 2', value = 2,
      HTML('<p><b>Sensitivities</b></p>'),
      tags$hr(),
      numericInput('range',
                   'Proportional perturbation of parameters in global sensitivity',
                   value = 0.1, step = 0.01),
      helpText('Press the button and wait until full color plots appear'),
      actionButton('sensitivities', 'Calculate sensitivities'),
      br(),
      HTML('If you pressed this button in another Tab <i>(Global/Local sensitivities)</i>, do not press it again at least you want to recalculate sensitivities.')
    ),
    
    conditionalPanel(
      condition = 'input.conditionedPanels == 3', value = 3,
      HTML('<p><b>Scenarios</b></p>'),
      
      tags$hr(),
      sliderInput('s.range',
                  'Sterilization range',
                  min = 0, max = 0.8, value = c(0, 0.2), step = 0.01),
      numericInput('s.intr', 'Number of intervals', value = 10),
      checkboxInput('s.fm', 'Sterilization of only females'),
      
      tags$hr(),
      sliderInput('ab.range',
                  'Abandonment range',
                  min = 0, max = 0.8, value = c(0, 0.2), step = 0.01),
      
      tags$hr(),
      sliderInput('ad.range',
                  'Adoption range',
                  min = 0, max = 0.8, value = c(0, 0.2), step = 0.01),
      
      tags$hr(),
      sliderInput('im.1',
                  'Immigration rate 1',
                  min = 0, max = 0.8, value = 0, step = 0.01),
      sliderInput('im.2',
                  'Immigration rate 2',
                  min = 0, max = 0.8, value = 0.2, step = 0.01),
      tags$hr(),
      selectInput('output_var2', strong('Choose a population'), 
                  choices = list(
                    'Intact females (f)',
                    'Sterilized females (fs)',
                    'Intact males (m)',
                    'Sterilized males (ms)',
                    'Intact animals (n)',
                    'Sterilized animals (ns)',
                    'Total population (N)'),
                  'Sterilized animals (ns)'
      ),
      
      tags$hr(),
      helpText('Press the button and wait until full color plots appear'),
      actionButton('create.scenarios', 'Create scenarios')
    )
    
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel('Introduction', value = 1,
               HTML(
                 '<p>With this App, you will run a mathematical model of population dynamics. We are preparing a peer-reviewed paper to explain technical details. See also the link at the bottom.</p>
                 
<p>Set the initial conditions and parameters. Define realistic values to avoid awkward results. The longer the simulation time, the longer the time to calculate and display results. The lesser the time step, the greater the output resoultion and the time to calculate and display results.</p>

<p><i>Abbreviations:</i> steril: sterilization; Prop: proportion.</p>

<b>Point estimates</b><br>
<p>In <i>Point estimates - plot</i> and <i>Point estimates - table</i> Tabs you will see the simulation results, which use values defined in the left side panel from <i>Introduction</i>, <i>Point estimates - plot, Point estimates - table</i> Tabs.</p>

<p>Each column in <i>Point estimates - table</i> Tab corresponds to a specific population. The meaning of header can be seen when clicking on the drop-down list at the bottom of the left side panel from the first three Tabs.</p>

<b>Sensitivities</b><br>
<p>In global sensitivity analysis, you will use the point estimates defined in the left side panel form the first three Tabs. However, each parameter will be allowed to take values from a range around the point estimate (parameters are perturbated). That range are definied in the first field of the left side panel from the <i>Global sensitivities</i> Tab. The model will be run 100 times, each time with a random set of parameter values. The output plot will have a reddish envelop representing all results from the 100 simulations. The bluish envelop represents the most common results (deviations from mean no greater than a standard deviation). The plot at the top shows results from simulations in which all parameters are perturbated. The plot at the bottom shows results from simulations in which one parameter is perturbated and the others remain fixed.</p>

<p>In local sensitivity analysis, parameters suffer small perturbations around its point estimates (nominal values). The plots of the norms (y axis labelled as L1 and L2) represent the importance of each parameter. The greater the bar, the greater the importance of the parameter.</p>

<p><i>After pressing the "Calculate sensitivities" button, sensitivities will be calculated and you must to wait. Remember that the model runs 100 times! Be patient and avoid screen resizing.</i></p>

<p>Plots will be displayed, then they vanish a little and then are plotted again. I will try to fix this behaviour in future versions.</p>

<b>Scenarios</b><br>
<p>In <i>Scenarios</i> Tab, You will use the point estimates defined in the left side panel from the first three Tabs, plus a set of values for the following parameters:
<ul>
<li>Sterilization rate.</li>
<li>Abandonment rate.</li>
<li>Adoption rate.</li>
<li>Immigration rate.</li>
</ul>
</p>
<p>Scenarios will correspond to all combinations given by those additional sets. For sterilization rate, define a range using the slider and the number of intervals. The greater the number of intervals, the greater the output resolution and the time to calculate and display results. For abandonment, adoption and immigration rates, just use the sliders.</p>

<p><i>After pressing the "Create scenarios" button, scenarios will be created and you must to wait. Remember that the model runs as many times as combination given by the sets of rates! Be patient and avoid screen resizing.</i></p>

<p>Plots will be displayed, then they vanish a little and then are plotted again. I will try to fix this behaviour in future versions.</p>

<p>
<b>Further information</b><br>
<ul>
<li>The display of some outputs might take a few seconds, be patient!</li>
<li>Reload the page to reset the fields.</li>
<li>Avoid screen resizing.</li>
<li>Working from command line, you will have more options and flexibility.</li>
<li>Tutorials with more detailed information can be found in <a href="https://github.com/oswaldosantos/capm">https://github.com/oswaldosantos/capm</a></li>
<li>If you find errors, have suggestions or any question, I will be glad to know it <a href="mailto:oswaldosant@gmail.com">oswaldosant@gmail.com</a></li>
</ul>
</P>
                 ')
      ),
      tabPanel('Point estimates - plot', value = 1,
               plotOutput('points_p', height = 600)),
      tabPanel('Point estimates - table', value = 1,
               dataTableOutput('points_t')),
      tabPanel('Global sensitivities', value = 2,
               plotOutput('globalall', height = 600),
               plotOutput('global', height = 600)),
      tabPanel('Local sensitivities', value = 2,
               plotOutput('local', height = 1200)),
      tabPanel('Scenarios', value = 3,
               plotOutput('scenarios', height = 600)),
      id = 'conditionedPanels'
    )
  )
))