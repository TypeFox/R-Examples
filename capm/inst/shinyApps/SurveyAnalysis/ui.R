shinyUI(fluidPage(
  
  titlePanel("Survey analysis"),
  
  sidebarPanel(
    
    radioButtons('design',  HTML('<b>Survey design</b>'),
                 c('Two-stage cluster sampling' = 'twostage',
                   'Simple (systematic) random sampling' = 'systematic',
                   'Stratified rendom sampling' = 'stratified')),
    tags$hr(),
    
    conditionalPanel(
      condition = 'input.design == "twostage"',
      helpText('Choose a csv file having PSU unique identifiers in the first column and PSU sizes in the second column.'),
      checkboxInput('header', 'Header', TRUE),
      fluidRow(
        column(6,
               radioButtons('sep.two1', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ';')
        ),
        column(6,
               radioButtons('quote.two1', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"')
        )
      ),
      fileInput('universe', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      helpText('Choose a csv file with sample data.'),
      checkboxInput('header', 'Header', TRUE),
      fluidRow(
        column(6,
               radioButtons('sep.two2', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ';')
        ),
        column(6,
               radioButtons('quote.two2', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"'))
      ),    
      fileInput('sample', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      tags$hr(),
      checkboxInput('examples.twostage', 'Use example files from the capm package instead of your own data.', F),
      
      tags$hr(),
      HTML('<b>Design</b>'),
      numericInput('psu.col',
                   'Column in sample data that contain PSU identifiers:',
                   value = NULL, min = 0),
      numericInput('ssu.col',
                   'Column in sample data that contain SSU identifiers:',
                   value = NULL, min = 0),
      numericInput('psu.2cdl',
                   'Number  of PSU included in the design:',
                   value = NULL, min = 0)
    ),
    
    conditionalPanel(
      condition = 'input.design == "systematic"',
      helpText('Choose a csv file with sample data.'),
      checkboxInput('header', 'Header', TRUE),
      fluidRow(
        column(6,
               radioButtons('sep.syst', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ';')
        ),
        column(6,
               radioButtons('quote.syst', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"')
        )
      ),    
      fileInput('sample', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      tags$hr(),
      checkboxInput('example.systematic', 'Use example file from the capm package instead of your own data.', F),
      
      tags$hr(),
      HTML('<b>Design</b>'),
      numericInput('N', 'Total of sampling units in the population:',
                   value = NA, min = 1)
    ),
    
    conditionalPanel(
      condition = 'input.design == "stratified"',
      helpText('Choose a csv file with sample data.'),
      checkboxInput('header', 'Header', TRUE),
      fluidRow(
        column(6,
               radioButtons('sep.strat', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ';')
        ),
        column(6,
               radioButtons('quote.strat', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"')
        )
      ),    
      fileInput('sample', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      tags$hr(),
      checkboxInput('example.stratified', 'Use example file from the capm package instead of your own data.', F),
      
      tags$hr(),
      HTML('<b>Design</b>'),
      textInput('strata.sizes', HTML('Column in sample data that contain strata sizes:'), value = ''),
      textInput('strata.membership', HTML('Column in sample data that contain strat membership:'), value = '')
    ),
    
    numericInput('conf.level', HTML('<b>Confidence level</b><br>'),
                 value = 0.95, min = 0, max = 1),
    
    tags$hr(),
    HTML('<b>Population pyramid</b><br><br>'),
    numericInput('age.col',
                 'Column in sample data containing age variable.',
                 value = NULL, min = 0),
    numericInput('sex.col',
                 'Column in sample data containing sex variable.',
                 value = NULL, min = 0),
    numericInput('cas.col',
                 'Column in sample data containing reproductive status variable.',
                 value = NULL, min = 0),
    
    tags$hr(),
    HTML('<p><b>Type of estimates</b><br> If you are using example files, copy and paste the following terms: <br> total, prop, mean, prop, prop, total, prop, prop, prop, prop, prop, prop, prop, prop'),
    textInput('variables', '',
              value = NULL),
    
    tags$hr(),
    HTML('<b>Press the button and wait</b><br>'),
    actionButton('calc.estimates', 'Get estimates')
  ),
  
  mainPanel(
    tabsetPanel(
      
      tabPanel(
        'Introduction', value = 1,
        HTML(
          '
<p>Here, you can upload data collected in a survey, make some descriptive analysis, and estimate population parameters.</p>

<p>Population parameters can be estimated according to one of the three choices listed on the top of the left side panel. After choosing the survey design and defining the required information, click on the <i>Selected sampling units</i> Tab.</p><br>

<b>Two-stage cluster sampling</b>
<p>Upload the two required csv files and make sure you choose the appropriate options in the left side panel to avoid error messages or awkward results. You can check if files are appropriatley imported, looking at the <i>Uploaded files</i> Tab. A summary of the sample data can be viewed in the <i>Sample data summary</i> Tab. Insted of using your own files, you can use example files from the capm package, checking the respective box on the left side panel.</p>

<p>Fill in the required fields to define the survey design. For the example file, use 2, 1 and 20 in the first three fields and define the desired confidence level for the estimates. The numbers 2 and 1 correspond to the required column positions in the sample data (in the <i>Uploaded files</i> Tab). Column names insted of colum positions can be used too. The number 20 is the number of PSU included in the example survey design.</p><br>

<b>Systematic sampling</b>
<p>Upload the required csv file and make sure you choose the appropriate options in the left side panel to avoid error messages or awkward results. You can check if file are appropriatley imported, looking at the <i>Uploaded files</i> Tab. A summary of the sample data can be viewed in the <i>Sample data summary</i> Tab. Alternatively, you can use the example file from the capm package, checking the respective box on the left side panel. To define the survey design, just specify the total of sampling units in the population and the desired confidence level for the estimates.</p><br>

<b>Stratified sampling</b>
<p>Upload the required csv file and make sure you choose the appropriate options in the left side panel to avoid error messages or awkward results. You can check if file are appropriatley imported, looking at the <i>Uploaded files</i> Tab. A summary of the sample data can be viewed in the <i>Sample data summary</i> Tab. Alternatively, you can use the example file from the capm package, checking the respective box on the left side panel.</p>

<p>Fill in the required fields to define the survey design. For the example file, use 16, and 15 in the first two fields and define the desired confidence level for the estimates. The numbers 16 and 15 correspond to the required column positions in the sample data (in the <i>Uploaded files</i> Tab). Column names insted of colum positions can be used too.</p><br>

<b>Sample data summary</b>
<p>To see a description of the meaning of each variable, go to RStudio and run <code>help(survey.data)</code>.</p><br>

<p>Population pyramids</p>
<p>If you checked the box to use example files, fill in the fields from top to bottom with, 5, 4 and 6 respectively (see this column positions in the output from <i>Survey data</i> Tab).</p><br>

<b>Estimates</b>
<p>Specify the type of estimate for each variable. To estimate a total, type <i>total</i>, to estimate a mean, type <i>mean</i> and to estimate a proportion, type <i>prop</i>. To get estimates of more than one variable, type the appropriate term for each variable and separate them by commas. The order of terms from left to right must be the order (left to right) of variables (columns) in the survey file.</p>

<p>From the command line in RStudio, it is possible to get estimates for specific subpopulations (i.e. by sex). To estimate the total for categorical variables such as "castrated" from the example file (<i>Survey data</i> Tab), the command line is the place to go.</p>

<b>Further information</b>
<p>
<ul>
<li>The display of some outputs might take a few seconds, be patient!</li>
<li>Reload the page to reset the fields</li>
<li>Working from command line, you will have more options and flexibility.</li>
<li>Tutorials with more detailed information can be found in <a href="https://github.com/oswaldosantos/capm">https://github.com/oswaldosantos/capm</a></li>
<li>If you find errors, have suggestions or any question, I will be glad to know it <a href="mailto:oswaldosant@gmail.com">oswaldosant@gmail.com</a></li>
</ul>
</p>
          ')),
      tabPanel('Uploaded files',
               h4(textOutput('file.title1')),
               dataTableOutput('universe'),
               tags$hr(),
               h4(textOutput('file.title2')),
               dataTableOutput('sample')),
          
      tabPanel('Sample data summary',
               verbatimTextOutput('summ.sample')),
      
      tabPanel('Population pyramid', value = 3,
               plotOutput('pyramid', height = 600)),
      
      tabPanel('Estimates', value = 4,
               tableOutput('variables'),
               tags$hr(),
               tableOutput('estimates')),
      
      id = "conditionedPanels"
    )
  )
))