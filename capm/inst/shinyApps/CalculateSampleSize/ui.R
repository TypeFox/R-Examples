shinyUI(fluidPage(
  
  titlePanel("Calculate sample size and composition"),
  
  sidebarPanel(
    radioButtons('design',  HTML('<b>Survey design</b>'),
                 c('Two-stage cluster sampling' = 'twostage',
                   'Simple (systematic) random sampling' = 'systematic',
                   'Stratified rendom sampling' = 'stratified')),
    tags$hr(),
    conditionalPanel(
      condition = "input.design == 'twostage'",
      HTML('<p>Choose a csv file having PSU unique identifiers in the first column and PSU sizes in the second column.</p>'),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ','),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),    
      fileInput('psu.ssu', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      HTML('<p>Choose a csv file having PSU unique identifiers in the first column and the totals observed in a pilot sample in the second column.</p>'),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Semicolon=';',
                     Tab='\t'),
                   ','),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),    
      fileInput('psu.x', '',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv')),
      
      tags$hr(),
      checkboxInput('examples', 'Instead of choosing your own csv files, use the example files from the capm package.', F),
      
      tags$hr(),
      numericInput('level', 'Confidence level',
                   value = 0.95, min = 0, max = 1),
      numericInput('error', 'Accepted error',
                   value = 0.1, min = 0, max = 1),
      numericInput('cost', 'Cost',
                   value = 4, min = 0),
      numericInput('min.ssu', 'Minimum number of SSU',
                   value = 15, min = 0)
    ),
    
    conditionalPanel(
      condition = "input.design == 'systematic'",
      numericInput('N', 'Total of sampling units in the population.',
                   value = NULL, min = 0, max = 1),
      numericInput('expected.mean', 'Expected mean.', value = NULL, min = 1),
      numericInput('expected.var', 'Expected variance.',
                   value = NULL, min = 1),
      numericInput('level', 'Confidence level',
                   value = 0.95, min = 0, max = 1),
      numericInput('error', 'Accepted error',
                   value = 0.1, min = 0, max = 1)
    ),
    
    conditionalPanel(
      condition = "input.design == 'stratified'",
      textInput('strata.names', HTML('Names of each strata:<br>Use "," to separate values (e.g. Urabn,Rural).'), value = ''),
      textInput('strata.N', HTML('Total of sampling units per strata:<br>Use "," to separate values (e.g. 100,50).'), value = ''),
      textInput('strata.mean',HTML('Expected mean per strata:<br>
                   Use "," to separate values.'), value = ''),
      textInput('strata.var', HTML('Expected variance per strata:<br>
                                      Use "," to separate values.'),
                value = NULL),
      numericInput('level', 'Confidence level',
                   value = 0.95, min = 0, max = 1),
      numericInput('error', 'Accepted error',
                   value = 0.1, min = 0, max = 1)
    )
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(
        'Introduction',
        HTML(
          '<p>The sample size can be calculated according to one of the three designs listed on the top of the left side panel. After choosing the survey design and defining the required information, click on the <i>Sample size</i> Tab.</p><br>

<p><b>Two-stage cluster sampling</b></p>
<p>Suppose that census tracks are PSU and the number of households in each census track represent primary sampling units (PSU) sizes. In the left side panel, you are asked to choose two csv files. The first file must have just two columns with that information. The second file must have one row per secondary sampling unit (SSU) sampled in a pilot study. The first column contains the PSU identifier to which the SSU belongs to. The second column contains the total observed in that SSU, for the variable of interest.</P>

<p>Make sure you choose the appropriate options (header, separator and quote), otherwise, you will get an error or an awkward result. You can also use example files from capm package cheking the respective box in the left side panel. In this case you do not need to choose any csv file.</p>

<p> The confidence level and the accepted error must be numbers between 0 and 1 inclusive. <i>Cost </i>is the ratio of the cost of sampling a PSU to the cost of sampling a SSU. <i>Minimum number of SSU </i> is the minimum number of SSU to be selected per PSU.</p><br>

<p><b>Simple random sampling</b></p>
<p>Specify the total number of sampling units (i.e. the number of households) in the population and the expected mean and variance (both values can be obtained from a pilot study or from a survey made in a similar context). The confidence level and the accepted error must be numbers between 0 and 1 inclusive.</p><br>

<p><b>Stratified random sampling</b></p>
<p>Define a name for each strata and specify the total number of sampling units (i.e. the number of households) in the population. Specify also the expected mean and variance for each strata (both values can be obtained from a pilot study or from a survey made in a similar context). The confidence level and the accepted error must be numbers between 0 and 1 inclusive.
</p><br>

<p>
<b>Further information</b><br>
<ul>
<li>Reload the page to reset the fields.</li>
<li>Working from command line, you will have more options and flexibility.</li>
<li>Tutorials with more detailed information can be found in <a href="https://github.com/oswaldosantos/capm">https://github.com/oswaldosantos/capm</a></li>
<li>If you find errors, have suggestions or any question, I will be glad to know it <a href="mailto:oswaldosant@gmail.com">oswaldosant@gmail.com</a></li>
</ul>
</P>')
      ),
      tabPanel('Sample size',
               tableOutput('size'),
               tags$hr(),
               h4(textOutput('file.title1')),
               dataTableOutput('universe'),
               tags$hr(),
               h4(textOutput('file.title2')),
               dataTableOutput('sample.data'))
    )
  )
))