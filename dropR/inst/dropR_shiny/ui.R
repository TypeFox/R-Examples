library(shiny)
shinyUI(fluidPage(
  tags$style(HTML('
                  svg text{
                    font-size:15px;
                  }
                  
.tab-content {
height: 600px !important;
}

#quest_cols.shiny-bound-input{
height:250px !important;
}

.nvd3 .nv-groups path.nv-line{
stroke-width: 2.5px !important;
}

.alertme #chisq_table td{
background-color: rgba(255,0,0,.3) !important;

}

                  
                  ')),
  tags$h2("DropR"),
  fluidRow(
    tabsetPanel(
      tabPanel("Start",
               tags$div(align = 'center',
                        h2('DropR'),
                        h3('Analyze and Visualize Dropout'),
                        # tags$div(id = 'uploadData',tags$a("Upload Data")),
                        br(),
                        tags$img(src = 'decrease.svg',width=80)
                 )
               ),
      tabPanel('How It Works',
               h3('Upload Data'),
               p('Go to the Input Data Tab and upload your dataset as a .csv file. Choose the right delimiter, decimal and text quote characters. Make sure to include your questions in a one column per question format, preferably in the order of occurence in the questionnaire. After Uploading the data choose all columns that contain questions as well as the column that holds the experimental condition (groups). This needs to be done within the same tab before continuing.'),
               h3('Identifying Drop Out'),
               p('The subsequent tab show a visualization of drop out by experimental condition. Also the total is shown.')
               ),
      tabPanel("Input",column(width = 3,
               fileInput('file1', 'Choose CSV File',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv')),
               tags$hr(),
               checkboxInput('header', 'Header', TRUE),
               radioButtons('sep', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ';'),
               radioButtons('dec', 'Decimal Delimiter used in .csv',
                            c(Comma=',',
                              Dot='.'),
                            ','),
               radioButtons('quote', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"')
               ),
               column(width=6,tableOutput('out_table')),
               column(width = 8,
                 column(width = 3,uiOutput('choose_condition')),
                 column(width = 3,uiOutput('choose_questions')) 
                 )
               
               # column(width = 8,textOutput('test')) # lil debugger
               ),
      tabPanel("Remaining Participants by Condition",
               h3('Remaining Participants by Condition'),
               column(width = 8,lineChartOutput("mychart")),
               column(width = 2,
                      HTML('<h5>χ<sup>2</sup> Test for differences among all conditions</h5>'),
                      uiOutput('choose_question'),
                      uiOutput('color_table'),
                      tableOutput('contingency')
                      ),
               br(clear = 'all'),
               h5('Odds Ratio Condition Matrix (All Combinations)'),
               column(width = 8,tableOutput('or_table'))
               ),
      tabPanel("Credits",
               column(width = 6,
                      p("Visualization in this app uses the NVD3 JavaScript Charting library. The shiny app makes use of the shiny-examples provided by Joe Cheng. DropR follows an idea of Ulf-Dietrich Reips and Matthias Bannert.")
                      # This JS is a bit of hack, just for switching tabs... 
                      # anyway... if the number of tabs is changed this needs to be adjusted... 
                      # Did not work with the reactive stuff ... 
#                       HTML("<script>$('#uploadData').click(function() {
#     				 tabs = $('.tabbable .nav.nav-tabs li')
# 					 	 tabs.each(function() {
# 							$(this).removeClass('active')
# 					 	 })
# 						 $(tabs[2]).addClass('active')
# 						
# 						 tabsContents = $('.tabbable .tab-content .tab-pane')
# 					 	 tabsContents.each(function() {
# 							$(this).removeClass('active')
# 					 	 })
# 						 $(tabsContents[2]).addClass('active')
# 
# 						$('#input').trigger('change').trigger('shown');
# 						 
# 					 })</script>
# 			")
               )
      )
    )
  )
))




# https://github.com/novus/nvd3/wiki/API-Documentation


# this remains interesting: 
#Shiny.onInputChange("mydata", number);


