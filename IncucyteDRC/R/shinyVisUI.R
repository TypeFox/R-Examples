#' Shiny Visualisation UI
#'
#' Creates a shiny UI for the interactive web app for the running the IncucyteDRC workflow
#'
#' @return A shiny UI
#' @export
#' @import shiny

shinyVisUI <- function() {

    navbarPage(title = "IncucyteDRC",
               tabPanel("Input",
                        sidebarLayout(
                            sidebarPanel(
                                fileInput('platemap_file', 'Choose platemap file',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values',
                                                    'text/plain',
                                                   '.Platemap',
                                                   '.txt')),
                                fileInput('data_file', 'Choose data file',
                                          accept=c('text/csv',
                                                   'text/comma-separated-values',
                                                   'text/plain',
                                                   '.txt')),
                                hr(),
                                conditionalPanel(
                                    condition = "input.cut_time_mode == false",
                                    sliderInput('cut_time_slider', 'Specify cut time', 1,300, 175)
                                    ),
                                conditionalPanel(
                                    condition = "input.cut_time_mode == true",
                                    sliderInput('baseline_time_slider', 'Specify baseline time', 1,100, 24),
                                    sliderInput('max_val_slider', 'Specify maximum value', 20,100, 80),
                                    sliderInput('no_doublings_slider', 'Specify # doublings', 1,6, 4, 0.1)
                                    ),
                                checkboxInput('cut_time_mode', 'Calculate cut time', value=FALSE),
                                checkboxInput('include_control_mode', 'Include control data in dose response', value=FALSE),


                                selectInput('group_columns_select', 'Select group columns',
                                            choices=c('growthcondition', 'celltype', 'passage', 'seedingdensity'),
                                            selected='growthcondition',
                                            multiple=TRUE,
                                            selectize=TRUE)
                                ),
                            uiOutput('mainpage_ui')

                        )),
               tabPanel("Data",
                        selectInput('data_format_select', 'Select data format',
                                    choices=c('Data Frame' = 'dataframe',
                                              'PRISM' = 'prism',
                                              'Dotmatics' = 'dotmatics'),
                                    selected='dataframe',
                                    multiple=FALSE,
                                    selectize=FALSE),
                        uiOutput('datapage_ui')
                        ),
               tabPanel('EC50',
                        uiOutput('ec50page_ui')),
               navbarMenu('More...',

                   tabPanel('Platemap',
                            uiOutput('platemap_ui')),
                   tabPanel('Cut time calculation',
                            uiOutput('cut_time_res_ui'))
                   ),
               tabPanel("Help",
                        selectInput('help_select', 'Select help vignette',
                                    choices=c('None' = 'none',
                                              'Overview Vignette' = 'overview',
                                              'Export from Incucyte Zoom' = 'zoom',
                                              'Video Tutorial' = 'video'),
                                    selected='none',
                                    multiple=FALSE,
                                    selectize=FALSE),
                        uiOutput('help_ui')))




}
