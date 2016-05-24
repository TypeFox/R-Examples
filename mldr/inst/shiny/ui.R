#
# This file belongs to mldrGUI, an EDA GUI for multilabel datasets developed on top of the mldr package.
# (c) 2015 - Francisco Charte Ojeda (fcharte@ugr.es), David Charte Luque (fdavidcl@outlook.com)
# See the package LICENSE file for license information
#

shinyUI(
  fluidPage(
    tags$head(
      tags$style("
        body {
          background: #F5F5FA;
        }
        .container-fluid {
          padding: 0;
          margin: 0;
        }

        header {
          background: #565667;
        }
          header h1 {
            margin: 0px;
            padding: 0.4em 0px;
            text-align: center;
            color: #F0F0F4;
            font-size: 1.6em;
            font-weight: 400;
          }

        .container-fluid .col-sm-8 {
          margin: 0 auto;
          float: none;
          width: auto;
          padding: 0;
          text-align: center;
        }

        #pages {
          text-align: center;
          background: #565667;
          margin-bottom: 1.2em;
        }
          #pages li {
            float: none;
            display: inline-block;
          }
            #pages li a, #pages li.active a {
              background: transparent;
              color: #fff;
              font-size: 1.1em;
              border-radius: 0;
              border-bottom: 0.2em solid transparent;
              outline: none;
            }
            #pages li.active a {
              border-bottom-color: #F5F5FA;
            }
            #pages li a:hover {
              background: rgba(255,255,255,0.2);
            }

        .tab-content {
          max-width: 80%;
          min-width: 66.6667%;
          display: inline-block;
          text-align: left;
        }
          .well {
            background: #FFF;
            color: #202020;
            border: 1px solid #EAEAF0;
            border-radius: 0;
            box-shadow: 0px 0px 3px -1px rgba(0, 0, 0, 0.2), 0px 3px 6px -3px rgba(0, 0, 0, 0.2);
            padding: 0.8em 1.2em;
            margin: 0px 0px 1.6em;
            line-height: 1.4rem;
            overflow: auto;
          }
            .well img {
              max-width: 100%;
              height: auto;
            }

            tr {
              border-left: .2em solid transparent;
            }
            tr.selected {
              border-left-color: rgba(30,90,180,0.7);
            }
        /** TODO
          Loading icon over .recalculating graphs
        **/
      "),
      tags$script(HTML("
        const ONECOL_WIDTH = 767;

        $(document).ready(function () {
          var graph = $('#labelLC').parent();
          var origOffsetY = $('.well').offset().top + $(window).scrollTop();
          var origDocHeight = $(document.body).height();
          var maxScroll = origDocHeight - $(graph).height();

          document.onscroll = function () {
            // Update maximum scroll distance if document height has changed
            if (origDocHeight != $(document.body).height() && $(window).scrollTop() < origOffsetY) {
              origDocHeight = $(document.body).height();
              maxScroll = origDocHeight - $(graph).height();

              console.log('docheight ' + origDocHeight + ' maxScroll ' + maxScroll);
            }

            if ($(window).scrollTop() >= origOffsetY && $(window).width() > ONECOL_WIDTH) {
              // Prevent graphic from making the document larger
              var curOffset = Math.min($(window).scrollTop(), maxScroll) - origOffsetY - 1;
              graph.css('margin-top', curOffset + 'px');
            } else {
              graph.css('margin-top', '0px');
            }
          };

          // Set a handler for the Quit button
          $(\"a[data-value='finish']\").attr('title', 'Quit mldr').click(function(e) {
            $('h1').html('You can close this tab now.')
            window.close();
          });
        });
      ")),
      tags$title("mldr"),
      tags$script(src='//cdn.datatables.net/1.10.2/js/jquery.dataTables.min.js',type='text/javascript'),
      tags$script(src='//cdn.datatables.net/tabletools/2.2.2/js/dataTables.tableTools.min.js',type='text/javascript'),
      tags$link(href='//cdn.datatables.net/tabletools/2.2.2/css/dataTables.tableTools.css',rel='stylesheet',type='text/css')
    ),
    tags$header(
      tags$h1(textOutput("title"))
    ),
    mainPanel(
      tabsetPanel(id = "pages", type = "pills", selected = "Main",
        tabPanel("Main",
          fluidRow(
            column(6,
              wellPanel(
                h3("Active MLD"),
                selectInput("mldrs", "Select a dataset", c()),
                hr(),
                h4("Load a dataset"),
                p("Select an ARFF and a XML file to load a MULAN dataset, or select only an ARFF file to load a MEKA dataset."),
                fileInput('arffname', 'Select the ARFF file'),
                fileInput('xmlname', 'Select the XML file'),
                actionButton("loadButton", "Load dataset")
              ),
              wellPanel(
                h3("How to use mldrGUI"),
                tags$small(paste0(
                  "mldrGUI is an EDA tool for multilabel datasets (MLDs).")), br(),
                tags$small(HTML(paste0("<ul><li>Use the controls above to select one of the MLDs included in the package, ",
                                  "or select an .arff and .xml file in your system to load any MLD.</li>"))),
                tags$small(HTML(paste0("<li>Once the MLD has been loaded, you will see its basic traits in this page.</li>"))),
                tags$small(HTML(paste0("<li>Use the tabs at the top of the page to explore other information, such as label",
                                  " distribution, frequency of labelsets, data about attributes or label concurrence information.</li>"))),
                tags$small(HTML(paste0("<li>Use the Quit button (<i class='fa fa-power-off'></i>) to close the application.</li></ul>"
                )))
              )
            ),
            column(6,
              fluidRow(
                wellPanel(
                  h3("Visual summary"),
                  downloadButton("saveAT", "Save plot"),
                  plotOutput("attributeByType", height = "auto"),
                  hr(),
                  downloadButton("saveCH", "Save plot"),
                  plotOutput("cardHistogram", height = "auto"),
                  hr(),
                  downloadButton("saveLH", "Save plot"),
                  plotOutput("labelHistogram", height = "auto"),
                  hr(),
                  downloadButton("saveLSH", "Save plot"),
                  plotOutput("labelsetHistogram", height = "auto")
                ),
                wellPanel(
                  h3("General summary"),
                  tableOutput("summaryGeneral")
                ),
                wellPanel(
                  h3("Label summary"),
                  tableOutput("summaryLabels")
                ),
                wellPanel(
                  h3("Labelset summary"),
                  tableOutput("summaryLabelsets")
                )
              )
            )
          )
        ),
        tabPanel("Labels", fluidPage(
          fluidRow(
            column(6, wellPanel(dataTableOutput("labels"))),
            column(6, wellPanel(
              downloadButton("saveLabels", "Save plot"),
              uiOutput("labelRange"),
              plotOutput("labelHC",height="auto"))
            )
          )
        )),
        tabPanel("Labelsets", fluidPage(
          fluidRow(
            column(6, wellPanel(dataTableOutput("labelsets"))),
            column(6, wellPanel(
              downloadButton("saveLabelsets", "Save plot"),
              plotOutput("labelsetHC",height="auto"))
            )
          )
        )),
        tabPanel("Attributes", fluidPage(
          wellPanel(dataTableOutput("attributes"))
        )),
        tabPanel("Concurrence", fluidPage(
          fluidRow(
            column(5,
                   wellPanel(
                     h3("Concurrence analysis"),
                     p("The SCUMBLE level for each label is shown in the table at the left."),
                     p("In the following table the minority labels most affected by SCUMBLE are shown."),
                     p("For each one of them, the names of the majority labels with interactions are provided."),
                     tableOutput("ConcurrenceAnalysis")
                   ),
                   wellPanel(
                     h3("Select the labels to plot"),
                     dataTableOutput("tblConcurrence"))
            ),
            column(7,
                   wellPanel(
                     downloadButton("saveConcurrence", "Save plot"),
                     plotOutput("labelLC",height="auto"))
            )
          )
        )),
        tabPanel("About", fluidPage(
          wellPanel(
            h3("About mldrGUI"),
            p(HTML('mldrGUI is an EDA GUI for multilabel datasets developed on top of the <a href="https://github.com/fcharte/mldr">mldr</a> package.')),
            p(HTML('&copy; 2015 &mdash; Francisco Charte Ojeda (fcharte@ugr.es), David Charte Luque (fdavidcl@outlook.com)')),
            p('See the package LICENSE file for license information'))
        )),
        tabPanel(HTML("<span class='quit-tab'><i class='fa fa-power-off'></i></span>"), value = "finish")
      )
    )
  )
)
