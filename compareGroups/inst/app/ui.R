shinyUI(fluidPage(theme = shinytheme("spacelab"),
  
    shinyjs::useShinyjs(), # Set up shinyjs
    
    headerPanel("",windowTitle="compareGroups | Explore and Summarise Epidemiological Data in R"),

    HTML('<style type="text/css"> #inputpanel{ max-height:800px; overflow:auto; margin-right:-3%} </style>'),
    HTML('<style type="text/css"> #outpanel {min-height:150px;} </style>'), 
    HTML('<style type="text/css"> #showload { background-color: rgb(250,250,250); border: 2px solid grey;} </style>'),
    HTML('<style type="text/css"> #selevars { width: 120px;} </style>'), 
    HTML('<style type="text/css"> #discvars { width: 120px;} </style>'),
    HTML('<style type="text/css"> #changeselevars { font-size: 15px; border: 2px solid grey;} </style>'),
    HTML('<style type="text/css"> #varPlot { width: 120px;} </style>'),
    HTML('<style type="text/css"> #maxvalues { width: 100px;} </style>'),
    HTML('<style type="text/css"> #exampledata { width: 200px;} </style>'),
    HTML('<style type="text/css"> #initial { height: 0px; width: 0px} </style>'),
    HTML('<style type="text/css"> #count { height: 0px; width: 0px} </style>'),
    HTML('<style type="text/css"> #bsModalhelpcg .modal-content { width:120%;height:120%;} </style>'),
    HTML('<style type="text/css"> #bsModalhelpcg .modal-header { background-color:rgb(68,110,155);} </style>'),
    

    fluidRow(

    ##################################################     
    ################ CONTROL PANEL ###################
    ##################################################
    
    column(4,

      div(id="inputpanel",  

          bsCollapse(id = "collapseInput", open = c("collapseLoad"),
            
            ## load data
            bsCollapsePanel(title="Step 1. Load data", value="collapseLoad", style = "primary",
              uiOutput("initial"),
              column(1,bsButton("infoLoad","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoLoadModal",HTML('<p> <strong>Step 1. Load data</strong></p>'), "infoLoad",uiOutput("helpload")),
              #radioButtons("exampledata", "", choices = c("Own data","REGICOR","PREDIMED","SNPS")),
              HTML('
              <div id="exampledata" class="form-group shiny-input-radiogroup shiny-input-container">
                <label class="control-label" for="exampledata"></label>
                <div class="shiny-options-group">
                <div class="radio">
                <label>
                <input type="radio" name="exampledata" value="Own data" checked="checked"/>
                <span>Upload your own data set from a file</span>
                </label>
                <div><br><font style="weight:bold">...or choose an example data set</font></p></div>
                </div>
                <div class="radio" style="margin-left:10px">
                <label>
                <input type="radio" name="exampledata" value="REGICOR"/>
                <span>REGICOR</span>
                </label>
                </div>
                <div class="radio" style="margin-left:10px">
                <label>
                <input type="radio" name="exampledata" value="PREDIMED"/>
                <span>PREDIMED</span>
                </label>
                </div>
                <div class="radio" style="margin-left:10px">
                <label>
                <input type="radio" name="exampledata" value="SNPS"/>
                <span>SNPS</span>
                </label>
                </div>
                </div>
                </div>'), 
              conditionalPanel(
                condition = "input.exampledata == 'Own data'",
                fileInput("files", "", accept=c('', 'sav', 'csv', 'dat', 'txt', 'xls', 'xlsx', 'mdb')),
                selectInput("datatype", "Data format", c('SPSS'='*.sav', 'TEXT'='*.txt','EXCEL'='*.xls','R'='*.rda'),'*.sav'),
                bsButton("encodingaction", "Encoding", size="extra-small", style="info"),
                radioButtons("encoding", "", c('default'='default','latin1'='latin1','utf8'='utf8'),'default',inline=TRUE),
                uiOutput("loadoptions")
              ),
              br(),
              actionButton("loadok","OK")
            ),
            
            ## vars list
            bsCollapsePanel(title="Step 2. Select variables", value="collapseSelect", style = "primary", 
              column(1,bsButton("infoSelect","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoSelectModal",HTML('<p> <strong>Step 2. Select variables</strong></p>'), "infoSelect",uiOutput("helpselect")),
              uiOutput("selevarslist")
            ),
            ## settings
            bsCollapsePanel(title="Step 3. Settings", value="collapseSettings", style = "primary", 
              column(1,bsButton("infoSettings","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoSettingsModal",HTML('<p> <strong>Step 3. Settings</strong></p>'), "infoSettings",
                tabsetPanel(
                  tabPanel("Type",uiOutput("helptype")),
                  tabPanel("Response",uiOutput("helpresponse")),
                  tabPanel("Hide",uiOutput("helphide")),
                  tabPanel("Subset",uiOutput("helpsubset")),
                  tabPanel("OR/HR", uiOutput("helpratio"))
                )
              ),
              conditionalPanel(
                condition = "input.initial",
                tabsetPanel(
                  ## methods
                  tabPanel("Type",
                    uiOutput("selemethod"),
                    uiOutput("selemethodNA")
                  ),
                  ## response
                  tabPanel("Response",
                    radioButtons("resptype", "", c("None","Group","Survival"),"None"),                 
                    actionButton("changeresp","","Update"),
                    uiOutput("response")
                  ),            
                  ## hide
                  tabPanel("Hide",
                    br(),
                    wellPanel(
                      fluidRow(
                        column(6,uiOutput("selehidevar")),
                        column(6,uiOutput("selehidecat"))
                      )
                    ),
                    textInput('hideno', "Hide 'no' category", ''),
                    actionButton("changehide","Update")
                  ),
                  ## subset
                  tabPanel("Subset", uiOutput("selevarsubset")),
                  ## OR/HR ratio for row-variables
                  tabPanel("OR/HR", uiOutput("ratio"))
                )
              )
            ),
            ## display
            bsCollapsePanel(title="Step 4. Display", value="collapseDisplay", style = "primary",
              column(1,bsButton("infoDisplay","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoDisplayModal",HTML('<p> <strong>Step 4. Display</strong></p>'), "infoDisplay",
                tabsetPanel(
                  tabPanel("Show",uiOutput("helpshow")),
                  tabPanel("Format",uiOutput("helpformat")),
                  tabPanel("Decimals",uiOutput("helpdecimals")),
                  tabPanel("Labels",uiOutput("helplabel"))
                )
              ),    
              conditionalPanel(
                condition = "input.initial",
                tabsetPanel(
                  ## show
                  tabPanel("Show",uiOutput("show")),
                  ## Format display
                  tabPanel("Format", uiOutput("format")),
                  ## number of decimals
                  tabPanel("Decimals", uiOutput("decimals")),
                  ## header labels
                  tabPanel("Labels", uiOutput("labels"))
                )
              )
            ),
            ## save
            bsCollapsePanel(title="Step 5. Save table", value="collapseSave", style = "primary",
              column(1,bsButton("infoSave","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoSaveModal",HTML('<p> <strong>Step 5. Save table</strong></p>'), "infoSave", uiOutput("helpsave")),
              conditionalPanel(
                condition = "input.initial",  
                selectInput("downloadtabletype", "Select format", choices = c("PDF","CSV","HTML","TXT","Word","Excel"),selectize=FALSE),
                conditionalPanel(
                  condition="input.downloadtabletype == 'PDF'",
                  wellPanel(
                    selectInput('sizepdf', 'Resize', c("tiny","scriptsize","footnotesize","small","normalsize","large","Large","LARGE","huge","Huge"),"normalsize", selectize=FALSE),
                    h4(""),
                    checkboxInput('landscape', 'Landscape', FALSE)
                  )
                ),
                conditionalPanel(        
                  condition="input.downloadtabletype == 'CSV'",
                  wellPanel(
                    radioButtons('sepcsv', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ',')
                  )
                ),
                downloadButton('actiondownloadtable', 'Download')
              )
            )
          )
      )
    ),
    
    #################################
    #######  RESULTS PANEL ##########
    #################################

    column(8,
           
      div(id="outpanel",

          navbarPage(title="", id="results", inverse=FALSE, 
            # About
            tabPanel(value="resultsAbout",HTML('<p title="About compareGroups project">ABOUT</p>'),
              # compareGroups
              uiOutput("helpabout"),
              br(),
              column(4,bsButton("helpcg","view examples",size="small",type="info"),offset=4),
              bsModal(id="bsModalhelpcg", title=HTML('<font style="color:white;">Examples</font>'), trigger="helpcg",size="large",
                fluidRow(
                  column(1,bsButton("dec",HTML('<p style="font-size:20px"><</p>'),style="primary",size="extra-small")),
                  column(10,uiOutput("helpModalContents")),
                  column(1,bsButton("inc",HTML('<p style="font-size:20px">></p>'),style="primary",size="extra-small"))
                )
              ),  
              br(),
              # WUI
              uiOutput("helpwui"),
              # Security
              uiOutput("helpsecurity")
            ),                     
            # Data
            navbarMenu("DATA",
              # summary         
              tabPanel(value="resultsSummary",title=HTML('<p title="Short summary from loaded data set">Summary</p>'),
                column(1,bsButton("infoSummary","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
                bsModal("infoSummaryModal",HTML('<p> <strong>Summary</strong></p>'), "infoSummary",uiOutput("helpsummary")),
                uiOutput("values")
              ),
              # VALUES (extended)
              tabPanel(value="resultsValues",title=HTML('<p title="Navigate thru the whole data set">Values</p>'),
                column(1,bsButton("infoExtended","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
                bsModal("infoExtendedModal",HTML('<p> <strong>Extended</strong></p>'), "infoExtended",uiOutput("helpvalues")),
                conditionalPanel(
                  condition="input.initial",
                  bsButton("valuextoptionsaction","View",style="info"),
                  wellPanel(id="valuextoptions",sliderInput("valueextsize", "Resize (%):", min=10, max=300, value=100, step=10))
                ),
                uiOutput("valuesext")
              )
            ),
            tabPanel(value="resultsTable",title=HTML('<p title="View descriptive table" style="font-size:140%;">TABLE</p>'),
              column(1,bsButton("infoTable","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoTableModal",HTML('<p> <strong>TABLE</strong></p>'), "infoTable",uiOutput("helptable")),
              conditionalPanel(
                condition = "input.initial",
                br(),
                bsButton("tableoptionsaction", "View (hide)", style="info"),
                wellPanel(id="tableoptions",
                  sliderInput("htmlsizerestab", "Resize:", min=0.5, max=2, value=1, step=0.1),
                  bsButton("tableinfo","Info"),
                  bsModal("tableinfoModal",title="Info table",trigger="tableinfo",
                    htmlOutput("sumtab")
                  )
                ),
                uiOutput("table")
              )        
            ),
            tabPanel(value="resultsPlot",title=HTML('<p title="Visualize data">PLOT</p>'),
              column(1,bsButton("infoPlot","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoPlotModal",HTML('<p> <strong>PLOT</strong></p>'), "infoPlot",uiOutput("helpplot")),
              uiOutput("uiplot")
            ),
            tabPanel(value="resultsSNPs",title=HTML('<p title="single nucleotide polymorphisms (SNPs) analyses">SNPs</p>'),
              column(1,bsButton("infoSNPs","",size="extra-small",style="info",icon=icon("info-circle")),offset=11),
              bsModal("infoSNPsModal",HTML('<p> <strong>SNPs</strong></p>'), "infoSNPs",uiOutput("helpsnps")),
              uiOutput("snps")
            )
          )
    
      )
    )
  )

))                               
