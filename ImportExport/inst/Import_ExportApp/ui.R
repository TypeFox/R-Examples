shinyUI(fluidPage(
  titlePanel("Import - Export"),
  HTML("<style type='text/css'> #vars { background-color: rgb(250,250,250);
       height:250px} </style>"),
  fluidPage(
    column(4,
      bsCollapse(id="collapse",multiple=T,open="Step 1) Format detection + import",
                 bsCollapsePanel("Step 1) Format detection + import",style="primary",
      fileInput("file",label="Select file"),
                 
                 checkboxInput("correctFormat",label="format well detected",value=T),
                 
                 conditionalPanel( condition="!input.correctFormat",selectInput("newformat","select import format",
                                                                      choices=c("TXT","XLS","XLSX","SPSS_.sav","MICROSOFT_OFFICE_ACCES_.mdb","STATA","CSV","R_.Rda","SAS_.sas7bdat","SAS_.xport"))),
                 
                 
                 
                 ########### Import arguments ################################
                 
                 conditionalPanel(condition= "output.format=='SPSS_.sav'",
                                  actionButton("import_options_spss",label="Show import options"),
                                  bsModal("modal_SPSS","Import options","import_options_spss",
                                          checkboxInput("button_spss_import","replace values with labels"),
                                          textInput("date_spss_import","Dates format",value="d-m-yy"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_spss","introduce manually the arguments inside spss_import()"),
                                          a("spss_import help", href="spss_import.html"))),
                 
                 conditionalPanel(condition= "output.format=='TXT' | output.format=='CSV'",
                                  actionButton("import_options_table",label="Show import options"),
                                  bsModal("modal_table","Import options","import_options_table",
                                          checkboxInput("button_txt_import","Header"),
                                          textInput("dec_txt_import","Dec",value="."),
                                          textInput("sep_txt_import","Field separator character (leave as false to perform automatic detection)",value=F),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_table","introduce manually the arguments inside table_import()"),
                                          a("table_import help", href="table_import.html"))),

                 
 
                conditionalPanel(condition= "output.format=='XLS' | output.format=='XLSX'",
                                 actionButton("import_options_excel",label="Show import options"),
                                 bsModal("modal_excel","Import options","import_options_excel",
                                         checkboxInput("button_excel_import","Header",value=T),
                                         numericInput("sheetIndex_excel_import","sheet Index",value=1),
                                         textInput("sheetName_excel_import","sheet name",value=NULL),
                                         helpText("---------------------------------------------------------------------------------------------------------------"),
                                         helpText("The use of the following command will disable selected options, including the file select."),
                                         textInput("arguments_excel","introduce manually the arguments inside read.xlsx()"))),

     
                 conditionalPanel(condition= "output.format=='MICROSOFT_OFFICE_ACCES_.mdb'",
                                  actionButton("import_options_acces",label="Show import options"),

                                  bsModal("modal_acces","Import options","import_options_acces",
                                          textInput("uid_access_import","UID"),
                                          textInput("pwd_access_import","PWD(password)"),
                                          uiOutput("listtables") ,
                                          textInput("tableNames_access_import","SQL query"),
                                          textInput("date_access_import","Dates format",value="d-m-yy"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_acces","introduce manually the arguments inside access_import()"),
                                          a("access_import help", href="access_import.html"))),
                 
                 conditionalPanel(condition= "output.format=='STATA'",
                                  actionButton("import_options_stata",label="Show import options"),
                                  bsModal("modal_stata","Import options","import_options_stata",
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_stata","introduce manually the arguments inside read_dta()"))),
                 
                 conditionalPanel(condition= "output.format=='SAS_.sas7bdat'",
                                  actionButton("import_options_SAS",label="Show import options"),
                                  bsModal("modal_SAS","Import options","import_options_SAS",
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_SAS","introduce manually the arguments inside read_sas()"))),
                 
                 conditionalPanel(condition= "output.format=='SAS_.xport'",
                                  actionButton("import_options_SASxport",label="Show import options"),
                                  bsModal("modal_SASxport","Import options","import_options_SASxport",
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options, including the file select."),
                                          textInput("arguments_SASxport","introduce manually the arguments inside read.xport() )"))),
                  
                 ######################################################################################
                 
                 
                  helpText(""),actionButton("goButton", "Import")   ),
                 
                 bsCollapsePanel("Step 2) Database checkup",style="primary",
                 
                 a("What's format corrector?", href="format_corrector.html") ,
                 checkboxInput("format_corrector","Run format_corrector",value=F),
                 
                 uiOutput("listvars") ,
                 checkboxInput("allvar","Select all",value=T)  ),

                 
      
      
      
      
      
                 bsCollapsePanel("Step 3) Export ",style="primary",
                                 
                 selectInput("outputFormat","Select export format",
                            choices=c("TXT","XLS","XLSX","SPSS_.sav","MICROSOFT_OFFICE_ACCES_.mdb","STATA","CSV","R_.Rda","SAS_.xport")),
                 
                 conditionalPanel(condition="input.outputFormat != 'MICROSOFT_OFFICE_ACCES_.mdb' & input.outputFormat != 'SPSS_.sav'",
                                  textInput("fileName","Output file name")),
      
                
                 

                 

                 ########### export arguments #######################################
                 conditionalPanel(condition="input.outputFormat=='TXT'",
                                  actionButton("export_options_txt",label="Show export options"),
                                  bsModal("modal_txt_export","Export options","export_options_txt",
                                          textInput("sep_txt_export","Field separator character",value=" "),
                                          textInput("dec_txt_export","Dec",value="."),
                                          textInput("NA_txt_export","NA writing value",value="NA"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_TXT_export","introduce manually the arguments inside write.table()"))),
                 
                 conditionalPanel(condition="input.outputFormat=='CSV'",
                                  actionButton("export_options_csv",label="Show export options"),
                                  bsModal("modal_csv_export","Export options","export_options_csv",
                                          textInput("dec_csv_export","dec",value="."),
                                          textInput("NA_csv_export","NA writing value",value="NA"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_csv_export","introduce manually the arguments inside write.csv()"))),
                 
                 conditionalPanel(condition="input.outputFormat=='SPSS_.sav'",
                                  actionButton("export_options_SPSS",label="Show export options"),
                                  bsModal("modal_SPSS_export","Export options","export_options_SPSS",
                                          textInput("final_spss_export","SPSS file name",value="example.sav"),
                                          textInput("data_spss_export","txt value file name",value="data.txt"),
                                          textInput("var.keep_spss","variables to save",value="ALL"),
                                          textInput("filerunsyntax","Location of runsyntax.exe or pspp.exe",value="C:/Program files/SPSS/runsyntx.exe"),
                                          checkboxInput("run_spss","Run spss and create the file",value=T),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_SPSS_export","introduce manually the arguments inside spss_export()"),
                                          a("spss_export help", href="spss_export.html"))),
                 
                 conditionalPanel(condition="input.outputFormat=='XLSX' | input.outputFormat=='XLS'",
                                  actionButton("export_options_excel",label="Show export options"),
                                  bsModal("modal_excel_export","Export options","export_options_excel",
                                          checkboxInput("varview_excel","save var_view in another sheet?"),
                                          textInput("tablename_excel","sheet name",value="sheet1"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_excel_export","introduce manually the arguments inside excel_export()"),
                                          a("excel_export help", href="excel_export.html"))),
                 
                 conditionalPanel(condition="input.outputFormat=='MICROSOFT_OFFICE_ACCES_.mdb'",
                                  actionButton("export_options_acces",label="Show export options"),
                                  bsModal("modal_acces_export","Export options","export_options_acces",
                                        
                                          textInput("uid_access_export","UID"),
                                          textInput("pwd_access_export","PWD(password)"),
                                          textInput("tableNames_access_export","table name"),
                                          checkboxInput("varview_access","save var_view in another sheet?"),
                                          helpText("The acces file to save the data will be specified after clicking 'add table'."),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_access_export","introduce manually the arguments inside access_export()"))),
                 

                 conditionalPanel(condition="input.outputFormat=='SAS_.xport'",
                                  actionButton("export_options_SAS",label="Show export options"),
                                  bsModal("modal_SAS_export","Export options","export_options_SAS",
                                          textInput("codefile","Codefile name"),
                                          helpText("---------------------------------------------------------------------------------------------------------------"),
                                          helpText("The use of the following command will disable selected options."),
                                          textInput("arguments_SAS_export","introduce manually the arguments inside write.foreign()"))),
                 ######################################################################################
                 
                  conditionalPanel(condition="input.outputFormat != 'MICROSOFT_OFFICE_ACCES_.mdb' & input.outputFormat != 'SPSS_.sav'",
                                   helpText(""),downloadButton("file","Download file")),
                 conditionalPanel(condition="input.outputFormat == 'MICROSOFT_OFFICE_ACCES_.mdb'",
                                  helpText(""),actionButton("access_export_button","Add table")),
                  conditionalPanel(condition="input.outputFormat == 'SPSS_.sav'",
                                   helpText(""),actionButton("SPSS_export_button","Download files"))
                 )   ) ),
    
    column(8,navbarPage("Tables",theme=shinytheme("flatly"),
      tabPanel("file",HTML(paste("<strong>Current format: </strong>",textOutput("format"),
               HTML("<p><strong>Selected variables:</strong>"),textOutput("variables"),
               HTML("<br><strong>Unselected variables:</strong>"),textOutput("unselected_variables"),
               HTML("<br><strong>Format_corrector implemented</strong>"),textOutput("format_corrector2") , sep=" "))),
              tabPanel( "var_view",tableOutput("var_view"))   ,
              # tabPanel("head",tableOutput("head") ),
              tabPanel("datatable",dataTableOutput("datatable")),
              tabPanel("summary",
                       radioButtons("summary_radiobuttons",label=NULL,
                                               choices = list("Display mean and SD" =0, "Display median and quartiles" = 1), 
                                               selected = 0,inline=T),
                       verbatimTextOutput("summary"),verbatimTextOutput("datesummary")),
              tabPanel("App guide",includeHTML("www/App-guide.html"))    
              
      )
      ) 
    
           
  
  )
  
))
