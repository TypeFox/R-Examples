
shinyUI(fluidPage(
  titlePanel(title="EffectLiteR"),
  
  sidebarLayout(
  sidebarPanel(
    tabsetPanel(
      ######### Data ############
      tabPanel('Data',
        img(src='effectliter_logo.png', align = "right"),
        br(),
        actionButton("newanalysis","Start a New Analysis"),
        uiOutput("reload"),
        hr(),        
        selectizeInput(inputId="exdata", label="Example Data", selected="",
                       choices= c("","nonortho","example01","example02lv",
                                  "example_multilevel"),
                       options = list(placeholder = 'choose example data'),
                       width='50%'),    
        # hr(),
        tryCatch(
          fileInput("file1", "Data File", 
                    accept=c(".csv", ".txt", ".sav", ".xpt", 
                              ".CSV", ".TXT", ".SAV", ".XPT",
                             ".DAT", ".dat"))
        ),
        helpText('Select either a .csv, .dat, .txt, .sav or a .xpt file to be uploaded. The corresponding R function (read.csv, read.table, read.spss, or read.xport) will be chosen automatically with the default settings for arguments. Some default arguments can be overwritten (see additional options below).'),
        br(),
        h5(strong("Additional Options to Read Data")),
        br(),
        selectizeInput(inputId="vallabels", 
                       label=h5("Use value labels (SPSS data)"), 
                       selected="default",
                       choices= c("default","yes","no"),
                       width='90%'),
        selectizeInput(inputId="header", 
                       label=h5("File contains variable names (csv, dat, and txt data)"), 
                       selected="default",
                       choices= c("default","yes","no"),
                       width='90%'),
        selectizeInput(inputId="sep", 
                       label=h5("Character separating columns (csv, dat, and txt data)"), 
                       selected="default",
                       choices= c("default","semicolon","white space"),
                       width='90%'),
        selectizeInput(inputId="dec", 
                       label=h5("Decimal character (csv, dat, and txt data)"), 
                       selected="default",
                       choices= c("default","decimal point","decimal comma"),
                       width='90%'),
        textInput(inputId="na.strings", 
                  label=h5("Missing value code (csv, dat, and txt data)"), 
                  value = "NA",
                  width='90%'),
        br()
      ),
      ########## Manifest Variables ############
      tabPanel('Manifest Variables',
               img(src='effectliter_logo.png', align = "right"),
               br(),
               conditionalPanel(
                  condition = "!input.latenty",
                  selectizeInput(
                    inputId="variabley", 
                    label="Dependent Variable Y", 
                    choices="",
                    options = list(placeholder = 'select dependent variable'),
                    width='60%')
               ),
               selectizeInput(
                 inputId="variablex", 
                 label="Treatment Variable X", 
                 choices="",
                 options = list(placeholder = 'select treatment variable'),
                 width='60%'),
               selectizeInput("variablek", "Categorical Covariates K", "", 
                   multiple=TRUE, selected="",
                   options = list(placeholder = 'select categorical covariates')),
               selectizeInput("variablez", "Continuous Covariates Z", "", 
                              multiple=TRUE, selected="",
                    options = list(placeholder = 'select continuous covariates')),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
br(),
br()
      ),
      ############ Latent Variables ###########
      tabPanel('Latent Variables',
        img(src='effectliter_logo.png', align = "right"),       
        br(),       
        checkboxInput("latenty", "Latent dependent variable", FALSE),        
        conditionalPanel(
          condition = "input.latenty",
          selectizeInput("indicatorsy", 
                        "Indicators of Latent Dependent Variable", 
                        "", multiple=TRUE,
                        options = list(placeholder = 'select indicators'))
        ),
        checkboxInput("latentz", "Add latent covariates", FALSE),
        conditionalPanel(
          condition = "input.latentz",
          ###
          ### Indicators for latent covariates
          column(8, p("Number of Latent Covariates")),
          column(4, numericInput("nlatentz", label=NULL, 0,
                       min = 0, max = 5, width='60%')),
          ### Cov 1
          conditionalPanel(
            condition = "input.nlatentz > 0",
            selectizeInput("indicatorsz1", 
                           "Indicators of Latent Covariate 1", 
                           "", multiple=TRUE, 
                           options = list(placeholder = 'select indicators'))
          ),
          ### Cov 2
          conditionalPanel(
            condition = "input.nlatentz > 1",
            selectizeInput("indicatorsz2", 
                           "Indicators of Latent Covariate 2", 
                           "", multiple=TRUE, 
                           options = list(placeholder = 'select indicators'))
          ),  
          ### Cov 3
          conditionalPanel(
            condition = "input.nlatentz > 2",
            selectizeInput("indicatorsz3", 
                           "Indicators of Latent Covariate 3", 
                           "", multiple=TRUE, 
                           options = list(placeholder = 'select indicators'))
          ),  
          ### Cov 4
          conditionalPanel(
            condition = "input.nlatentz > 3",
            selectizeInput("indicatorsz4", 
                           "Indicators of Latent Covariate 4", 
                           "", multiple=TRUE, 
                           options = list(placeholder = 'select indicators'))
          ),  
          ### Cov 5
          conditionalPanel(
            condition = "input.nlatentz > 4",
            selectizeInput("indicatorsz5", 
                           "Indicators of Latent Covariate 5", 
                           "", multiple=TRUE, 
                           options = list(placeholder = 'select indicators'))
          ),
          helpText("Latent covariates will be added to the list of continuous manifest covariates (if specified).")
        ),
        hr(),
        p("Additional Options"),
        ###
        ### Names of latent variables
        checkboxInput("latnames", "Names of latent variables", FALSE),        
        conditionalPanel(
          condition = "input.latnames",
          ### Dependent Variable
          conditionalPanel(
            condition = "input.latenty",
            textInput("name.etay", "Latent Dependent Variable", 
                      value = "eta",
                      width='60%')
          ),
          ### Cov 1
          conditionalPanel(
            condition = "input.nlatentz > 0",
            textInput("name.etaz1", "Latent Covariate 1", 
                      value = "xi1",
                      width='60%')
          ),
          ### Cov 2
          conditionalPanel(
            condition = "input.nlatentz > 1",
            textInput("name.etaz2", "Latent Covariate 2", 
                      value = "xi2",
                      width='60%')
          ),  
          ### Cov 3
          conditionalPanel(
            condition = "input.nlatentz > 2",
            textInput("name.etaz3", "Latent Covariate 3", 
                      value = "xi3",
                      width='60%')
          ),  
          ### Cov 4
          conditionalPanel(
            condition = "input.nlatentz > 3",
            textInput("name.etaz4", "Latent Covariate 4", 
                      value = "xi4",
                      width='60%')
          ),  
          ### Cov 5
          conditionalPanel(
            condition = "input.nlatentz > 4",
            textInput("name.etaz5", "Latent Covariate 5", 
                      value = "xi5",
                      width='60%')
          )
        ),
        ###
        ### Customized Measurement Models
        checkboxInput("custommeasmodels", "Customize measurement models", FALSE),      
        conditionalPanel(
          condition = "input.custommeasmodels",
          ### Dependent Variable
          conditionalPanel(
            condition = "input.latenty",
            selectizeInput(
              inputId="mm.etay", 
              label="Measurement Model for Latent Dependent Variable", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          ),
          ### Cov 1
          conditionalPanel(
            condition = "input.nlatentz > 0",
            selectizeInput(
              inputId="mm.etaz1", 
              label="Measurement Model for Latent Covariate 1", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          ),
          ### Cov 2
          conditionalPanel(
            condition = "input.nlatentz > 1",
            selectizeInput(
              inputId="mm.etaz2", 
              label="Measurement Model for Latent Covariate 2", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          ),  
          ### Cov 3
          conditionalPanel(
            condition = "input.nlatentz > 2",
            selectizeInput(
              inputId="mm.etaz3", 
              label="Measurement Model for Latent Covariate 3", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          ),  
          ### Cov 4
          conditionalPanel(
            condition = "input.nlatentz > 3",
            selectizeInput(
              inputId="mm.etaz4", 
              label="Measurement Model for Latent Covariate 4", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          ),  
          ### Cov 5
          conditionalPanel(
            condition = "input.nlatentz > 4",
            selectizeInput(
              inputId="mm.etaz5", 
              label="Measurement Model for Latent Covariate 5", 
              choices = c("default" = "default",
                          "equivalent measures" = "parallel",
                          "essentially equivalent measures" = "tau-equi",
                          "congeneric measures" = "tau-cong"),
              width='100%')
          )
        ),
        helpText("Names of latent variables and customize measurement models are optional. If not specified, EffectLiteR picks default names and tries to guess a reasonable measurement model: Congeneric for latent variables with three or more indicators, essentially tau-equivalent for latent variables with less than three indicators and for latent variables with cross-loadings (e.g., method factors), and parallel for single-indicator latent variables.")
      ),
      ########## Options #############
      tabPanel('Options',
        img(src='effectliter_logo.png', align = "right"),
        br(),
        selectInput("control", "Control Group", "", width='60%'),
        radioButtons("missing", "Missing Data", 
                     choices=c("listwise","fiml"), 
                     selected = "listwise"),
        radioButtons("fixed.cell", "Sampling Model", 
                     choices=c("stochastic"="stochastic",
                               "fixed cell sizes"="fixed",
                               "fixed cell sizes and fixed means of Z"="fixed+e"), 
                     selected = "stochastic"),
        radioButtons("se", "Standard Errors", 
                     choices=c("standard","boot","first.order",
                               "robust.sem","robust.huber.white"),
                     selected = "standard"),
        conditionalPanel(
          condition = "input.se == 'boot'",
          numericInput("bootstrap", "Number of bootstrap draws", 
                      min=1, max=5000, value=5)
        ),
        br(),
        h5("Additional Options"),
        checkboxInput("homoscedasticity", "Homoscedastic residual variances", 
                      value=FALSE),
        br(),
        h5("Add to lavaan Syntax"),
        tags$textarea(id="add.syntax", rows=5, cols=40, ""),
        helpText('This option is for advanced users. The text will be appended to the lavaan syntax generated by EffectLiteR. It can for example be used to add additional (in-)equality constraints or to compute user-defined conditional effects.')
        
    ),
########## Interactions #############
    tabPanel('Interactions',
         img(src='effectliter_logo.png', align = "right"),
         br(),
         radioButtons("interactions", "Interactions", 
                      choices=c("Full model"="all",
                                "Only two-way interactions"="2-way",
                                "Only X:K interactions"="X:K",
                                "Only X:Z interactions"="X:Z",                                
                                "No treatment*covariate interactions"="none"), 
                      selected = "all",
                      width='60%'),
         br(),
         helpText('These are some pre-defined sets of constraints on interactions. You can consult the regression model in the main output to see which coefficients are fixed to zero. In addition, you can use the additional options input to specify customized constraints.')      
      ),
########## Propensity Scores #############
      tabPanel('Propensity Scores',
         img(src='effectliter_logo.png', align = "right"),
         br(),
         checkboxInput(
           inputId="propscoreformula", 
           label="Specify formula for propensity score model yourself", 
           value=FALSE,
           width='60%'),
         conditionalPanel(
           condition = "!input.propscoreformula",
           br(),
           selectizeInput("propscore", "Covariates in Propensity Score Model", "",
                          multiple=TRUE, selected="",
                          options = list(placeholder = 'select covariates'))
         ),
         conditionalPanel(
           condition = "input.propscoreformula",
           br(),           
           textInput("prop.formula", "Specify formula", 
                     value = "")
         ),
         br(),
         helpText('Propensity scores are predicted probabilities from a multinomial regression of the treatment variable on the above selected covariates. The logit transformed propensity score(s) are included as continuous covariates in the EffectLiteR analysis. For more flexibility, you can specify the R formula yourself.')
      ),

########## Complex Survey #############
tabPanel('Complex Survey',
         img(src='effectliter_logo.png', align = "right"),
         br(),
         selectizeInput("ids", "Cluster Variable", "",
                        multiple=FALSE, selected="",
                        options = list(placeholder = 'select cluster ID'),
                        width='60%'),
         selectizeInput("weights", "Sampling Weights", "",
                        multiple=FALSE, selected="",
                        options = list(placeholder = 'select sampling weights'),
                        width='60%'),
         helpText('The lavaan.survey package is called to account for cluster variables and sampling weights. This is an experimental feature.'),
         helpText('Note: Only use weights if you know what you are doing. For example, some conditional treatment effects may require different weights than average effects.')
         
      )
  )),
  
  mainPanel(
    tabsetPanel(
      ######### Data Table ##########
      tabPanel('Data', dataTableOutput("mytable1")),
      
      ######### EffectLiteR ##########
      tabPanel("EffectLiteR", verbatimTextOutput("summary")),
      
      ######### lavaan Syntax ##########
      tabPanel("lavaan Syntax",
               br(),
               downloadLink('downloadLavData', 'Download Data for lavaan'),
               br(),
               br(),
               verbatimTextOutput("elrcall"),
               verbatimTextOutput("lavcall"),
               verbatimTextOutput("lavsyntax")),
      
      ######### lavaan Results ##########
      tabPanel("lavaan Results", verbatimTextOutput("lavresults")),
      
      ######### Conditional Effects I ##########
      tabPanel('Conditional Effects I', 
               verbatimTextOutput("helptextcondeffects"),
               downloadLink('downloadConditionalEffects', 'Download Conditional Effects Data'),
               br(),
               br(),
               dataTableOutput("condeffs")),
      
      ######### Conditional Effects II ##########
      tabPanel("Conditional Effects II", 
          verbatimTextOutput("helptextcondeffects2"),
          br(),
          br(),
          column(3, wellPanel(
            h5("Values of Covariates"),
            uiOutput("ui")
          )) ,
          column(9, wellPanel(
            h5("Conditional Effects"),
            verbatimTextOutput("outputcondeffect2")
          ))
               
      ),
      
      
      ######### Plot 1 ##########
      tabPanel("Plot 1", 
               verbatimTextOutput("helptextplot1"),
               plotOutput("plot1")
      ),
      
      ######### Plot 2 ##########
      tabPanel("Plot 2", 
        verbatimTextOutput("helptextplot2"),
        selectizeInput(
          inputId="zselect", 
          label="Continuous Covariate", 
          choices="",
          options = list(placeholder = 'select continuous covariate'),
          width='30%'),
      plotOutput("plot2")),
      
      ######### Plot 3 ##########
      tabPanel("Plot 3",
          verbatimTextOutput("helptextplot3"),
          column(3, wellPanel(
            selectInput("gxselect", 
                        "Effect function", 
                        "g1(K,Z)", 
                        multiple=FALSE,
                        selectize=TRUE),
            selectInput("zselect2", 
                        "Regressor", 
                        "", 
                        multiple=FALSE,
                        selectize=TRUE),
            selectInput("zselect3", 
                        "Colour variable", 
                        "", 
                        multiple=FALSE,
                        selectize=TRUE)
          )) ,
          column(9, 
                 plotOutput("plot3"))
          )
      
    )    
  ))

))