library(shiny)
# runApp(".", launch.browser = TRUE)

# Define UI for STSRE application
shinyUI(pageWithSidebar(
  
  # STSRE title
  headerPanel(
    img(src="logo.png", height=600, width=600),
    tags$head(tags$style("body {background-color: black;}"))),
  
  # Sidebar with tabs
  sidebarPanel(
    tabsetPanel(
      tabPanel(
        headerPanel(img(src="b0.png", height=35, width=35)),
        wellPanel(
          p(strong("The Sea Turtle Sex Ratio Estimator (seaturtleSRE or STSRE) is a new tool to easily estimate sea turtle hatchling sex ratio (Primary Sex Ratio, from here on ‘SR’) produced in nests from various beaches, based on their incubation temperature.")),
          tags$hr(),
          p("STSRE has been first developped by Maria Sousa Martins and is now co-developped by her and ",a("Marc Girondot", 
                                                    href="http://max2.ese.u-psud.fr/epc/conservation/Girondot/Publications/Marc.html",
                                                    target="_blank")),
          tags$hr(),
          p("STSRE will allow users from all over the world to easily calculate SR values for specific work (Individual Use), and also help to build (and continuously complete and update) a global database of theoretical sea turtle SRs (Global Use). Ultimately, this will help scientists to evaluate the vulnerability found in sea turtle populations. Thus, STSRE's main objective is to help prioritize conservation areas in relation to climate change impacts, particularly global warming.")
          ),
        
        value=0),
      
      
      
      tabPanel(
        headerPanel(img(src="b1.png", height=35, width=35)),
        wellPanel(
          p(strong("SPECIES")),
          radioButtons("species", "Which species are you working
                       with?",
                       list("Caretta caretta"="Cc", 
                            "Chelonia mydas"="Cm", 
                            "Dermochelys coriacea"="Dc",
                            "Eretmochelys imbricata"="Ei",
                            "Lepidochelys kempii"="Lk",
                            "Lepidochelys olivacea"="Lo",
                            "Natator depressus"="Nd")),
          numericInput("hatchling.mean", "Mean size of hatchlings:",
                       value=39, min=20, max=60, step=0.01),
          numericInput("hatchling.sd", "Standard deviation of size of hatchlings:",
                       value=3, min=0.1, max=6, step=0.01)
        ),
        value=1),
      
      
      
      tabPanel(
        headerPanel(img(src="b2.png", height=35, width=35)),
        
        wellPanel(
          p(strong("NESTING AREA")),
          radioButtons("country_type", "Where do your turtles nest?",
                       list("List", "Coordinates")),
          
          conditionalPanel(
            condition="input.country_type=='List'",
            uiOutput("countrySelect")),
          
          conditionalPanel(
            condition="input.country_type=='Coordinates'",
            numericInput("longitude", "Longitude:",
                         value=00, min=-180, max=180, step=1),
            numericInput("latitude", "Latitude:",
                         value=00, min=-90, max=90, step=1)),
          
          uiOutput("rmuSelect"),
          
          helpText("If you don't know what is your turtle's RMU, please check",
                   a(em("Wallace et al., 2010."), 
                     href="http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0015465#pone-0015465-g007",
                     target="_blank")),
          
          sliderInput(inputId="zoom", 
                      label="Zoom:", 
                      min=0.1, max=18, value=18, step=1),
          tags$hr()),
        
        value=2),
      
      
      
      tabPanel(
        headerPanel(img(src="b3.png", height=35, width=15)),
        wellPanel(
          p(strong("FIELD NEST THERMAL DATA")),
          fileInput("file1", "Please input a Incubation Temperature Profile", accept=c('text/csv', 
                    'text/comma-separated-values,text/plain', '.csv')),
          uiOutput("results"),
          helpText("Note: The input file must be as .csv. The first column must be the date-time and the second be the temperature in °C"),
          tags$hr(),
          checkboxInput('header', 'Header', TRUE),
          radioButtons('sep', 'Separator',
                       c(Comma=',',
                         Semicolon=';',
                         Tab='\t'),
                       'Semicolon'),
          radioButtons('quote', 'Quote',
                       c(None='',
                         'Double Quote'='"',
                         'Single Quote'="'"),
                       'Double Quote'),
          radioButtons('dateformat', 'Date and Time format',
                       c("DD/MM/YYYY HH:MM:SS"='%d/%m/%Y %H:%M:%S',
                         "DD/MM/YY HH:MM:SS"='%d/%m/%y %H:%M:%S',
                         "YYYY-MM-DD HH:MM:SS"='%Y-%m-%d %H:%M:%S',
                         "YY-MM-DD HH:MM:SS"='%y-%m-%d %H:%M:%S'),
                       'DD/MM/YYYY HH:MM:SS')
          
        ),      
        value=3),
      
      
      
      tabPanel(
        headerPanel(img(src="b4.png", height=35, width=35)),
        wellPanel(
          p(strong("CONSTANT TEMPERATURE INCUBATION DATA")),
          p(strong(em("Caretta caretta"))),
          mapply(function(ref) {checkboxInput(inputId=paste0("c",ref), 
                                              label=paste0("RMU=",as.character(STSRE_TSD$RMU[which(STSRE_TSD$Ref==ref)[1]])," in ", as.character(STSRE_TSD$Reference[which(STSRE_TSD$Ref==ref)[1]])) , 
                                              value=FALSE)}, 
                 ref=levels(as.factor(as.character(STSRE_TSD$Ref[which(STSRE_TSD$Sp=="Cc")]))), SIMPLIFY=FALSE),
          p(strong(em("Chelonia mydas"))),
          mapply(function(ref) {checkboxInput(inputId=paste0("c",ref), 
                                              label=paste0("RMU=",as.character(STSRE_TSD$RMU[which(STSRE_TSD$Ref==ref)[1]])," in ", as.character(STSRE_TSD$Reference[which(STSRE_TSD$Ref==ref)[1]])) , 
                                              value=FALSE)}, 
                 ref=levels(as.factor(as.character(STSRE_TSD$Ref[which(STSRE_TSD$Sp=="Cm")]))), SIMPLIFY=FALSE),
          p(strong(em("Dermochelys coriacea"))),
          mapply(function(ref) {checkboxInput(inputId=paste0("c",ref), 
                                              label=paste0("RMU=",as.character(STSRE_TSD$RMU[which(STSRE_TSD$Ref==ref)[1]])," in ", as.character(STSRE_TSD$Reference[which(STSRE_TSD$Ref==ref)[1]])) , 
                                              value=FALSE)}, 
                 ref=levels(as.factor(as.character(STSRE_TSD$Ref[which(STSRE_TSD$Sp=="Dc")]))), SIMPLIFY=FALSE),
          p(strong(em("Eretmochelys imbricata"))),
          mapply(function(ref) {checkboxInput(inputId=paste0("c",ref), 
                                              label=paste0("RMU=",as.character(STSRE_TSD$RMU[which(STSRE_TSD$Ref==ref)[1]])," in ", as.character(STSRE_TSD$Reference[which(STSRE_TSD$Ref==ref)[1]])) , 
                                              value=FALSE)}, 
                 ref=levels(as.factor(as.character(STSRE_TSD$Ref[which(STSRE_TSD$Sp=="Ei")]))), SIMPLIFY=FALSE),
          p(strong(em("Lepidochelys olivacea"))),
          mapply(function(ref) {checkboxInput(inputId=paste0("c",ref), 
                                              label=paste0("RMU=",as.character(STSRE_TSD$RMU[which(STSRE_TSD$Ref==ref)[1]])," in ", as.character(STSRE_TSD$Reference[which(STSRE_TSD$Ref==ref)[1]])) , 
                                              value=FALSE)}, 
                 ref=levels(as.factor(as.character(STSRE_TSD$Ref[which(STSRE_TSD$Sp=="Lo")]))), SIMPLIFY=FALSE)
        ),        
        value=4),
      
      
      
      tabPanel(
        headerPanel(img(src="b5.png", height=35, width=35)),
        wellPanel(
          p(strong("RESULTS"))),
        
        value=5),
      
      id="conditionedPanels")),
  
  
  # Mainbar with tabs
  mainPanel(
    conditionalPanel(
      condition="input.conditionedPanels==0",
      helpText("STSRE"),
      img(src="00.jpg", height=600, width=600)),
    
    
    
    conditionalPanel(
      condition="input.conditionedPanels==1",
      helpText("Species"),
      
      conditionalPanel(
        condition="input.species=='Cc'",
        img(src="Cc.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Cm'",
        img(src="Cm.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Dc'",
        img(src="Dc.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Ei'",
        img(src="Ei.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Lk'",
        img(src="Lk.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Lo'",
        img(src="Lo.jpg", height=400, width=400)),
      
      conditionalPanel(
        condition="input.species=='Nd'",
        img(src="Nd.jpg", height=400, width=400))),
    
    
    
    conditionalPanel(
      condition="input.conditionedPanels==2",
      helpText("Nesting Area"),
      
      p(strong("Countries:")),
      verbatimTextOutput(outputId="NestingAreaMap"),
      helpText("The names of the countries reported are the names of the countries closest to the nearest border"),
      
      plotOutput(outputId="MAP", width="800px", height="600px", clickId="click")
      ),
    
    conditionalPanel(
      condition="input.conditionedPanels==3",
      verbatimTextOutput(outputId="resultsInfo"),
      plotOutput(outputId="resultPlot", width="400px", height="400px")
    ),
    
    conditionalPanel(
      condition="input.conditionedPanels==4",
      verbatimTextOutput(outputId="TSDInfo"),
      plotOutput(outputId="TSDPlot", width="600px", height="400px")
    )
    
    )))