library(shiny)
shrinkLabel <- "Shrinkage Type:"
shrinkChoices <- list("Hard" = "Hard", "Soft" = "Soft", "Garrote" = "Garrote")
# Define UI
sideWidth = 3
mainWidth = 12 - sideWidth
shinyUI(navbarPage(
  # Application title
  "mWaveD Interactive Demo",
#   "Multichannel wavelet deconvolution with additive long memory errors",
  # First Table Panel
  tabPanel("Multichannel Signal",
    sidebarLayout(
      sidebarPanel(
        selectInput("signalShow", "Show Signal Details:",
                   list("Underlying signal" = 1,
                        "Blurred signal" = 2,
                        "Noisy & blurred signal" = 3), selected=3),
        selectInput("sig", "Underlying Signal:",
                  list("Blocks" = "blocks",
                       "Bumps" = "bumps",
                       "Cusp" = "cusp",
                       "Doppler" = "doppler",
                       "HeaviSine" = "heavisine",
                       "LIDAR" = "lidar")),
        selectInput("blur", "Type of blur:",
                  list("No blur (direct)" = "direct",
                       "Smooth (Gamma)" = "smooth",
                       "Box Car" = "box.car")),
        conditionalPanel(condition="input.blur=='smooth'",
                       sliderInput("gshape", "Shape parameter range:",
                                   min = 0.25, max = 1.5, value = c(0.5,0.7)),
                       sliderInput("gscale", "Scale parameter range:",
                                   min = 0.25, max = 0.5, value = c(0.25,0.5))),
        sliderInput("m","m = Number of channels",1,14,1,locale='au'),
        sliderInput("alpha", "alpha = dependence parameter:",0.1,1, value = c(0.85,0.99)),
        sliderInput("SNR","SNR = Signal to Noise Ratio (dB):",5,30,value=c(25,30),locale='au'),
        sliderInput("J","n = 2^J : # obs per channel:",8,14,12,locale='au'), 
        verbatimTextOutput("summarySignal"),
        verbatimTextOutput("signalCalls"), width = sideWidth),
      mainPanel(
        plotOutput("reactiveSignal", height = "auto"))), icon = icon("signal")),
  # Second Tab Panel
  tabPanel("Resolution selection",
      sidebarLayout(
        sidebarPanel(
          selectInput("resolution", "Finest resolution selection method:",
                      list("Smooth stopping rule" = "smooth",
                           "Blockwise variance" = "block")),
        verbatimTextOutput("summaryResolution"),
        verbatimTextOutput("resolutionCalls"), width = sideWidth),
      mainPanel(
        plotOutput("resolutionPlot", height = "auto"), width = mainWidth)), icon = icon("check-square-o")),
  # Third tab panel
  tabPanel("Re-constructed Signal",
    sidebarLayout(
      sidebarPanel(
        selectInput("shrinkage1", label = shrinkLabel, choices = shrinkChoices, selected = 'hard'),
        radioButtons("etaChoose", "Smoothing parameter:",
                     c("Default" = "FALSE",
                       "Custom" = "TRUE")),
        conditionalPanel(condition="input.etaChoose=='TRUE'",
                         uiOutput("etaSlider")),
        selectInput("degree", label = 'Degree of Meyer Wavelet polynomial:', choices = 0:4, selected = 3),
        verbatimTextOutput("summaryWVD"),
        conditionalPanel(condition="input.m > 1",
          radioButtons("wvdshow", "Show alternative estimates:",
            list("Using only best channel" = 1,
                  "Using naive mean" = 2,
                  "No alternatives" = 3), selected = 3)),
        verbatimTextOutput("mWaveDCalls"), width = sideWidth),
      mainPanel(
        plotOutput("wvdPlot", height = "auto"), width = mainWidth)), icon = icon("pencil-square-o")),
  # Fourth tab panel
  tabPanel("Multiresolution Analysis",
    sidebarLayout(
      sidebarPanel(
        selectInput("shrinkage2", label = shrinkLabel, choices = shrinkChoices, selected = "Hard"),
        verbatimTextOutput("summaryMRA"),
        verbatimTextOutput("mraCalls"), width = sideWidth),
      mainPanel(plotOutput("multiPlot", height = "auto"), width = mainWidth)), icon = icon("search-plus"))
))