library(shiny)
library(shinyjs)
library(rglwidget)

shinyUI(
  navbarPage('AFM Image Analysis',
                   tabPanel('File',
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput('choose_inputtype'),
                                uiOutput('choose_type'),
                                uiOutput('choose_dataset'),
                                tags$hr(),
                                htmlOutput("displayIn3DFileButton"),
                                #actionButton('displayIn3DFileButton', label = 'Display 3D model'),
                                downloadButton('saveRdataFileButton', label = 'Export calculations')
                              ),
                              mainPanel(
                                uiOutput('imageInformationsUI'),
                                tableOutput('basicInfoFileTable'),
                                uiOutput('roughnessUI'),
                                tableOutput('roughnessesFileTable')
                                # ,tags$script(src="AFM.js")
                              )
                            )),
                   tabPanel('PSD',
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput('imageNamePSD'),
                                tags$hr(),
                                sliderInput('breaksSliderPSD', label = 'Breaks in PSD 2D to calculate PSD 1D',
                                            min = 1, max = 7, value = 5, step=1, ticks=FALSE),
                                actionButton('RoughnessByLengthScaleButton', label = 'Calculate'),
                                tags$hr(),
                                uiOutput('downloadPSDPSDButton'),
                                uiOutput('downloadRoughnessVsLengthscalePSDButton')
                              )
                              ,
                              mainPanel(
                                uiOutput('plotPSDUI'),
                                plotOutput('plotPSD'),
                                uiOutput('plotPSDRvsLUI'),
                                plotOutput('plotPSDRvsL')
                              )
                            )),
                   
                   navbarMenu('Variance',
                              tabPanel('Checks',
                                       sidebarLayout(
                                         sidebarPanel(
                                           uiOutput('imageNameCheck'),
                                           tags$hr(),
                                           sliderInput('sampleIsotropyVarianceCheckSlider', label = 'Sample to calculate directional variograms',
                                                       min = 1, max = 100, value = 100, step=1),
                                           actionButton('checkNormalityIsotropyCheckButton',label='Check normality and isotropy')
                                         )
                                         ,
                                         mainPanel(
                                           uiOutput('normalityVarianceCheckUI'),
                                           imageOutput('normalityIsotropyVarianceCheckImage'),
                                           uiOutput('isotropyVarianceCheckUI'),
                                           plotOutput('directionalVariogramsVarianceCheckImage')
                                         )
                                       )
                              ),
                              tabPanel('Models',
                                       sidebarLayout(
                                         sidebarPanel(
                                           uiOutput('imageNameVarianceModels'),
                                           tags$hr(),
                                           sliderInput('sampleVariogramModelsSlider', label = 'Sample to calculate directional variograms',
                                                       min = 1, max = 100, value = 100, step=1),
                                           sliderInput('sampleFitVarianceModelsSlider', label = 'Sample to fit models',
                                                       min = 0, max = 4, value = 3.43, step=0.01),
                                           sliderInput('sampleValidateVarianceModelsSlider', label = 'Sample to validate models',
                                                       min = 1, max = 100, value = 100),
                                           actionButton('fitVariogramVarianceModelsButton',label='Fit variogram models')
                                         )
                                         ,
                                         mainPanel(
                                           uiOutput('bestmodeltableVarianceModelsUI'),
                                           tableOutput('bestmodeltableVarianceModelsPlot'),
                                           imageOutput('allmodelsModelImage')
                                         )
                                       )
                              )
                   ),
                   tabPanel('Fractal',
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput('imageNameFractal'),
                                tags$hr(),
                                actionButton('calculateFractalDimensionsButton',label='Calculate')
                              )
                              ,
                              mainPanel(
                                uiOutput('fractalDimensionsFractalUI'),
                                tableOutput('fractalDimensionsFractalTable'),
                                imageOutput('fractalDimensionsFractalPlots_fd2d_isotropic'),
                                imageOutput('fractalDimensionsFractalPlots_fd2d_squareincr'),
                                imageOutput('fractalDimensionsFractalPlots_fd2d_filter1')
                              )
                            )),
#              tabPanel('Networks',
#                       sidebarLayout(
#                         sidebarPanel(
#                           uiOutput('imageNameNetworks'),
#                           registerSceneChange(),
#                           tags$hr(),
#                           sliderInput('heightNetworksslider', label = 'Height multiplier',
#                                       min = 0.1, max = 10, value = 1, step=0.1),
#                           sliderInput('filterNetworksslider', label = 'Filter',
#                                       min = 0.1, max = 10, value = c(1,10), step=0.1),
#                           actionButton('checkFilterNetworksButton', label = 'Check filter'),
#                           actionButton('calculateNetworksNetworksButton', label = 'Calculate networks')
#                         ),
#                         mainPanel(
#                           uiOutput('panelNetworksUI'),
#                           plotOutput("distNetworksPlot"),
#                           plotOutput("newImageNetworksPlot"),
#                           plotOutput("skeletonImageNetworksPlot")
#                         )
#                       )),
                   tabPanel('3D',
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput('imageName3D'),
                                registerSceneChange(),
                                tags$hr(),
                                sliderInput('height3Dslider', label = 'Height multiplier',
                                            min = 0.1, max = 10, value = 1, step=0.1),
                                actionButton('displayIn3D3DButton', label = 'Display 3D image'),
                                downloadButton('snapshot3DButton', label = 'Snapshot'),
                                tags$hr(),
                                actionButton('calculate3DModel3DButton', label = 'Calculate 3D model for printing'),
                                downloadButton('export3DModel3DButton',label='Export model for 3D printing')
                              ),
                              mainPanel(
                                uiOutput('panel3DUI')
                                ,rglwidgetOutput('thewidget', width = "100%", height = 600)
                              )
                            )),                   
                   tabPanel('Reports',
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput('imageNameReports'),
                                tags$hr(),                                
                                downloadButton('generateCheckReport', label = 'Download check report'),
                                tags$hr(),
                                downloadButton('generateReport', label = 'Download full report')
                              ),
                              mainPanel(
                                tableOutput('alreadyCalculatedPlot')
                              )
                            )),
                   tabPanel('About',
                            mainPanel(
                              tags$iframe(
                                seamless="seamless",
                                src="http://www.openbmap.com/AFM/index_afmapp.html",style="width: 400px; height: 400px")
                            )
                   ),
             tags$head(tags$script(src="google-analytics.js")),
             useShinyjs()
        ))

