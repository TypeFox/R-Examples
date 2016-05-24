#' Interactive Plotting for Longtiudinal Functional Data Analysis using FPCA
#'
#' Produces an interactive plot illustrating longitudinal functioanl data analysis (Park and Staicu, 2015).
#' 
#' @param obj lfpca object to be plotted.
#' @param xlab x axis label
#' @param ylab y axis label
#' @param title plot title
#' @param ... additional arguments passed to plotting functions
#' 
#' @author So Young Park \email{spark13@@ncsu.edu}, Ana-Maria Staicu \email{astaicu@@ncsu.edu}
#' @export
#' 
#' @references Park, S.Y. and Staicu, A.M. (2015). Longitudinal functional data analysis. Stat 4 212-226.
#' 
#' @seealso \code{\link{plot_shiny}}; \code{fpca.lfda} in the refund package for estimation method. 
#' @import shiny
#' @import ggplot2
#' @import lme4
#' @importFrom grDevices rainbow
#' @importFrom utils sessionInfo
#'  
plot_shiny.lfpca <- function(obj, xlab = "", ylab="", title = "", ...){
  
  ##################################
  # for exploratory analysis
  ##################################
  x <- obj
  
  ### NULLify global values called in ggplot
  funArg = y = curveID = MEAN = COV = lambda = i = NULL
  
  result <- x
  myDat <- result$obsData
  n <- length(unique(myDat$i))
  numFuncArg <- length(myDat$funcArg)
  numCurves <- length(myDat$i)
  vecData <- data.frame(curveID = rep(1:numCurves, each=numFuncArg),
                         i=rep(myDat$i, each = numFuncArg), 
                         j=rep(myDat$j,each = numFuncArg), 
                         y=as.vector(t(myDat$y)), 
                         funArg = rep(myDat$funcArg, numCurves), 
                         longTimes = rep(myDat$Tij, each = numFuncArg))
  overallMean <- as.vector(colMeans(myDat$y) )
  mgin <- diff(range(vecData$y))*0.01
  my.Ylim <- c(range(vecData$y)[1] - mgin, range(vecData$y)[2] + mgin)
  
  ##################################
  # for longitudinal fpca 
  ##################################
  
  
  result$TT <- seq(min(result$visitTime), max(result$visitTime), length.out=nrow(result$bivariateSmoothMeanFunc))
  
  b<-list()
  for(k in 1:result$mFPCA.npc){
    a <- do.call(rbind, lapply(result$sFPCA.xiHat.bySubj, function(i) i[,k]))
    b[[k]] <- a
  }
  result$sFPCA.xiHat <- b
  result$mFPCA.pve <- cumsum(result$mFPCA.scree.eval)/sum(result$mFPCA.scree.eval)
  
  
  
  shinyApp(
    
    ui = navbarPage(title = strong(style = "color: #ACD6FF; padding: 0px 0px 10px 10px; opacity: 0.95; ", "LFPCA Plot"), windowTitle = "PlotInteractiveLFDA", 
                    collapsible = FALSE, id = "nav",
                    inverse = TRUE, header = NULL,
                    
                    tabPanel("Exploratory Plots", icon = icon("arrow-circle-right", lib = "font-awesome"), 
                             tabsetPanel(
                               #################################
                               ## Tab 1: Data
                               #################################   
                               
                               tabPanel("Data", icon = icon("archive", lib = "font-awesome"),
                                        column(3,
                                               checkboxInput("allDat", h4("Show all observations in background"), FALSE),
                                               hr(),
                                               checkboxInput("overallMean", h4("Show overall mean"), FALSE),
                                               hr(),
                                               selectInput('subj', h4('Select subject ID'), choices = unique(myDat$i), selected=1),
                                               hr(),
                                               helpText("Plot shows observed curves for the subject selected above."),
                                               
                                               br(), downloadButton('downloadPlotObsPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotObs", "Download Plot as Object", class = "plot-download")
                                        ),
                                        
                                        column(9, h4("Plot of Longitudinal Functional Data"),
                                               plotOutput('Obs')
                                        )
                               ),
                               
                               #################################
                               ## Tab 2: Data by subject
                               #################################   
                               
                               tabPanel("Data by Subject", icon = icon("caret-square-o-right", lib="font-awesome"),
                                        column(3,
                                               selectInput('subj1', h4('Subject ID'), unique(myDat$i)),
                                               hr(),
                                               checkboxInput("allDatbySubj", h4("Show all observations from this subject"), TRUE), 
                                               hr(),
                                               checkboxInput("meanBySubj", h4("Show subject mean"), FALSE),
                                               hr(),
                                               checkboxInput("overallMean1", h4("Show overall mean"), FALSE),
                                               hr(),
                                               uiOutput("slider"),
                                               #                                      sliderInput("longiTime", h4("actual time of visits (scaled)"),
                                               #                                                  min = min(myDat$Tij), max = max(myDat$Tij),
                                               #                                                  value = min(myDat$Tij), animate=animationOptions(interval=400)),
                                               hr(),
                                               helpText("Grey lines indicate observed curves for the subject selected above. 
                                                        Solid blue line indicates subject-specific mean curve and black dashed line indicates overall mean."),
                                               
                                               #hr(),
                                               helpText("Play \'actual time of visits\' to see when each curve is observed (indicated with red)."),
                                               br(), downloadButton('downloadPlotObsbySubjPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotObsbySubj", "Download Plot as Object", class = "plot-download")
                                               
                                               
                                               ),
                                        
                                        column(9, h4("Plot of Longitudinal Functional Observations by Subject"),
                                               plotOutput("bySubj")
                                        )
                               ),
                               
                               ##############################################
                               ## Tab 3: Dist of Longitudinal Times
                               ##############################################   
                               
                               tabPanel("Distribution of Longitudinal Times", icon = icon("bar-chart", lib = "font-awesome"),
                                        column(3,
                                               selectInput('subj3', h4('Subject ID'), unique(myDat$i)),
                                               hr(),
                                               checkboxInput("allDatbySubj2", h4("Show all observations from this subject"), TRUE),
                                               hr(),
                                               helpText("Empty circles indicate sampling time points of all curves. 
                                                        Red dots indicate sampling time points of curves observed for the subject selected above. 
                                                        Histogram shows the distribution of sampling time points for all subject."),
                                               br(), downloadButton('downloadPlotLongiTimePDF', 'Download Plot as PDF'),
                                               br(), br(), downloadButton("downloadPlotLongiTime", "Download Plot as Object", class = "plot-download"),
                                               br(), br(), downloadButton('downloadPlotLongiHistPDF', 'Download Plot as PDF (Histogram)'), 
                                               br(), br(), downloadButton("downloadPlotLongiHist", "Download Plot as Object (Histogram)", class = "plot-download")
                                               
                                               ),
                                        
                                        column(9, h4("Sampling Distribution of Longitudinal Times"),
                                               plotOutput("lonTime"),
                                               plotOutput("lonTimeHist")
                                        )
                                        )
                             )
                             
                             
                    ),
                    
                    
                    tabPanel("Estimation Plots", icon = icon("arrow-circle-right", lib = "font-awesome"),
                             tabsetPanel(
                               
                               tabPanel("Mean Surface", icon = icon("check-square", lib = "font-awesome"),
                                        column(3, 
                                               helpText("Plot shows estimated bivariate smooth mean function; x-axis represents visit time (T), y-axis represents functional argument (s), and color represents the estimated mean function. "),
                                               br(), downloadButton('downloadPlotmuhatPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotmuhat", "Download Plot as Object", class = "plot-download")
                                               
                                        ),
                                        
                                        column(9, h4("Bivariate Mean Surface"),
                                               plotOutput('bivMean')
                                        )
                               ),
                               
                               # estimated marginal COVARIANCE
                               tabPanel("Marginal Covariance", icon = icon("star", lib="font-awesome"),
                                        column(3, 
                                               helpText("Plot shows estimated marginal covariance."),
                                               
                                               br(), downloadButton('downloadPlotMarginalCOVPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotMarginalCOV", "Download Plot as Object", class = "plot-download")
                                               
                                        ),
                                        
                                        column(9, h4("Estimated Marginal Covariance"),
                                               plotOutput('MarginalCOV')
                                        )            
                               ),
                               
                               # estimated marginal eigenfunctions
                               tabPanel("Marginal Eigenfunctions", icon=icon("eyedropper", lib="font-awesome"),
                                        column(3,
                                               selectInput("PCchoiceMEIG", label = h4("Select FPC"), choices = seq_len(result$mFPCA.npc), selected = 1),
                                               hr(),
                                               helpText("Plot shows estimated marginal eigenfunction for the FPC selected above."),
                                               br(), downloadButton('downloadPlotMEIGPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotMEIG", "Download Plot as Object", class = "plot-download")
                                               
                                               
                                               
                                        ),
                                        
                                        column(9, h4("Estimated Marginal Eigenfunction"),
                                               plotOutput('MEIG')
                                        ) 
                                        
                                        
                                        
                                        
                                        
                               ),
                               
                               tabPanel("Mean +/- mFPCs", icon = icon("stats", lib = "glyphicon"),
                                        column(3, 
                                               selectInput("PCchoice", label = h4("Select FPC"), choices = seq_len(result$mFPCA.npc), selected = 1),
                                               hr(),
                                               helpText("Solid black line indicates population mean. For the FPC selected above, blue and red lines
                                                        indicate the population mean +/- the FPC times 2 SDs of the associated marginal score distribution (square root of the corresponding marignal eigenvalue)."),
                                               br(), downloadButton('downloadPlotLinComPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotLinCom", "Download Plot as Object", class = "plot-download")
                                               
                                               
                                               ),
                                        
                                        column(9, h4("Mean +/- marginal FPCs"),
                                               plotOutput('PCplot')
                                        )
                                        
                               ),
                               
                               tabPanel("Scree Plot", icon = icon("medkit"),
                                        column(3, 
                                               helpText("Scree plots are displayed to the right. The first panel shows the plot of eigenvalues, and 
                                                        the second plot shows the cumulative percent variance explained."),
                                               br(), downloadButton('downloadPlotScreePDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotScree", "Download Plot as Object", class = "plot-download")
                                               
                                               ),
                                        column(9, h4("Scree Plots (marginal FPCA)"), 
                                               plotOutput('Scree')
                                        )     
                               ),
                               
                               # covariance of longitudinal dynamics (by marginal PCs)
                               tabPanel("Covariance of Longitudinal Dynamics", icon = icon("star", lib="font-awesome"),
                                        column(3,  
                                               selectInput("PCchoiceLDCOV", label = h4("Select FPC"), choices = seq_len(result$mFPCA.npc), selected = 1),
                                               hr(),
                                               helpText("Plot shows estimated covariance of longitudinal dynamics for the marginal FPC selected above."),
                                               br(), downloadButton('downloadPlotLDCOVkPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotLDCOVk", "Download Plot as Object", class = "plot-download")
                                               
                                        ),       
                                        column(9, h4("Estimated Covariance of Longitudinal Dynamics"), plotOutput('LDCOVk'))
                               ),
                               
                               tabPanel("Basis Coefficients", icon=icon("user", lib = "font-awesome"),
                                        column(3, hr(),
                                               checkboxInput("allSubject", h4("Show all subjects in background"), TRUE),
                                               hr(),
                                               selectInput('subject', h4('Select Subject ID'), unique(result$i)),
                                               hr(),
                                               selectInput("PCchoice1", label = h4("Select FPC"), choices = seq_len(result$mFPCA.npc), selected = 1),
                                               hr(),
                                               helpText("Grey lines indicate predicted basis coefficient functions of longitudinal time (T). 
                                                        Red line highlights a predicted curve for the marginal FPC and subject selected above."),
                                               br(), downloadButton('downloadPlotScorePDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotScore", "Download Plot as Object", class = "plot-download")
                                               
                                               ),
                                        column(9, h4("Predicted Basis Coefficients Plot"),
                                               plotOutput('xiHat_T')
                                        )
                               ),
                               
                               tabPanel("Yhat", icon=icon("dot-circle-o", lib = "font-awesome"),
                                        column(3, hr(),
                                               sliderInput("longiTime1", "Longitudinal Time",
                                                           min = min(result$TT), max = max(result$TT), step = 2*(result$TT[2]-result$TT[1]),
                                                           value = min(result$TT), animate=animationOptions(interval=800)),
                                               hr(),
                                               selectInput('subject1', h4('Select subject ID'), unique(result$i)),
                                               hr(),
                                               helpText("Plot shows predicted trajectory of the subject selected above. Play \'Longituidnal Time\' to see longitudinal change of the predicted trajectory. "),
                                               br(), downloadButton('downloadPlotPredictedYPDF', 'Download Plot as PDF'), br(), br(),
                                               downloadButton("downloadPlotPredicted", "Download Plot as Object", class = "plot-download")
                                               
                                        ),
                                        column(9, h4("Estimated individual Trajectory"),
                                               plotOutput('Yhat')
                                        )
                               )
                               
                               
                             )
                    )         
                    
                    

    ),
    
    server <- function(input, output) {
      
      #################################
      #################################
      ###   Exploratory Plots
      #################################
      #################################
      
      #################################
      ## Code for Tab 1: Data
      #################################   
      
      plotInputObs <- reactive({
        showAllDat <- input$allDat
        showOverallMean <- input$overallMean
        subjectChoice <- input$subj
        
        selectedData <- vecData[which(vecData$i==subjectChoice),]
        subjNumVisits <- max(unique(selectedData$j))
        myColors <- rainbow(max(vecData$j))[1:subjNumVisits]
        
        p <- ggplot(data=as.data.frame(vecData[which(vecData$i==1),]), aes(x=funArg, y=y, group=curveID)) + geom_line(linetype="blank") + theme_bw() + ylim(my.Ylim)
        
        if(showAllDat) p <- p + geom_line(aes(x=funArg, y=y, group=curveID), colour= "#b6b6b6", data=as.data.frame(vecData))
        if(showOverallMean) p <- p + geom_line(mapping=aes(x=x,y=y, group=rep(1,length(x))), data=data.frame(x=myDat$funcArg, y=overallMean), colour="black", linetype=2, size=1)
        
        p <- p + geom_line(aes(x=funArg, y=y, group=curveID), colour= rep(myColors, each=numFuncArg), data=as.data.frame(selectedData))
        p <- p + theme(legend.position="none") + xlab("s") + ylab("y") + 
          ggtitle(paste("Subject ", subjectChoice, "; Total of ", max(selectedData$j), " visits", sep=""))
        
        p
      })
      output$Obs <- renderPlot({
        print(plotInputObs())
      })
      
      output$downloadPlotObsPDF <- savePDF("all_observations.pdf", plotInputObs())  
      output$downloadPlotObs <- savePlot("all_observations.Rdata", plotInputObs())  
      
      #################################
      ## Code for Tab 2: Data by Subject
      ################################# 
      output$slider <- renderUI({
        
        subjectChoice <- input$subj1
        vecSelectedData <- vecData[which(vecData$i==subjectChoice),]
        TijSubj <- vecSelectedData$longTimes
        
        args <- list(inputId="longiTime", label=h4("actual time of visits (scaled)"), 
                     value = min(unique(TijSubj)), 
                     animate=animationOptions(interval=2000), 
                     ticks = c(unique(TijSubj)),
                     round = -2)
        args$min = 1
        args$max = length(args$ticks)
        if (sessionInfo()$otherPkgs$shiny$Version>="0.11") {
          
          ticks <- paste0(round(args$ticks, digits=2), collapse=',')
          args$ticks <- TRUE
          html  <- do.call('sliderInput', args)
          html$children[[2]]$attribs[['data-values']] <- ticks;
          
        }else {
          html  <- do.call('sliderInput', args)     
        }
        
        html
      })
      
      plotInputDataBySubj <- reactive({
        
        subjectChoice <- input$subj1
        vecSelectedData <- vecData[which(vecData$i==subjectChoice),]
        indSubj <- which(myDat$i==input$subj1)
        obsBySubj <- myDat$y[indSubj,]
        TijSubj <- vecSelectedData$longTimes
        nrep <- length(unique(TijSubj))
        
        ind <- input$longiTime
        if(input$allDatbySubj){
          p <- ggplot(data=vecSelectedData, aes(x=funArg, y=y, group=curveID)) + geom_line(colour="#b6b6b6")  + theme_bw() + ylim(my.Ylim)
#           if(any(c((ind==0), (ind==nrep+1), (length(ind) < 1)))){
#             p <- p + xlab("s") + ylab("y")+ggtitle(paste("Subject ", subjectChoice, "; Total of ", nrep, " visits", sep=""))
#           }else{
            p <- p + geom_line(data=vecSelectedData[which(vecSelectedData$j==(ind+1)), ], aes(x=funArg, y=y, group=curveID), colour="red")
            p <- p + xlab("s") + ylab("y") + 
              ggtitle(paste("Subject ", subjectChoice, "; Total of ", nrep, " visits; ", "Actual time of ", ind, "-th visit (red) = ",
                            round(unique(TijSubj)[ind+1], digits=2), sep=""))
          # }
          
          if(input$meanBySubj) p <- p + geom_line(mapping=aes(x=x,y=y, group=rep(1,length(myDat$funcArg))), data=data.frame(x=myDat$funcArg, y=colMeans(obsBySubj)), colour="blue", linetype=1, size=1)
          if(input$overallMean1) p <- p + geom_line(mapping=aes(x=x,y=y, group=rep(1,length(myDat$funcArg))), data=data.frame(x=myDat$funcArg, y=overallMean), colour="black", linetype=2, size=1)
          
          
          
        }else{
          p <- ggplot(data=vecSelectedData, aes(x=funArg, y=y, group=curveID)) + geom_line(linetype="blank")  + theme_bw() + ylim(my.Ylim)
#           if(any(c((ind==0), (ind==nrep+1), (length(ind) < 1)))){
#             p <- p + xlab("s") + ylab("y")+ggtitle(paste("Subject ", subjectChoice, "; Total of ", nrep, " visits", sep=""))
#           }else{
            p <- p + geom_line(data=vecSelectedData[which(vecSelectedData$j==(ind+1)), ], aes(x=funArg, y=y, group=curveID), colour="red")
            p <- p + xlab("s") + ylab("y") + 
              ggtitle(paste("Subject ", subjectChoice, "; Total of ", nrep, " visits; ", "Actual time of ", ind, "-th visit (red) = ",
                            round(unique(TijSubj)[ind+1], digits=2), sep=""))
          # } 
          if(input$meanBySubj) p <- p + geom_line(mapping=aes(x=x,y=y, group=rep(1,length(myDat$funcArg))), data=data.frame(x=myDat$funcArg, y=colMeans(obsBySubj)), colour="blue", linetype=1, size=1)
          if(input$overallMean1) p <- p + geom_line(mapping=aes(x=x,y=y, group=rep(1,length(myDat$funcArg))), data=data.frame(x=myDat$funcArg, y=overallMean), colour="black", linetype=2, size=1)
          
        }
        p
      })
      
      output$bySubj <- renderPlot({
        print(plotInputDataBySubj())
      })
      
      output$downloadPlotObsbySubjPDF <- savePDF("obs_by_subj.pdf", plotInputDataBySubj())  
      output$downloadPlotObsbySubj <- savePlot("obs_by_subj.Rdata", plotInputDataBySubj())  
      
      ####################################################
      ## Code for Tab 3: Dist of Longitudinal Times
      ####################################################    
      
      plotInputLongTime <- reactive({
        
        subjectChoice1 <- as.numeric(input$subj3)
        index <- which(as.vector(myDat$i)==subjectChoice1)
        TijSubj1 <- as.vector(myDat$Tij[index])
        
        if(input$allDatbySubj2){
          p <- ggplot(mapping=aes(x=myDat$Tij, y=rep(1,length(myDat$Tij))))+geom_line() + geom_point(size=4,shape=21, fill="white")+theme_bw()
          p <- p + xlab("Longitudinal Times") + ylab("") + theme(axis.text.y=element_blank()) +
            ggtitle(paste("Longitudinal Times of Subject ", subjectChoice1, sep="")) + 
            geom_point(mapping=aes(x= x, y=y), data=data.frame(x=TijSubj1, y=rep(1, length(TijSubj1))), size=4, shape=21, fill="red")
        }else{
          p <- ggplot(mapping=aes(x=myDat$Tij, y=rep(1,length(myDat$Tij))))+theme_bw()+geom_line(linetype=0)
          p <- p + xlab("Longitudinal Times") + ylab("") + theme(axis.text.y=element_blank()) +
            ggtitle(paste("Longitudinal Times of Subject ", subjectChoice1, sep=""))
          p <- p + geom_line(mapping=aes(x= x, y = y), data=data.frame(x=TijSubj1, y=rep(1, length(TijSubj1))), linetype=1)
          p <- p + geom_point(mapping=aes(x= x, y = y), data=data.frame(x=TijSubj1, y=rep(1, length(TijSubj1))),
                              size=4, shape=21, fill="red")
        }
        p
      })
      
      plotInputLongTimeHist <- reactive({
        p <- ggplot(mapping=aes(x=myDat$Tij)) + theme_bw() 
        p <- p + geom_histogram(binwidth=0.05,  col="black", fill="#b6b6b6", alpha = 1) 
        p <- p + xlab("Longitudinal Times") + ylab("Frequency") + ggtitle("Histogram of Longitudinal Times")
        p  
      })
      output$lonTime <- renderPlot({
        print(plotInputLongTime())
      })
      output$lonTimeHist <- renderPlot({
        print(plotInputLongTimeHist())
      })
      
      output$downloadPlotLongiTimePDF <- savePDF("sampling_pts.pdf", plotInputLongTime())  
      output$downloadPlotLongiHistPDF <- savePDF("histogram.pdf", plotInputLongTimeHist())  
      
      output$downloadPlotLongiTime <- savePlot("sampling_pts.Rdata", plotInputLongTime())  
      output$downloadPlotLongiHist <- savePlot("histogram.Rdata", plotInputLongTimeHist())  
      
      #################################
      #################################
      ###   Estimation Plots
      #################################
      #################################
      
      ##########################################
      ## Code for Tab 1: Bivariate Mean Surface 
      ##########################################
      plotBivMean <- reactive({
        p <- ggplot(data = data.frame(x=rep(result$funcArg, each = length(result$TT)), y= rep(result$TT, length(result$funcArg)), 
                                      MEAN=as.vector(result$bivariateSmoothMeanFunc)), 
                    aes(x=x,y=y,fill=MEAN)) + geom_tile() + theme_bw() + xlab("Longitudinal Time (T)") + ylab("s") 
        p <- p  + scale_fill_gradient(low="red", high="white")
        p
      })
      
      output$bivMean <- renderPlot(
        print(plotBivMean())
      )
      
      output$downloadPlotmuhatPDF <- savePDF("biv_mean_surface.pdf", plotBivMean())  
      output$downloadPlotmuhat <- savePlot("biv_mean_surface.Rdata", plotBivMean())  
      
      
      ##############################################################
      ## Code for Tab 2.0.0: estimated marginal Covariance plot
      ##############################################################    
      
      plotInputMCOV <- reactive({
        
        p <- ggplot(data = data.frame(x=rep(result$funcArg, length(result$funcArg)), y= rep(result$funcArg, each=length(result$funcArg)), 
                                      COV=as.vector(result$mFPCA.covar)), 
                    aes(x=x,y=y,fill=COV)) + geom_tile() + theme_bw() + xlab("s") + ylab("s") 
        p <- p  + scale_fill_gradient(low="red", high="white")
        p
        
      })
      
      output$MarginalCOV <- renderPlot(
        print(plotInputMCOV())
      )
      
      output$downloadPlotMarginalCOVPDF <- savePDF("marginal_cov.pdf", plotInputMCOV())  
      output$downloadPlotMarginalCOV <- savePlot("marginal_cov.Rdata", plotInputMCOV())  
      
      ##############################################################
      ## Code for Tab 2.0.1: estimated marginal eigenfunctions
      ##############################################################  
      
      plotInputMEIG <- reactive({
        PCchoiceMEIG <- as.numeric(input$PCchoiceMEIG)
        p <- ggplot(data=data.frame(x=result$funcArg, y=result$mFPCA.efunctions[,PCchoiceMEIG]), aes(x=x, y=y)) + geom_line(lwd=1, colour="red") + theme_bw()
        p <- p + geom_line(data=data.frame(x=result$funcArg, y=rep(0, length(result$funcArg))), aes(x=x,y=y), colour="grey", lwd=1, linetype="dashed")
        p <- p + xlab("s") + ylab("")
        p
      })
      
      output$MEIG <- renderPlot(
        print(plotInputMEIG())
      )
      output$downloadPlotMEIGPDF <- savePDF("marginal_eigenfns.pdf", plotInputMEIG())  
      output$downloadPlotMEIG <- savePlot("marginal_eigenfns.Rdata", plotInputMEIG())  
      
      ##########################################
      ## Code for Tab 2.1: marginal FPCs plot
      ##########################################    
      
      plotInputPC <- reactive({
        
        PCchoice = as.numeric(input$PCchoice)
        
        p <- ggplot(data=data.frame(x=result$funcArg, y=colMeans(result$bivariateSmoothMeanFunc) ), aes(x=x, y=y)) + geom_line(lwd=1) + theme_bw()
        p <- p + geom_point(data = data.frame(x=result$funcArg, y=(colMeans(result$bivariateSmoothMeanFunc) + 2*sqrt(result$mFPCA.evalues[PCchoice])*result$mFPCA.efunctions[,PCchoice])),
                            color = "blue", size = 4, shape = '+')
        p <- p + geom_point(data = data.frame(x=result$funcArg, y=(colMeans(result$bivariateSmoothMeanFunc)  - 2*sqrt(result$mFPCA.evalues[PCchoice])*result$mFPCA.efunctions[,PCchoice])),
                            color = "red", size = 4, shape = '-')
        p <- p + ylim(range(result$fitted.values.all)) + xlab("s") + ylab("")+
          ggtitle(bquote(""~psi[.(PCchoice)](s) ~ "," ~.(100*round(result$mFPCA.pve[(PCchoice)],3)) ~ "% Variance explained"))   
        p
      })
      
      output$PCplot <- renderPlot(
        print(plotInputPC())  
      )
      
      output$downloadPlotLinComPDF <- savePDF("marginal_LinCom.pdf", plotInputPC())  
      output$downloadPlotLinCom <- savePlot("marginal_LinCom.Rdata", plotInputPC())  
      
      
      #################################
      ## Code for Tab 3: scree plot
      #################################
      
      plotInputScree <- reactive({
        scree <- data.frame(k = rep(1:(result$mFPCA.npc+1), 2), 
                            lambda = c(result$mFPCA.scree.eval[1:(result$mFPCA.npc+1)], 
                                       result$mFPCA.pve[1:(result$mFPCA.npc+1)]),
                            type = rep(c("Eigenvalue", "Percent Variance Explained"), each = (result$mFPCA.npc+1)))
        
        p <- ggplot(data=scree, aes(x=k, y=lambda))+geom_line(linetype=1, lwd=1.5, color="black")+
          geom_point(size = 4, color = "black")+ theme_bw() + xlab("Principal Component") + ylab("") +
          facet_wrap(~type, scales = "free_y") + ylim(0, NA) 
        
        p
      })
      
      output$Scree <- renderPlot(
        print(plotInputScree())
      )
      
      output$downloadPlotScree <- downloadHandler(
        filename = function(){'screeplots_mFPCA.png' },
        content = function(file) {
          ggsave(file,plotInputScree())
        }
      )
      
      output$downloadPlotScreePDF <- savePDF("scree_plot.pdf", plotInputScree())  
      output$downloadPlotScree  <- savePlot("scree_plot.Rdata", plotInputScree())  
      
      ######################################################
      ## Code for Tab 4.0.0: covariance of longitudinal dynamics
      ######################################################
      
      plotInputLDCOVk <- reactive({
        PCchoice_LDCOVk = as.numeric(input$PCchoiceLDCOV)
        selectedCOV <- result$sFPCA.longDynCov.k[[PCchoice_LDCOVk]]
        p <- ggplot(data = data.frame(x=rep(result$TT, length(result$TT)), y= rep(result$TT, each=length(result$TT)), COV=as.vector(selectedCOV)), 
                    aes(x=x,y=y,fill=COV)) + geom_tile() + theme_bw() + xlab("Longitudinal Time (T)") + ylab("Longitudinal Time (T)") 
        p <- p + scale_fill_gradient(low="red", high="white")
        p
        
      })
      
      output$LDCOVk <- renderPlot({
        print(plotInputLDCOVk())
      })
      
      output$downloadPlotLDCOVk <- downloadHandler(
        filename = function(){ paste('LDynamics_k', input$PCchoiceLDCOV, '.png', sep="") },
        content = function(file) {
          ggsave(file,plotInputLDCOVk())
        }
      )
      output$downloadPlotLDCOVkPDF <- savePDF("dynamics_cov.pdf", plotInputLDCOVk())  
      output$downloadPlotLDCOVk  <- savePlot("dynamics_cov.Rdata", plotInputLDCOVk())  
      
      #####################################################
      ## Code for Tab 4.1: basis coefficients plot
      #####################################################
      
      plotInputxiHat <- reactive({
        
        n <- length(unique(result$i))
        subjChoice <- input$subject
        kChoice <- as.numeric(input$PCchoice1)
        vecDat <- data.frame(y=as.vector(result$sFPCA.xiHat[[kChoice]]),
                             i = rep(unique(result$i), length(result$TT)),
                             x = rep(result$TT, each=n))
        mylim <- c(range(vecDat$y)[1] - diff(range(vecDat$y))*0.1, range(vecDat$y)[2]+diff(range(vecDat$y))*0.1)
        if(input$allSubject){
          p <- ggplot(data=vecDat, aes(x=x, y=y, group=i))  + geom_line(colour="#b6b6b6") + theme_bw() + ylim(mylim)
          p <- p+geom_line(data=vecDat[which(vecDat$i==subjChoice),], aes(x=x,y=y,group=i), colour="red", lwd=1)
          p <- p+xlab("Longitudinal Time (T)") + ylab(bquote(hat(xi)(T))) + 
            ggtitle(paste("Predicted Basis Coefficients for k =", kChoice,sep=""))
        }else{
          p <- ggplot(data=vecDat, aes(x=x, y=y, group=i))  + geom_line(linetype="blank") + theme_bw()
          p <- p+geom_line(data=vecDat[which(vecDat$i==subjChoice),], aes(x=x,y=y,group=i), colour="red", lwd=1)
          p <- p+xlab("Longitudinal Time (T)") + ylab(bquote(hat(xi)(T))) + 
            ggtitle(paste("Predicted Basis Coefficients for k =", kChoice,sep=""))
        }
        p 
      })
      
      output$xiHat_T <- renderPlot(
        print(plotInputxiHat())
      )
      
      output$downloadPlotScorePDF <- savePDF("score.pdf", plotInputxiHat())  
      output$downloadPlotScore  <- savePlot("score.Rdata", plotInputxiHat())  
      
      #################################
      ## Code for Tab 5: individual trajectory
      #################################
      
      plotInputYhat <- reactive({
        tempDat <- result$fitted.values.all[[as.numeric(input$subject1)]]
        v <- which(abs(result$TT-input$longiTime1)== min(abs(result$TT-input$longiTime1)))
        fixedT = result$TT[v]
        
        p <- ggplot(data.frame(x=result$funcArg, y=tempDat[v,]), aes(x=x,y=y))+geom_line(colour="red", lwd=1) + theme_bw()
        p <- p+xlab("s")+ylab(bquote(hat(Y)))+ggtitle(paste("Subject i = ", input$subject1))+ ylim(range(result$fitted.values.all))
        p
      })
      output$Yhat <- renderPlot({
        print(plotInputYhat())
      })
      output$downloadPlotYhat <- downloadHandler(
        filename = function(){paste('Yhat_Subject', input$subject1,"T",fixedT,'.png') },
        content = function(file) {
          ggsave(file,plotInputYhat())
        }
      )
      
      output$downloadPlotPredictedYPDF <- savePDF("predictedY.pdf", plotInputYhat())  
      output$downloadPlotPredictedY  <- savePlot("predictedY.Rdata", plotInputYhat())  
      
    }
    
    )
  
  
  
  
  
}  