#' Interactive Plotting for Functional Principal Component Analysis
#'
#' Produces an interactive plot illustrating a functional principal component
#' analysis.
#'
#' @param obj fpca object to be plotted.
#' @param xlab x axis label
#' @param ylab y axis label
#' @param title plot title
#' @param ... additional arguments passed to plotting functions
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu},
#' Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#'
#' @seealso \code{\link{plot_shiny}}
#' @import shiny
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @export
#'
plot_shiny.fpca = function(obj, xlab = "", ylab="", title = "", ...) {

  fpca.obj <- obj

  ### NULLify global values called in ggplot
  V1 = V2 = V3 = V4 = k = lambda = value = subj = time = NULL

  ################################
  ## code for processing tabs
  ################################

  ## Tab 2: scree plot

  scree <- data.frame(k = rep(1:fpca.obj$npc, 2),
                      lambda = c(fpca.obj$evalues, cumsum(fpca.obj$evalues)/ sum(fpca.obj$evalues)),
                      type = rep(c("Eigenvalue", "Percent Variance Explained"), each = fpca.obj$npc))

  ## Tab 3: linear combination of PCs

  varpercent = lapply(fpca.obj$evalues, function(i){100*round(i/sum(fpca.obj$evalues),3)}) # calculates percent variance explained
  calls <- as.list(rep(NA, fpca.obj$npc))
  PCs <- rep(NA, fpca.obj$npc)
  for(i in 1:fpca.obj$npc){

    PCnum = paste("PC", i, sep="")

    calls[[i]] =  eval(call("sliderInput", inputId= PCnum, label = paste(PCnum, ": ", varpercent[[i]],  "% Variance", sep=""),
                            min = -2, max = 2, step = .1, value = 0, post = " SD", animate = animationOptions(interval=400, loop=T)))

    PCs[i] = PCnum
  }

  #################################
  ## App
  #################################

  shinyApp(

  #################################
  ## UI
  #################################

    ui = navbarPage(title = strong(style = "color: #ACD6FF; padding: 0px 0px 10px 10px; opacity: 0.95; ", "FPCA Plot"), windowTitle = "refund.shiny",
                    collapsible = FALSE, id = "nav",
                    inverse = TRUE, header = NULL,
                    tabPanel("Mean +/- FPCs", icon = icon("stats", lib = "glyphicon"),
                             column(3,
                                    helpText("Solid black line indicates population mean. For the FPC selected below, blue and red lines
                                             indicate the population mean +/- the FPC times 2 SDs of the associated score distribution."), hr(),
                                    selectInput("PCchoice", label = ("Select FPC"), choices = 1:fpca.obj$npc, selected = 1), hr(),
                                    downloadButton('downloadPDFMuPC', "Download Plot as PDF"), br(), br(),
                                    downloadButton("downloadPlotMuPC", "Download Plot as Object", class = "plot-download")
                                    ),
                             column(9, h4("Mean and FPCs"), plotOutput('muPCplot') )
                            ),
                    tabPanel("Scree Plot", icon = icon("medkit"),
                             column(3,
                                    helpText("Scree plots; the left panel shows the plot of eigenvalues and
                                             the right panel shows the cumulative percent variance explained."), hr(),
                                   downloadButton("downloadPDFScree", "Download Plot as PDF"), br(), br(),
                                   downloadButton("downloadPlotScree", "Download Plot as Object", class = "plot-download")
                                    ),
                             column(9, h4("Scree Plots"),  plotOutput('Scree') )
                            ),
                    tabPanel("Linear Combinations", icon = icon("line-chart"),
                             withMathJax(),
                             column(3,
                                    helpText("Plot shows the linear combination of mean and FPCs with the scores specified using the sliders below."), hr(),
                                    eval(calls), hr(),
                                    downloadButton("downloadPDFLinCom", "Download Plot as PDF"), br(), br(),
                                    downloadButton("downloadPlotLinCom", "Download Plot as Object", class = "plot-download")
                                    ),
                             column(9, h4("Linear Combination of Mean and FPCs"),
                                      plotOutput('LinCom')
                                    )
                             ),
                    tabPanel("Subject Fits",  icon = icon("user"),
                             column(3,
                                    helpText("Plot shows observed data and fitted values for the subject selected below"),
                                    selectInput("subject", label = ("Select Subject"), choices = 1:dim(fpca.obj$Yhat)[1], selected =1), hr(),
                                    downloadButton("downloadPDFSubject", "Download Plot as PDF"), br(), br(),
                                    downloadButton("downloadPlotSubject", "Download Plot as Object", class = "plot-download")
                                    ),
                             column(9, h4("Fitted and Observed Values for Selected Subject"),
                                      plotOutput("Subject")
                                    )
                             ),
                    tabPanel("Score Scatterplot",icon = icon("binoculars"),
                             fluidRow(
                             column(3, helpText("Use the drop down menus to select FPCs for the X and Y axis. Plot shows observed score
                                             scatterplot for selected FPCs; click and drag on the scatterplot to select subjects."), hr(),
                                    selectInput("PCX", label = ("Select X-axis FPC"), choices = 1:fpca.obj$npc, selected = 1),
                                    selectInput("PCY", label = ("Select Y-axis FPC"), choices = 1:fpca.obj$npc, selected = 2)
                                    ),
                             column(9, h4("Score Scatterplot for Selected FPCs"),
                                      plotOutput("ScorePlot",
                                                 brush=brushOpts(
                                                   id = "ScorePlot_brush",
                                                   resetOnNew = TRUE)
                                                 )
                                    )),
                             fluidRow(
                               column(3, helpText("Black curves are fitted values for all subjects. Blue curves correspond to subjects
                                                  selected in the graph above. If no points are selected, the mean curve is shown.")
                                      ),
                               column(9,
                                      plotOutput("ScorePlot2")
                                      )
                               )
                             )
                    ),

    #################################
    ## Server
    #################################

    server = function(input, output){

      mu = as.data.frame(cbind(1:length(fpca.obj$mu), fpca.obj$mu))
      efunctions = fpca.obj$efunctions; sqrt.evalues = diag(sqrt(fpca.obj$evalues))
      scaled_efunctions = efunctions %*% sqrt.evalues

      plotDefaults = list(theme_bw(), xlab(xlab), ylab(ylab), ylim(c(range(fpca.obj$Yhat)[1], range(fpca.obj$Yhat)[2])),
                          scale_x_continuous(breaks = seq(0, length(fpca.obj$mu)-1, length=6), labels = paste0(c(0, 0.2, 0.4, 0.6, 0.8, 1))) )

      #################################
      ## Code for mu PC plot
      #################################

      plotInputMuPC <- reactive({
        PCchoice = as.numeric(input$PCchoice)
        scaled_efuncs = scaled_efunctions[,PCchoice]

        p1 <- ggplot(mu, aes(x = V1, y = V2)) + geom_line(lwd=1) + plotDefaults +
          geom_point(data = as.data.frame(cbind(1:length(fpca.obj$mu), fpca.obj$mu + 2*scaled_efuncs)), color = "blue", size = 4, shape = '+')+
          geom_point(data = as.data.frame(cbind(1:length(fpca.obj$mu), fpca.obj$mu - 2*scaled_efuncs)), color = "indianred", size = 4, shape = "-")+
          ggtitle(bquote(psi[.(input$PCchoice)]~(t) ~ "," ~.(100*round(fpca.obj$evalues[as.numeric(input$PCchoice)]/sum(fpca.obj$evalues),3)) ~ "% Variance"))
      })

      output$muPCplot <- renderPlot(
        print(plotInputMuPC())
      )

      output$downloadPDFMuPC <- savePDF("muPC.pdf", plotInputMuPC())
      output$downloadPlotMuPC <- savePlot("muPC.RData", plotInputMuPC())

      #################################
      ## Code for scree plot
      #################################

      plotInputScree <- reactive({
        p2 <-screeplots <- ggplot(scree, aes(x=k, y=lambda))+geom_line(linetype=1, lwd=1.5, color="black")+
          geom_point(size = 4, color = "black")+ theme_bw() + xlab("Principal Component") + ylab("") +
          facet_wrap(~type, scales = "free_y") + ylim(0, NA)
      })

      output$Scree <- renderPlot(
        print(plotInputScree())
      )

      output$downloadPDFScree <- savePDF("scree.pdf", plotInputScree())
      output$downloadPlotScree <- savePlot("scree.RData", plotInputScree())


      #################################
      ## Code for linear combinations
      #################################

      plotInputLinCom <- reactive({
        PCweights = rep(NA, length(PCs)); for(i in 1:length(PCs)){PCweights[i] = input[[PCs[i]]]}
        df = as.data.frame(cbind(1:length(fpca.obj$mu), as.matrix(fpca.obj$mu)+efunctions %*% sqrt.evalues %*% PCweights ))

        p3 <- ggplot(mu, aes(x=V1, y=V2))+geom_line(lwd=0.75, aes(color= "mu"))+  plotDefaults + theme(legend.key = element_blank()) +
          geom_line(data = df, lwd = 1.5, aes(color = "subject")) + xlab(xlab) + ylab(ylab) + ggtitle(title)+
          scale_color_manual("Line Legend", values = c(mu = "gray", subject = "cornflowerblue"), guide = FALSE)
      })

      output$LinCom <- renderPlot( print(plotInputLinCom()) )

      output$downloadPDFLinCom <- savePDF("LinCom.pdf", plotInputLinCom())
      output$downloadPlotLinCom <- savePlot("LinCom.RData", plotInputLinCom())

      #################################
      ## Code for subject plots
      #################################

      plotInputSubject <- reactive({
        subjectnum = as.numeric(input$subject)
        df = as.data.frame(cbind(1:length(fpca.obj$mu), fpca.obj$mu, fpca.obj$Yhat[subjectnum,], fpca.obj$Y[subjectnum,]))

        p4 <- ggplot(data = df, aes(x=V1,y=V2)) + geom_line(lwd=0.5, color = "gray") + plotDefaults +
          geom_line(data = df, aes(y=V3), size=1, color = "cornflowerblue") +
          geom_point(data = df, aes(y=V4), color = "blue", alpha = 1/3)
      })

      output$Subject <- renderPlot( print(plotInputSubject()) )

      output$downloadPDFSubject <- savePDF("subject.pdf", plotInputSubject())
      output$downloadPlotSubject <- savePlot("subject.RData", plotInputSubject())

      #################################
      ## Code for score plots
      #################################
      scoredata = as.data.frame(cbind(fpca.obj$scores, fpca.obj$Yhat))
      colnames(scoredata) = c(paste0("PC", 1:fpca.obj$npc), paste0("subj", 1:dim(fpca.obj$Yhat)[2]))

      ## get PCs selected for X and Y axis
      PCX <- reactive({ paste0("PC", input$PCX) })
      PCY <- reactive({ paste0("PC", input$PCY) })


      ## Tab 5 Plot
      output$ScorePlot <- renderPlot({
        ggplot(scoredata, aes_string(x = PCX(), y = PCY()))+geom_point(color = "blue", alpha = 1/5, size = 3)+theme_bw()+
          xlab(paste("Scores for FPC", input$PCX))+ylab(paste("Scores for FPC", input$PCY))
      })

      ### second score plot
      Yhat.all.m = melt(fpca.obj$Yhat)
      colnames(Yhat.all.m) = c("subj", "time", "value")
      baseplot = ggplot(Yhat.all.m, aes(x=time, y=value, group = subj)) + geom_line(alpha = 1/5, color="black") + plotDefaults

      output$ScorePlot2 <- renderPlot({

        brush <- input$ScorePlot_brush
        if(!is.null(brush)){
          points = brushedPoints(scoredata, input$ScorePlot_brush, xvar=PCX(), yvar = PCY())
          Yhat.m = melt(as.matrix(points[,-c(1:fpca.obj$npc)]))

        }else{
          Yhat.m = as.data.frame(cbind(1, 1:length(fpca.obj$mu), fpca.obj$mu))
        }

        colnames(Yhat.m) <- c("subj", "time", "value")
        baseplot+geom_line(data= Yhat.m, aes(x=as.numeric(time), y=value, group = subj), color="cornflowerblue")

      })

    } ## end server
  )
}

