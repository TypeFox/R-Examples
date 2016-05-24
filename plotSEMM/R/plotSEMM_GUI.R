#' PlotSEMM GUI
#'
#' Graphical user interface with the shiny package. Supports manual input as well as importing
#' from precomputed Mplus files. An online tutorial and additional materials can be found at
#' \url{http://www.yorku.ca/pek/index_files/appendices.htm}
#'
#' @param ... additional arguments passed to \code{shiny::runApp}, such as
#'   \code{launch.browser = TRUE}
#'
#' @aliases plotSEMM_GUI
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com} and Jolynn Pek
#' @keywords shiny GUI
#' @export plotSEMM_GUI
#' @references
#'
#' Bauer, D.J. (2005). A semiparametric approach to modeling nonlinear relations among latent variables.
#' \emph{Structural Equation Modeling: A Multidisciplinary Journal, 12}(4), 513-535.
#'
#' Pek, J. & Chalmers, R. P. (2015). Diagnosing Nonlinearity With Confidence Envelopes for a
#' Semiparametric Approach to Modeling Bivariate Nonlinear Relations Among Latent Variables.
#' \emph{Structural Equation Modeling, 22}, 288-293.
#'
#' Pek, J., Chalmers, R. P., Kok B. E., & Losardo, D. (2015). Visualizing Confidence Bands for
#' Semiparametrically Estimated Nonlinear Relations among Latent Variables.
#' \emph{Journal of Educational and Behavioral Statistics, 40}, 402-423.
#'
#' Pek, J., Losardo, D., & Bauer, D. J. (2011). Confidence intervals for a semiparametric
#' approach to modeling nonlinear relations among latent variables. \emph{Structural Equation
#' Modeling, 18}, 537-553.
#'
#' Pek, J., Sterba, S. K., Kok, B. E., & Bauer, D. J. (2009). Estimating and visualizing non-linear
#' relations among latent variables: A semiparametric approach.
#' \emph{Multivariate Behavioral Research, 44}, 407-436.
#'
#' @examples
#' \dontrun{
#' plotSEMM_GUI()
#' plotSEMM_GUI(launch.browser=TRUE) #if using RStudio, will launch system browser default
#' }
plotSEMM_GUI <- function(...){
    runApp(plotSEMM_GUI.internal(), ...)
}

plotSEMM_GUI.internal <- function(){

    #custom function to put user text input boxes side-by-side
    numberInputRow<-function (inputId, label, value = "")
    {
        div(style="display:inline-block",
            tags$label(label, `for` = inputId),
            tags$input(id = inputId, type = "number", value = value,class="input-small"))
    }

        ret <- list(

            #--------------------------------------------------------------------
            ui = pageWithSidebar(

                # Application title
                headerPanel("PlotSEMM"),

                sidebarPanel(
                    h5 ('The accompanying Online Appendix for this web application is located at:
                        http://www.yorku.ca/pek/index_files/appendices.htm'),
                    h4('User Input:'),
                    #h5('Note: The \'Mplus Files\' option requires a single
                       #Mplus output file (.out) is in the directory.'),

                    selectInput(inputId="method", label="Type of Input:",
                                choices=c("Mplus Files"="Mplusfile", "Manual Input"="Manually", " "=" "), selected=" "),

                    conditionalPanel(condition = "input.method != ' '",
                                     selectInput(inputId="plottype",label="Type of Output:",
                                                 choices=c("Contour"="contour", 'Probability'='probability',
                                                           "Confidence Bands (Mplus Input Only)"="ci",
                                                           " "=" "),
                                                 selected=" ")
                    ),

                    conditionalPanel(condition = "input.method == 'Mplusfile'",
                                     textInput(inputId='Mpath', label='Directory containing Mplus files:',
                                               value=getwd())
                    ),

                    ###
                    # tech extras for plotting features
                    conditionalPanel(condition = "input.method != ' '",
                                     checkboxInput(inputId='tech_extras',
                                                   label='Override various plotting defaults.',
                                                   value=FALSE)
                    ),

                    conditionalPanel(condition = "input.tech_extras == true",

                                     selectInput(inputId="legend_location", label="Legend location",
                                                 choices=c("default"="default", "bottomright"="bottomright", "bottom"="bottom",
                                                           "bottomleft"="bottomleft", "left"="left", "topleft"="topleft",
                                                           "top"="top", "topright"="topright", "right"="right", "center"="center",
                                                           "none"="none"),
                                                 selected="default"),

                                     textInput(inputId='xlab', label='X-axis label:',
                                               value='Latent Predictor'),

                                     textInput(inputId='ylab', label='Y-axis label:',
                                               value='Latent Outcome'),

                                     numericInput('npoints', 'Number of points to evaluate', 50,
                                                  min = 25, max = 1000),

                                     checkboxInput(inputId='class_info',
                                                   label='Show class specific distributions, regression lines,
                                                   and mixing probabilities for Contour plot.',
                                                   value=TRUE),

                                     checkboxInput(inputId='include_title',
                                                   label='Allow plot title to indicate "No Line was Found within the Confidence Envelope(s)". ',
                                                   value=TRUE),

                                     numericInput('cex', 'Font size (cex):', 1.5,
                                                  min = 1e-2, max = 20),

                                     checkboxInput(inputId='save_data',
                                                   label='Save the data used to plot the graphics to working
                                                   directory in R.',
                                                   value=FALSE),

                                     textInput(inputId='save_filename', label='Saved data file name:',
                                               value='plotSEMM_saved_data.Rdata')

                    ),

                    ###

                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     selectInput(inputId='CI', label='Confidence Level:',
                                                 choices=c("95%", "90%"), selected="95%")
                    ),

                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     numberInputRow(inputId="citable_value",
                                                    label=HTML("Conditional Value of Latent Predictor (&eta;<sub>1</sub>) for Confidence Interval(s) about E[&eta;<sub>2</sub>|&eta;<sub>1</sub>]:"))
                    ),

                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     checkboxInput(inputId='plot_deltaci',
                                                          label='Delta Method Wald-type Confidence Intervals.',
                                                          value=TRUE)
                    ),

                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     checkboxInput(inputId='plot_bsci',
                                                          label='Parametric Bootstrap Confidence Intervals.',
                                                          value=FALSE)
                    ),

                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     checkboxInput(inputId='plot_deltace',
                                                          label='Delta Method Wald-type Confidence Envelope.',
                                                          value=TRUE)
                    ),



                    conditionalPanel(condition = "input.plottype == 'ci'",
                                     checkboxInput(inputId='boot',
                                                          label='Parametric Bootstrap Confidence Envelope.',
                                                          value=FALSE)
                    ),

		     conditionalPanel(condition = "input.plottype == 'ci'",
                                     checkboxInput(inputId='linesearch',
                                               label='Run Line Finding Algorithm.',
                                               value=TRUE)
                    ),

                    #Manual input
                    conditionalPanel(condition = "input.method == 'Manually'",
                                     sliderInput(inputId = "nclass",
                                                 label = "Number of latent classes:",
                                                 min = 1, max = 7, value = 2, step = 1),

                                     hr(),
                                     h5 ('\nNote: Numbers must contain values on both sides of the decimal.
                                         \nE.g., .520 must be input as 0.520'),

                                     hr(),
                                     h6('Class 1:'),
                                     numberInputRow(inputId="pi1", label=HTML("&pi;")),
                                     numberInputRow(inputId="alpha1.1", label=HTML("&alpha;<sub>1</sub>")),
                                     numberInputRow(inputId="alpha2.1", label=HTML("&alpha;<sub>2</sub>")),
                                     numberInputRow(inputId="beta21.1", label=HTML("&beta;<sub>21</sub>")),
                                     numberInputRow(inputId="psi11.1", label=HTML("&psi;<sub>11</sub>")),
                                     numberInputRow(inputId="psi22.1", label=HTML("&psi;<sub>22</sub>")),

                                     conditionalPanel(condition = "input.nclass > 1",
                                                      hr(),
                                                      h6('Class 2:'),
                                                      numberInputRow(inputId="pi2", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.2", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.2", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.2", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.2", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.2", label=HTML("&psi;<sub>22</sub>"))
                                                      ),

                                     conditionalPanel(condition = "input.nclass > 2",
                                                      hr(),
                                                      h6('Class 3:'),
                                                      numberInputRow(inputId="pi3", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.3", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.3", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.3", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.3", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.3", label=HTML("&psi;<sub>22</sub>"))
                                     ),

                                     conditionalPanel(condition = "input.nclass > 3",
                                                      hr(),
                                                      h6('Class 4:'),
                                                      numberInputRow(inputId="pi4", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.4", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.4", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.4", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.4", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.4", label=HTML("&psi;<sub>22</sub>"))
                                     ),

                                     conditionalPanel(condition = "input.nclass > 4",
                                                      hr(),
                                                      h6('Class 5:'),
                                                      numberInputRow(inputId="pi5", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.5", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.5", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.5", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.5", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.5", label=HTML("&psi;<sub>22</sub>"))
                                     ),

                                     conditionalPanel(condition = "input.nclass > 5",
                                                      hr(),
                                                      h6('Class 6:'),
                                                      numberInputRow(inputId="pi6", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.6", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.6", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.6", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.6", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.6", label=HTML("&psi;<sub>22</sub>"))
                                     ),

                                     conditionalPanel(condition = "input.nclass > 6",
                                                      hr(),
                                                      h6('Class 7:'),
                                                      numberInputRow(inputId="pi7", label=HTML("&pi;")),
                                                      numberInputRow(inputId="alpha1.7", label=HTML("&alpha;<sub>1</sub>")),
                                                      numberInputRow(inputId="alpha2.7", label=HTML("&alpha;<sub>2</sub>")),
                                                      numberInputRow(inputId="beta21.7", label=HTML("&beta;<sub>21</sub>")),
                                                      numberInputRow(inputId="psi11.7", label=HTML("&psi;<sub>11</sub>")),
                                                      numberInputRow(inputId="psi22.7", label=HTML("&psi;<sub>22</sub>"))
                                     )

                    ),

                    hr(),

                    submitButton(text = "Update")

                    ), #end sidebarPanel

                #------------------------------------------------------------------------

                mainPanel(
                    plotOutput(outputId = "main_plot", height = '800px', width = "100%")
                    )

                ), #end pageWithSidebar

            #--------------------------------------------------------------------

            server = function(input, output) {

                #preamble function here to grab input data
                GUI_setup <- function(input){

                    nclass <- input$nclass
                    ret <- NULL
                    original_dir <- getwd()
                    on.exit(setwd(original_dir))
                    if(input$method == 'Manually'){

                        pi <- alpha1 <- alpha2 <- beta21 <- psi11 <- psi22 <- numeric(nclass)
                        if(is.na(input$pi1))
                            return(NULL)
                        Names <- names(input)
                        pi.names <- paste0('pi', 1:7)
                        alpha1.names <- paste0('alpha1.', 1:7)
                        alpha2.names <- paste0('alpha2.', 1:7)
                        beta21.names <- paste0('beta21.', 1:7)
                        psi11.names <- paste0('psi11.', 1:7)
                        psi22.names <- paste0('psi22.', 1:7)
                        for(i in 1L:nclass){
                            pi[i] <- input[[Names[which(Names %in% pi.names)[i]]]]
                            alpha1[i] <- input[[Names[which(Names %in% alpha1.names)[i]]]]
                            alpha2[i] <- input[[Names[which(Names %in% alpha2.names)[i]]]]
                            beta21[i] <- input[[Names[which(Names %in% beta21.names)[i]]]]
                            psi11[i] <- input[[Names[which(Names %in% psi11.names)[i]]]]
                            psi22[i] <- input[[Names[which(Names %in% psi22.names)[i]]]]
                        }
                        test <- c(pi, alpha1, alpha2, beta21, psi11, psi22)
                        if(any(is.na(test)))
                            stop('Must include all input values for each class.')
                        if(abs(sum(pi) - 1) > 1e-7)
                            stop('Class probabilities (pi) should sum to 1.')
                        if(any(c(psi11, psi22) < 0))
                            stop('Negative variances supplied. Please fix.')
                        ret <- plotSEMM_setup(pi, alpha1, alpha2, beta21, psi11, psi22, points=input$npoints)
                    } else if(input$method == 'Mplusfile'){

                        if(!is.null(input$Mpath)){
                            setwd(input$Mpath)
                            files <- dir()
                            file <- files[grepl("*\\.out$", files)]
                            if(!length(file))
                                return(NULL)
                            if(length(file) > 1L)
                                stop('Multiple .out files in specifed directory')
                            read <- suppressWarnings(MplusAutomation::readModels(file, recursive=TRUE))
                        } else {
                            return(NULL)
                        }
                        if(input$plottype != 'ci'){
                            ovars <- strsplit(toupper(read$input$variable$usevariables), split = ' ')[[1L]]
                            pi <- read$class_counts$modelEstimated$proportion
                            pars <- read$parameters[[1L]]
                            pars <- pars[!(pars$param %in% ovars), ] #latents only
                            tmp <- min(which(grepl("*\\.ON$", pars$paramHeader)))
                            ON <- pars$paramHeader[tmp]
                            DV <- strsplit(ON, '.ON')[[1L]]
                            IV <- pars$param[tmp]
                            pars <- pars[pars$param %in% c(IV,DV), ]
                            alpha1 <- pars$est[pars$paramHeader == 'Means']
                            alpha2 <- pars$est[pars$paramHeader == 'Intercepts']
                            beta21 <- pars$est[pars$paramHeader == ON]
                            psi11 <- pars$est[pars$paramHeader == 'Variances']
                            psi22 <- pars$est[pars$paramHeader == 'Residual.Variances']
                            if(any(c(psi11, psi22) < 0))
                                stop('Negative variances supplied. Please fix.')
                            ret <- plotSEMM_setup(pi, alpha1, alpha2, beta21, psi11, psi22, points=input$npoints)
                        } else {
                            boot <- input$boot
                            setup <- read.plotSEMM_wACOV(read)
                            ret <- plotSEMM_setup2(setup, boot=read, boot.CE=boot, boot.CI=input$plot_bsci,
                                                   alpha = ifelse(input$CI == "95%", .025, .05),
                                                   points=input$npoints, fixed_value=input$citable_value)
                            if(!is.na(input$citable_value)){
                                customtmp <- ret[nrow(ret), ]
                                ret <- ret[-nrow(ret), ]
                            }
                            if(input$linesearch){
                                lines <- .Call('linear', ret$delta_CElo, ret$delta_CEhi, ret$Eta1)
                                line <- which(rowSums(t(ret$delta_CElo<= t(lines)) &
                                                              t(ret$delta_CEhi >= t(lines))) == ncol(lines))
                                if(length(line)){
                                    line <- min(line)
                                    attr(ret, "search") <- rbind(ret$Eta1, lines[line,])
                                }
                                if(boot){
                                    lines <- .Call('linear', ret$bs_CElo, ret$bs_CEhi, ret$Eta1)
                                    line <- which(rowSums(t(ret$bs_CElo <= t(lines)) &
                                                              t(ret$bs_CEhi >= t(lines))) == ncol(lines))
                                    if(length(line)){
                                        line <- max(line)
                                        attr(ret, "search.bs") <- rbind(ret$Eta1, lines[line,])
                                    }
                                }
                            }
                            if(!is.na(input$citable_value)){
                                ret <- rbind(ret, customtmp)
                            }
                        }
                    }
                    return(ret)
                }

                #-----------------------------------------------------------------

                output$main_plot <- renderPlot({
                    ret <- GUI_setup(input)
                    if(input$save_data) save(ret, file = input$save_filename)
                    if(!is.null(ret)){
                        plottype <- input$plottype
                        if(plottype == 'contour') plotSEMM_contour(ret, input=input,
                                                                   cex=input$cex)
                        if(plottype == 'probability') plotSEMM_probability(ret, input=input,
                                                                           cex=input$cex)
                        if(plottype == 'ci')
                            plotSEMM_ci(ret, linesearch=input$linesearch,
                                        deltaci=input$plot_deltaci,
                                        bsci=input$plot_bsci,
                                        deltace=input$plot_deltace,
                                        ninty_five=input$CI == "95%",
                                        input=input,
                                        include_title=input$include_title,
                                        use_fixed_value=!is.na(input$citable_value),
                                        cex=input$cex)
                    } else examplePlot()
                })

            } #end server function

        ) #end list collector

        return(ret)
}

