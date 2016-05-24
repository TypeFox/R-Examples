library(shiny)
library(shinyAce)
library(shinyBS)
library(meta)
library(metafor)
#library(metamisc)
library(MAd)
library(MAc)
library(quantreg)
library(ggplot2)
library(compute.es)
library(SCMA)
library(SCRT)
shinyUI(navbarPage(title=div(img(src="http://kylehamilton.com/wp-content/uploads/2015/04/mavis1a-e1430059956805.png")), windowTitle="MAVIS v1.1.1",
#shinyUI(navbarPage("MAVIS: Meta Analysis Via Shiny v1.0.4", windowTitle="MAVIS v1.0.4",
# If I want to add a theme here is how to do it
# shinyUI(navbarPage(theme = shinytheme("flatly"),"MAVIS: Meta Analysis Via Shiny v1.0.3",
      tabPanel("Main", icon = icon("stats", lib = "glyphicon"),
               sidebarPanel(
                 
                 radioButtons("type", strong("Data Analysis and Input Options:"),
                              list("Mean Differences (n, M, SD)" = "mdms",
                                   "Mean Differences (n, Effect size d)" = "mdes",
                                   "Correlations (n, r)" = "cor",
                                   "Dichotomous Models"="or"
                                   ),
                 ),
                 bsTooltip("type", "Make sure to review the options on the toolbar before running your analysis",
                           "right", trigger = "click", options = list(container = "body")),
                 helpText("See Dichotomous Model Options, default is set to log odds ratio"),
                 checkboxInput("moderator", label = ("The data contains a categorical moderator (subgroup) variable."), value = T),
                 br(),
                 submitButton("Update View"),
                 helpText("Click here to update your results, you need to do this after you change the data, model, or any of the settings"),
                 br(),
                 actionButton("quit", "Quit", icon("sign-out")),
                 helpText("Press Quit to exit the application")
               
                 
               ),
               br(),

               p('Note: Input values must be separated by tabs. Copy and paste from Excel.'),
               p("Your data needs to have exactly the same header (variable names) in the first row."),
               p("For examples of how this data should look click on the Input Examples tab"),

               aceEditor("text", value="Study\tModerator\tN1\tM1\tSD1\tN2\tM2\tSD2\nStudy 01\tJH\t30\t51.57\t16.5\t30\t72.97\t13.23\nStudy 02\tUNI\t23\t75.09\t23.01\t24\t81.63\t14.42\nStudy 03\tSH\t83\t30.08\t14.29\t81\t35.38\t16.13\nStudy 04\tSH\t21\t2.95\t1.28\t21\t3.48\t0.68\nStudy 05\tSH\t15\t53.8\t17.4\t15\t60.47\t17.37\nStudy 06\tSH\t7\t15.7\t4.1\t7\t27.3\t4.1\nStudy 07\tSH\t28\t27.9\t9.57\t28\t33.2\t15.65\nStudy 08\tUNI\t40\t17.53\t8.87\t40\t19.23\t9.55\nStudy 09\tUNI\t18\t11.86\t13.24\t17\t29.92\t16.67\nStudy 10\tUNI\t21\t29.76\t16\t25\t27.98\t16.52\nStudy 11\tUNI\t26\t8.23\t3.59\t26\t9.65\t2.99\nStudy 12\tUNI\t49\t13.71\t4.07\t48\t16\t3.47\nStudy 13\tUNI\t27\t2.8\t1.7\t27\t5.9\t1.4\nStudy 14\tSH\t41\t10.05\t2.52\t34\t11.03\t1.78\nStudy 15\tUNI\t58\t3.62\t1.79\t57\t4.26\t1.61\nStudy 16\tSH\t60\t7.36\t2.8\t63\t8.82\t2.5\nStudy 17\tUNI\t15\t5.93\t3.55\t15\t12.27\t4.95\nStudy 18\tJH\t37\t13.68\t3.68\t142\t17.53\t4.34\nStudy 19\tJH\t27\t3.3\t2.3\t54\t12.98\t7.67\nStudy 20\tJH\t35\t5.49\t3.88\t39\t12.36\t7.68\nStudy 21\tJH\t32\t5.81\t3.14\t34\t12.44\t5.66\nStudy 22\tJH\t62\t17.84\t4.09\t60\t18.18\t4.09\nStudy 23\tSH\t39\t8.77\t5\t39\t13.72\t5.32\nStudy 24\tSH\t213\t59.8\t15.3\t39\t79.8\t9.5\nStudy 25\tUNI\t34\t14.32\t2.79\t42\t16\t2.05\nStudy 26\tUNI\t77\t70.85\t11.74\t56\t78.17\t9.94\nStudy 27\tUNI\t28\t80.83\t22.47\t28\t85.06\t23.52\nStudy 28\tUNI\t33\t25.38\t4.71\t36\t25.02\t3.36\nStudy 29\tUNI\t66\t0.45\t0.29\t66\t0.93\t0.59",
                         mode="r", theme="mono_industrial"),

               br(),

               h3("Effect size and sampling variance"),

               verbatimTextOutput("data.out"),

               br(),

               h3("Fixed effects model"),
               verbatimTextOutput("fe.out"),

               br(),

               h3("Random effects model"),
               verbatimTextOutput("re.out"),

               p('[Criteria for checking heterogeneity]',
                 br(),
                 br(),
                 'I^2 (How much effect sizes across studies differ)', br(),
                 "25-50: Little different", br(),
                 '50-75: Quite different', br(),
                 '75-100: Considerably different', br(),
                 br(),
                 'Test for Heterogeneity: p-val < .05 (not homogeneous)', br(),
                 br(),
                 'H (sqrt(H^2)) > 1: There is unexplained heterogeneity.'
               ),

               br(),
               br(),


               h3("Forest plot (Fixed effects model)"),
               downloadButton('downloadfePlot', 'Download the plot as pdf'),
               plotOutput("fePlot", height = "550px"),



               br(),

               h3("Forest plot (Random effects model)"),
               downloadButton('downloadrePlot', 'Download the plot as pdf'),
               plotOutput("rePlot", height = "550px"),

               br(),

               h3("Funnel plot (Fixed effects model)"),
               downloadButton('downloadFunFixPlot', 'Download the plot as pdf'),
               plotOutput("FunFixPlot"),
               p('Open circles (if any) on the right side show missing NULL studies estimated with the trim-and-fill method, added in the funnel plot.'),

               br(),

               h3("Funnel plot (Random effects model)"),
               downloadButton('downloadFunRandPlot', 'Download the plot as pdf'),
               plotOutput("FunRandPlot"),
               p('Open circles (if any) on the right side show missing NULL studies estimated with the trim-and-fill method, added in the funnel plot.'),
               br(),
               br(),
               h3("Publication Bias"),
               verbatimTextOutput("asy.out"), # regression tests for funnel plot asymmetry
               p('Fail-safe N is the number of nonsignificant studies necessary to make the result nonsignificant. "When the fail-safe N is high, that is interpreted to mean that even a large number of nonsignificant studies may not influence the statistical significance of meta-analytic results too greatly."',
                 a('(Oswald & Plonsky, 2010, p. 92)', href='http://dx.doi.org/10.1017/S0267190510000115', target="_blank"), '.'),
               br(),

               br(),

               # Display this only if "moderator" is checked
               conditionalPanel(condition = "input.moderator == true",
                                h3("Moderator (subgroup) analysis"),
                                verbatimTextOutput("modAnalysis.out")
               ),

               br(),

               # Display this only if "moderator" is checked
               conditionalPanel(condition = "input.moderator == true",
                                h4("Categorical moderator graph (Fixed effects model)"),
                                plotOutput("ModFixGraph")
               ),

               br(),

               # Display this only if "moderator" is checked
               conditionalPanel(condition = "input.moderator == true",
                                h4("Categorical moderator graph (Random effects model)"),
                                plotOutput("ModRandGraph")
               ),

               br(),
               br(),

               strong('R session and package information'),
               verbatimTextOutput("info.out")
      ),



      tabPanel("Input Examples", icon = icon("table", lib = "font-awesome"),

               p('Note: Input values must be separated by tabs. Copy and paste from Excel.'),

               p(HTML("<b><div style='background-color:#FADDF2;border:1px solid black;'>Your data needs to have exactly the same header (variable names) in the first row.</div></b>")),

               br(),

               p(strong("Mean Differences (n, M, SD)")),
               aceEditor("text1", value="Study\tModerator\tN1\tM1\tSD1\tN2\tM2\tSD2\nStudy 01\tJH\t30\t51.57\t16.5\t30\t72.97\t13.23\nStudy 02\tUNI\t23\t75.09\t23.01\t24\t81.63\t14.42\nStudy 03\tSH\t83\t30.08\t14.29\t81\t35.38\t16.13\nStudy 04\tSH\t21\t2.95\t1.28\t21\t3.48\t0.68\nStudy 05\tSH\t15\t53.8\t17.4\t15\t60.47\t17.37\nStudy 06\tSH\t7\t15.7\t4.1\t7\t27.3\t4.1\nStudy 07\tSH\t28\t27.9\t9.57\t28\t33.2\t15.65\nStudy 08\tUNI\t40\t17.53\t8.87\t40\t19.23\t9.55\nStudy 09\tUNI\t18\t11.86\t13.24\t17\t29.92\t16.67\nStudy 10\tUNI\t21\t29.76\t16\t25\t27.98\t16.52\nStudy 11\tUNI\t26\t8.23\t3.59\t26\t9.65\t2.99\nStudy 12\tUNI\t49\t13.71\t4.07\t48\t16\t3.47\nStudy 13\tUNI\t27\t2.8\t1.7\t27\t5.9\t1.4\nStudy 14\tSH\t41\t10.05\t2.52\t34\t11.03\t1.78\nStudy 15\tUNI\t58\t3.62\t1.79\t57\t4.26\t1.61\nStudy 16\tSH\t60\t7.36\t2.8\t63\t8.82\t2.5\nStudy 17\tUNI\t15\t5.93\t3.55\t15\t12.27\t4.95\nStudy 18\tJH\t37\t13.68\t3.68\t142\t17.53\t4.34\nStudy 19\tJH\t27\t3.3\t2.3\t54\t12.98\t7.67\nStudy 20\tJH\t35\t5.49\t3.88\t39\t12.36\t7.68\nStudy 21\tJH\t32\t5.81\t3.14\t34\t12.44\t5.66\nStudy 22\tJH\t62\t17.84\t4.09\t60\t18.18\t4.09\nStudy 23\tSH\t39\t8.77\t5\t39\t13.72\t5.32\nStudy 24\tSH\t213\t59.8\t15.3\t39\t79.8\t9.5\nStudy 25\tUNI\t34\t14.32\t2.79\t42\t16\t2.05\nStudy 26\tUNI\t77\t70.85\t11.74\t56\t78.17\t9.94\nStudy 27\tUNI\t28\t80.83\t22.47\t28\t85.06\t23.52\nStudy 28\tUNI\t33\t25.38\t4.71\t36\t25.02\t3.36\nStudy 29\tUNI\t66\t0.45\t0.29\t66\t0.93\t0.59", mode="r", theme="monokai"),


               br(),
               p(strong("Mean Differences (n, Effect size d)")),
               aceEditor("text2", value="Study\tModerator\tN1\tN2\td\nStudy 01\tJH\t30\t30\t-1.431\nStudy 02\tUNI\t23\t24\t-0.3423\nStudy 03\tSH\t83\t81\t-0.3481\nStudy 04\tSH\t21\t21\t-0.5171\nStudy 05\tSH\t15\t15\t-0.3837\nStudy 06\tSH\t7\t7\t-2.8293\nStudy 07\tSH\t28\t28\t-0.4086\nStudy 08\tUNI\t40\t40\t-0.1845\nStudy 09\tUNI\t18\t17\t-1.2039\nStudy 10\tUNI\t21\t25\t0.1093\nStudy 11\tUNI\t26\t26\t-0.4298\nStudy 12\tUNI\t49\t48\t-0.605\nStudy 13\tUNI\t27\t27\t-1.9907\nStudy 14\tSH\t41\t34\t-0.4422\nStudy 15\tUNI\t58\t57\t-0.3758\nStudy 16\tSH\t60\t63\t-0.5508\nStudy 17\tUNI\t15\t15\t-1.4719\nStudy 18\tJH\t37\t142\t-0.9136\nStudy 19\tJH\t27\t54\t-1.5079\nStudy 20\tJH\t35\t39\t-1.111\nStudy 21\tJH\t32\t34\t-1.4368\nStudy 22\tJH\t62\t60\t-0.0831\nStudy 23\tSH\t39\t39\t-0.9588\nStudy 24\tSH\t213\t39\t-1.3729\nStudy 25\tUNI\t34\t42\t-0.6976\nStudy 26\tUNI\t77\t56\t-0.6642\nStudy 27\tUNI\t28\t28\t-0.1839\nStudy 28\tUNI\t33\t36\t0.0886\nStudy 29\tUNI\t66\t66\t-1.0326", mode="r", theme="monokai"),


               br(),
               p(strong("Correlations (n, r)")),
               aceEditor("text3", value="Study\tN\tr\tModerator\nIzumi (2000)\t175\t0.78\tcollege\nYu (2009)\t53\t0.38\tJS high\nThuy (1996)\t250\t0.69\tcollege\nOckey (2002)\t90\t0.89\tcollege\nAraru (2005)\t86\t0.52\tJS high\nWee (1997)\t182\t0.98\tcollege\nOzoda (2007)\t591\t0.91\tcollege\nHala (2004)\t30\t0.95\tcollege\nTapio (2008)\t37\t0.47\tJS high\nAndarani (2008)\t107\t0.84\tcollege\nDavis (1999)\t74\t0.99\tcollege\nPlonsky (2002)\t217\t0.86\tcollege\nGassel (1993)\t203\t0.99\tcollege",mode="r", theme="monokai"),

               br(),
               p(strong("Dichotomous (upoz, uneg, NU, kpoz, kneg, NK)")),
               aceEditor("text4", value="Study\tModerator\tupoz\tuneg\tNU\tkpoz\tkneg\tNK\nStudy 01\tra\t4\t119\t123\t11\t128\t139\nStudy 02\tra\t6\t300\t306\t29\t274\t303\nStudy 03\tra\t3\t228\t331\t11\t209\t220\nStudy 04\tsys\t17\t1699\t1716\t65\t1600\t1665\nStudy 05\tsys\t5\t2493\t2498\t3\t2338\t2341\nStudy 06\tra\t29\t7470\t7499\t45\t7232\t7277", mode="r", theme="monokai"),
               
               br()

      ),

navbarMenu("Model Options and Settings", icon = icon("cog", lib = "font-awesome"),
#                       tabPanel("Bayesian Model Options", icon = icon("tasks", lib = "font-awesome"),
#                                
#                      strong('Bayesian Analysis Options'),
#                      selectInput("bayoption1", label = h3("Run Bayesian Analysis"), 
#                                  choices = list("No" = FALSE, "Yes" = TRUE)
#                                
#                       )),
           tabPanel("Correlation model options", icon = icon("line-chart", lib = "font-awesome"),
                    
                    radioButtons("cormeasures", strong("Correlation model measures"),
                                 c("raw correlation coefficient" = "COR",
                                   "raw correlation coefficient corrected for its slight negative bias" = "UCOR",
                                   "Fisher’s r-to-z transformed correlation coefficient" = "ZCOR"
                                 ), selected = "ZCOR"),
                    p(h6('Fisher’s r-to-z transformed correlation coefficient is the default estimator for the metafor package')),
                    
                    verbatimTextOutput('cormeasures.out')
                    
           ),
           tabPanel("Dichotomous Model Options", icon = icon("ellipsis-v", lib = "font-awesome"),
                    
                    radioButtons("dichotomousoptions", strong("Measure Selection"),
                                 c("log relative risk" = "RR",
                                   "log odds ratio" = "OR",
                                   "risk difference" = "RD",
                                   "arcsine square-root transformed risk difference (Rücker et al., 2009)." = "AS",
                                   "log odds ratio estimated with Peto’s method (Yusuf et al., 1985)." = "PETO",
                                   "probit transformed risk difference as an estimate of the standardized mean difference." = "PBIT",
                                   "transformed odds ratio as an estimate of the standardized mean difference (normal distributions)." = "OR2DN",
                                   "transformed odds ratio as an estimate of the standardized mean difference (logistic distributions)." = "OR2DL"
                                 ), selected = "OR"),
                    p(h6('logs odds ratio is the default option and is the one you should use for the example provided in the Input Examples tab.')),
                    
                    verbatimTextOutput('dichotomousoptions.out')
                    
           ),
           tabPanel("Random-effects model estimators", icon = icon("random", lib = "glyphicon"),
                    
                    radioButtons("model", strong("Random-effects model estimators"),
                                 c("DerSimonian-Laird" = "DL",
                                   "Hedges" = "HE",
                                   "Hunter-Schmidt" = "HS",
                                   "Sidik-Jonkman" = "SJ",
                                   "Maximum-likelihood" = "ML",
                                   "Restricted maximum-likelihood" = "REML",
                                   "Empirical Bayes (Paule-Mandel)" = "EB",
                                   "Generalized Q-statistic" = "GENQ"
                                 ), selected = "REML"),
                    p(h6('Restricted maximum-likelihood is the default estimator for the metafor package')),
                    
                    checkboxInput("khadjust", label = "Knapp & Hartung Adjustment", value = FALSE),
                    p(h6('Knapp & Hartung Adjustment is turned off by the default in the metafor package')),
                    p("The Knapp and Hartung (2003) method is an adjustment to the standard errors of the estimated
coefficients, which helps to account for the uncertainty in the estimate of the amount of
(residual) heterogeneity and leads to different reference distributions."),
                    h3("References"),
                    p("Knapp, G., & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. Statistics in Medicine, 22, 2693–2710.")
                    
                    
           )),

 navbarMenu("Publication Bias", icon = icon("book", lib = "font-awesome"),
      tabPanel("Trim and Fill Options", icon = icon("chevron-right", lib = "font-awesome"),

               radioButtons("trimfillopt", strong("Trim and Fill Estimator"),
                            c("L0" = "L0",
                              "R0" = "R0",
                              "Q0" = "Q0"
                            ), selected = "L0"),
               p(h6('Three different estimators for the number of missing studies were proposed by Duval and Tweedie (2000a, 2000b; see also Duval, 2005). The default estimator for the metafor package is L0')),
               
               verbatimTextOutput('trimfillopt.out'),
               h3("References"),
               p("Duval, S. J., & Tweedie, R. L. (2000a). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics, 56, 455–463."),
               p("Duval, S. J., & Tweedie, R. L. (2000b). A nonparametric trim and fill method of accounting for publication bias in meta-analysis. Journal of the American Statistical Association, 95, 89–98."),
               p("Duval, S. J. (2005). The trim and fill method. In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.) Publication bias in meta-analysis: Prevention, assessment, and adjustments (pp. 127–144). Chichester, England: Wiley.")
               

      ),

tabPanel("Tests for Funnel Plot Asymmetry", icon = icon("chevron-right", lib = "font-awesome"),
verticalLayout(
  
  wellPanel(
    fluidRow(
      column(4,
             p(strong("Regression Test Options")),
                            radioButtons("regtestpredictor", strong("Predictor"),
                                         c("standard error" = "sei",
                                           "sampling variance" = "vi",
                                           "sample size" = "ni",
                                           "inverse of the sample size" = "ninv"
                                         ), selected = "sei"),
                            radioButtons("regtestmodeltype", strong("Model Selection"),
                                          c("Weighted Regression with a Multiplicative Dispersion Term" = "lm",
                                            "Meta-analytic Models" = "rma"
                                          ), selected = "lm"),
             
             
             p(br())
      ),
      column(4,
             strong('Funnel Plot Options'),
             checkboxInput("contourenhancedbox", "Contour enhanced plots", FALSE),
             helpText("Check this box if you would like to have your funnel plots contour enhanced see (Peters et al., 2008)"),
             checkboxInput("regtestfullmodel", "Results from the fitted model", FALSE),
             helpText("Check this box if you would like to see the full results from the fitted model")
             
      )
      
    )),
  p("For more information about the different methods of detecting publication bias in a meta-analysis see (Jin, Zhou, & He, 2015)"),
  h3("References"),
  p("Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. British Medical Journal, 315, 629--634."),
  p("Jin, Zhi-Chao, Zhou, Xiao-Hua & He, Jia (2015). Statistical methods for dealing with publication bias in meta-analysis. Statistics in Medicine, 34, 343-360."),
  p("Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2008). Contour-enhanced meta-analysis funnel plots help distinguish publication bias from other causes of asymmetry. Journal of Clinical Epidemiology, 61(10), 991–-996."),
  p("Sterne, J. A. C., & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: Guidelines on choice of axis. Journal of Clinical Epidemiology, 54(10), 1046--1055."),
  br()
  
)

),
      
      tabPanel("File Drawer Analysis", icon = icon("chevron-right", lib = "font-awesome"),
               
               radioButtons("filedraweranalysis", strong("File Drawer Analysis"),
                            c("Rosenthal" = "Rosenthal",
                              "Orwin" = "Orwin",
                              "Rosenberg" = "Rosenberg"
                            ), selected = "Rosenthal"),
               p(h6('Method for running file drawer analysis. The default in the metafor package is Rosenthal')),
               p('The Rosenthal method (sometimes called a ‘file drawer analysis’) calculates the number of studies
averaging null results that would have to be added to the given set of observed outcomes to reduce
the combined significance level (p-value) to a target alpha level (e.g., .05). The calculation is based
on Stouffer’s method to combine p-values and is described in Rosenthal (1979).'),
p('The Orwin method calculates the number of studies averaging null results that would have to be
added to the given set of observed outcomes to reduce the (unweighted) average effect size to a
target (unweighted) average effect size. The method is described in Orwin (1983).'),
p('The Rosenberg method calculates the number of studies averaging null results that would have to be
added to the given set of observed outcomes to reduce significance level (p-value) of the (weighted)
average effect size (based on a fixed-effects model) to a target alpha level (e.g., .05). The method is
described in Rosenberg (2005).'),
               
               verbatimTextOutput('filedraweranalysis.out'),
h3("References"),
p("Rosenthal, R. (1979). The file drawer problem and tolerance for null results. Psychological Bulletin, 86, 638--641."),
p("Orwin, R. G. (1983). A fail-safe N for effect size in meta-analysis. Journal of Educational Statistics, 8, 157--159."),
p("Rosenberg, M. S. (2005). The file-drawer problem revisited: A general weighted method for calculating fail-safe numbers in meta-analysis. Evolution, 59, 464--468.")

               
      )),
navbarMenu("Effect Size Calculator", icon = icon("calculator", lib = "font-awesome"),
          tabPanel("One Study with Means, SDs, Ns", icon = icon("chevron-right", lib = "font-awesome"),
         
                   verticalLayout(
                     
                     wellPanel(
                       fluidRow(
                         column(3,
                      p(strong("Group 1:")),
                     
                     numericInput("nx", " Sample size (n)", 21),
                     
                     numericInput("mx", " Mean", 61.33),
                     
                     numericInput("sdx", " SD", 16.43),
                     
                     p(br())
                         ),
                     column(4, offset = 1,
                     p(strong("Group 2:")),
                     
                     numericInput("ny", " Sample size (n)", 24),
                     
                     numericInput("my", " Mean", 59.79),
                     
                     numericInput("sdy", " SD", 18.50),
                     
                     p(br())
                       ),
                     column(4,
                     strong('Option:'),
                     
                     
                     checkboxInput("varequal", "t-test with equal variances assumed", FALSE),
                     
                     
                     checkboxInput("vartest", "Show test for equality of variances", FALSE),
                     helpText("Click here to update your results."),
                     submitButton("Update View")
                     )
                     
                   )),
                   
                     h3("Checking the input data"),
                     tableOutput("values"),
                                
                    br(),
                                
                     h3("Mean of the differences and 95% CI"),
                     verbatimTextOutput("difference.out"),
                                
                     br(),
                                
                     h3("t-test"),
                     verbatimTextOutput("ttest.out"),
                     h3(""),
                     verbatimTextOutput("vartest.out"),
                                
                     br(),
                                
                     h3("Effect size indices"),
                     verbatimTextOutput("es.out"),                   
                     br()
                                
                       )
         
          ),
          tabPanel("ANCOVA F-statistic to Effect Size", icon = icon("chevron-right", lib = "font-awesome"),
                   verticalLayout(
                     
                     wellPanel(
                       fluidRow(
                         column(3,
                                p(strong("ANCOVA F-statistic to Effect Size")),
                                
                                numericInput("ancovaf", " F value from ANCOVA", 21),
                                
                                numericInput("ancovafn1", " Treatment group sample size", 50),
                                
                                numericInput("ancovafn2", " Comparison group sample size", 50),
                                
                                p(br())
                         ),
                         column(4, offset = 1,
                                
                                numericInput("anovafcovar", " Covariate outcome correlation or multiple correlation", 0.24),
                                
                                numericInput("anovafcovarnum", " Number of covariates", 3),
                                
                                numericInput("sdy", " SD", 18.50),
                                helpText("Click here to update your results"),
                                submitButton("Update View"),
                                p(br())
                         )

                         
                       )),
                     
                     
                     h3("Effect size indices"),
                     verbatimTextOutput("ancovaf.out"),
                     p(br())
                     
                   )
                   
          ),
          tabPanel("Mean Values from ANCOVA F-statistic to Effect Size", icon = icon("chevron-right", lib = "font-awesome"),
                   verticalLayout(
                     
                     wellPanel(
                       fluidRow(
                         column(3,
                                p(strong("Mean Values from ANCOVA F-statistic to Effect Size")),
                                
                                numericInput("ancovamean1", " Adjusted mean of treatment group from ANCOVA", 21.7),
                                
                                numericInput("ancovamean2", " Adjusted mean of comparison group from ANCOVA", 33.5),
                                
                                numericInput("ancovameansd", " Adjusted standard deviation", 50),
                                
                                p(br())
                         ),
                         column(4, offset = 1,
                                
                                numericInput("ancovameann1", " Treatment group sample size", 142),
                                
                                numericInput("ancovameann2", " Comparison group sample size", 133),
                                
                                numericInput("ancovameancovar", " Covariate outcome correlation or multiple correlation", 0.24),

                                numericInput("ancovameancovarnumber", " Number of covariate", 3),
                                
                                
                                helpText("Click here to update your results"),
                                submitButton("Update View"),
                                p(br())
                         )
                         
                         
                       )),
                     
                     
                     h3("Effect size indices"),
                     verbatimTextOutput("ancovamean.out"),
                     p(br())
                     
                   )
                   
          ),
tabPanel("Chi-Squared Statistic to Effect Size", icon = icon("chevron-right", lib = "font-awesome"),
         verticalLayout(
           
           wellPanel(
             fluidRow(
               column(3,
                      p(strong("Chi-Squared Statistic to Effect Size")),
                      
                      numericInput("chisquaredstat", " Chi squared statistic from primary study.", 5.3),
                      
                      numericInput("chisquaredn1", " Sample size in primary study.", 50),
                      
                      p(br())
               ),
               column(4, offset = 1,
                      helpText("Click here to update your results"),
                      submitButton("Update View"),
                      p(br())
               )
               
               
             )),
           
           
           h3("Effect size indices"),
           verbatimTextOutput("chisquared.out"),
           p(br())
           
         )
         
),
tabPanel("Outcome Measures for Individual Groups", icon = icon("chevron-right", lib = "font-awesome"),
         
         verticalLayout(
           
           wellPanel(
             fluidRow(
               column(3,
                      p(strong("Dichotomous Variables")),
                      
                      numericInput("xi", "Frequencies of the event of interest", 6),
                      
                      numericInput("ni", "Sample size", 323),
                      
                      #numericInput("n1i", "Total", 121),
                      
                      p(br())
               ),
               column(4,
                      radioButtons("divari1", strong("Measure Selection"),
                                   c("raw proportion" = "PR",
                                     "log transformed proportion" = "PLN",
                                     "logit transformed proportion (i.e., log odds)" = "PLO",
                                     "arcsine square-root transformed proportion (i.e., the angular transformation)" = "PAS",
                                     "Freeman-Tukey double arcsine transformed proportion (Freeman & Tukey, 1950)." = "PFT"
                                   ), selected = "PR"),
                      submitButton("Update View")
               )
               
             )),
           
           
           h3("Effect Size Estimates and Corresponding Sampling Variances"),
           verbatimTextOutput("divari1.out"),
           
           br()
           
         )
         
),
tabPanel("Outcome Measures for Two-Group Comparisons", icon = icon("chevron-right", lib = "font-awesome"),
         
         verticalLayout(
           
           wellPanel(
             fluidRow(
               column(3,
                      p(strong("Group 1:")),
                      
                      numericInput("ai", "Outcome 1", 100),
                      
                      numericInput("bi", "Outcome 2", 21),
                      
                      #numericInput("n1i", "Total", 121),
                      
                      p(br())
               ),
               column(4, offset = 1,
                      p(strong("Group 2:")),
                      
                      numericInput("ci", "Outcome 1", 120),
                      
                      numericInput("di", "Outcome 2", 67),
                      
                      #numericInput("n2i", "Total", 187),
                      
                      p(br())
               ),
               column(4,
                      radioButtons("twoxtwovalue", strong("Measure Selection"),
                                   c("log relative risk" = "RR",
                                     "log odds ratio" = "OR",
                                     "risk difference" = "RD",
                                     "arcsine square-root transformed risk difference (Rücker et al., 2009)." = "AS",
                                     "log odds ratio estimated with Peto’s method (Yusuf et al., 1985)." = "PETO"
                                   ), selected = "OR"),
                      submitButton("Update View")
               )
               
             )),
           
           
           h3("Effect Size Estimates and Corresponding Sampling Variances"),
           verbatimTextOutput("twobytwogroups.out"),
           
           br()
           
         )
         
),
tabPanel("p-value to Effect Size", icon = icon("chevron-right", lib = "font-awesome"),
         verticalLayout(
           
           wellPanel(
             fluidRow(
               column(3,
                      p(strong("p-value to Effect Size")),
                      
                      numericInput("pvaluenum", " p-value.", 0.01),
                      numericInput("pvaluen1", " Sample size of treatment group.", 50),                      
                      numericInput("pvaluen2", " Sample size of comparison group.", 50),
                      
                      radioButtons("pvaluetail", strong("One or two-tailed p-value."),
                                   c("One tail" = "one",
                                     "Two tail" = "two"
                                   ), selected = "two"),
                      
                      p(br())
               ),
               column(4, offset = 1,
                      helpText("Click here to update your results"),
                      submitButton("Update View"),
                      p(br())
               )
               
               
             )),
           
           
           h3("Effect size indices"),
           verbatimTextOutput("pvaluees.out"),
           p(br()),
           
           br()
           
         )
         
),
         tabPanel("Single Case Designs", icon = icon("chevron-right", lib = "font-awesome"),
                  verticalLayout(
                    
                    wellPanel(
                      fluidRow(
                        column(3,
                               p(strong("Single Case Design Type")),
                               
                               radioButtons("SCDtype", strong("Type of Single Case Design"),
                                            c("AB" = "AB",
                                              "ABA" = "ABA",
                                              "ABAB" = "ABAB",
                                              "Completely Random Design" = "CRD",
                                              "Randomized Block Design" = "RBD",
                                              "Alternating Treatments Design" = "ATD",
                                              "Multiple-baseline AB design" = "MBD"
                                            ), selected = "AB"),
                               radioButtons("SCDes", strong("Effect Size"),
                                            c("Standardized Mean Difference" = "SMD",
                                              "Pooled Standardized Mean Difference" = "SMDpool",
                                              "Percentage of Nonoverlapping Data (Positive)" = "PND+",
                                              "Percentage of Nonoverlapping Data (Negative)" = "PND-",
                                              "Percentage of Data Points Exceeding the Median (Positive)" = "PND+",
                                              "Percentage of Data Points Exceeding the Median (Negative)" = "PND-"
                                            ), selected = "SMD"),
                               helpText("Click here to update your results"),
                               bsAlert("alert"),
                               submitButton("Update View"),
                               
                               p(br())
                        ),
                        p(strong("Single Case Design Data Entry")),
                        p("The left column should contain the condition labels and the right column should contain the obtained scores"),
                        aceEditor("SCDdata", value="A, 9.523465\nA, 12.371462\nA, 13.265618\nA, 10.182837\nA, 10.987079\nA, 8.161392\nA, 10.655287\nA, 9.563863\nA, 9.381336\nA, 8.822936\nA, 10.227932\nA, 11.961484\nA, 9.425201\nA, 12.199128\nB, 16.212489\nB, 17.657583\nB, 18.45166\nB, 16.645105\nB, 14.618445\nB, 15.769643\nB, 16.017145\nB, 14.000921\nB, 17.081538\nB, 14.06722\nB, 20.423526\nB, 14.123096\nB, 16.728538", mode="r", theme="terminal"),
                        p("Below is your computed effect size, unless you've selected either Percentage of Nonoverlapping Data or Percentage of Data Points Exceeding the Median in which case the number below is the percentage."),
                        verbatimTextOutput('SCDES.out'),
                        p(br()),
                        h3("References"),
                        p("Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467--478."),
                        p("Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477--485.")
                        
                        
                        
                      )),
   
                    p(br())
                    
                  )
         
)),
navbarMenu("About MAVIS", icon = icon("dot-circle-o", lib = "font-awesome"),
           tabPanel("About MAVIS", icon = icon("bar-chart-o", lib = "font-awesome"),
                    
                    HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2015/04/mavis3.png" alt="" style="float: top; margin-right:5px" /></div>'),
                    br(),
                    strong('About MAVIS'),
                    p("MAVIS was designed from the beginning to help users run a meta-analysis as effortlessly as possible. 
                      The software accomplishes this by leveraging the R programming language for data analysis and the Shiny 
                      package from RStudio to power the user interface and server software. These two things combined give 
                      MAVIS a positive user experience with an easy to use interface along with the power of R to provide 
                      the best possible user experience."),
                    br(),
                    strong("MAVIS Version 1.1.1"),
                    p("Last Updated July 20th 2015"),
                    p("Number of monthly downloads from CRAN"),
                    img(src = "http://cranlogs.r-pkg.org/badges/MAVIS", seamless=NA),
                    
                    
                    br()
                    
           ),
           tabPanel("Authors and Contributors", icon = icon("users", lib = "font-awesome"),
                    
                    strong('Acknowledgments'),
                    
                    p('William Kyle Hamilton would like to thank the ',
                      a("Health Communications and Interventions Lab at the University of California, Merced", href="http://cameronhcilab.com/", target="_blank"),
                      'for their comments and beta testing efforts on this application ', 'as well as',
                      a("Kathleen Coburn", href="http://psychology.ucmerced.edu/content/kathleen-coburn", target="_blank"),
                      'for her feedback and evaluation of the statistical methods related to this project.'),
                    
                    p('Atsushi Mizumoto would like to thank',
                      a("Dr. Luke Plonsky", href="http://oak.ucc.nau.edu/ldp3/", target="_blank"), 'and',
                      a("Dr. Yo In'nami", href="https://sites.google.com/site/yoinnami/", target="_blank"),
                      'for their support and feedback to create this web application.'),
                    
                    br(),
                    
                    strong('Authors'),
                    
                    HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/kyle80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                    p(a("William Kyle Hamilton - University of California, Merced", href="http://www.kylehamilton.com", target="_blank")),
                    p("William Kyle Hamilton maintains this application and has authored new features."),
                    
                    br(),
                    HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/atsushi80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                    p(a("Atsushi Mizumoto, PhD - Kansai University", href="http://mizumot.com", target="_blank"),br(),
                      p("Atsushi Mizumoto wrote the first version of this application; this application is a fork of the original which can be found", a("here.", href="https://github.com/mizumot/meta", target="_blank"))
                      
                      
                    ),

                    br(),
                    
                    strong('Contributors and Translators'),
                    
                    HTML('<div style="clear: left;"><img src="http://oi59.tinypic.com/2mnrcci.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                    p(a("Burak Aydin, PhD - Recep Tayyip Erdoğan University", href="http://akademisyen.erdogan.edu.tr/akademisyen.php?uyeid=827a0e170c32e5ce6e7b31ebda784148", target="_blank"),br(),
                      p("Burak Aydin is working on a Turkish version of MAVIS and contributed the dichotomous data entry feature.")
                    ),
                    
                    br(),
                    br(),
                    HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2015/04/katie80.png" alt="" style="float: left; margin-right:5px" /></div>'),
                    p(a("Kathleen Coburn - University of California, Merced", href="http://psychology.ucmerced.edu/content/kathleen-coburn", target="_blank"),br(),
                      p("Kathleen Coburn contributed technical advice on how to run a meta-analysis as well as information on publication bias.")
           ),
           br()
           ),
           tabPanel("Bug Reports", icon = icon("bug", lib = "font-awesome"),
                    
                    strong('Bug Reports'),

                    p("If you discover a problem with MAVIS please submit it to the project GitHub page", 
                      a("https://github.com/kylehamilton/MAVIS/issues", href="https://github.com/kylehamilton/MAVIS/issues", target="_blank"),br()),

                    p("MAVIS is an Open Source project, you are more than welcome to submit patches or features and help the project grow."),
                    
                    
                    br()
                    
           ),
           tabPanel("Feedback", icon = icon("comments", lib = "font-awesome"),
                    
                    strong('Feedback about MAVIS'),
                    
                    p("Feedback about your MAVIS experience is always welcome and highly encouraged!"),
                    p("Feel free to contact the project maintainer with any questions, user experiences, uses of MAVIS, or
                       feature requests at kyle.hamilton@gmail.com"),
                    
                    br()
                    
           ),
           
           tabPanel("License", icon = icon("legal", lib = "font-awesome"),
                    
                    strong('License'),
                    
                    p("MAVIS: Meta Analysis via Shiny"),
                    p(" Copyright 2015  William Kyle Hamilton and Atsushi Mizumoto"),

                    p(" This program is free software you can redistribute it and or modify
                      it under the terms of the GNU General Public License as published by
                      the Free Software Foundation either version 3 of the License or
                      at your option any later version."),

                    p("This program is distributed in the hope that it will be useful,
                      but WITHOUT ANY WARRANTY; without even the implied warranty of
                      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                      GNU General Public License for more details."),

                    p("You should have received a copy of the GNU General Public License
                      along with this program.  If not, see", a("http://www.gnu.org/licenses/gpl.html", href="http://www.gnu.org/licenses/gpl.html", target="_blank"),br()),
                    img(src = "http://www.gnu.org/graphics/gplv3-127x51.png", seamless=NA),
                    
                    
                    br(),
                    
                    strong('Futher Infomation'),
                    p("If you would like to learn more about the GNU General Public License and what it means tl'dr legal has a simple explaination which can be found here", a("https://www.tldrlegal.com/l/gpl-3.0", href="https://www.tldrlegal.com/l/gpl-3.0", target="_blank"),br()),
                    

                    
                    br()
                    
           ),
           
           tabPanel("Support", icon = icon("chevron-right", lib = "font-awesome"),
                    
                    
                    strong('Support'),
                    
                    p("If you're having problems with MAVIS feel free to refer to our GitHub wiki or the documentation available on CRAN."),
                    a("CRAN page for MAVIS", href="http://cran.r-project.org/web/packages/MAVIS/index.html", target="_blank"),
                    br(),
                    a("GitHub Wiki page for MAVIS", href="https://github.com/kylehamilton/MAVIS/wiki", target="_blank"),
                    br(),
                    p("As always you are more than welcome to contact the project maintainer at kyle.hamilton@gmail.com"),
                    br()
                    
                    )),

#This is just so I can get ui.R to run, I'll fix this later
tabPanel(" ",
         h5(" ")
),

p(br())

)
)


