library("shiny")
# install.packages("shinyBS")
library("shinyBS")


shinyUI(fluidPage(theme = "bootstrap.css",

  titlePanel("The Vevea and Hedges Weight-Function Model for Publication Bias"),
  sidebarLayout(
    sidebarPanel(
      "Choose a .csv or .txt file:", br(),
      fileInput('file1', ' ',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv','.txt')),
      #tags$hr(),
      "Does your data file contain a header?",
      checkboxInput('header','Yes', TRUE),
      "Are your data separated by commas, semicolons, tabs, or spaces?",
      radioButtons('sep', ' ',
                   c(Commas=',',
                     Semicolons=';',
                     Tabs='\t',
                     Spaces=' '),
                   ','),
      "For columns of your data file including text, should quotes be included?",
      radioButtons('quote', ' ',
                   c('No'='',
                     'Double Quotes'='"',
                     'Single Quotes'="'"),
                   ''),p(),
      #conditionalPanel(
      #  condition = "!is.null(filedata())",
#       #bootstrapPage(div(
#         #class="container-fluid",
#         #div(class="row-fluid",
#             #div(class="span6",
#                 h2("Source Code"),
#                 aceEditor("code", mode="r", value="x <- c(1, 2, 3); x"),
#                 actionButton("eval", "Evaluate"),
#             #),
#             #div(class="span6",
#                 #h2("Output"),
#                 #verbatimTextOutput("output"),
#             #),
#         #)
#       #)),
#       strong("Please enter the column numbers of the variables containing your effect sizes and sampling variances, in that order, separated by commas. Example: '1,2'"),p(),
#       textInput('effects', ' ', value="0"),p(),
      uiOutput("selecteffects"),p(),
      uiOutput("selectvariances"),p(),
      #strong("How many moderators would you like to include?"),p(),
      uiOutput("selectmods"),p(),
      #numericInput('npred', label = " ", value=0,
      #             min = 0, max = Inf),p(),
      #strong("Please enter the column numbers of these moderators, separated by commas. Example: '3,4'"),p(),
      #textInput('moderators', ' ', value="0"),p(),
      #strong("Select at least one p-value cutpoint to include in your model. To include an interval not provided, type it in and press enter."),
#p(),
      uiOutput("selectsteps"),p(),
      strong("Estimate the Vevea and Woods (2005) model?"),p(),
      checkboxInput("woods",
              label = "Yes", value = FALSE),p(),
      conditionalPanel(
        strong("Enter a pre-specified weight for each of the p-value intervals."),p(),
        condition = "input.woods",
        uiOutput("presetweights")
      ),p(),
      uiOutput("selectdigits")
#       selectizeInput(
#         'steps', ' ', choices=c(0.001,
#                                                             0.005,
#                                                             0.010,
#                                                             0.020,
#                                                             0.025,
#                                                             0.050,
#                                                             0.100,
#                                                             0.200,
#                                                             0.250,
#                                                             0.300,
#                                                             0.350,
#                                                             0.500,
#                                                             0.600,
#                                                             0.650,
#                                                             0.700,
#                                                             0.750,
#                                                             0.800,
#                                                             0.900,
#                                                             0.950,
#                                                             0.975,
#                                                             0.980,
#                                                             0.990,
#                                                             0.995,
#                                                             0.999),
#         multiple=TRUE,
#         selected=c(0.050), options=list(create=TRUE,openOnFocus=TRUE))
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("About", icon=icon("book", lib="font-awesome"),
                    #strong('What is the Vevea and Hedges Weight-Function Model?'), p("This model was introduced by Jack L. Vevea and Larry V. Hedges in 1995. It was ..."), strong('How to Use the Model'), p("Use the panel on the left to read in your data file. While doing so, you can click on the tab labeled 'Data File' to verify that it has been uploaded correctly. You can read in your data as either a .txt file or a .csv file. You can also specify whether your data file is comma-, tab-, or space-delimited."), p("Once you are satisfied that your data was read correctly, scroll down and enter the numbers of the columns holding your effect sizes and sampling variances, in that order. The columns are numbered from left to right. For instance, if your effects are in the first column and your variances in the second, you should enter '1, 2' (without quotes.) Make sure that the first number represents the column of effects and the second the column of variances."), p("Next, you can specify the model. If you would like to include moderators, enter the number of moderators you'd like to include (i.e., if one moderator, enter '1' and so on). Then enter the column numbers of these moderators; if they are in columns three and four, enter '3, 4' (again, without quotes)."), p("Lastly, check the boxes corresponding to the p-value intervals you want to estimate. If you make no changes, the default model is estimated with two intervals (p-values less than 0.05, and p-values greater than 0.05 but less than 1.00 -- in other words, a distinction between significant vs. non-significant p-values). Note that you can assess one-tailed or two-tailed publication bias by checking the appropriate boxes; for more information, see Vevea & Hedges (1995) or Vevea & Woods (2005)."), p("Now that you've entered the data and specified your model, you can click on the 'Unadjusted Random-Effects Model' and 'Adjusted Random-Effects Model' tabs to view your results. There are several components of interest. The unadjusted model will present you with a variance component, intercept, and any moderator estimates from a random- or mixed-effects model that has not been adjusted for publication bias. The adjusted model will present you with adjusted estimates, in addition to weights that correspond to your specified p-value intervals. For instance, if you specified the default intervals of 0.05 and 1.00, the first interval will correspond to p-values less than 0.05, and the second to p-values between 0.05 and 1.00. The adjusted weight-function model always fixes the first weight to one for estimation purposes. This means your output will only contain one weight, which can be interpreted in relation to the first weight. In other words, if you specified the above intervals and obtained a weight of 0.50, your results indicate that effects with p-values greater than 0.05 are half as likely to survive the selection process as significant effects."), p("Keep in mind that all the weights in your output should be interpreted relative to the first weight, which is fixed to one -- i.e., if a weight is estimated at 1.76, p-values in that interval are 1.76 times as likely to survive the selection process as p-values in your first interval. If there are no observed effects or less than three observed effects in one of your p-value intervals, the program will print a warning. This is because the weights may be biased if there are not enough effect sizes in each interval. If this is the case, try re-specifying the p-value intervals."), p("The last part of the 'Adjusted Model' output that is of interest is the likelihood ratio test. This test compares the adjusted model to the unadjusted model, with degrees of freedom corresponding to the number of intervals estimated, and indicates whether the adjusted model fits your data significantly better than the unadjusted model -- that is, whether significant information is gained by estimating the selection model and adjusting the mean and moderator estimates."), strong('Concluding Thoughts'), p("Keep in mind that, as with all other assessments of publication bias, the results this model produces should be considered primarily as a sensitivity analysis. It may be a good idea to run this model multiple times on your data, specifying different p-value intervals each time, and to assess the magnitude of the adjusted mean estimate each time. For instance, if the unadjusted model produces a mean effect of d = 0.50 and, using the adjusted model, you obtain means ranging from d = 0.01 to d = 0.85, this is a sign that your results may not be robust to publication bias. Again, this, like all other assessments, is a sensitivity analysis; it will not allow you to firmly say that publication bias is or is not present, but it will allow you to say, 'If publication bias were present, my data would likely be robust to its effects.'"),
                    p("To use this model, follow the steps in the panel on the left to read in your data file. You can click the tab labeled 'Data File' to ensure your data were read in correctly. Use the drop-down menus to select the columns containing your effect sizes and sampling variances. Note that your data file must contain a column of effect sizes and their corresponding sampling variances for estimation to work. Lastly, select any moderator variables and p-value cutpoints to include in your model. The tabs labeled 'Funnel Plot' and 'Model Results' will then display a funnel plot of your effect sizes and the model output, respectively. For more information on use of the model, we refer you to Vevea and Hedges (1995) and Hedges and Vevea (1996)."),p("We advise users to keep in mind that, as with all other assessments of publication bias, the use of this model should be primarily considered a sensitivity analysis. It is a good idea to run this model multiple times on your data, specifying different p-value cutpoints each time, and to assess the magnitude of the adjusted mean estimate. For instance, if the unadjusted model produces a mean effect of d = 0.50 and, using the adjusted model, you obtain means ranging from d = -0.01 to d = 0.85, this is a sign that publication bias may be present in your data."), p("If you have questions or comments about using or interpreting the model, or if you have a bug or issue to report about this Shiny app, please feel free to send an email to Dr. Jack Vevea or Kathleen Coburn."), br(), strong('Authors:'), p(a("Dr. Jack L. Vevea", href="http://faculty.ucmerced.edu/jvevea/", target="_blank"), "(model and application)"),p(a("Kathleen Coburn", href="http://www.katiecoburn.weebly.com", target="_blank"),"(application)"),br(), strong('References:'), p(a("Hedges, L. V. & Vevea, J. L. (1996). Estimating effect size under publication bias: Small sample properties and robustness of a random effects selection model. Journal of Educational and Behavioral Statistics, 21, 299-333.", href="http://bayes-1.ucmerced.edu/jvevea/pdfs/Hedges%20&%20Vevea%201996%20JEBS.pdf", target="_blank")), p(a("Vevea, J. L. & Hedges, L. V. (1995). A general linear model for estimating effect size in the presence of publication bias. Psychometrika, 60(3), 419-435.", href="http://bayes-1.ucmerced.edu/jvevea/pdfs/Vevea%20&%20Hedges%201995.pdf", target="_blank"                                                                                                                                                                                                                                                                                                                                                          )),p(a("Vevea, J. L. & Woods, C. M. (2005). Publication bias in research synthesis: Sensitivity analysis using a priori weight functions. Psychological Methods, 10, 428.", href="http://bayes-1.ucmerced.edu/jvevea/pdfs/vevea%20and%20woods.pdf", target="_blank"))),
                  tabPanel("Data File",icon=icon("list", lib="font-awesome"), tableOutput("contents")),
                  #tabPanel("Sampling Variance Computation", radioButtons("ES", label="What is your effect size metric?", choices=list("Odds Ratio" = 1, "Correlation" = 2, "Standardized Mean Difference" = 3, "Risk Ratio" = 4), selected=2),p(),c("Please enter the column number of the variable containing your effect sizes:"),p(), textInput("effectsb", ' ', value="0"),tableOutput("effectsc")),
                  tabPanel("Funnel Plot",icon=icon("desktop", lib="font-awesome"),checkboxInput('flip','Plot effect sizes on x-axis', FALSE),p(),checkboxInput('estimates','Show unadjusted mean estimate (in red)',FALSE),p(),checkboxInput('estimates2','Show adjusted mean estimate (in blue)',FALSE),p(),checkboxInput('contour','Add contour lines to funnel plot at p-value cutpoints', FALSE), p(), sliderInput(inputId = "height", label = "Plot Height (px):", min = 0, max = 400, step = 1, value = 400),p(), sliderInput(inputId = "width", label = "Plot Width (%):", min = 0, max = 100, step = 1, value = 100),p(), uiOutput("funnelplot2"),downloadButton('downloadfunnelplot','Download the plot as a .pdf')),
                  #tabPanel("Unadjusted Random-Effects Model", tableOutput("unadjustedweightfunct")),
                 # tabPanel("Adjusted Random-Effects Model", textOutput("errormessage"),tableOutput("adjustedweightfunct"), strong("Likelihood ratio test comparing this model to its unadjusted version:"), p(),tableOutput("likelihoodratio")),
                  tabPanel("Model Results", icon=icon("desktop", lib="font-awesome"), strong("Unadjusted Model:"), p(),tableOutput("unadjustedweightfunction"),p(),strong("Adjusted Model:"),imageOutput("questionmark",width=17,height=17,inline=TRUE), bsTooltip("questionmark", "Remember: All weights are interpreted relative to the first weight, which is fixed to 1. A weight of 0.50 indicates that p-values in that interval are half as likely to survive selection as those in the first interval.", placement="right"), p(),tableOutput("adjustedweightfunction"),p(), strong("Likelihood Ratio Test:"),imageOutput("questionmark2",width=17,height=17,inline=TRUE),bsTooltip("questionmark2", "This test compares the adjusted and unadjusted models.", placement="right"),p(),tableOutput("likelihoodratio"),p(),strong("Effect Sizes Per Interval:"),p(),tableOutput("samplesizes"),p(),tableOutput("numberofeffects"))


  )
))))
