numberInput<-function (inputId,label, value = "",...){
  div(class="form-group",style="display:block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "number", value = value,class="input-small",...))
}


shinyServer(function(input, output, session){

  ############################################################################################
  # reset app on close
  session$onSessionEnded(function() {
    stopApp()
  })

  simons=reactiveValues()
  custom=reactiveValues()


  # update the simons two stage designs
  update.simon=observe({

    simons$optimal=freq_binom_one_simons_twostage(input$p0,input$p1,input$alpha,input$power,input$prior.a,input$prior.b,round=TRUE,method="optimal")
    simons$minimax=freq_binom_one_simons_twostage(input$p0,input$p1,input$alpha,input$power,input$prior.a,input$prior.b,round=TRUE,method="minmax")

  })

  # Update the custom two stage design
  custom.design=observe({
    custom$design=properties_binom_one(failure=c(input$s1,input$r-1),success=c(input$r1,input$r),reviews=c(input$n1,input$n),input$p0,input$p1,input$prior.a,input$prior.b,round=TRUE)
  })

  output$optimal=renderText({
    paste(h4("Maximum number of successes to stop trial at interim for futility:",strong(simons$optimal@failure[1]),"with",strong(simons$optimal@reviews[1]),"patients"),
          h4("Minimum number of successes to declare successful trial at final analysis:", strong(simons$optimal@success[2]), "with",strong(simons$optimal@reviews[2]),"patients"),br(),
          strong("Frequentist properties"),br(),
          "Expected sample size under H<sub>0</sub>=",strong(simons$optimal@exp.p0),br(),br(),
          "Type 1 error (&alpha;) = <strong>P</strong>(data&ge;r|p=p<sub>0</sub>, n) =",strong(simons$optimal@alpha),br(),
          "Power (&beta;) = <strong>P</strong>(data&ge;r|p=p<sub>1</sub>, n) =",strong(simons$optimal@power),br(),br(),
          strong("Bayesian properties"),br(),
          "&zeta; = <strong>P</strong>(p&lt;p<sub>1</sub>|data,prior) =",strong(min(simons$optimal@zeta,na.rm =TRUE)),br(),
          "&eta; = <strong>P</strong>(p&gt;p<sub>0</sub>|data,prior) =",strong(min(simons$optimal@eta,na.rm =TRUE)),br())
  })

  output$minimax=renderText({
    paste(
      h4("Maximum number of successes to stop trial at interim for futility:",strong(simons$minimax@failure[1]),"with",strong(simons$minimax@reviews[1]),"patients"),
      h4("Number of successes to declare successful trial at final analysis:", strong(simons$minimax@success[2]), "with",strong(simons$minimax@reviews[2]),"patients"),br(),
      strong("Frequentist properties"),br(),
      "Expected sample size under H<sub>0</sub>=",strong(simons$minimax@exp.p0),br(),br(),
      "Type 1 error (&alpha;) = <strong>P</strong>(data&ge;r|p=p<sub>0</sub>, n) =",strong(simons$minimax@alpha),br(),
      "Power (&beta;) = <strong>P</strong>(data&ge;r|p=p<sub>1</sub>, n) =",strong(simons$minimax@power),br(),br(),
      strong("Bayesian properties"),br(),
      "&zeta; = <strong>P</strong>(p&lt;p<sub>1</sub>|data,prior) =",strong(min(simons$minimax@zeta[2],na.rm =TRUE)),br(),
      "&eta; = <strong>P</strong>(p&gt;p<sub>0</sub>|data,prior) =",strong(min(simons$minimax@eta[2],na.rm =TRUE)),br())
  })

  output$custom=renderText({
    paste(


    h4("Number of patients required at interim analysis",strong(custom$design@reviews[1])),
    h4("Maximum number of successes to stop trial at interim for futility:",strong(ifelse(custom$design@failure[1]<0,"Not stopping for futility",custom$design@failure[1]))),
    h4("Minimum number of successes to stop trial at interim for efficacy:",strong(ifelse(custom$design@success[1]>custom$design@reviews[1],"Not stopping for efficacy",custom$design@success[1]))),
    h4("Number of patients required at final analysis",strong(custom$design@reviews[2])),
    h4("Number of successes to declare successful trial at final analysis:", strong(custom$design@success[2])),br(),
    strong("Frequentist properties"),br(),
    "Expected sample size under H<sub>0</sub>=",strong(custom$design@exp.p0),br(),
    "Expected sample size under H<sub>1</sub>=",strong(custom$design@exp.p1),br(),br(),
    "Type 1 error (&alpha;) = <strong>P</strong>(data&ge;r|p=p<sub>0</sub>, n) =",strong(custom$design@alpha),br(),
    "Power (&beta;) = <strong>P</strong>(data&ge;r|p=p<sub>1</sub>, n) =",strong(custom$design@power),br(),br(),
    strong("Bayesian properties"),br(),
    "&zeta; = <strong>P</strong>(p&lt;p<sub>1</sub>|data,prior) =",strong(min(custom$design@zeta,na.rm =TRUE)),br(),
    "&eta; = <strong>P</strong>(p&gt;p<sub>0</sub>|data,prior) =",strong(min(custom$design@eta,na.rm =TRUE)),br())
  })

  output$bayes.prior.ss=renderText({
    paste("Effective sample size of prior n<sub>prior</sub>", strong(input$prior.a+input$prior.b))
  })



})
