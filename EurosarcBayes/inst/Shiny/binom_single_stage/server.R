shinyServer(function(input, output, session){

  ############################################################################################
  # reset app on close
  session$onSessionEnded(function() {
    stopApp()
  })

  model.prop=reactiveValues()

  update=observe({

    eta=input$eta
    zeta=input$zeta
    p0=input$p0
    p1=input$p1
    prior.a=input$prior.a
    prior.b=input$prior.b

    model.prop$bayes=bayes_binom_one_postprob_onestage(p0,p1,eta,zeta,prior.a,prior.b,round=TRUE)
    model.prop$freq=freq_binom_one_onestage(p0,p1,1-eta,zeta,prior.a,prior.b,round=TRUE)

    cat("\nBayesian model call: bayes_binom_one_postprob_onestage(",p0,",",p1,",",eta,",",zeta,",",prior.a,",",prior.b,",round=TRUE)")
    cat("\nFrequentist model call: freq_binom_one_onestage(",p0,",",p1,",",1-eta,",",zeta,",",prior.a,",",prior.b,",round=TRUE)\n")
    })


    output$eff.sample.size=renderText({
      paste("Effective sample size of prior n<sub>prior</sub> =", strong(input$prior.a+input$prior.b))
  })

  output$sample.size=renderText({
   paste(h4("Required sample size n =",strong(model.prop$bayes@reviews)))
  })

 output$sample.success=renderText({
   paste(h4("Number of successes observed on trial to recommend treatment r =",strong(model.prop$bayes@success)))
 })

 output$sample.eta=renderText({
   paste("&eta; = <strong>P</strong>(p&gt;p<sub>0</sub>|data,prior) =", strong(model.prop$bayes@eta))
 })

 output$sample.zeta=renderText({
   HTML(paste("&zeta; = <strong>P</strong>(p&lt;p<sub>1</sub>|data,prior) =", strong(model.prop$bayes@zeta)))
 })

 output$sample.alpha=renderText({
   paste("Type 1 error (&alpha;) = <strong>P</strong>(data&ge;r|p=p<sub>0</sub>, n) =", strong(model.prop$bayes@alpha))
 })

 output$sample.power=renderText({
   paste("Power (&beta;) = <strong>P</strong>(data&ge;r|p=p<sub>1</sub>, n) =", strong(model.prop$bayes@power))
 })

 output$f.sample.size=renderText({
   paste(h4("Required sample size n =",strong(model.prop$freq@reviews)))
 })

 output$f.sample.success=renderText({
   paste(h4("Number of successes observed on trial to recommend treatment r =", strong(model.prop$freq@success)))
 })

 output$f.sample.alpha=renderText({
   paste("Type 1 error (&alpha;) = <strong>P</strong>(data&ge;r|p=p<sub>0</sub>, n) =", strong(model.prop$freq@alpha))
 })

 output$f.sample.power=renderText({
   paste("Power (&beta;) = <strong>P</strong>(data&ge;r|p=p<sub>1</sub>, n) =", strong(model.prop$freq@power))
 })

 output$f.sample.eta=renderText({
   HTML(paste("&eta; = <strong>P</strong>(p&gt;p<sub>0</sub>|data,prior) =", strong(model.prop$freq@eta)))
 })

 output$f.sample.zeta=renderText({
   HTML(paste("&zeta; = <strong>P</strong>(p&lt;p<sub>1</sub>|data,prior) =", strong(model.prop$freq@zeta)))
 })

 output$priorplot=renderPlot({
   input$prior.a
   input$prior.b
   x=0:1000/1000
   y=dbeta(x,input$prior.a,input$prior.b)
   par(cex=1.3,mar=c(5,4,2,2))
   plot(x,y,type="l",col=4,xlab="Probability",ylab="Probability density",xaxs="i",yaxs="i",ylim=c(0,range(pretty(y))[2]+0.1),xlim=c(0,1))
   polygon(c(0,x,1),c(0,y,0),col=adjustcolor("blue",alpha.f=0.25),border=4)
   abline(v=model.prop$bayes@success/model.prop$bayes@reviews,col=3)
   abline(v=model.prop$freq@success/model.prop$freq@reviews,col=2)
    legend("topright",col=c(3,2),lty=c(1,1),legend=c("Bayesian","Frequentist"),cex=1.75)
 },height=600)


})
