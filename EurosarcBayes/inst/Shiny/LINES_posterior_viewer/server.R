library(shiny)
x=0:1000/1000
cex.main=1.5
cex.axis=1.5

# Define server logic required to draw a histogram
shinyServer(function(input, output, session){
  ############################################################################################
  # reset app on close
  session$onSessionEnded(function() {


    stopApp()
  })

  prior.data=reactive({

    r0=dbeta(x,input$rpriora,input$rpriorb)
    t0=dbeta(x,input$tpriora,input$tpriorb)

    return(list(r0=r0,t0=t0))
  })

  post.data=reactive({

    rdata=dbeta(x,input$rdataa,input$rdatab)
    tdata=dbeta(x,input$tdataa,input$tdatab)

    rposteriora=input$rpriora+input$rdataa
    rposteriorb=input$rpriorb+input$rdatab
    tposteriora=input$tpriora+input$tdataa
    tposteriorb=input$tpriorb+input$tdatab

    r1=dbeta(x,rposteriora,rposteriorb)
    t1=dbeta(x,tposteriora,tposteriorb)

    return(list(rdata=rdata,tdata=tdata,r1=r1,t1=t1,rposteriora=rposteriora,rposteriorb=rposteriorb,tposteriora=tposteriora,tposteriorb=tposteriorb))
  })



  output$prior.graph <- renderPlot({
    par(mfrow=c(1,2))

    # set graph limits
    prior.r=max(prior.data()$r0)
    data.r=max(post.data()$rdata)
    prior.t=max(prior.data()$t0)
    data.t=max(post.data()$tdata)
    maxr=max(3,ifelse(is.na(prior.r)==F & prior.r<10^6,prior.r,3),ifelse(is.na(data.r)==F & data.r<10^6,data.r,3))
    maxt=max(3,ifelse(is.na(prior.t)==F & prior.t<10^6,prior.t,3),ifelse(is.na(data.t)==F & data.t<10^6,data.t,3))
    ylimr=c(0,maxr+0.6)
    ylimt=c(0,maxt+0.6)
    #####################################################################################################

    plot(x,prior.data()$r0,type="l",col=adjustcolor("blue",alpha.f=0.3),xaxs="i",yaxs="i",ylim=ylimr,
         xlab="Probability of response", ylab="Density",cex.lab=cex.axis)
    lines(x,post.data()$rdata,col="blue")
    title(main="Response",cex.main=cex.main)

    plot(x,prior.data()$t0,type="l",col=adjustcolor("red",alpha.f=0.3),xaxs="i",yaxs="i",ylim=ylimt,
         xlab="Probability of toxicity", ylab="Density",cex.lab=cex.axis)
    lines(x,post.data()$tdata,col="red")
    title(main="Toxicity",cex.main=cex.main)
  },width = "auto", height = "auto")

  output$posterior.graph <- renderPlot({

    par(mfrow=c(1,2))

    plot(x,post.data()$r1,type="l",xaxs="i",yaxs="i",col=1,ylim=c(0,ceiling(max(post.data()$r1)+0.6)),
         xlab="Probability of response", ylab="Density",cex.lab=cex.axis)
    title(main="Response",cex.main=cex.main)
    abline(v=c(input$rlower,input$rupper),h=0)

    if(input$resp.endpoint=="Futility"){
      polygon(c(0,x[1+0:(input$rupper*1000)],x[1+(input$rupper*1000)]),c(0,post.data()$r1[1+0:(input$rupper*1000)],0),col=adjustcolor("red",alpha.f=0.5),border="red")
      text(input$rupper+0.2,max(post.data()$r1)+0.3,paste0("P(R<",input$rupper,") = ",round(pbeta(input$rupper,post.data()$rposteriora,post.data()$rposteriorb),3)),cex=2,col="red")
    } else if(input$resp.endpoint=="Efficacy") {
      polygon(c(x[1+(input$rlower*1000)],x[1+(input$rlower*1000):1000],1),c(0,post.data()$r1[1+(input$rlower*1000):1000],0),col=adjustcolor("blue",alpha.f=0.5),border="blue")
      text(input$rupper+0.2,max(post.data()$r1)+0.3,paste0("P(R>",input$rlower,") = ",round(1-pbeta(input$rlower,post.data()$rposteriora,post.data()$rposteriorb),3)),cex=2,col="blue")
    }

    plot(x,post.data()$t1,type="l",xaxs="i",yaxs="i",col=1,ylim=c(0,ceiling(max(post.data()$t1)+0.6)),
         xlab="Probability of toxicity", ylab="Density",cex.lab=cex.axis)
    title(main="Toxicity",cex.main=cex.main)
    abline(v=c(input$tlower,input$tupper),h=0)

    if(input$tox.endpoint=="Toxicity"){
      polygon(c(x[1+(input$tlower*1000)],x[1+(input$tlower*1000):1000],1),c(0,post.data()$t1[1+(input$tlower*1000):1000],0),col=adjustcolor("red",alpha.f=0.5),border="red")
      text(input$tupper+0.2,max(post.data()$t1)+0.3,paste0("P(T>",input$tlower,") = ",round(1-pbeta(input$tlower,post.data()$tposteriora,post.data()$tposteriorb),3)),cex=2,col="red")
    } else if(input$tox.endpoint=="No Toxicity")  {
      polygon(c(0,x[1+0:(input$tupper*1000)],x[1+(input$tupper*1000)]),c(0,post.data()$t1[1+0:(input$tupper*1000)],0),col=adjustcolor("blue",alpha.f=0.5),border="blue")
      text(input$tupper+0.2,max(post.data()$t1)+0.3,paste0("P(T<",input$tupper,") = ",round(pbeta(input$tupper,post.data()$tposteriora,post.data()$tposteriorb),3)),cex=2,col="blue")
    }

  },width = "auto", height = "auto")
})

"\u2119"
