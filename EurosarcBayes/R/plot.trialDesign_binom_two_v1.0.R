

.plot.trialDesign_binom_two=function(x , R_min=NULL, R_max=NULL, T_min=NULL, T_max=NULL, ...){

  BayesTwoEndpoint=x
  decision=BayesTwoEndpoint@decision
  reviews=BayesTwoEndpoint@reviews

  # Empty plot for lines graph
  par(mfrow=(c(2,ceiling((length(reviews)+1)/2))),mai=0.17*c(5, 5.4, 4, 1))

  if(length(names(BayesTwoEndpoint@graph))>0){
    par(mfrow=(c(2,ceiling((length(reviews)+1)/2))),mai=0.17*c(5, 5.4, 4, 1))
    BayesTwoEndpoint@graph$fun.graph(BayesTwoEndpoint@graph)
  } else {
    mar=par()$mar
    par(mar=c(0,0,0,0))
    plot(0,0,xlim=c(0,1),ylim=c(0,1),col=0,axes=FALSE,xlab="",ylab="")

    polygon(c(0.05,0.05,0.25,0.25),c(0.7,0.9,0.9,0.7),col="red")
    polygon(c(0.05,0.05,0.25,0.25),c(0.4,0.6,0.6,0.4),col="orange")
    polygon(c(0.05,0.05,0.25,0.25),c(0.1,0.3,0.3,0.1),col="seagreen4")
    text(0.3,0.8,"Stop research for futility/safety",adj = 0,cex=1.5)
    text(0.3,0.5,"Continue recruiting",adj = 0,cex=1.5)
    text(0.3,0.2,"Recommend further research",adj = 0,cex=1.5)
    par(mar=mar)
  }
  for(i in 1:length(reviews)){

    image(decision[[i]],zlim=c(1,3),xlab="Number of responses",ylab="Number of toxic events",col=c("red","orange","seagreen4"),main=paste("Analysis at",reviews[i],"patients"),xaxs = "i", yaxs = "i", cex.lab = 1.5, cex.axis = 1.5,cex.main=1.5,axes=FALSE)

    axis(1,at=(0:reviews[i])/reviews[i],labels=0:reviews[i],cex.axis=1.5)
    axis(2,at=(0:reviews[i])/reviews[i],labels=0:reviews[i],cex.axis=1.5)
    box()
    abline(v=c(R_min,R_max))
    abline(h=c(T_min,T_max))
  }

}

setMethod(f="plot",signature=c(x="trialDesign_binom_two",y=NULL),definition=.plot.trialDesign_binom_two)

# ended
