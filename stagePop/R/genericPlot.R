#' Generic Plotting Function
#'
#' Plot the state variables over time (sums over all strains if there are multiple strains in a species) 
#'
#' @param output Model output, from \code{\link{popModel}}
#' @param numSpecies Number of species in the model
#' @param numStages Number of life stages
#' @param speciesNames species names (vector of strings). Default is NULL.
#' @param stageNames stage names (a list of vectors if there is more than one species).
#' @param saveFig Choose to save the figure (TRUE or FALSE). Default is FALSE. 
#' @param figType Figure format can be 'eps', 'tiff' or 'png'. Default is 'eps' 
#' @param figName filepath to save figure to. Default is 'stagePopFig'
#'
#' @seealso \code{\link{popModel}}
#' @export

genericPlot=function(output,numSpecies,numStages,speciesNames,stageNames,saveFig=FALSE,figType='eps',figName='stagePopFig'){

# print('enter genericPlot')
  if (!is.list(stageNames)){stageNames=list(stageNames)}
 
  #set size of labels
  cexlabsize=2.0
  cexaxissize=2.0
  cexlegendsize=1.75
  cexmainsize=1.75
  
  if (numSpecies>1 & max(numStages)>1){
      wlen=5*numSpecies;hlen=5
      cexlegendsize=1.25
  }else{
      wlen=7;hlen=7}

  dev.new(bg="white",horizontal=FALSE,onefile = FALSE, paper = 'special',width=wlen,height=hlen)

  time=output[,1]
  ct=2

  par(mar=c(5,5,2,2))
  
  if (max(numStages)==1){  #no stage structure in any species
    cols=rainbow(numSpecies)
    data=output[,ct:(ct+numSpecies-1)]
    ct=numSpecies+ct
    plot(c(min(time),max(time)),c(min(data),1.1*max(data)),type='n',xlab='time',ylab='density',cex.lab=cexlabsize,cex.axis=cexaxissize,cex.main=cexmainsize)
    if (numSpecies>1){
      for (j in seq(1,numSpecies)){
        lines(time,data[,j],col=cols[j],lwd=2)
        legend('topright',speciesNames[1:numSpecies],lty=1,col=cols,cex=cexlegendsize,lwd=2,bty='n')}
    }else{lines(time,data,col=cols)}
    

  }else{

    if (numSpecies==1){     #only one species but it has stage structure
      cols=rainbow(numStages)
      data=output[,ct:(ct+numStages-1)]
      ct=numStages+ct
      plot(c(min(time),max(time)),c(min(data),1.1*max(data)),type='n',xlab='time',ylab='density',cex.lab=cexlabsize,cex.axis=cexaxissize,cex.main=cexmainsize)
      if (numStages>1){
        for (j in seq(1,numStages)){
          lines(time,data[,j],col=cols[j],lwd=2)
            legend('topright',stageNames[[1]],lty=1,col=cols,cex=cexlegendsize,lwd=2,bty='n')
        }
      }else{lines(time,data,col=cols)}
      
    }else{ #multiple species with stage structure
    
      par(mfrow=c(1,numSpecies))
      ct=2
      for (i in seq(1,numSpecies)){
        cols=rainbow(numStages[i])
        data=output[,ct:(ct+numStages[i]-1)]
        if (is.finite(max(data))){
          titleText=speciesNames[i]
          plot(c(min(time),max(time)),c(min(data),1.1*max(data)),type='n',main=titleText,xlab='time',ylab='density',cex.lab=cexlabsize,cex.axis=cexaxissize,cex.main=cexmainsize)
          if (numStages[i]>1){
            for (j in seq(1,numStages[i])){
              lines(time,data[,j],col=cols[j],lwd=2)}
            legend('topright',stageNames[[i]],lty=1,col=cols,cex=cexlegendsize,lwd=2,bty='n')
          }else{lines(time,data,col=cols,lwd=2)}
          ct=numStages[i]+ct
        }
      }
    }
  }

  if (saveFig){
    if (figType=='pdf'){dev.copy2pdf(file=paste(figName,'.pdf',sep=""))}
    if (figType=='eps'){dev.copy2eps(file=paste(figName,'.eps',sep=""))}
    if (figType=='png'){dev.print(png,filename=paste(figName,'.png',sep=""),res=100,width=wlen,height=hlen,units='in')}
    if (figType=='tiff'){dev.print(tiff,filename=paste(figName,'.tiff',sep=""),res=100,width=wlen,height=hlen,units='in')}
    #dev.off()
  }

  
}

