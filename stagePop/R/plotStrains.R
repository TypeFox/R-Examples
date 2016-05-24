plotStrains=function(out,numSpecies,numStages,numStrains,speciesNames,stageNames,figName='strainFig',figType='eps',saveFig=FALSE){
  
  #set size of labels
  cexlabsize=2.0
  cexaxissize=2.0
  cexlegendsize=1.0
  cexmainsize=1.75
  
  if (!is.list(stageNames)){stageNames=list(stageNames)}
  
  for (species in seq(1,numSpecies)){

    if (numStrains[species]>1){
      cols=rainbow(numStrains[species])
      
      if (numStages[species]>1){wlen=5*numStages[species];hlen=5}else{wlen=7;hlen=7}
      dev.new(bg="white",horizontal=FALSE,onefile = FALSE, paper = 'special',width=wlen,height=hlen)
      par(mfrow=c(1,numStages[species]))

      for (stage in seq(1,numStages[species])){
        str1=paste(speciesNames[species],'.',stageNames[[species]][stage],sep='')
        titlestr=paste(speciesNames[species],' (',stageNames[[species]][stage],')',sep='')
        plot(c(min(out[,'time']),1.1*max(out[,'time'])),c(min(0,min(out[,grepl(paste('^',str1,sep=''),colnames(out))])),max(out[,grepl(paste('^',str1,sep=''),colnames(out))])),type='n',main=titlestr,xlab='time',ylab='density',cex.lab=cexlabsize,cex.axis=cexaxissize,cex.main=cexmainsize)
#        plot(c(min(out[,'time']),1.1*max(out[,'time'])),c(min(0,min(out[,grepl(paste('^',str1,sep=''),colnames(out))])),max(out[,grepl(paste('^',str1,sep=''),colnames(out))])),type='n',main=titlestr,xlab='time',ylab='density',cex.lab=2,cex.axis=2,cex.main=1.75)

        for (i in seq(1,numStrains[species])){
          str=paste(str1,'.strain',i,sep='')
          lines(out[,'time'],out[,str],col=cols[i],lwd=2)
        }
        legend('topleft',paste('strain',seq(1,numStrains[species])),col=cols,lty=1,lwd=2,bty='n',cex=cexlegendsize)
      }
      if (saveFig){
        if (numSpecies>1){figName=paste(figName,speciesNames[species],sep=''); print(paste('Fig saved to',figName))}
        if (figType=='pdf'){dev.copy2pdf(file=paste(figName,'.pdf',sep=""))}
        if (figType=='eps'){dev.copy2eps(file=paste(figName,'.eps',sep=""))}
        if (figType=='png'){dev.print(png,filename=paste(figName,'.png',sep=""),res=100,width=wlen,height=hlen,units='in')}
        if (figType=='tiff'){dev.print(tiff,filename=paste(figName,'.tiff',sep=""),res=100,width=wlen,height=hlen,units='in')}
      }
    }
  }
}
