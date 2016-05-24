
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

simulation.plot <- function(problem=NULL,simdf=NULL,files=NULL){
  
  # Get the number of different trts
  trts <- unique(simdf$TRT)
  
  # Create storage list for plots
  plotlist <- list()
  newplot <- 1
  
  # Plotting loop
  for (i in 1:length(trts)) {
    itrt <- trts[i]
    if (i<10){
      figName <- sprintf('0%d.sim.Dose.%s.pdf',i,itrt)
    } else {
      figName <- sprintf('%d.sim.Dose.%s.pdf',i,itrt)
    }
    
    # Extract data from itrt level
    tmpdf <- simdf[which(simdf$TRT==itrt),]
    tmpdf$DVID <- factor(tmpdf$DVID,levels=unique(tmpdf$DVID))
    
    # Get the number of ids and cmts in tmpdf
    ids <- unique(tmpdf$ID)
    cmts <- unique(tmpdf$DVID)
    
    # Open device
    if (length(cmts)==1){
      pdf(file=figName,
          width=6.5,
          height=6.5)
    } else {
      pdf(file=figName,
          width=6.5,
          height=9)
    }
    
    trellis.par.set(superpose.symbol=list(col=c("blue","black"),
                                          pch=c(1,3)),
                    superpose.line=list(col=c("blue","black")))
    
    for (iid in ids){
      # Extract data from iid level
      iddf <- tmpdf[which(tmpdf$ID==iid),]
      
      if (problem$method=='subject'){
        ititle <- paste('Subject', iid, '- Treatment', i)
      } else {
        ititle <- paste('Population - Treatment', i)
      }
      
      # Create plots
      simplot <- xyplot(SIM+OBS~TIME|DVID,
                        data=iddf,
                        as.table=TRUE,
                        distribute.type=TRUE,
                        type=c('l','p'),
                        scales=list(x=list(relation='free'),
                                    y=list(relation='free')),
                        strip=strip.custom(var.name='Output',
                                           style=1,
                                           bg=0,
                                           sep=' ',
                                           strip.name=c(TRUE,TRUE)),
                        xlab='Time',
                        ylab='Observations / Predictions',
                        main=ititle,
                        layout=get.layout(length(cmts)))

    print(simplot)
    plotlist[[newplot]] <- simplot
    newplot <- newplot + 1
    }
    dev.off()
  }
  
  # Display plots in iteractive mode
  if (interactive()){
    dev.new()
    trellis.par.set(superpose.symbol=list(col=c("blue","black"),
                                          pch=c(1,3)),
                    superpose.line=list(col=c("blue","black")))
    print(plotlist[[1]])
    par(ask=TRUE)
    
    if (length(plotlist)>=2){
      for (i in 2:length(plotlist)){
        print(plotlist[[i]])
      }
    }
    par(ask=FALSE)
  }
  
  # Print a message to screen
  cat(sprintf('\n%s%s\n%s\n','Diagnostic figures have been created and saved in the ',
              'working directory.','You may open and edit them at your convenience.'))
  
}
