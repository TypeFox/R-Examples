ModelSelectionVisualization <- function(outdir){
                                        # Loading the model performance tables
betamultiplier <- NULL
AICc <- NULL
AUC.Test <- NULL
variables <- NULL

    ModelPerformances <- read.table(paste(outdir,"/ModelPerformance.txt",sep=""), header=TRUE, stringsAsFactors=FALSE)
    ModelPerformances$AIC[ModelPerformances$AIC=="x"]=NA
    ModelPerformances$AICc[ModelPerformances$AICc=="x"]=NA
    ModelPerformances$AICc <- as.numeric(ModelPerformances$AICc)
    ModelPerformances$BIC[ModelPerformances$BIC=="x"]=NA
    ModelPerformances$AUC.Test <- as.numeric(ModelPerformances$AUC.Test)
    
    
                                        # Select the lowest AICc values 
    ModelPerformances.bestsets <- ModelPerformances[
                                                    which(ModelPerformances$AICc==min(ModelPerformances$AICc,na.rm=TRUE))
                                                   ,]
    
    ModelPerformances <- ModelPerformances[
                                           -which(ModelPerformances$AICc==min(ModelPerformances$AICc,na.rm=TRUE))
                                          ,]
                                    
    
    ModelPerformances.added <- ModelPerformances.bestsets # Just to add this in the end so that the legends don't appear red
    ModelPerformances.added$AUC.Test <- ModelPerformances.added$AUC.Test+5
    
    png(filename = paste(outdir,"/ModelSelectionAICc_MarkedMinAICc.png",sep=""), height=120, width=250, units="mm", pointsize=12, res=600)
    
    print(
        ggplot2::ggplot(ModelPerformances,ggplot2::aes(betamultiplier,AICc,size=log(variables),colour=variables))+
            ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red")+
                ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red",shape=1,size=8)+
                    ggplot2::geom_point()+
                        ggplot2::ylim(min(ModelPerformances$AICc)-100,max(ModelPerformances$AICc)+100)+
                            ggplot2::guides(size=FALSE)+ # Exclude legend on size
                                ggplot2::theme(axis.text.y=ggplot2::element_text(size=14),
                                      axis.text.x = ggplot2::element_text(size=14),
                                      axis.title.x=ggplot2::element_text(size=15),
                                      axis.title.y=ggplot2::element_text(size=15,vjust=1),
                                      legend.text=ggplot2::element_text(size=13))+
                                          ggplot2::theme(legend.title = ggplot2::element_text(size = 13))+
                                                  ggplot2::labs(colour="Variables")
        
    )
    dev.off()
    
    
    png(filename = paste(outdir,"/ModelSelectionAUCTest_MarkedMinAICc.png",sep=""), height=120, width=250, units="mm", pointsize=12, res=600)
    
    print(
        ggplot2::ggplot(ModelPerformances,ggplot2::aes(betamultiplier,AUC.Test,size=log(variables),colour=variables))+
            ggplot2::geom_point()+
                ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red")+
                    ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red",shape=1,size=8)+
                        ggplot2::geom_point(data=ModelPerformances.added)+    
                            ggplot2::ylim(min(ModelPerformances$AUC.Test)-0.1,max(ModelPerformances$AUC.Test)+0.1)+
                                ggplot2::guides(size=FALSE)+ # Exclude legend on size
                                    ggplot2::theme(axis.text.y=ggplot2::element_text(size=14),
                                          axis.text.x = ggplot2::element_text(size=14),
                                          axis.title.x= ggplot2::element_text(size=15),
                                          axis.title.y= ggplot2::element_text(size=15,vjust=1),
                                          legend.text=ggplot2::element_text(size=13))+
                                              ggplot2::theme(legend.title = ggplot2::element_text(size = 13))+
                                                      ggplot2::labs(colour="Variables")
        
    )
    dev.off()
    
    ModelPerformances <- read.table(paste(outdir,"/ModelPerformance.txt",sep=""), header=TRUE, stringsAsFactors=FALSE)
    ModelPerformances$AIC[ModelPerformances$AIC=="x"]=NA
    ModelPerformances$AICc[ModelPerformances$AICc=="x"]=NA
    ModelPerformances$AICc <- as.numeric(ModelPerformances$AICc)
    ModelPerformances$BIC[ModelPerformances$BIC=="x"]=NA
    ModelPerformances$AUC.Test <- as.numeric(ModelPerformances$AUC.Test)
    

             # Select the highest AUC.Test values
    ModelPerformances.bestsets <- ModelPerformances[
                                                    which(ModelPerformances$AUC.Test==max(ModelPerformances$AUC.Test,na.rm=TRUE))
                                                   ,]
    
    ModelPerformances <- ModelPerformances[
                                           -which(ModelPerformances$AUC.Test==max(ModelPerformances$AUC.Test,na.rm=TRUE))
                                          ,]
                                    
    
    ModelPerformances.added <- ModelPerformances.bestsets # Just to add this in the end so that the legends don't appear red
    ModelPerformances.added$AUC.Test <- ModelPerformances.added$AUC.Test+5
    
    png(filename = paste(outdir,"/ModelSelectionAICc_MarkedMaxAUCTest.png",sep=""), height=120, width=250, units="mm", pointsize=12, res=600)
    
    print(
        ggplot2::ggplot(ModelPerformances,ggplot2::aes(betamultiplier,AICc,size=log(variables),colour=variables))+
            ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red")+
                ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red",shape=1,size=8)+
                    ggplot2::geom_point()+
                        ggplot2::ylim(min(ModelPerformances$AICc)-100,max(ModelPerformances$AICc)+100)+ 
                            ggplot2::guides(size=FALSE)+ # Exclude legend on size
                                ggplot2::theme(axis.text.y=ggplot2::element_text(size=14),
                                      axis.text.x = ggplot2::element_text(size=14),
                                      axis.title.x=ggplot2::element_text(size=15),
                                      axis.title.y=ggplot2::element_text(size=15,vjust=1),
                                      legend.text=ggplot2::element_text(size=13))+
                                          ggplot2::theme(legend.title = ggplot2::element_text(size = 13))+
                                                  ggplot2::labs(colour="Variables")
        
    )
    dev.off()
    
    
    png(filename = paste(outdir,"/ModelSelectionAUCTest_MarkedMaxAUCTest.png",sep=""), height=120, width=250, units="mm", pointsize=12, res=600)
    
    print(
        ggplot2::ggplot(ModelPerformances,ggplot2::aes(betamultiplier,AUC.Test,size=log(variables),colour=variables))+
            ggplot2::geom_point()+
                ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red")+
                    ggplot2::geom_point(data=ModelPerformances.bestsets,colour="red",shape=1,size=8)+
                        ggplot2::geom_point(data=ModelPerformances.added)+    
                            ggplot2::ylim(min(ModelPerformances$AUC.Test)-0.1,max(ModelPerformances$AUC.Test)+0.1)+
                                ggplot2::guides(size=FALSE)+ # Exclude legend on size
                                    ggplot2::theme(axis.text.y=ggplot2::element_text(size=14),
                                          axis.text.x = ggplot2::element_text(size=14),
                                          axis.title.x=ggplot2::element_text(size=15),
                                          axis.title.y=ggplot2::element_text(size=15,vjust=1),
                                          legend.text=ggplot2::element_text(size=13))+
                                              ggplot2::theme(legend.title = ggplot2::element_text(size = 13))+
                                                      ggplot2::labs(colour="Variables")
        
    )
    dev.off()
    
}
