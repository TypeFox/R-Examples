ModelWithBestPerformance <- function(outdir){

    # This function extracts the best model based on AUC.Test values and based on AIC.c values

    ModelPerformances <- read.table(paste(outdir,"/ModelPerformance.txt",sep=""),header=TRUE,stringsAsFactors=FALSE)
    VariableSelection <- read.table(paste(outdir,"/VariableSelectionProcess.txt",sep=""),header=FALSE,stringsAsFactors=FALSE)
    row.names(VariableSelection) <- VariableSelection[,1]
    VariableSelection <- VariableSelection[,-1]


    p <- as.numeric(as.vector(ModelPerformances$AUC.Test))
                                        # Select the model of highest AUC.Test
    pbestAUC.Test <-  ModelPerformances[which(ModelPerformances$AUC.Test==max(p)),]
    
    
    p <- as.numeric(as.vector(ModelPerformances$AICc))
                                        # Select the model of highest AUC.Test
    p <- as.numeric(as.vector(p[p!="x"]))
    pbestAICc <-  ModelPerformances[which(ModelPerformances$AICc==min(p)),]
    

    VariableSelectionSteps.AUCTest <- VariableSelection[,as.numeric(VariableSelection[3,])==as.numeric(pbestAUC.Test$betamultiplier)]
    VariableSelectionSteps.AICc <- VariableSelection[,as.numeric(VariableSelection[3,])==as.numeric(pbestAICc$betamultiplier)]

    write.table(pbestAUC.Test, file = paste(outdir,"/ModelWithMaxAUCTest.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE)
    write.table(pbestAICc, file = paste(outdir,"/ModelWithMinAICc.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = TRUE)

    write.table(VariableSelectionSteps.AUCTest, file = paste(outdir,"/VariableSelectionMaxAUCTest.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = FALSE)
    write.table(VariableSelectionSteps.AICc, file = paste(outdir,"/VariableSelectionMinAICc.txt",sep=""), append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                col.names = FALSE)

    cat("\n","==========================================================================","\n")
    cat(" Model with the maximum AUC.Test value specified in:","\n")
    cat(" ModelWithMaxAUCTest.txt","\n")
    cat("\n","Relevant variable selection steps associated with this model are saved in:","\n")
    cat(" VariableSelectionMaxAUCTest.txt","\n")
    cat("\n","==========================================================================","\n")
    cat("\n","\n")
    cat("\n","==========================================================================","\n")
    cat(" Model with the minimum AICc value specified in:","\n")
    cat(" ModelWithMinAICc.txt","\n")
    cat("\n","Relevant variable selection steps associated with this model are saved in:","\n")
    cat(" VariableSelectionMinAICc.txt","\n")
    cat("\n","==========================================================================","\n")
    
}
