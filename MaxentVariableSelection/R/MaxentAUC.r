MaxentAUC <-
function(outdir){
    Resulttable <- read.csv(paste(outdir,"/maxentResults.csv",sep=""),sep=",",header=TRUE)
    # Select those columns where the column name contains the word "Test AUC"
    TestAUCs <- Resulttable[,(colnames(Resulttable) %in% grep("Test.AUC",colnames(Resulttable),value=T))]
    # The last one of these values is the average value
    TestAUC.average <- TestAUCs[length(TestAUCs)]
        # Select those columns where the column name contains the word "Training AUC"
    TrainAUCs <- Resulttable[,(colnames(Resulttable) %in% grep("Training.AUC",colnames(Resulttable),value=T))]
        # The last one of these values is the average value
    TrainAUC.average <- TrainAUCs[length(TrainAUCs)]
    return(c(TestAUC.average,TrainAUC.average))
    }
