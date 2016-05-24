Summarize <- 
function(GraphQuant){
## THIS FUNCTION IS USED TO PRINT THE SUMMARY OF THE METHOD 
## AND PLACE IN THE OPENING TAB FOR ANY GIVEN METHOD.
Summary <- "\n";

## ADD MODEL NAME
Summary <- paste(Summary,  "Model: ", GraphQuant$Name,"\n");

## NUMBER OF SAMPLES
Summary <- paste(Summary, "Number of Samples:",  as.character(nrow(GraphQuant$Coords)), "\n \n");
## INFORMATION FROM TRUE CLUSTER MEMBERSHIP
if (!is.null(GraphQuant$TrueMem)){
        ##NUMBER OF TRUE CLUSTERS
        Summary <- paste(Summary, "True Number of Clusters: ",  as.character(length(table(GraphQuant$TrueMem))), "\n \n");

        ##COUNTS IN EACH CLUSTER
        TrueComp <- as.matrix(table(GraphQuant$TrueMem));
        Summary <- paste(Summary, "Composition of true clusters \n \n");
        Summary <- paste(Summary,  "Cluster No: \t \t");
        for (i in 1:length(table(GraphQuant$TrueMem))){
                    Summary <- paste(Summary, "\t", as.character(i));
         }
        Summary <- paste(Summary, "\n Counts: \t \t");
        for (i in 1:length(table(GraphQuant$TrueMem))){
                    Summary <- paste(Summary,  "\t",  as.character(TrueComp[i]));
         }
}

## INFORMATION FROM ESTIMATED CLUSTER MEMBERSHIP
## NUMBER OF ESTIMATED CLUSTERS
Summary <- paste(Summary, "\n \n Estimated Number of Clusters: ",  as.character(length(table(GraphQuant$ClusMem))), "\n \n");

##COUNTS IN EACH CLUSTER
        EstComp <- as.matrix(table(GraphQuant$ClusMem));
        Summary <- paste(Summary, "Composition of true clusters \n \n");
        Summary <- paste(Summary,  "Cluster No: \t \t");
        for (i in 1:length(table(GraphQuant$ClusMem))){
                    Summary <- paste(Summary, "\t", as.character(i));
         }
        Summary <- paste(Summary, " \n Counts: \t \t");
        for (i in 1:length(table(GraphQuant$ClusMem))){
                    Summary <- paste(Summary,  "\t",  as.character(EstComp[i]));
         }

## MISCLASSIFICATION ERROR
if (!is.null(GraphQuant$TrueMem)){
        MCE <- 100*MCError(GraphQuant$ClusMem, GraphQuant$TrueMem);
        Summary <- paste(Summary, "\n \n Error of Misclassification: \t \t \t ", as.character(MCE), "% \n");
}

if (is.null(GraphQuant$TrueMem)){   
        Summary <- paste(Summary, "\n \n Misclassification error not calculated as true cluster membership is not provided.");
}
return(Summary);
}