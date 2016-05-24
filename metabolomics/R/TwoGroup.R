TwoGroup <- function(inputdata, alternative="two.sided", paired=FALSE, 
    padjmethod="BH", saveoutput=FALSE, outputname="results", ...)
{
    # Perform t-Tests
    ttests <- TTests(inputdata=inputdata, alternative=alternative, 
        paired=paired, padjmethod=padjmethod, saveoutput=FALSE, ...
    )
    # Perform fold change comparison
    fc <- FoldChange(inputdata=inputdata, paired=paired, saveoutput=FALSE, 
        plot.hist=FALSE, outputname="fc.results"
    )
    # Generate output
    results <- cbind(ttests, fc[, 2],  fc[, 2] / ttests[, 1])
    colnames(results) <- c("t-statistic", "p-value", 
        paste(padjmethod, "Adjusted p-value"), "fold change", "standard error"
    )
    
    if (saveoutput)
        write.csv(results, paste(c(outputname, ".csv"), collapse=""))
    output<-c()
    output$output<-results
    
    return(structure(output, class="results"))
}
