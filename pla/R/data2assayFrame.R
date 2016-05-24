data2assayFrame <-
function (dataframe,
          dr = 1.5,
          Z = log((1/dr)^(max(Dilution) - Dilution)),
          design = "lsd") 
{
    OK <- .checkPlaFrame(dataframe, design = design)
    if (OK) {
        Dilution   <- as.numeric(unlist(dataframe["Dilution"]))
        names(Z)   <- "Z"
        Sample     <- as.numeric(unlist(dataframe["Sample"]))
        Sample     <- 1 + max(Sample) - Sample
        Step       <- as.numeric(unlist(dataframe["Dilution"]))
        Step       <- 1 + max(Step) - Step
        SampleStep <- paste0(Sample, ":", Step)
        if (design == "lsd") {
            Replicate <- unlist(dataframe["Row"])
            names(Replicate) <- "Replicate"
            dataframe <- cbind(dataframe, Replicate)
        } else
            Replicate <- unlist(dataframe["Replicate"])
        Order      <- order(Sample, Step, Replicate)
     ## Response   <- unlist(dataframe["Response"])
        J          <- 1:length(Replicate)
        data       <- cbind(dataframe, SampleStep, Z, J)
        data       <- data[Order, ]
        I          <- 1:length(Replicate)
        data       <- cbind(data, I)
        return(data)
    }
}
