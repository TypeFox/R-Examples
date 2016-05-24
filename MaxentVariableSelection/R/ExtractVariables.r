ExtractVariables <- function(uncorrelated.variables,occurrencesites,backgroundsites){ 
    
                                        # Function to create new csv
                                        # files that contain only
                                        # variables that are not
                                        # correlated to the most
                                        # important variable and that
                                        # have a model contribution
                                        # above
                                        # contributionthreshold
Original.occurrencetable <- read.csv(occurrencesites,header=TRUE)

Original.backgroundtable <- read.csv(backgroundsites,header=TRUE)

# Identify the column position of the uncorrelated variables in the
# original tables

Extracted.occurrencetable <-
Original.occurrencetable[,colnames(Original.occurrencetable)%in%c(uncorrelated.variables,"species","longitude","latitude")]

Extracted.backgroundtable <-
Original.backgroundtable[,colnames(Original.backgroundtable)%in%c(uncorrelated.variables,"species","longitude","latitude")]

write.table(Extracted.occurrencetable,file=occurrencesites, append = FALSE, quote = FALSE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = "double")

write.table(Extracted.backgroundtable,file=backgroundsites, append = FALSE, quote = FALSE, sep = ",",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = "double")

}
