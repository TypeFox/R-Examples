                                        # Function used within the calc.r function
Calculate.Hj.values <- function(allelefrequency2.1table){

  # splitting this table into separate tables for each population
  allelefrequency3 <- split(allelefrequency2.1table,allelefrequency2.1table$population)

#  Herein I have to use another function sapply
Hj.one.population <- sapply(allelefrequency3,Calculate.Hj.values.1pop)
actual.locus <- as.character((allelefrequency2.1table$locus)[1])
Hj.values <- cbind(as.character(actual.locus),names(allelefrequency3),as.numeric(as.vector(Hj.one.population)))
  # The Hj-values are combined with a column of the names of the actual population
                                                                            # and a column of the names of the actual locus.
Hj.values
}
