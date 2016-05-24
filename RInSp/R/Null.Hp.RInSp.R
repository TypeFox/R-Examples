Null.Hp.RInSp = function(dataset, prop = "sum"){
  #
  # The procedure let you build an integer/count dataset of class RInSp:
  #  - by fixing the number of items/resources per individual and the population diet;
  #  - by providing a specific dataset of class RIS and and option of "sum"
  #    or "average"} for the population diet.
  # It builds a table with a specified number of individuals and for each, a resources
  # vector is created using a multinomial distribution with probability equal to the
  # specified population diet and number of items.
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
if (mode(dataset) == "list") 
     { 
       cat("\n For big data set this can take a long time! \n")
      # when input is a RIS object from the input.RIS procedure
       if (class(dataset) != "RInSp") stop("The input must be an object of class RInSp.") 
       if (dataset$data.type != "integer") stop("Only integer data type is accepted.")
       if (mode(prop) %in% c("numeric", "character") == FALSE) stop("Wrong data type for diet proportions.")
       if (mode(prop) == "numeric")
         {
          if (length(prop) <2) stop("The diet proportion vector must be bigger then one resource.")
          popdiet = abs(prop) / sum(abs(prop))
          } else
         {
          if (prop %in% c("sum", "average") == TRUE) 
            { 
              if (prop == "sum") {
                  popdiet = apply(dataset$resources, 2, sum) / sum(dataset$resources)} else {
                  popdiet = apply(dataset$proportion, 2, mean)}
             }
            else 
            {
             stop("Wrong proportion option.")
            }
         }
    NInds = dataset$num.individuals
    if (NInds <= 1) stop("The input dataset is to short.")
    rows.counts = apply(dataset$resources, 1, sum)
    } else {
    cat("\n For big data set this can take a long time! \n")
    # When dataset is a vector of rows.counts, prop must be a vector type
    if (mode(dataset) != "numeric") stop("The input dataset must be a numeric vector.")
     rows.counts = as.integer(abs(dataset))
     NInds = length(rows.counts)
     if (NInds <= 1) stop("The input dataset is to short.")
     if (mode(prop) == "numeric")
        {
         if (length(prop) <2) stop("The diet proportion vector must be bigger then one resource.")
         popdiet = abs(prop) / sum(abs(prop))
        } else
         {
         stop("Wrong proportion option.")
         }
    }

exit = 0
count2= 0
count = 0 # variable for progress-bar
progress = 0 # variable for progress-bar
pb = txtProgressBar(min=0, max= 80, char= "+", style = 3)
while (exit == 0) {
  Nulldata = data.frame()
  for(i in 1: NInds){
    Nulldata = rbind(Nulldata, t(rmultinom(1, rows.counts[i], popdiet)))}
  rows = sum(apply(Nulldata, 1, sum) > 0)
  cols = sum(apply(Nulldata, 2, sum) > 0)
  count = count + 1
  count2 = count2 +1
  if (count == 5) {
    progress = progress + 1
    setTxtProgressBar(pb, progress)
    count = 0}
  if ((rows  == length(rows.counts)) & (cols == length(popdiet))) exit = 1
  }
close(pb)
Nulldata.RInSp = import.RInSp(Nulldata)
cat("\n Results reached after", count2, "iterations \n")
return(Nulldata.RInSp)
}

