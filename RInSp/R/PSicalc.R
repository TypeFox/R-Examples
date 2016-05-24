PSicalc = function(dataset, pop.diet = "sum", exclude = FALSE, replicates=999){
  #
  # Measure of individual specialization, based on the average pairwise overlap of
  # the niche distribution of individuals and the population after Bolnick et al. (2003)
  
  #
  # Author: Nicola ZACCARELLI, Giorgio MANCINELLI, Dan BOLNICK
  # E-mail: nicola.zaccarelli@gmail.com,
  #         giorgio.mancinelli@unisalento.it
  #         danbolnick@mail.texas.edu
  #
  # Version: 1.0
  # Date: 10/11/2012
  #
# some checking 
if (class(dataset) != "RInSp") stop("The input must be an object of class RSI")
if (dataset$data.type != "integer") stop("Input data type must be integer.")
if (pop.diet %in% c("sum", "average") == FALSE) stop("The specified population diet type is wrong.")
if (exclude %in% c("TRUE", "T", "FALSE", "F") == FALSE) stop("The specified exclusion option is wrong.")
if (!is.double(replicates)) stop("The specified replicates option is wrong.") else replicates = abs(as.integer(replicates))
cat("\n If your dataset is big, this can take time. Please be patient. \n")
if (exclude %in% c("TRUE", "T")) 
 {
     cat("\n When individuals are excluded no Monte Carlo resampling is allowed. \n")
     pb = txtProgressBar(min=0, max=2*dataset$num.individuals,  char= "+", style = 3) # start textbar
     Y = matrix(0, dataset$num.individuals, 1)
     for (i in 1:dataset$num.individuals)
         {
           setTxtProgressBar(pb, i)
           indkeep = c(1:dataset$num.individuals) != i
           datatmp = subset(dataset$resources, subset=indkeep)
           datatmppro = subset(dataset$proportions, subset=indkeep)
           if (pop.diet == "average") dietpop = apply(datatmppro, 2, sum) / (dataset$num.individuals - 1) else dietpop = apply(datatmp, 2, sum) / sum(datatmp)
           if (i == 1) newpopdiet = dietpop else newpopdiet = rbind(newpopdiet, dietpop)
           Y[i, 1] = sum(datatmp)
         }
     tmp = abs(dataset$proportions - newpopdiet)
     PSi = 1 - 0.5*apply(tmp, 1, sum)
     IS = sum(PSi) / dataset$num.individuals
     plessthanq = dataset$proportions * (dataset$proportions <= newpopdiet)
     qlessthanp = newpopdiet * (dataset$proportions > newpopdiet)
     t1 = matrix(apply(plessthanq*(1 - plessthanq), 1, sum), dataset$num.individuals, 1)
     t2 = matrix(0, dataset$num.individuals, 1)
     t4 = matrix(0, dataset$num.individuals, 1)
     for (i in 1:dataset$num.individuals)
        { 
         setTxtProgressBar(pb, i + dataset$num.individuals)
         for (j in 1:(dataset$num.prey - 1)) 
             {
             for (k in (j + 1):dataset$num.prey)
                 {
                  t2[i, 1] = t2[i, 1] + 2*plessthanq[i, j] * plessthanq[i, k]
                  t4[i, 1] = t4[i, 1] + 2*qlessthanp[i, j] * qlessthanp[i, k]
                 }
             }
         }
     t3 = matrix(apply(qlessthanp*(1 - qlessthanp), 1, sum), dataset$num.individuals, 1)
     # Alternative formula with the correct Y term
     idiet = matrix(apply(dataset$resources,1, sum), dataset$num.individuals,1)
     varPSi = (t1 - t2)/idiet + (t3 - t4)/Y
     ris = list(PSi= matrix(PSi, dataset$num.individuals,1), IS= IS, population.diet = newpopdiet, num.individuals= dataset$num.individuals, VarPSi = varPSi, parameter = 0)
     close(pb) # close textbar  
} else 
     {
     if (pop.diet == "sum") diet.pop = 0 else diet.pop = 1
     # coerce vectors to be double to assure correct transfer to C code
     if (!is.double(dataset$resources)) dataset$resources = matrix(as.double(dataset$resources), dataset$num.individuals, dataset$num.prey)
     ris2 = .Call("PSicalc", dataset$resources, as.vector(diet.pop), as.vector(replicates), PACKAGE="RInSp")
     NRows = dataset$num.individuals
     PSi = ris2[1:NRows, 1]
     varPSi = ris2[(NRows+1):(2*NRows), 1]
     IS = ris2[(2*NRows + 1),1]
     IS.sim = apply(ris2[1:NRows, 1:(replicates +1)], 2, mean)
     IS.sim = matrix(IS.sim, length(IS.sim), 1)
     colnames(IS.sim) = "IS"
     cum.distr = ecdf(IS.sim)
     ISpvalue = cum.distr(IS.sim[1])
     if (pop.diet == "average") dietpop = apply(dataset$proportions, 2, sum) / (dataset$num.individuals) else dietpop = apply(dataset$resources, 2, sum) / sum(dataset$resources)
     ris= list(PSi = PSi, IS = IS, PSi.montecarlo = ris2[1:NRows, 1:(replicates +1)], Var.montecarlo = ris2[(NRows +1):(2*NRows), 1:(replicates+1)], VarPSi = varPSi, population.diet = dietpop, IS.pvalue = ISpvalue, montecarlo = IS.sim, num.individuals = dataset$num.individuals, parameter = 1)
  }
class(ris) = "RInSp"
cat("\n The value of PSi and its variance are: \n")
tmp = cbind(PSi, varPSi)
attributes(tmp)$dimnames[[2]] = c("PSi", "VarPSi")
print(tmp, digits=4)
cat("\n The value of IS is: ", IS)
if (exclude %in% c("TRUE", "T"))
  { cat("\n") }
 else 
  { cat("\n The p-value for P(ISsim => ISobs) is: ", ISpvalue, "\n") }
return(ris)}

