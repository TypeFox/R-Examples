### Routine written by:
### Richard Nielsen <nielsen.rich@gmail.com>

##################################################################################
## Imbens and Rubin procedure ####################################################

  ## DROP IN DEVIANCES FUNCTION ##################################################
    ## modified code from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2008q3/001175.html
  drop.deviance.test <- function(full.model, reduced.model){
    L1 <- logLik(full.model)
    L0 <- logLik(reduced.model)
    d0 <- deviance(reduced.model)
    d1 <- deviance(full.model)
    LR <- d0 - d1
    df <- attr(L1, "df") - attr(L0, "df")
    return(list("d1"=d1, "d0"=d0, "df"=df, "LR"=LR, "pvalue" = pchisq(LR, 1, lower.tail = FALSE)))
  }
  ################################################################################


pscoreSelect <- function(formula, data, C.L=2*(pnorm(-1,0,1)), C.Q=0.1){

    ## extract the treatment and possible control variables

  treatvarname <- unlist(strsplit(as.character(formula)[2], split="+", fixed = T))
  covariates <- unlist(strsplit(as.character(formula)[3], split="+", fixed = T))
    ## take the leading and trailing spaces out of the variable names
  covariates <- gsub(" ", "", covariates, fixed = T)
    ## get rid of any accidental variables that are just spaces
  if(length(which(covariates==""))>0){
    covariates <- covariates[-c(which(covariates==""))]
  }
 	
    ## there is a dependency below on the data being called "dat"
  dat <- data

    ## Enumerate all possible combinations of the main effects

    ## k is the number of variables
  k <- length(covariates)
  pb <- txtProgressBar(min = 0, max = k, initial = 1, style = 3)

  base.formula <- (paste(treatvarname, "~ 1"))
  current.formula <- base.formula
  current.model <- (glm(as.formula(noquote(current.formula)), family="binomial", data=dat))

  cat("\nTesting main effects\n")
  for(i in 1:length(covariates)){
    if(i>i) 
     setTxtProgressBar(pb, i)	
    test.formula <- (paste(current.formula, "+", covariates[i]))
    test.model <- (glm(as.formula(noquote(test.formula)), family="binomial", data=dat))
    if(drop.deviance.test(test.model,current.model)$pvalue < C.L){
      current.formula <- test.formula
      current.model <- test.model
    }
  }
  current.formula    
  summary(current.model)
  close(pb)

    ## How to specify every possible interaction:
  
    ## Specify the first order terms in the final model
  final.covariates <- unlist(strsplit(as.character(as.formula(noquote(current.formula)))[3], split="+", fixed = T))[-1]
    ## sub out the spaces
  final.covariates <- gsub(" ", "", final.covariates, fixed = T)
    ## specify the second k
  k2 <- length(final.covariates)
    ## specify the first stage data (to exclude squares of any factors)
  final.data <- subset(data, select=final.covariates)


    ## This elaborates all of the combinations
    ## the possible interactions with other variables
    ## if there is more than one variable left
  if(k2>=2){
    interactions <- t(combn(k2,2))
      ## This makes sure I'm including the right interactions
    interaction.names <- matrix(final.covariates[interactions],
                            nrow(interactions),ncol(interactions))
    int.names <- paste(interaction.names[,1],":", interaction.names[,2], sep="")
    } else {
    int.names <- c()
  }
    ## the possible squared terms
  squaredterms <- t(rbind(seq(k2),seq(k2)))

  squaredterm.names <- matrix(final.covariates[squaredterms],
                            nrow(squaredterms),ncol(squaredterms))

    ## First, see which variables are factors
    ## Then, don't use these because they don't "square"
  factor.holder <- rep(NA,ncol(final.data))
  for(i in 1:ncol(final.data)){
    factor.holder[i] <- is.factor(final.data[,i])
  }
  factor.holder
  sqr.names <- paste("I(",squaredterm.names[factor.holder!=T,1],"*",
                   squaredterm.names[factor.holder!=T,1],")", sep="")
  int.sqr.names <- c(int.names, sqr.names)


    ## start testing the interactions
    ## "current.formula" is the final formula from the previous section
  current.model <- (glm(as.formula(noquote(current.formula)), family="binomial", data=dat))

  cat("\nTesting interactions\n")

  pb <- txtProgressBar(min = 1, max = length(int.sqr.names), initial = 1, style = 3)

  for(i in 1:length(int.sqr.names)){
    if(i>1)
     setTxtProgressBar(pb, i)	
    test.formula <- (paste(current.formula, "+", int.sqr.names[i]))
    test.model <- (glm(as.formula(noquote(test.formula)), family="binomial", data=dat))
    if(drop.deviance.test(test.model,current.model)$pvalue < C.L){
      current.formula <- test.formula
      current.model <- test.model
    }
  }
  close(pb)

  return(as.formula(noquote(current.formula)))
   
}


