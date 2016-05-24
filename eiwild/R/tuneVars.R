#' Generate Variance parameters for proposal density
#' 
#' @description
#' \code{tuneVars} generates variance parameters dependent on the acceptance ratio of
#' drawn parameters in the Metropolis algorithm used in \code{indAggEi}. The target 
#' acceptance rate will be given by parameter \code{accRat}
#' 
#' @param form formula in this format 
#' cbind(column_1,column_2, ...,column_c)~cbind(row_1,row_2,...,row_r))
#' @param aggr data.frame with aggregate data. One district per line and one column giving one ID per district. (see Details)
#' @param indi data.frame with individual data. One district per line and one column giving one ID per district. (see Details)
#' If no individual data are present it defaults to NULL
#' @param IDCols vector of length 2 (or 1) giving the columnnames or numbers of ID column
#' @param whichPriori character string defining the hyperpriori. default="gamma"
#' @param prioriPars vector giving the parameters of the hyperpriori
#' 
#' @param accRat vector with two elements describing the wished range of the acceptance ratios
#' @param minProp numeric between 0 and 1 describing the percentage of parameters 
#' to have the wished acceptance ratios (\code{accRat}). \code{maxiter} will be the maximum iteration
#' @param maxiter numeric how many times the algorithm should run maximum. 
#' If NULL tuning will run until minProp is reached
#' @param sample the sample size to be drawn each tuning run.
#' @param verbose an integer specifying whether the progress of the sampler is printed to the screen (defaults to 0).
#' If verbose is greater than 0, the iteration number is printed to the screen every verboseth iteration
#' @param verboseTune logical if tuning iteration should be printed (default=TRUE)
#' @param improv numeric vector with 2 elements c(a,b). 
#'  standard deviation will be calculated with the last b percentages of parameters to have the wished acceptance ratio.
#'  If standard deviation is lower than a, than tuning is finished. Default is NULL
#' 
#' @param betaVars array of dimensions (rows, columns, districts) giving variance of proposal density for betavalues
#' @param alphaVars matrix of dimensions (rows, columns) giving variance of proposal density for alphavalues.
#' @param startValsAlpha matrix with dimension=c(rows,columns) giving the starting values for alpha.
#' If \code{NULL} random numbers of rdirichlet with prioriPars will be drawn
#' @param startValsBeta array with dimension=c(rows,columns,districts) giving the starting values of beta
#' If \code{NULL} random multinomial numbers with startValsAlpha or prioriPars will be draws
#' @param seed Default is NULL. Can be given the "seed" attribute of an eiwild-object to reproduce an eiwild-object

#' 
#' @details
#' \code{indi} is a districts x [(r*c)+1] data.frame containing one district per line. 
#' One column gives the ID of the districts. This will we connected to the ID column in the \var{aggr}-data.frame.
#' The rest of one line in \var{indi} is every row beside the nex row.
#' For example a 2x3 ecological Inference problem with \var{formula}
#' \code{cbind(col1,col2,col3) ~ cbind(row1,row2)}
#' will have the  row format :
#' \code{[ID, row1.col1, row1.col2, row1.col3, row2.col1, row2.col2, row2.col3]}
#' 
#' It is important that the formula names correspond to the exact column number in the \var{indi}-data.frame.
#' 
#' The \var{aggr} data.frame can have more columns than the names given in \var{formula} as long as the colnames exist
#' 
#' Priorities for finishing of tuning are as follows:
#' If \var{improv} isn't specified: \var{minProp} and \var{maxiter} are checked.
#' If \var{improv} is specified: 1) \var{improv} is checked, 2) \var{minProp} and \var{maxiter} are checked.
#' 
#' @return A list containing matrices of variance parameters for the proposal densities
#' 
#' @seealso
#' \code{\link[eiwild]{convertEiData}}, \code{\link[eiwild]{runMBayes}}, \code{\link[coda]{mcmc}}
#' \code{\link[eiwild]{tuneVars}}, \code{\link[eiwild]{indAggEi}}
#' 
#' @examples
#' \dontrun{
#' data(topleveldat)
#' out1 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),sample=10000, verbose=11000)
#' out2 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=10000, verbose=11000,
#'                  maxiter=NULL, improv=c(0.01,5))
#' out3 <- tuneVars(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"), sample=10000, verbose=11000,
#'                  maxiter=NULL, accRat=c(0.45,0.55), improv=c(0.01,5))
#' str(out3)
#' out4 <- indAggEi(form=form, aggr=aggr, indi=indi, IDCols=c("ID","ID"),
#'                  betaVars=out1$betaVars, alphaVars=out1$alphaVars,
#'                  sample=10000,thinning=1,burnin=100, verbose=1000)
#' out4
#' }
#' 
#' 
#' @export
#' @useDynLib eiwild

tuneVars <- function(form, aggr, indi=NULL, 
                     IDCols=c("ID"), whichPriori="gamma",prioriPars=list(shape=4, rate=2),
                     accRat=c(0.4,0.6), minProp=0.7, maxiter=20,
                     sample=10000, verbose=10000, verboseTune=TRUE,
                     improv=NULL,
                     betaVars=NULL, alphaVars=NULL, 
                     startValsAlpha=NULL, startValsBeta=NULL,
                     seed=NULL){
   
   conv <- convertEiData(form, aggr, indi, IDCols)
   
   ##checks
   if(any(accRat<0)||any(accRat>1))
      stop("\"accRat\" has to be within [0,1]", call.=FALSE)
   if(minProp<0||minProp>1)
      stop("\"minProp\" has to be within [0,1]", call.=FALSE)
   if(is.null(maxiter)){
      warning("\"maxiter\" is NULL. This means tuning will run indefinitely until \"minProp\" is reached.",
              immediate. = TRUE, call.=FALSE)
      maxiter<-9999999999999999
   }

   r <- ncol(conv$rowdf)
   c <- ncol(conv$coldf)
   p <- nrow(conv$rowdf)
#    ##generate fake accRat matrices for while loop
#    alphaAcc <- rep(-1, r*c)
#    betaAcc <- rep(-1, r*(c-1)*p)

   if(is.null(improv)){
      lastProps <- rep(0,5)
   } else {
      if(improv[2]<2) stop("2nd element of \"improv\" has to be greater than 1!",call.=FALSE)
      if(improv[1]<=0) stop("1st element of \"improv\" has to be positve!", call.=FALSE)
      lastProps <- rep(0,improv[2])
   }
   
   ii <- 1
   while(lastProps[1]<minProp && ii<=maxiter){
      if(verboseTune==TRUE) cat("\n","Tuning iteration: ",ii,"\n")
      draws <- runMBayes(convList=conv,
                         whichPriori=whichPriori, prioriPars=prioriPars,
                         startValsAlpha, startValsBeta,
                         betaVars, alphaVars,
                         sample, burnin=0, thinning=1,
                         verbose, seed=seed)
      
      alphaVars <- draws$alphaVars
      alphaAcc <- draws$alphaAcc
      for(rr in 1:r) ## anpassung der standardabweichungen
         for(cc in 1:c)
            alphaVars[rr,cc] <- adaptVar(alphaAcc[rr+r*(cc-1)],
                                         alphaVars[rr,cc], accRat)
      
      betaVars <- draws$betaVars
      betaAcc <- draws$betaAcc
      for(pp in 1:p) ## anpassung der standardabweichungen
         for(rr in 1:r)
            for(cc in 1:(c-1))
               betaVars[rr,cc,pp] <- adaptVar(betaAcc[rr+r*(cc-1)+r*(c-1)*(pp-1)],
                                              betaVars[rr,cc,pp], accRat)
      
      # abbruch, wenn standard deviation der letzten improv[2]-Elemente kleiner als improv[1] ist
      if(is.null(improv)){
         lastProps <- c(calcProp(alphaAcc, betaAcc, accRat), lastProps[1:4])
         cat("Perc of Param in Range: ", lastProps[1],"\n")
      } else{
         lastProps <- c(calcProp(alphaAcc, betaAcc, accRat), lastProps[1:(improv[2]-1)])
         cat("Perc of Param in Range: ", lastProps[1]," -- sd: ",sd(lastProps), "\n")
         if(sd(lastProps)<=improv[1])
            break
      }
      
      if(ii==1)
         seed <- attributes(draws)$seed
      ii <- ii+1
   }
   cat(lastProps[1],"\n")
   
   ret <- list(betaVars=betaVars, alphaVars=alphaVars,
               betaAcc=betaAcc, alphaAcc=alphaAcc)
   attr(ret, "seed")<-seed
   return(ret)
}


# function for returning percentage of acceptance ratios in range of accRat
calcProp <- function(alphaAcc, betaAcc, accRat){
   
   alphaProp<-length(which(alphaAcc>=accRat[1]&alphaAcc<=accRat[2]))
   betaProp <-length(which(betaAcc>=accRat[1]&betaAcc<=accRat[2])) 
   prop <- (alphaProp+betaProp)/length(c(alphaAcc,betaAcc))
#    cat(prop,"\n\n")
   return(prop)
}

# function for adapting variance
adaptVar <- function(acc, sdd, accRat){
   
   if(acc<accRat[1]*0.75) # bei zu kleiner acc wird Bereich verkleinert, damit mehr um den Mittelwert gezogen wird
      return(sdd*0.8)
   if(acc<accRat[1])
      return(sdd*0.95)   
   if(acc>accRat[2]*1.25) # bei zu großer acc wird Bereich vergrößert, um weiter entfernte Werte zu erreichen
      return(sdd*1.2)
   if(acc>accRat[2])
      return(sdd*1.05)
   if(accRat[1]<=acc & acc<=accRat[2])
      return(sdd*1)
}
