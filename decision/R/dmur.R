#' Statistical Decision Analysis: with the help of this package we construct the payoff table and the opportunity loss table and calculate various statistical parameters like EMVs, MEMV, EPPI, EVPI, EOL which help in taking effective decisions.
#'@param dAmount A data vector to store values
#'@param prob A data vector to store the probabilities
#'@param sPrice A numeric value to store selling price
#'@param pPrice A numeric value to store preparation/ purchase price
#'@examples
#'dmur(c(50, 100, 150, 200), c(0.2,0.4,0.3,0.1), 20, 10)
#'
#' @references
#' 1. Khalili, k., Damghani, M. T., Taghavifard, R., Tavakkoli M., (2009) Decision making under uncertain and risky situations, Enterprise risk management symposium monograph society of Actuaries- Schaumburg, Illinois, vol. 15.
#'
#' 2. Marakas, G. M. (2006) Decision support system- In the 21st century, 2nd edition, Pearson education, New Delhi.
#'
#' 3. Turban, E., Aronson, J. E. and Liang, T.P. (2006) Decision support system and intelligent systems, Pearson education, New Delhi.

#'
#'@return MEMV, bestValue, EPPI, EVPI, EOL
#'@export

dmur<- function(dAmount, prob, sPrice, pPrice)
{
  p<-dAmount
  c<-p
  pt<-{}
  SP <- sPrice; PC <-pPrice; Pro <-(SP-PC); # SP: selling price, PC: preparation cost, Pro: profit

  for (i in 1:length(p)) {
    for (j in 1:length(p)){
      if (c[j]<p[i])
      {
        pt<-c(pt,((c[j]-p[i])*PC + c[j]*Pro))
      }
      else if (c[j]>p[i]) {
        pt<-c(pt, p[i]*Pro)

      }
      else {
        pt<-c(pt, p[i]*Pro)
      }
    }}

  pfM = matrix(pt, nrow = length(p)) # construct 4 by 4 matrix
  rownames(pfM) = c(p)
  colnames(pfM) = c(p)
  pfM <- cbind(pfM, prob) # calculation of probability and entering into table
  cat("The payoff table\n");print(pfM); cat("\n")

  #Expected monitory value (EMV)
  EMV<-apply(pfM[,(length(p)+1)]*pfM[,1:length(p)], 2, sum)
  MEMV<-max(EMV)
  bestValue<- p[which(EMV == MEMV)]

  cat("The expected monitory values\n"); print(apply(pfM[,(length(p)+1)]*pfM[,1:length(p)], 2, sum)); cat("\n")
  cat("The maximum (optimum) expected monitory value (MEMV) is: ", MEMV);cat("\n")
  cat("The best choice is: ", bestValue)

  #Expected profit with perfect information (EPPI)
  EPPI<- sum(apply(pfM[,1:length(p)],2,max)*pfM[,(length(p)+1)]);cat("\n")
  cat("\nThe expected profit with perfect information (EPPI) is: ", EPPI); cat("\n")

  #Expected value of the perfect information (EVPI)
  EVPI<-EPPI- MEMV;cat("\n")
  cat("The expected value of the perfect information (EVPI) is: ", EVPI);

  # Expected opportunity loss (EOL)
  # Construction of opportunity loss table
  pt1<-{}
  for (i in 1:length(p)) {
    for (j in 1:length(p)){
      if (c[j]<p[i])
      {
        pt1<-c(pt1,((p[i]- c[j])*PC))
      }
      else if (c[j]>p[i]) {
        pt1<-c(pt1, (c[j]-p[i])*Pro)

      }
      else {
        pt1<-c(pt1, 0)
      }
    }}

  pfM1 = matrix(pt1, nrow = length(p)) # 4 by 4 matrix
  rownames(pfM1) = c(p) # Enter row names
  colnames(pfM1) = c(p) # Enter column names
  pfM1 <- cbind(pfM1, prob)
  cat("\n\nThe opportunity loss table \n"); print(pfM1)

  OLT<-apply(pfM1[,(length(p)+1)]*pfM1[,1:length(p)], 2, sum)
  EOL<-min(OLT)
  cat ("\nThe expected opportunity loss (EOL): ", EOL); cat("\n")
  result_list <- list("MEMV"= MEMV,"bestValue"= bestValue, "EPPI" = EPPI, "EVPI" = EVPI, "EOL" = EOL)
  return(result_list)
}
