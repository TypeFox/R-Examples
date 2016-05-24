empCreditScoring <- function(scores, classes, p0=0.55, p1=0.1, ROI=0.2644){
  # This software comes with absolutely no warranty. Use at your own risk.
  #
  # Adapted from:
  # Verbraken T., Bravo C., Weber, R., and Baesens, B. 2014. Development and 
  # application of consumer credit scoring models using profit-based 
  # classification measures. European Journal of Operational Research.
  # 238 (2): 505-513.
  #
  # Estimates the EMP for credit risk scoring, considering constant ROI and
  # a bimodal LGD function with point masses p0 and p1 for no loss and total 
  # loss, respectively.
  #
  #
  # Arguments:
  #   scores: A vector of probability scores.
  #   classes: A vector of true class labels.
  #   p0:Percentage of cases on the first point mass of the LGD distribution (complete recovery).
  #   p1: Percentage of cases on the second point mass of the LGD distribution (complete loss).
  #   ROI: Constant ROI per granted loan. A percentage.
  # Value:
  #   An EMP object with two components.
  #     EMP: The Expected Maximum Profit of the ROC curve at EMPfrac cutoff.
  #     EMPfrac: The percentage of cases that should be excluded, that is, 
  #     the percentual cutoff at EMP profit.
  
  roc = .empRocInfo(scores, classes)
  
  alpha <- 1-p0-p1
  
  lambda <- c(0,(roc$pi1*ROI/roc$pi0)*diff(roc$F1)/diff(roc$F0))
  lambda <- c(lambda[lambda<1],1)
  
  lambdaii <- head(lambda,n=-1)
  lambdaie <- tail(lambda,n=-1)
  F0 <- roc$F0[1:length(lambdaii)]
  F1 <- roc$F1[1:length(lambdaii)]
  
  EMPC <- sum(alpha*(lambdaie-lambdaii)*(roc$pi0*F0*(lambdaie+lambdaii)/2 - ROI*F1*roc$pi1)) + 
    (roc$pi0*tail(F0,n=1)-ROI*roc$pi1*tail(F1,n=1))*p1
  EMPCfrac <- sum(alpha*(lambdaie-lambdaii)*(roc$pi0*F0+roc$pi1*F1)) +
    p1*(roc$pi0*tail(F0,n=1) + roc$pi1*tail(F1,n=1))
  
  list(EMPC=EMPC, EMPCfrac=EMPCfrac)
}
