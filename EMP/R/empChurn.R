empChurn <- function(scores, classes, alpha = 6, beta = 14, clv = 200, d = 10, f = 1) {
  # This software comes with absolutely no warranty. Use at your own risk.
  #
  # Adapted from:
  # Verbraken, T., Wouter, V. and Baesens, B. (2013). A Novel Profit Maximizing 
  # Metric for Measuring Classification Performance of Customer Churn Prediction
  # Models. Knowledge and Data Engineering, IEEE Transactions on. 25 (5): 
  # 961-973.
  # Available Online: http://ieeexplore.ieee.org/iel5/69/6486492/06165289.pdf
  #
  # Estimates the EMP for customer churn prediction, considering constant CLV
  # and a given cost of contact f and retention offer d.
  #
  #
  # Arguments:
  #   scores: A vector of probability scores.
  #   classes: A vector of true class labels.
  #   p0: Percentage of cases on the first point mass of the LGD distribution 
  #   (complete recovery).
  #   p1: Percentage of cases on the second point mass of the LGD distribution 
  #   (complete loss).
  #   ROI: Constant ROI per granted loan. A percentage.
  # Value:
  #   An EMP object with two components.
  #     EMP: The Expected Maximum Profit of the ROC curve at EMPfrac cutoff.}
  #     EMPfrac: The percentage of cases that should be excluded, that is, 
  #     the percentual cutoff at EMP profit.  
  
  roc = .empRocInfo(scores, classes)
  
  egamma <- alpha / (alpha+beta);
  delta  <- d/clv
  phi    <- f/clv
  
  B <- function(x,a,b){ pbeta(x,a,b)*beta(a,b) } # incomplete beta function
  
  gamma <- c(0,(roc$pi1*(delta+phi)*diff(roc$F1)+roc$pi0*phi*diff(roc$F0))/(roc$pi0*(1-delta)*diff(roc$F0)))
  gamma <- c(gamma[gamma<1],1)
  
  indE <- max(which((gamma<egamma)==T))
  MP = clv*((egamma*(1-delta)-phi)*roc$pi0*roc$F0[indE]-(delta+phi)*roc$pi1*roc$F1[indE]);
  MPfrac = roc$pi0*roc$F0[indE] + roc$pi1*roc$F1[indE];
  
  gammaii <- head(gamma, n=-1)
  gammaie <- tail(gamma, n=-1)
  F0 <- roc$F0[1:length(gammaii)]
  F1 <- roc$F1[1:length(gammaii)]
  
  contr0 <- ( clv*(1-delta)*roc$pi0*F0)                       * (B(gammaie,alpha+1,beta)-B(gammaii,alpha+1,beta))/B(1,alpha,beta)
  contr1 <- (-clv*(phi*roc$pi0*F0 + (delta+phi)*roc$pi1*F1))  * (B(gammaie,alpha  ,beta)-B(gammaii,alpha,  beta))/B(1,alpha,beta)
  EMP <- sum(contr0+contr1)
  EMPfrac <- c(t(((B(gammaie,alpha,beta)-B(gammaii,alpha,beta))/B(1,alpha,beta))) %*% (roc$pi0*F0+roc$pi1*F1))
  
  list(MP=MP, MPfrac=MPfrac, EMP=EMP, EMPfrac=EMPfrac)
}

