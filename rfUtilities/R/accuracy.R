#' @title Accuracy
#' @description Classification accuracy measures for pcc, kappa, users accuracy, producers accuracy
#'    
#' @param x   vector of predicted data or table/matrix contengency table
#' @param y   vector of observed data, if x is not table/matrix contengency table 
#'
#' @return A list class object with the following components:
#' \itemize{ 
#' \item   PCC                    percent correctly classified (accuracy)
#' \item   users.accuracy         The users accuracy  
#' \item   producers.accuracy     The producers accuracy
#' \item   kappa                  Cohen's Kappa (chance corrected accuracy)
#' \item   sensitivity            Sensitivity
#' \item   specificity            Specificity 
#' \item   plr                    Positive Likelihood Ratio   
#' \item   nlr                    Negative Likelihood Ratio  
#' \item   typeI.error            Type I error
#' \item   typeII.error           Type II error
#' \item   gain                   Information gain
#' \item   f.score                F-score
#' \item   auc                    Area Under the ROC Curve
#' \item   confusion              A confusion matrix 
#'  }
#' 
#' @note sensitivity = true positives / ( true positives + false positives ) 
#' @note specificity = true negatives / ( true negatives + false positives ) 
#' @note Type I error = 1 - specificity
#' @note Type II error = 1 - sensitivity
#' @note Positive Likelihood Ratio  = sensitivity / (1 - specificity) 
#' @note Negative Likelihood Ratio  = (1 - sensitivity) / specificity
#' @note gain  = sensitivity / ( (true positives + true negatives) / n )
#' @note recall = true positives / (true positives + false negatives)
#' @note auc = (tpr - fpr + 1) / 2
#' @note F-Score = 2 * (precision * recall) / (precision + recall) 
#' @note Hanssen-Kuiper skill score (aka true score statistic) = [(tp * tn) - (fp * fn)] / [(tp + fn) + (fp + tn)]
#' @note The true skill score has an expected -1 to +1, with 0 representing no discrimination.   
#' 
#' @author Jeffrey S. Evans    <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Cohen, J. (1960) A coefficient of agreement for nominal scales. Educational and Psychological Measurement 20 (1):37-46
#' Cohen, J. (1968) Weighted kappa: Nominal scale agreement with provision for scaled disagreement or partial credit. Psychological Bulletin 70 (4):213-220  
#' Powers, D.M.W., (2011). Evaluation: From Precision, Recall and F-Measure to ROC, Informedness, Markedness & Correlation. Journal of Machine Learning Technologies 2(1):37-63.
#'
#' @examples 
#'  # Two classes (vector)
#'  observed <- sample(c(rep("Pres",50),rep("Abs",50)), 100, replace=TRUE )
#'  accuracy(observed[sample(1:length(observed))], observed)
#'
#'  # Two classes (contingency table)
#' accuracy(cbind(c(15,11), c(2,123)))
#'
#'  # Multiple classes
#'  accuracy(iris[sample(1:150),]$Species, iris$Species)
#'
#' @export 
accuracy <- function (x, y) {
  if(inherits(x, c("table", "matrix"))) {
    t.xy <- x 
  } else {
    t.xy <- table(x, y)
  }
    if(inherits(t.xy, "matrix")) { 
      rownames(t.xy) <- c("0","1")
	  colnames(t.xy) <- c("0","1")
    }
	tc <- match(colnames(t.xy), rownames(t.xy))
	mtc <- matrix(ncol = ncol(t.xy), nrow = length(tc[tc == "NA"]), 0)
        nrn <- colnames(t.xy)[is.na(tc) == TRUE]
          rownames(mtc) <- nrn
            t1 <- rbind(t.xy, mtc)
              tr <- match(rownames(t1), colnames(t1))
                mtr <- matrix(nrow = nrow(t1), ncol = length(tr[tr == "NA"]), 0)
                  ncn <- rownames(t1)[is.na(tr) == TRUE]
                colnames(mtr) <- ncn
              t2 <- cbind(t1, mtr)
            sr <- sort(rownames(t2))
          mr <- match(sr, rownames(t2))
        t3 <- t(t2[mr, ])
      sc <- sort(rownames(t3))
    mc <- match(sc, rownames(t3))
    t4 <- t(t3[mc, ])
      agree <- diag(t4)
        prod1 <- apply(t4, 1, sum)
          prod2 <- agree/prod1
             user1 <- apply(t4, 2, sum)
            user2 <- agree/user1
          N <- sum(t4)
        k1 <- sum(agree)
      k2 <- sum(prod1 * user1)
    khat <- abs(((N * k1) - k2)/(N^2 - k2))
	
    if( dim(t.xy)[1] == 2 ) {	
      n = sum(t.xy)
      prevalence <- t.xy[1] / n 
      precision <- t.xy[1] / (t.xy[1] + t.xy[3])  
      tpr <- t.xy[1] / ( t.xy[1] +  t.xy[2] )    
      tnr <- t.xy[4] / ( t.xy[4] +  t.xy[3] )    
      type1.error <- 1 - tnr                     
      type2.error <- 1 - tpr                     
      plr <- tpr / (1 - tnr)                     
      nlr <- (1 - tpr) / tnr                     
      f.score <- 2 * (precision * tpr) / (precision + tpr)
	  true.skill <- ( (t.xy[1] * t.xy[4]) - (t.xy[3] * t.xy[2]) ) / 
	                  ( (t.xy[1] + t.xy[2]) * (t.xy[3] + t.xy[4]) )
	  gain <- precision / ( (t.xy[1] + t.xy[4]) / n ) 
	  mcc <- ( t.xy[1] * t.xy[4] - t.xy[3] * t.xy[2] ) /
              sqrt( (t.xy[1] + t.xy[3]) * 
			    (t.xy[1] + t.xy[2]) * 
			    (t.xy[4] + t.xy[3]) * 
			    (t.xy[4] + t.xy[2]) )
      auc <- (tpr - ( t.xy[3] / ( t.xy[2] +  t.xy[3] ) ) + 1) / 2
	  acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100,
                 auc = auc, 	
	             users.accuracy = round(user2 * 100, 1),  
	             producers.accuracy = round(prod2 * 100, 1),
                 kappa = round(khat, 4),
				 true.skill = true.skill, 
				 sensitivity = tpr,
				 specificity = tnr,
				 plr = plr,
				 nlr = nlr,
				 typeI.error = type1.error,
				 typeII.error = type2.error,
				 f.score = f.score,
				 gain = gain,
                 matthews = mcc, 				 
                 confusion = t.xy)	
        } else {
	  acc <- list( PCC = (sum(diag(t4))/sum(t4)) * 100, 
	             users.accuracy = round(user2 * 100, 1),  
	             producers.accuracy = round(prod2 * 100, 1),
                 kappa = round(khat, 4),
                 confusion = t.xy)
	}			 
	class(acc) <- c("accuracy", "list") 			   
          return( acc )
}
