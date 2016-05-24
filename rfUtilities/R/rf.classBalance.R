#' @title Random Forest Class Balance (Zero Inflation Correction) Model
#' @description Implements Evans & Cushman (2008) Random Forests class-balance (zero inflation) modeling approach. 
#'                                                                                                                                                                                                    
#' @param ydata Response variable using index (i.e., [,2] or [,"SPP"] )                        
#' @param xdata Independent variables using index (i.e., [,3:14] or [3:ncol(data)] )
#' @param p p-value of covariance convergence (do not recommend changing)
#' @param cbf        Scaling factor to test if problem is imbalanced, default is size of majority class * 3
#' @param sf         Majority subsampling factor. If sf=1 then random sample would be perfectly balanced with smallest class [s|0=n|1] whereas; sf=2 provides [s|0=(n|1*2)]
#' @param ...        Additional arguments passed to randomForest
#'  
#' @return A list class object with the following components:
#'  @return   model Final Combined Random Forests ensemble
#'  @return   oob.error Median out-of-bag error
#'  @return   confusion Confusion matrix (summed across models)
#'  @return   pcc Percent correctly classified
#'
#' @note
#' This approach runs independent Random Forest models using random subsets of the majority class until covariance convergences on full data. The final model is obtained by combining independent ensembles.  
#'
#' @author Jeffrey S. Evans   <jeffrey_evans<at>tnc.org>
#'
#' @references
#' Evans, J.S. and S.A. Cushman (2009) Gradient Modeling of Conifer Species Using Random Forest. Landscape Ecology 5:673-683.
#' @references
#' Evans J.S., M.A. Murphy, Z.A. Holden, S.A. Cushman (2011). Modeling species distribution and change using Random Forests CH.8 in Predictive Modeling in Landscape Ecology eds Drew, CA, Huettmann F, Wiersma Y. Springer 
#'                                                                       
#' @examples 
#' require(randomForest)
#' data(iris)
#'   iris$Species <- as.character(iris$Species)
#'     iris$Species <- ifelse(iris$Species == "setosa", "virginica", iris$Species)
#'       iris$Species <- as.factor(iris$Species)	
#' 	
#' # Percent of "virginica" observations
#' length( iris$Species[iris$Species == "virginica"] ) / dim(iris)[1]*100
#' 	
#' rf.classBalance( ydata=iris[,"Species"], xdata=iris[,1:4], cbf=1 )
#'
#' @export
rf.classBalance <- function (ydata, xdata, p=0.005, cbf=3, sf=2, ...) 
 {
  if (  class(ydata) != "factor" ) { ydata <- as.factor(ydata) }
  CompCov <- function(m1, m2, pVal=p) {
       k = 2
        p = 2
         n1 = dim(m1)[1]
          n2 = dim(m2)[1] 
           n = n1 + n2
            s1 <- crossprod(m1[1:dim(m1)[1]])
             s2 <- crossprod(m2[1:dim(m2)[1]])
              c1 = (1/(n1-1)) * s1
              c2 = (1/(n2-1)) * s2
             c3 = (s1+s2)/(n-k)
            d = det(c3)
            d1 = det(c1)
           d2 = det(c2) 
          m = ( (n - k) * log(d) ) - ( (n1 - 1) * log(d1) + (n2 - 1) * log(d2) )
         h = 1 - ((2 * p * p + 3 * p - 1) / (6 * (p + 1) * (k - 1)) * 
		         (1 / (n1 - 1) + 1 / (n2 - 1) + 1 / (n - k)))
        chi = round(abs(m * h),digits=6)
        dfree = p * (p + 1) * (k - 1) / 2
        print( paste("EQUIVALENCE p", chi, sep=": ") )
      if ( (chi <= pVal ) == TRUE & (i > 2) |  (i > 20)  == TRUE ) { 
        ( "TRUE" )
      } else {
        ( "FALSE" ) 
      }
    }  		 
    y <- ydata
    x <- xdata  		 		 
    class.ct <- table(y)
    maj.class <- names(class.ct)[which.max(class.ct)]; maj.idx <- which.max(class.ct) 
    min.class <- names(class.ct)[which.min(class.ct)]; min.idx <- which.min(class.ct)  
      if ( ( class.ct[maj.idx] <= class.ct[min.idx] * cbf ) == TRUE) 
        stop("CLASSES ARE BALANCED!")  	 
        tmp.data <- data.frame(y, x)
	    majority <- tmp.data[tmp.data[,"y"] == maj.class ,]       
        minority <- tmp.data[tmp.data[,"y"] == min.class ,]    
	    all.cov <- stats::cov(majority[,names(x)])     
	    test <- as.data.frame(array(0, dim=c( 0, dim(tmp.data)[2] )))
        names(test) <- names(majority) 
      if ( !is.na(match("rf.model",ls()))) rm(rf.model)
        n <- dim(minority)[1] * sf                 
    i=0; converge = c("FALSE")  
      while (converge != "TRUE" )
       {
       i=i+1
        ns <- sample(1:nrow(majority), n) 
        class.sample <- majority[ns, ]
        mdata <- rbind(minority, class.sample)   
          if (  class(mdata[,1]) != "factor" ) { mdata[,1] <- as.factor(mdata[,1]) }
          if ( !is.na(match("rf.model",ls()))) {               
            rf.fit <- randomForest::randomForest(x=mdata[,2:ncol(mdata)], y=mdata[,1], ...)                           
            rf.model <- randomForest::combine(rf.fit, rf.model)           
            OOB <- ( OOB + stats::median(rf.fit$err.rate[,1]) ) 
            CM <- (CM + rf.fit$confusion)                 
          } else {
            rf.model <- randomForest::randomForest(x=mdata[,2:ncol(mdata)], y=mdata[,1], ...)  
            OOB <- stats::median(rf.model$err.rate[,1]) 
            CM <- rf.model$confusion                                    
          }
        test <- rbind(test, class.sample)    
        test.cov <- stats::cov( test[,names(x)] )
        converge <- CompCov(all.cov, test.cov)  
    }
        OOB <- OOB / i
      CM[,3] <- CM[,3] / i
    PCC <- (sum(diag(CM))/sum(CM[1:dim(CM)[1],1:dim(CM)[1]])) * 100
  list( model=rf.model, OOB.error=OOB, confusion=CM, pcc=PCC )
}
