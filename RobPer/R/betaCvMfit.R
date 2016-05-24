betaCvMfit <- function(data,CvM=TRUE, rob=TRUE) {
    # Robustly fits a Beta distribution to data using Cramer-von-Mises (CvM)
    # distance minimization.
	
    if(any(is.na(data))) data <- data[-which(is.na(data))]  # delete NA-observations from data
    data[which(data<0)]<-0                                  # set negative observations to zero
    if(rob) {
        x. <- median(data)
        s.. <- mad(data)^2
    } else {
        x. <- mean(data)
        s.. <- var(data)
    }
    if(s..!=0) {
        shape1  <- max(0.00001,-(x.*(-x.+x.^2+s..))/s..)  # moment estimators for the beta distribution forced to be positive
        shape2  <- max(0.00001,(shape1-shape1*x.)/x.)
        erg<- c(shape1, shape2)
        CvMbeta <- function(x, vals) {
            nn<- length(vals)
            return(sum((suppressWarnings(pbeta(vals, shape1=x[1], shape2=x[2]))-(rank(vals)-0.5)/nn)^2)/nn)
            #+ 1/(12*nn^2)
        }
        # Minimization criterion to minimize the Cramer-von-Mises (CvM) distance between 
        # a Beta distribution with parameters (x[1], x[2]) and a data sample "vals".  
        # Based on R-Code provided by Brenton C. Clarke to minimize the CvM distance
        # between a Gamma distribution and a data sample.
        # See Clarke, McKinnon and Riley (2012): A fast robust method for fitting gamma
        # distributions. Statistical Papers, 53 Nr. 4, 1001-1014.
        if(CvM) { ff<-function(x) CvMbeta(x, vals=data)
            shapes  <- optim(c(shape1, shape2), ff)$par  # minimize CvM distance
            erg<- shapes
        }
        return(erg)
    } else {return(c(NA, NA))} # if variance of the data is zero, no moment estimators can be calculated and no parameters are determined.
}
