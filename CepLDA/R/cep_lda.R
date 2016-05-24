#' @name cep.lda
#' @aliases cep.lda
#' @title Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability
#' 
#' @description
#' The main program. 
#' 
#' @usage
#' cep.lda(y,x,xNew,L,mcep,nw,k,cv,tol)
#' 
#' @param y n-vector indicating group membership of training time series.
#' @param x \emph{N} by \emph{n} matrix containing \emph{n} training time series each with length \emph{N}.
#' @param xNew \emph{N} by \emph{nNew} matrix, containing nNew time series whose memberships are predicted.
#' @param L Number of cepstral coefficients used in the lda.  If FALSE, cross-validation is used for the data driven selection of L.  Default is FALSE.
#' @param mcep Maximum number of cepstral coefficients considerd. Default is set to 10.
#' @param nw Width of tapers used in multitaper spectral estimation. Default is set to 4.
#' @param k Number of tapers used in multitaper spectral estimation. Default is set to 7.
#' @param cv If TRUE, returns results (classes and posterior probabilities) for leave-one-out cross-validation. Note that, if the prior is estimated, the proportions in the whole dataset are used. As with the standard lda function, 
#'        if used, prediction on a test data set cannot be done and weight functions are not produced (simular to the predict.lda). Default is FALSE.
#' @param tol Tolerance to decide if a matrix is singular; it will reject variables and linear combinations of unit-variance variables whose variance is less than \eqn{tol^2}.
#' @export
#' @importFrom stats formula
#' @importFrom stats predict
#' @return List with 5 elements
#' \item{C.lda}{lda output on the cepstral scale. Similar to output of \code{lda(MASS)} function.}
#' \item{cep.data}{Data frame containing cepstral coefficients and group information from training data.}
#' \item{Lopt}{Number of cepstral coefficients used.}
#' \item{lspec}{Estimated log-spectral weight functions.}
#' \item{predict}{Results of classification.  If external data xNew is supplied, these data are classified.  If not, biased classification of the training data x is returned.  For unbiased leave-out-one
#'      cross-validated classification of training data, use cv=TRUE.}
#' @references Krafty, RT (2016) Discriminant Analysis of Time Series in the Presence of Within-Group Spectral Variability.  \emph{Journal of Time Series Analysis} 
#' @author Zeda Li \email{<zeda.li@@temple.edu>}; Robert Krafty \email{<rkrafty@@pitt.edu>}
#' @seealso 
#'  \code{\link{predict.ceplda}}, \code{\link{plot.ceplda}}, \code{\link{print.ceplda}}, \code{\link{Lopt.get}}
#' @import MASS astsa class multitaper
#' @examples 
#' ## Simulate training data
#' nj = 50  #number of series in training data
#' N = 500  #length of time series
#' traindata1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' traindata2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' traindata3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' train <- cbind(traindata1$X,traindata2$X,traindata3$X)
#' y <- c(rep(1,nj),rep(2,nj),rep(3,nj))
#' 
#' ## Fit the discriminant analysis
#' fit <- cep.lda(y,train)
#' fit  #displays group means and cepstral weight functions
#' 
#' ## Discriminant plot
#' plot(fit)
#' 
#' ## Plot log-spectral weights
#' par(mfrow=c(1,2))
#' plot(fit$lspec$frq, fit$lspec$dsc[,1],type='l',xlab="frequency", ylab="log-spectral weights")
#' plot(fit$lspec$frq, fit$lspec$dsc[,2],type='l',xlab="frequency", ylab="log-spectral weights")
#' 
#' ## Bias classification of training data
#' mean(fit$predict$class == y) #classifictaion rate
#' table(y,fit$predict$class)
#' 
#' ## Fit the discriminant analysis while classifing training data via cross-validation
#' fit.cv <- cep.lda(y,train, cv=TRUE)
#' mean(fit.cv$predict$class == y) #classifictaion rate
#' table(y,fit.cv$predict$class)
#' 
#' ## Simulate test data
#' testdata1 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.01,.7),r.phi2=c(-.12,-.06),r.sig2=c(.3,3))
#' testdata2 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.5,1.2),r.phi2=c(-.36,-.25),r.sig2=c(.3,3))
#' testdata3 <- r.cond.ar2(N=N,nj=nj,r.phi1=c(.9,1.5),r.phi2=c(-.56,-.75),r.sig2=c(.3,3))
#' test <- cbind(testdata1$X,testdata2$X,testdata3$X)
#' yTest <- c(rep(1,nj),rep(2,nj),rep(3,nj))
#' 
#' ## Fit discriminant analysis and classify new data
#' fit.pre <- cep.lda(y,train,test)
#' mean(fit.pre$predict$class == y)
#' table(yTest,fit.pre$predict$class)





cep.lda <- function(
        y,              # an n-vector indicating the group a time series belongs to
        x,              # N by n matrix, containing n time series
        xNew=NULL,      # N by nNew matrix containing time series
        L=FALSE,        # Number of cepstral coefficients.  CV is used to select if false.
        mcep=10,        # the maximum number of cepstral coefficient considerd
        nw=4,           # width of the taper
        k=7,            # number of tapers
        cv=FALSE,       # If doing cross-validation
        tol=1.0e-4
        ){
        if(is.matrix(x)==FALSE)
                stop("\n x must be a matrix or cannot be a signle time series")
        if(ncol(x) < 2)
                stop("\n Must have at least two time series observations")
        if(length(y) != dim(x)[2])
                stop("\n Length of y and numbers of time series in x must agree")
        if(is.null(xNew) == 0 & cv != 0)
                stop("\n Cannot do precition of new data and leave-out-one prediction of training data.")
        
        #################################################
        ## get random spectra and cepstral coefficeints
        D.hat0 <- cep.get(y,x,nw,k)

        #############################################################
        ## Cross-validation to get the optimal number of coefficients
        if(isTRUE(L) ){
                Lopt<-L
        }else{
          Lopt <- Lopt.get(D.hat0,mcep)
        }


        ################################################
        ## Run the discriminant analysis
        b <- as.formula(paste("y ~ ",paste(colnames(D.hat0[,1:(Lopt+1)]), collapse="+"),sep = ""))
        if(isTRUE(cv)){
                C.lda <- lda(b, data=D.hat0, CV=TRUE,tol=tol)
                pre = list(class=C.lda$class, posterior=C.lda$posterior)
                cep.lda <- list(C.lda=C.lda, cep.data=D.hat0, Lopt=Lopt, lspec = NULL, predict=pre)
        } else{
                C.lda <- lda(b, data=D.hat0, CV=FALSE,tol=tol)
                # Compute Log-spectral Weights
                Q = min(Lopt,length(unique(y))-1)
                frq <- seq(from=0, to=.5, by=1/(dim(D.hat0)[2]-1))
                dsc <- matrix(NA, nrow=length(frq), ncol=min(Lopt,length(unique(y))-1))
                              for(q in 1:Q){
                                dsc[,q] = recon( C.lda$scaling[,q],frq)
                              }
                lspec <- list(dsc=dsc, frq=frq)               
                # Do predictions on new data or test data (if no new data is supplied)
                if(is.null(xNew)==0){
                  pre <- predict(C.lda, cep.get(array(0,dim=dim(xNew)[2]),xNew), prior=C.lda$prior)
                }else{
                  pre <- predict(C.lda, D.hat0, prior=C.lda$prior)
                }
                #
                cep.lda <- list(C.lda=C.lda, cep.data=D.hat0, Lopt=Lopt, lspec=lspec, predict=pre)
                

        }
        class(cep.lda) <- "ceplda"
        return(cep.lda)
}