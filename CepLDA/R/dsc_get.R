#' @importFrom stats predict
#' @importFrom graphics plot
#' @importFrom graphics par

dsc.get <- function(object, plot1=FALSE, plot2=FALSE){
        
        if(!inherits(object,"ceplda"))
                stop("Object must be of class 'ceplda'")
        
        test <- predict(object$X.lda,object$cep.data)
        frq <- seq(from=0, to=.5, by=1/(dim(object$cep.data)[2]-1))
        dsc1 <- recon( object$X.lda$scaling[,1],frq)
        dsc2 <- recon( object$X.lda$scaling[,2],frq)
        if(plot1==TRUE){
                type = as.character(object$cep.data$y)
                plot(test$x[,1],test$x[,2],pch=type, ylab="2nd Discriminant", xlab="1st Discriminant")
        }
        if(plot2==TRUE){
                par(mfrow=c(1,2))
                plot(frq,dsc1,type="l")
                plot(frq,dsc2,type="l")
                par(mfrow=c(1,1))
        }
        return(list(frq=frq, weight1=dsc1, weight2=dsc2, discriminat=test$x))            
}
