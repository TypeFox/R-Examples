#' compare_AIC compares the AIC of several outputs obtained with the same data.
#' @title Compares the AIC of several outputs
#' @author Marc Girondot
#' @return A list with DeltaAIC and Akaike weight for the models.
#' @param ... Successive results to be compared as lists.
#' @param factor.value The $value of the list object is multiplied by factor.value to calculate AIC.
#' @description This function is used to compares the AIC of several outputs obtained with the same data but with different set of parameters.\cr
#' The parameters must be lists with $aic or $AIC or $value and $par elements or if AIC(element) is defined.\cr
#' if \code{$value} and \code{$par} are present in the object, the AIC is calculated as \code{2*factor.value*value+2*length(par)}. If \code{$value} is -log(likeihood), then factor.value must be 1 and if \code{$value} is log(likeihood), then factor.value must be -1.\cr
#' If several objects are within the same list, their AIC are summed.\cr
#' For example, compare_AIC(g1=list(group), g2=list(separe1, separe2)) can be used to compare a single model onto two different sets of data against each set of data fited with its own set of parameters.\cr
#' Take a look at \code{ICtab} in package \code{bbmle} which is similar.
#' @examples
#' library("HelpersMG")
#' # Here two different models are fitted
#' x <- 1:30
#' y <- rnorm(30, 10, 2)+log(x)
#' plot(x, y)
#' d <- data.frame(x=x, y=y)
#' m1 <- lm(y ~ x, data=d)
#' m2 <- lm(y ~ log(x), data=d)
#' compare_AIC(linear=m1, log=m2)
#' # Here test if two datasets can be modeled with a single model
#' x2 <- 1:30
#' y2 <- rnorm(30, 15, 2)+log(x2)
#' plot(x, y, ylim=c(5, 25))
#' plot_add(x2, y2, col="red")
#' d2 <- data.frame(x=x2, y=y2)
#' m1_2 <- lm(y ~ x, data=d2)
#' x_grouped <- c(x, x2)
#' y_grouped <- c(y, y2)
#' d_grouped <- data.frame(x=x_grouped, y=y_grouped)
#' m1_grouped <- lm(y ~ x, data=d_grouped)
#' compare_AIC(separate=list(m1, m1_2), grouped=m1_grouped)
#' @export


compare_AIC <- function(..., factor.value=1) {
result <- list(...)

if (is.list(result) & length(result)==1) result <- unlist(result, recursive=FALSE)

if (!is.null(result)) {
	if (!is.list(result) || (is.null(names(result))) || (any(names(result)==""))) {
		print("The results must be included within a list with names; see example.")
		return(invisible())
	} else {
	out<-NULL
	l<-length(result)
		if (l<2) {
			print("A least two results must be provided.")
			return(invisible())
		} else {
			aic<-NULL
			name<-names(result)
			for (i in 1:l) {

        encours <- result[i]
        t <- (class(try(AIC(encours[[1]]), silent=TRUE))=="try-error")
        if (t & is.null(encours[[1]]$aic) & is.null(encours[[1]]$AIC) & is.null(encours[[1]]$value))
          encours <- encours[[1]]
        
 
# je n'ai pas de AIC ou rien pour le calculer
          sumAIC <- 0
            AICencours <- NULL
            for (j in 1:length(encours)) {
              encours2 <- encours[[j]]
              t <- (class(try(AIC(encours2), silent=TRUE))=="try-error")
              if (!t) AICencours <- AIC(encours2)
              if (!is.null(encours2$AIC)) AICencours <- encours2$AIC
              if (!is.null(encours2$aic)) AICencours <- encours2$aic
              if (!is.null(encours2$value) & !is.null(encours2$par))
                AICencours <- 2*factor.value*encours2$value+2*(length(encours2$par))
            if (is.null(AICencours))
              {
                print(paste("Object", name[i], "has not the required format"))
                return(invisible())
              }
            sumAIC <- sumAIC + AICencours
            }
			aic <- c(aic, sumAIC)	
			
			}
			
			bestaic<-min(aic)
			ser<-which(aic==bestaic)
			deltaaic<-aic-bestaic
			aw<-exp(-0.5*deltaaic)
			saw=sum(aw)
			aw<-aw/saw
			
			out<-data.frame(cbind(AIC=aic, DeltaAIC=deltaaic, Akaike_weight=aw), row.names=name)
			print(paste("The lowest AIC (",sprintf("%.3f", bestaic) ,") is for series ", name[ser], " with Akaike weight=", sprintf("%.3f", aw[ser]), sep=""))
			
			return(out)
		}
	}
}
}
