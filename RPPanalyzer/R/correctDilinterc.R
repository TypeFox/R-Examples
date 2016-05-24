## adapted script for dilution series intersect correction
## (based on Daniel Kaschek's version, June 2012)

correctDilinterc <- function(dilseries, arraydesc, timeseries, exportNo) {
  
    intercepts<-getIntercepts(dilseries, arraydesc)
    correctedData<-getSignals(timeseries, intercepts, arraydesc, exportNo)
    return(correctedData)
}

    as.my <- function(v) {
		as.numeric(as.character(v))
	}

    ## getIntercepts ##
    ##
    ## Purpose: Derive intercepts of dilution series in dependence of dilSeriesID (column in sampledescription.txt)
    ##          and slide/pad/incubationRun/spottingRun number (colnames of arraydescription). 
    ##          A smoothing spline is used to extrapolate to 0.
    ##          Nonparametric bootstrap is used to estimate uncertainty of the intercept estimate.
    ##
    ## Input: dilseries = dilution series (sample type "control")
    ##        arraydesc = element "arraydescription" of RPPanalyzers data list
    ##
    ## Output: data.frame with columns for dilseries_id and slide/pad/incubationRun/spottingRun number
    ##         and antibody, estimated intercept and estimated error of intercept. 
    ##         Additionally, each single fit is plotted to "getIntercepts_Output.pdf".
    
    getIntercepts <- function(dilseries, arraydesc) {
      result <- c()
      dilSeriesID<-NULL
      pdf("getIntercepts_Output.pdf")
      # select dilseries_id and pad_slide_incubRun_spotRun
      for(A1 in unique(dilseries$dilSeriesID)) {
        subseries1 <- subset(dilseries, dilSeriesID==A1)
        for(A2 in unique(colnames(subseries1)[grep(colnames(arraydesc)[1],colnames(dilseries)):ncol(dilseries)])) {
            lseries <- subseries1[,c("concentration",A2)]
            # menge vs. intensity
            y <- as.my(lseries[,A2])
            x <- as.my(lseries$concentration)
            # intercept for interpolation on y-axis (bootstrap --> variation)
            intercept <- c()
            for(k in 1:30) {
              ssp<-NULL
              while(is.null(ssp)){ 
                mysample <- sample(1:length(y), length(y), replace=TRUE)
                df <- 3
                ssp<-tryCatch(smooth.spline(x[mysample],y[mysample], df=df), error=function(e) NULL)
              }
              intercept <- c(intercept, predict(ssp, x=0)$y)
            }
            errIntercept <- sd(intercept)
            intercept <- mean(intercept)
            # plot limits
            t <- seq(0, max(x), len=100)
            xrange <- c(0, max(x))
            yrange <- c(0, max(c(y, predict(ssp, t)$y)))
            plot(x, y, ylim = yrange, xlim = xrange, xlab="concentration", ylab="signal intensity [a.u.]")
            matplot(t, predict(ssp, t)$y, type="l", add=TRUE)
            arrows(0, intercept, 0, intercept-errIntercept, length=max(x)/50, angle=90, lwd=2)
            arrows(0, intercept, 0, intercept+errIntercept, length=max(x)/50, angle=90, lwd=2)
            title(paste(A1, arraydesc["target",A2], arraydesc["AB_ID",A2]),
                        sub=paste("pad-slide-incubationRun-spottingRun:", A2))
            
            result <- rbind(result,	data.frame(dilseries_id=A1, pad_slide_incubRun_spotRun=A2, intercept=intercept, 
                                               errIntercept=errIntercept, ab=arraydesc["AB_ID",A2], 
                                               slide=arraydesc["slide",A2]))
          }
        }
      dev.off()
      return(result)
    }
     ## getSignals ##
    ## 
    ## Purpose: Update the timeseries signal to (FG - intercept).
    ## 
    ## Input: timeseries = timeseries raw data
    ##        intercepts = intercepts from getIntercepts()
    ##        arraydesc = element "arraydescription" of RPPanalyzers data list
    ##        exportNo = variable for analyzeIntercepts ("export")
    ##
    ## Output: A new data.frame which is the original timeseries with updated signal intensity
    
    getSignals <- function(timeseries, intercepts, arraydesc, exportNo) {
      
      dilSeriesID<-NULL
      anaIntercept <- analyzeIntercepts(intercepts, test="F", export=exportNo)
      fit <- attr(anaIntercept, "fit")
      result <- timeseries
      for(A1 in unique(timeseries$dilSeriesID)) {
        subseries1 <- subset(timeseries, dilSeriesID==A1)
        for(A2 in unique(colnames(subseries1)[grep(colnames(arraydesc)[1],colnames(timeseries)):ncol(timeseries)])) {
            newdata <- data.frame(ab=arraydesc["AB_ID",A2], slide=arraydesc["slide",A2], 
                                  dilseries_id=A1)
            intercept <- as.numeric(predict(fit, newdata=newdata))
            signals <- as.my(subseries1[,A2])-intercept
            result[which(result$dilSeriesID==A1),A2]<-signals
        }
      }
      return(result)	
    }
 
    ## analyzeIntercepts ##
    ## 
    ## Purpose: Analysis of Variances for nested models.
    ## 
    ## Input: intercepts = result from getIntercepts()
    ##        export = integer which of the linear fits should be exported to the attribute of the result
    ##        test = test parameter for ANOVA (see documentation "anova")
    ##
    ## Output: anova object, one of the fits (in attribute "fit"), pdf with anova RSS barplot
    
    analyzeIntercepts <- function(intercepts, test="F", export) {
      fit <- list()
      errIntercept <- intercepts[,"errIntercept"]
      fit[[1]] <- lm(intercept ~ 1, data=intercepts, weights=1/errIntercept^2)
      fit[[2]] <- lm(intercept ~ ab, data=intercepts, weights=1/errIntercept^2)
      fit[[3]] <- lm(intercept ~ ab + slide, data=intercepts, weights=1/errIntercept^2)
      if(length(unique(intercepts$dilseries_id))>1){
        fit[[4]] <- lm(intercept ~ ab + slide + dilseries_id, data=intercepts, weights=1/errIntercept^2)
        ano <- anova(fit[[1]], fit[[2]], fit[[3]], fit[[4]], test=test)
        rss <- ano$RSS
        names(rss) <- c("const.", "+ ab", "+ slide", "+ dilseries_id")
      }else{
        ano <- anova(fit[[1]], fit[[2]], fit[[3]], test=test)
        rss <- ano$RSS
        names(rss) <- c("const.", "+ ab", "+ slide")
      }
      pdf("anovaIntercepts_Output.pdf")
      barplot(rss, ylab="weighted RSS")
      dev.off()
      
      attr(ano, "fit") <- fit[[export]]
      
      return(ano)
    }
 


