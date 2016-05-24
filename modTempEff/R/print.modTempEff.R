`print.modTempEff` <-
function(x, digits = max(3, getOption("digits") - 3), ...){
      n<-length(x$fitted.values)
      dev<-x$dev
      tot.edf<-sum(x$edf)
      bic<-x$aic-2*tot.edf+log(n)*tot.edf
      ubre<- (x$dev + 2*tot.edf -n)/n
      edf.cold.tot<-sum(x$edf.cold)
      edf.heat.tot<-sum(x$edf.heat)
      cat("Model Summary", paste("(n = ", n,"):",sep=""), "\n")
      cat("AIC =",x$aic,"  BIC =",bic,"  ubre =",round(ubre,5),"  dev =" ,x$dev, "\n")
      if(length(x$ToTheat)<=0) {
        return(invisible(x))
        }
      cat("Degrees of freedom:\n")
      xx<-matrix(,2,4)
      colnames(xx)<-c("Model", "Cold", "Heat", "Seasonality")
      rownames(xx)<-c("edf","rank")
      xx[1,1]<-tot.edf; xx[2,1]<-length(x$coef)
      xx["edf","Cold"]<-edf.cold.tot; xx["rank","Cold"]<- x$rank.cold
      xx["edf","Heat"]<-edf.heat.tot; xx["rank","Heat"]<- x$rank.heat
      if(!is.null(x$edf.seas)){
        xx["edf","Seasonality"]<-sum(x$edf.seas)
        xx["rank","Seasonality"]<-x$rank.seas}
      #xx[2,]<-as.integer(xx[2,])
      print(xx,digits=digits,...)
      }

