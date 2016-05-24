`summary.modTempEff` <-
function(object, spar=TRUE, digits = max(3, getOption("digits") - 3), ...){
      if(length(object$ToTheat)<=0) {
        stop("the model does not include csdl coefficients")}
      n<-length(object$fitted.values)
      dev<-object$dev
      tot.edf<-sum(object$edf)
      bic<-object$aic-2*tot.edf+log(n)*tot.edf
      ubre<- (object$dev + 2*tot.edf -n)/n
      edf.cold.tot<-sum(object$edf.cold)
      edf.heat.tot<-sum(object$edf.heat)
      tf <- terms(object$formula, specials = c("csdl","seas","dl"))
      id.csdl<-attr(tf,"specials")$csdl 
      nomeTemp<-all.vars(object$formula)[id.csdl]
      is.ridge<-"omega.Cold"%in%rownames(object$sp.mio)|"omega.Heat"%in%rownames(object$sp.mio)
    	if(is.null(object$edf.seas)) { seasP<-NA 
          } else {
          seasP<-paste("edf =", round(sum(object$edf.seas),2), "(rank =",object$rank.seas, "; log(spar) =", 
              paste(round(log(object$sp.mio["lambda.seas",]),3),")",sep=""))
          object$sp.mio<-object$sp.mio[-match("lambda.seas",rownames(object$sp.mio),0),]
          }
      #print(object$call)
      cat("Model:\t")
      print(object$call) #formula is better??
#      cat("Temperature:", nomeTemp, "\t")
#	   if(!is.ridge) {cat("ridge penalty: ", NA, "\n")} else {
#      cat("ridge penalty: ", "cold=", ridgeC, "heat=", ridgeH,"\n")
#          }
      cat("\nSeasonality (smooth): ", seasP,"\n")
      cat("\nFit summary", paste("(model edf = ",round(sum(object$edf),2), "; n = ",n,"):",sep=""),"\n")
      cat("AIC =",object$aic," BIC =",bic," ubre =",round(ubre,5)," dev =" ,object$dev, "\n")
      xx<-matrix(,2,6)
      rownames(xx)<-c("Cold","Heat")
      colnames(xx)<-c("Est","SE.freq","SE.bayes","rank","edf","L")
      xx["Cold",1:2]<-object$ToTcold[1:2]; xx["Cold",3]<-object$ToTcold.bayes[2]
      xx["Heat",1:2]<-object$ToTheat[1:2]; xx["Heat",3]<-object$ToTheat.bayes[2]
      xx["Cold","rank"]<-object$rank.cold
      xx["Heat","rank"]<-object$rank.heat
      xx["Cold","edf"]<-edf.cold.tot; xx["Heat","edf"]<-edf.heat.tot
      xx["Cold","L"]<-length(object$betaCold)-1
      xx["Heat","L"]<-length(object$betaHeat)-1
#      cat("\nNet effects of ", paste(nomeTemp,":",sep=""), "\n")
      cat("\nNet effects of ", paste(nomeTemp, " (based on edf = ", round(edf.cold.tot+edf.heat.tot+
        length(object$delta),2), "):",sep=""), "\n")
      print(xx,digits=digits)
      if(spar){
      cat("\nlog(spar) for smooth DL curves: \n")
      print(drop(t(log(object$sp.mio))),digits=digits)
        }
      if(length(object$delta)<=0){
        cat("\nFIXED Threshold at:",object$psi,"\n")
        } else {
      cat("\nThreshold: \n")
      psi<-object$psi
      rownames(psi)<-if(nrow(object$psi)==1) {""} else {c("psi1","psi2")}
      print(psi,digits=3)
      cat("V variable(s):\n")
      print(object$Tdelta)
        }
      }

