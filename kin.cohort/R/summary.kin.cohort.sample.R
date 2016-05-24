summary.kin.cohort.sample<- function(object,...){

if (!inherits(object,"kin.cohort.sample"))
  stop("Please provide a family data.frame created by kc.simul")

cat("Simulated population for kin-cohort analysis\n\n")

#c("famid", "rel", "age", "gender", "agecancer", "cancer", "carrier", "real.carrier"

probands <- object$rel==0
npro <- sum(probands)
cancer   <- object$cancer[probands]==1
carrier  <- object$carrier[probands]==1


cat("Probands: ", npro,"\n")
cat("Affected probands: ",sum(cancer), "(",formatC(sum(cancer )/npro*100,1,format="f"),"%)\n")
cat("Carrier  probands: ",sum(carrier), "(",formatC(sum(carrier)/npro*100,1,format="f"),"%)\n")
#
# depends Table
cat("\nAffected probands by carrier status\n")

T<-table(carrier,cancer)
P<-formatC(prop.table(T)*100,1,format="f")

cat("Noncarriers: ", T[1,2], "(",P[1,2]," %)\n")
cat("Carriers:    ", T[2,2], "(",P[2,2]," %)\n")

nrel <- sum(!probands)
cancer   <- object$cancer[!probands]==1
carrier  <- object$real.carrier[!probands]==1

cat("\n\nRelatives: ", nrel,"\n")
cat("Affected relatives: ",sum(cancer), "(",formatC(sum(cancer )/nrel*100,1,format="f"),"%)\n")
cat("Carrier  relatives: ",sum(carrier), "(",formatC(sum(carrier)/nrel*100,1,format="f"),"%) \n")
#
# depends Table
cat("\nAffected relatives by carrier status\n")

T<-table(carrier,cancer)
P<-formatC(prop.table(T)*100,1,format="f")

cat("Noncarriers: ", T[1,2], "(",P[1,2]," %)\n")
cat("Carriers:    ", T[2,2], "(",P[2,2]," %)\n")
}


