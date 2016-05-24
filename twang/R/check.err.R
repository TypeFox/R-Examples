check.err<-function(cov.table, stage, alerts.stack, estimand, ess.ctrl, ess.treat)
{

if(estimand == "ATT")
{
   ind  <- (cov.table$tx.sd < .0001) | (cov.table$std.ef.sz > 500)
   prob <- cov.table$std.eff.sz[ind]
   if(length(prob)>0)
   {
      sink(alerts.stack, append=TRUE)
      cat("\n problematic standard deviations in stage ",stage,"\n\n")
      print(cov.table[which(ind),c("tx.sd","std.eff.sz")])
      cat("\n\n\n")
      sink()
   }
}


if(estimand == "ATE")
{
   sd.p = ((ess.treat * cov.table$tx.sd) + (ess.ctrl * cov.table$ct.sd))/(ess.treat + ess.ctrl)
   ind  <- (sd.p < .0001) | (cov.table$std.ef.sz > 500)
   prob <- cov.table$std.eff.sz[ind]
   if(length(prob)>0)
   {
      sink(alerts.stack, append=TRUE)
      cat("\n problematic standard deviations in stage ",stage,"\n\n")
      print(cov.table[which(ind),c("sd.p","std.eff.sz")])
      cat("\n\n\n")
      sink()
   }
}

}

