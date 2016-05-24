print.gains <- function(x,digits=2,...) 
{  
  if (x$percents==FALSE) {

   if ( (x$conf %in% c("normal","t","boot")) & x$optimal==TRUE)
   {
     cat("Depth                            Cume   Cume Pct                 Optimal Optimal     Mean       CI for   \n")
     cat(" of           Cume     Mean      Mean   of Total    Lift   Cume    Lift    Cume     Model      Mean Resp \n")
     cat("File     N      N      Resp      Resp      Resp    Index   Lift   Index    Lift     Score     Lower Upper\n")
     cat("---------------------------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(x$mean.resp[i],width=10,digits=digits,format="f"),
           formatC(x$cume.mean.resp[i],width=10,digits=digits,format="f"),
           formatC(100*x$cume.pct.of.total[i],width=10,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$opt.lift[i],width=8,digits=0,format="f"),
           formatC(x$opt.cume.lift[i],width=8,digits=0,format="f"),
           formatC(x$mean.prediction[i],width=10,digits=digits,format="f"),
           formatC(x$conf.lower[i],format="f",width=10,digits=digits),
           formatC(x$conf.upper[i],format="f",width=6,digits=digits),"\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")) & x$optimal!=TRUE)
   {
     cat("Depth                            Cume   Cume Pct                     Mean       CI for   \n")
     cat(" of           Cume     Mean      Mean   of Total    Lift   Cume     Model      Mean Resp \n")
     cat("File     N      N      Resp      Resp      Resp    Index   Lift     Score     Lower Upper\n")
     cat("-----------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(x$mean.resp[i],width=10,digits=digits,format="f"),
           formatC(x$cume.mean.resp[i],width=10,digits=digits,format="f"),
           formatC(100*x$cume.pct.of.total[i],width=10,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$mean.prediction[i],width=10,digits=digits,format="f"),
           formatC(x$conf.lower[i],format="f",width=10,digits=digits),
           formatC(x$conf.upper[i],format="f",width=6,digits=digits),"\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")==FALSE) & x$optimal==TRUE)
   {
     cat("Depth                            Cume   Cume Pct                 Optimal Optimal     Mean\n")
     cat(" of           Cume     Mean      Mean   of Total    Lift   Cume    Lift    Cume     Model\n")
     cat("File     N      N      Resp      Resp      Resp    Index   Lift   Index    Lift     Score\n")
     cat("-----------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(x$mean.resp[i],width=10,digits=digits,format="f"),
           formatC(x$cume.mean.resp[i],width=10,digits=digits,format="f"),
           formatC(100*x$cume.pct.of.total[i],width=10,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$opt.lift[i],width=8,digits=0,format="f"),
           formatC(x$opt.cume.lift[i],width=8,digits=0,format="f"),
           formatC(x$mean.prediction[i],width=10,digits=digits,format="f"),"\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")==FALSE) & x$optimal!=TRUE)
   {
     cat("Depth                            Cume   Cume Pct                     Mean\n")
     cat(" of           Cume     Mean      Mean   of Total    Lift   Cume     Model\n")
     cat("File     N      N      Resp      Resp      Resp    Index   Lift     Score\n")
     cat("-------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(x$mean.resp[i],width=10,digits=digits,format="f"),
           formatC(x$cume.mean.resp[i],width=10,digits=digits,format="f"),
           formatC(100*x$cume.pct.of.total[i],width=10,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$mean.prediction[i],width=10,digits=digits,format="f"),"\n",sep="")
     }
   }
  } else {  #x$percents==TRUE

   if ( (x$conf %in% c("normal","t","boot")) & x$optimal==TRUE)
   {
     cat("Depth                            Cume   Cume Pct                 Optimal Optimal     Mean        CI for     \n")
     cat(" of           Cume     Resp      Resp   of Total    Lift   Cume    Lift    Cume     Model       Resp Rate   \n")
     cat("File     N      N      Rate      Rate      Resp    Index   Lift   Index    Lift     Score     Lower   Upper \n")
     cat("------------------------------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(100*x$mean.resp[i],width=10,digits=digits,format="f"),"%",
           formatC(100*x$cume.mean.resp[i],width=9,digits=digits,format="f"),"%",
           formatC(100*x$cume.pct.of.total[i],width=9,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$opt.lift[i],width=8,digits=0,format="f"),
           formatC(x$opt.cume.lift[i],width=8,digits=0,format="f"),
           formatC(100*x$mean.prediction[i],width=10,digits=digits,format="f"),"%",
           formatC(100*max(x$conf.lower[i],0),format="f",width=8,digits=digits),"%",
           formatC(100*min(x$conf.upper[i],1),format="f",width=8,digits=digits),"%\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")) & x$optimal!=TRUE)
   {
     cat("Depth                            Cume   Cume Pct                     Mean        CI for     \n")
     cat(" of           Cume     Resp      Resp   of Total    Lift   Cume     Model       Resp Rate   \n")
     cat("File     N      N      Rate      Rate      Resp    Index   Lift     Score     Lower   Upper \n")
     cat("--------------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(100*x$mean.resp[i],width=10,digits=digits,format="f"),"%",
           formatC(100*x$cume.mean.resp[i],width=9,digits=digits,format="f"),"%",
           formatC(100*x$cume.pct.of.total[i],width=9,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(100*x$mean.prediction[i],width=10,digits=digits,format="f"),"%",
           formatC(100*max(x$conf.lower[i],0),format="f",width=8,digits=digits),"%",
           formatC(100*min(x$conf.upper[i],1),format="f",width=8,digits=digits),"%\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")==FALSE) & x$optimal==TRUE)
   {
     cat("Depth                            Cume   Cume Pct                 Optimal Optimal     Mean \n")
     cat(" of           Cume     Resp      Resp   of Total    Lift   Cume    Lift    Cume     Model \n")
     cat("File     N      N      Rate      Rate      Resp    Index   Lift   Index    Lift     Score \n")
     cat("------------------------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(100*x$mean.resp[i],width=10,digits=digits,format="f"),"%",
           formatC(100*x$cume.mean.resp[i],width=9,digits=digits,format="f"),"%",
           formatC(100*x$cume.pct.of.total[i],width=9,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(x$opt.lift[i],width=8,digits=0,format="f"),
           formatC(x$opt.cume.lift[i],width=8,digits=0,format="f"),
           formatC(100*x$mean.prediction[i],width=10,digits=digits,format="f"),"%\n",sep="")
     }
   } else if ( (x$conf %in% c("normal","t","boot")==FALSE) & x$optimal!=TRUE)
   {
     cat("Depth                            Cume   Cume Pct                     Mean \n")
     cat(" of           Cume     Resp      Resp   of Total    Lift   Cume     Model \n")
     cat("File     N      N      Rate      Rate      Resp    Index   Lift     Score \n")
     cat("--------------------------------------------------------------------------\n")
     for (i in 1:x$num.groups) {
       cat(format(x$depth[i],width=4),format(x$obs[i],width=6),format(x$cume.obs[i],width=7),
           formatC(100*x$mean.resp[i],width=10,digits=digits,format="f"),"%",
           formatC(100*x$cume.mean.resp[i],width=9,digits=digits,format="f"),"%",
           formatC(100*x$cume.pct.of.total[i],width=9,digits=1,format="f"),"%",
           formatC(x$lift[i],width=8,digits=0,format="f"),
           formatC(x$cume.lift[i],width=7,digits=0,format="f"),
           formatC(100*x$mean.prediction[i],width=10,digits=digits,format="f"),"%\n",sep="")
     }
   }
  }
}
