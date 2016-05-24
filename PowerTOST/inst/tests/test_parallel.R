# Author: dlabes
# -----------------------------------------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz <- function(alpha=0.05, power, CV, GMR, theta1, logscale=TRUE, 
                    design="2x2", method="exact")
{
# 'cartesian' product
  data <- merge(CVs, GMRs)
  names(data) <- c("CV","GMR")
  tbl  <- data.frame()
  for (j in seq_along(power))
  {
    data$n <- 1
    data$power <- power[j]
    for (i in seq_along(data$n)) {
      data$n[i] <- sampleN.TOST(alpha=alpha, CV=data[i,"CV"], 
                                theta0=data$GMR[i], logscale=logscale,
                                targetpower=power[j], theta1=theta1,
                                design=design, method=method, 
                                print=FALSE)[,"Sample size"]
    }
    data2 <- reshape(data, v.names="n", idvar=c("power","CV"), timevar="GMR", 
                     direction="wide")
    names(data2) <- gsub("n.","R",names(data2))
    names(data2)[names(data2)=="R1"] <- "R1.0"
    names(data2)[names(data2)=="R0"] <- "R0.0"
    cat("Power",power[j],"\n")
    print(data2[,-2],row.names=FALSE)
    cat("\n")
    tbl <- rbind(tbl, data2)
  }
  return(invisible(tbl))
}
cat(paste('S.A.Julious\n',
'"Tutorial in Biostatistics\n',
'Sample sizes for clinical trials with Normal data"\n',
'Statistics in Medicine, Vol. 23, 1921-1986 (2004)\n\n',sep=""))
CVs  <- seq(from=0.3,  to=0.85, by=0.05)
GMRs <- seq(from=0.95, to=1.05,  by=0.05)
power <- c(0.9)
cat("Table VIII parallel group, BEL 0.8-1.25, multiplicative, nct\n")
cat("Column level of BE =10%. PowerTOST gives ntotal!\n")
tSJ.VIII.10 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.9,
                       method="noncent", design="parallel")

CVs  <- seq(from=0.3,  to=0.85, by=0.05)
GMRs <- seq(from=0.85, to=1.2,  by=0.05)
power <- c(0.9)
cat("Table VIII parallel group, BEL 0.8-1.25, multiplicative, nct\n")
cat("Column level of BE =20%. PowerTOST gives ntotal!\n")
tSJ.VIII.20 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.8,
                       method="noncent", design="parallel")
                   
cat(paste('Chow and Wang\n',
'"On Sample Size Calculation in Bioequivalence Trials"\n',
'J. Pharmacokin. Biopharm. Vol. 28(2), 155-169 (2001)\n\n',sep=""))
CVs  <- seq(from=0.1,  to=0.4, by=0.02)
CVs  <- se2CV(CVs)
#GMRs <- seq(from=0, to=0.15, by=0.05)
#GMRs <- exp(GMRs)
GMRs <- c(1, 1.05, 1.11, 1.16)
power <- c(0.8, 0.9)
cat("Table I parallel group, BEL 0.8-1.25, log-transf. data, exact\n")
cat("Their s% translates into CV=se2CV(s%).\n")
cat("Their theta' is log(theta0), i.e. translates to theta0=exp(theta')!\n")
tCW.I <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.8,
                 method="exact", design="parallel")

CVs  <- seq(from=0.1,  to=0.4, by=0.02)
GMRs <- seq(from=0.0, to=0.15, by=0.05)
power <- c(0.8, 0.9)
cat("Table III parallel group, BEL 0.8-1.20, raw data, exact\n")
tCW.III <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=FALSE, theta1=-0.2,
                   method="exact", design="parallel")


