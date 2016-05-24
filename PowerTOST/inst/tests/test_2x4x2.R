# Sample size tables for Baalam's design
# Author: dlabes
# -----------------------------------------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz <- function(alpha=0.05, power, CV, GMR, theta1, logscale=TRUE, 
                    design="2x2", method="exact")
{
# 'cartesian' product
  data <- merge(CV, GMR)
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

# Chow and Liu use a model with carry-over, i.e. a degree of freedom
# which is one less than in PowerTOST

cat(paste('Chow S.C., Liu J.P.\n',
          '"Design and Analysis of Bioavailability\n',
          'and Bioequvalence Studies", Third edition\n',
          'CRC Press, Chapman & Hall, Boca Raton (2009)\n\n',sep=""))
cat("Table 9.6.1 Baalams, BEL 0.8-1.20, additive model,\n")
cat("approximate via shifted central t-distribution\n")
CVs   <- seq(from=0.1,  to=0.4, by=0.02)
GMRs  <- seq(from=0.0, to=0.15, by=0.05)
power <- c(0.8, 0.9)
t9.6.1 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=FALSE, theta1=-0.2,
                  method="shifted", design="2x4x2")

cat("Table 9.6.5 Baalams, BEL 0.8-1.25, multiplicative model model,\n");
cat("approximate via shifted central t-distribution\n")
cat("Attention! Liu, Chow's CV is se!\n")
se    <- seq(from=0.1,  to=0.4, by=0.02)
CVs   <- se2CV(se)
GMRs  <- seq(from=0.85, to=1.2, by=0.05)
power <- c(0.8, 0.9)
t9.6.5 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.8,
                  method="shift", design="2x4x2")

