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
  #shift power to first column
  tbl <- tbl[,c(2,1,3:ncol(tbl))]
  return(invisible(tbl))
} 

# Hauschke D., Steinijans V. and Pigeot I.
# "Bioequivalence studies in Drug Development"
# John Wiley & Sons, Chichester (2007)
cat(paste('Hauschke D., Steinijans V. and Pigeot I.\n',
'"Bioequivalence studies in Drug Development"\n',
'John Wiley & Sons, Chichester (2007)\n\n',sep=""))
# Table 5.1 in
# Sample size for 2x2 design, multiplicative model
# theta1=0.8, theta2=1.25(1/theta1)
CVs  <- seq(from=0.1,  to=0.4, by=0.025)
GMRs <- seq(from=0.85, to=1.2, by=0.05)
power <- c(0.8, 0.9)
cat("Table 5.1 2x2, BEL 0.8-1.25, multiplicative model, exact\n")
t5.1 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.8,
                method="exact", design="2x2")
# You can compare the result to the data object ct5.1
# for validation purposes
            
# Table 5.2
# Sample size for 2x2 design, multiplicative model
# theta1=0.75, theta2=1.3333(1/theta1)
CVs  <- seq(from=0.15,  to=0.45, by=0.025)
GMRs <- seq(from=0.80, to=1.25, by=0.05)
power <- c(0.8, 0.9)
cat("Table 5.2 2x2, BEL 0.75-1.3333, multiplicative model, exact\n")
t5.2 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.75,
                method="exact", design="2x2")

# Table 5.3
# Sample size for 2x2 design, multiplicative model
# theta1=0.9, theta2=1.1111 (1/theta1)
CVs   <- seq(from=0.05,  to=0.25, by=0.025)
GMRs  <- seq(from=0.925, to=1.075, by=0.025)
power <- c(0.8, 0.9)
cat("Table 5.3 2x2, BEL 0.9-1.1111, multiplicative model, exact\n")
t5.3 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=TRUE, theta1=0.9,
                method="exact", design="2x2")

cat(paste('Chow S.C., Liu J.P.\n',
          '"Design and Analysis of Bioavailability\n',
          'and Bioequvalence Studies", Third edition\n',
          'CRC Press, Chapman & Hall, Boca Raton (2009)\n\n',sep=""))
cat("Table 5.4.1 2x2, BEL 0.8-1.20, additive model, exact\n")
CVs   <- seq(from=0.1,  to=0.4, by=0.02)
GMRs  <- seq(from=0.0, to=0.15, by=0.05)
power <- c(0.8, 0.9)
t5.4.1 <- sampsiz(power=power, CV=CVs, GMR=GMRs, logscale=FALSE, theta1=-0.2,
                  method="exact", design="2x2")
