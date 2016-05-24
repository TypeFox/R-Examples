# Author: dlabes
# -----------------------------------------------------------------------------

library(PowerTOST)
# function to create sample size tables
sampsiz2 <- function(alpha=0.05, power, CV, CVb, GMR, theta1, design="2x2")
{
# 'cartesian' product
  data1 <- merge(CV, CVb, stringsAsFactors=FALSE)
  data  <- merge(CV, GMR, stringsAsFactors=FALSE)
  data  <- merge(data, data1, by="x", stringsAsFactors=FALSE)
  names(data) <- c("CV","GMR","CVb")
  data <- data[ ,c(1,3,2)]
  data <- data[order(data$CV,data$CVb,data$GMR),]
  tbl  <- data.frame()
  for (j in seq_along(power))
  {
    data$n <- 1
    data$power <- power[j]
    for (i in seq_along(data$n)) {
      data$n[i] <- sampleN.RatioF(alpha=alpha, CV=data[i,"CV"], 
                                  CVb=data[i,"CVb"], 
                                  theta1=theta1, theta0=data$GMR[i], 
                                  targetpower=power[j], design=design, 
                                  print=FALSE)[,"Sample size"]
    }
    #print(data)
    data2 <- reshape(data, v.names="n", idvar=c("power","CV","CVb"), timevar="GMR", 
                     direction="wide")
    names(data2) <- gsub("n.","R",names(data2))
    names(data2)[names(data2)=="R1"] <- "R1.0"
    names(data2)[names(data2)=="R0"] <- "R0.0"
    cat("Power",power[j],"\n")
    print(data2[,-3],row.names=FALSE)
    cat("\n")
    tbl <- rbind(tbl, data2)
  }
  return(invisible(tbl))
} 

# Hauschke D., Steinijans V. and Pigeot I.
# "Bioequivalence studies in Drug Development"
# John Wiley & Sons, Chichester (2007)
cat(paste('Hauschke D., Steinijans V. and Pigeot I.\n',
'"Bioequivalence studies in Drug Development"\n',
'John Wiley & Sons, Chichester (2007)\n\n',sep=""))

CVs  <- seq(from=0.1,  to=0.4, by=0.05)
CVbs <- ""
GMRs <- seq(from=0.85, to=1.2, by=0.05)
power <- c(0.8, 0.9)
cat("Table 10.1 parallel, BEL 0.8-1.25, multiplicative model\n")
t10.1 <- sampsiz2(alpha=0.025, power=power, CV=CVs, CVb=CVbs, GMR=GMRs, 
                  theta1=0.8, design="parallel")

CVs  <- seq(from=0.1,  to=0.4, by=0.05)
CVbs <- seq(from=0.2,  to=1, by=0.1)
GMRs <- seq(from=0.85, to=1.2, by=0.05)
power <- c(0.8)
cat("Table 10.3a 2x2, BEL 0.8-1.25, multiplicative model\n")
t10.3a <- sampsiz2(alpha=0.025, power=power, CV=CVs, CVb=CVbs, GMR=GMRs, 
                   theta1=0.8, design="2x2")
