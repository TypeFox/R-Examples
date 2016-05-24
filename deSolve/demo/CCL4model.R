### Functions to facilitate fitting the CCl4 inhalation model

initparms <- function(...) {
  arglist <- list(...)
  Pm <- numeric(36)
  ## The Changeable parameters are ones that can be modified on input
  Changeable <- c("BW", "QP", "QC", "VFC", "VLC", "VMC", "QFC", "QLC", "QMC",
                  "PLA", "PFA", "PMA", "PTA", "PB", "MW", "VMAX", "KM",
                  "CONC", "KL", "RATS", "VCHC")
  ## Computed parameters are strictly functions of the Changeable ones.
  Computed <- c("VCH", "AI0", "PL", "PF", "PT", "PM", "VTC", "VT", "VF",
                "VL", "VM", "QF", "QL", "QM", "QT")
  names(Pm) <- c(Changeable,
                 Computed
                 )

### Physiological parameters
  Pm["BW"] <- 0.182                     # Body weight (kg)
  Pm["QP"] <- 4.0                  # Alveolar ventilation rate (hr^-1)
  Pm["QC"] <- 4.0                       # Cardiac output (hr^-1)
  Pm["VFC"] <- 0.08                 # Fraction fat tissue (kg/(kg/BW))
  Pm["VLC"] <- 0.04               # Fraction liver tissue (kg/(kg/BW))
  Pm["VMC"] <- 0.74           # Fraction of muscle tissue (kg/(kg/BW))
  Pm["QFC"] <- 0.05         # Fractional blood flow to fat ((hr^-1)/QC
  Pm["QLC"] <- 0.15      # Fractional blood flow to liver ((hr^-1)/QC)
  Pm["QMC"] <- 0.32     # Fractional blood flow to muscle ((hr^-1)/QC)

  ## Chemical specific parameters for chemical
  Pm["PLA"] <- 16.17                 # Liver/air partition coefficient
  Pm["PFA"] <- 281.48                  # Fat/air partition coefficient
  Pm["PMA"] <- 13.3                 # Muscle/air partition coefficient
  Pm["PTA"] <- 16.17               # Viscera/air partition coefficient
  Pm["PB"] <- 5.487                  # Blood/air partition coefficient
  Pm["MW"] <- 153.8                     # Molecular weight (g/mol)
  Pm["VMAX"] <- 0.11          # Maximum velocity of metabolism (mg/hr)
  Pm["KM"] <- 1.3                   # Michaelis-Menten constant (mg/l)

  ## Parameters for simulated experiment
  Pm["CONC"] <- 1000                    # Inhaled concentration
  Pm["KL"] <- 0.02                  # Loss rate from empty chamber /hr
  Pm["RATS"] <- 1.0               # Number of rats enclosed in chamber
  Pm["VCHC"] <- 3.8                     # Volume of closed chamber (l)

  ## Now, change anything from the argument list
  ## First, delete anything in arglist that is not in Changeable
  whichdel <- which(! names(arglist) %in% Changeable)
  if (length(whichdel)) {
    warning(paste("Parameters", paste(names(arglist)[whichdel], collapse=", "),
                  "are not in this model\n"))
  }
  arglist[whichdel] <- NULL
  ## Is there anything else
  if (length(arglist)) {
    Pm[names(arglist)] <- as.vector(unlist(arglist))
  }
    
  
  ## Computed parameter values

  Pm["VCH"] <- Pm["VCHC"] - Pm["RATS"]*Pm["BW"] # Net chamber volume
  Pm["AI0"] <- Pm["CONC"]*Pm["VCH"]*Pm["MW"]/24450 # Initial amt. in chamber (mg)
  Pm[c("PL", "PF", "PT", "PM")] <- Pm[c("PLA", "PFA", "PTA", "PMA")]/Pm["PB"]

  ## Fraction viscera (kg/(kg BW))
  Pm["VTC"] <- 0.91 - sum(Pm[c("VLC", "VFC", "VMC")])

  Pm[c("VT", "VF", "VL", "VM")] <- Pm[c("VTC", "VFC", "VLC", "VMC")]*Pm["BW"]

  Pm[c("QF", "QL", "QM")] <- Pm[c("QFC", "QLC", "QMC")]*Pm["QC"]

  Pm["QT"] <- Pm["QC"] - sum(Pm[c("QF", "QL", "QM")])
  Pm
}

### We don't actually use these functions (though they work)
### They exist because cclmodel.orig is easier to read than ccl4modelG
### The model function also computes some values that are of interest in
### checking the model and for calculating a dose metric:
### the amount metabolized (AM)
### the area under the concentration-time curve in the liver (CLT)
### and the mass balance (MASS), which should be constant if everything
### worked right.

## State variable, y, assignments.

## CI   CM  CT  CF  CL     
## AI  AAM  AT  AF  AL  CLT  AM
## 1    2   3   4   5    6   7

initstate.orig <- function(Pm) {
  y <- rep(0, 7)
  names(y) <- c("AI", "AAM", "AT", "AF", "AL", "CLT", "AM")
  y["AI"] <- Pm["AI0"]
  y
}

parms <- initparms()

ccl4model.orig <- with(as.list(parms), function(t, y, parms) {
  conc <- y[c("AI", "AAM", "AT", "AF", "AL")]/c(VCH, VM, VT, VF, VL)
  ## Vconc[1] is conc in mixed venous blood
  Vconc <- c(0, conc[2:5]/parms[c("PM", "PT", "PF", "PL")]) # '0' is a placeholder
  Vconc[1] <- sum(Vconc[2:5]*c(QM, QT, QF, QL))/QC
  ## CA is conc in arterial blood
  CA <- (QC * Vconc[1] + QP * conc[1])/
    (QC + QP/PB)

  ## Exhaled chemical
  CX <- CA/PB
  ## return the derivatives and other computed items
  list(c(RATS*QP*(CX - conc[1]) - KL*y["AI"],
         QM*(CA - Vconc[2]),
         QT*(CA - Vconc[3]),
         QF*(CA - Vconc[4]),
         QL*(CA - Vconc[5]) -
           (RAM <- VMAX*Vconc[5]/(KM + Vconc[5])),
         conc[5],
         RAM),
       c(DOSE = as.vector(AI0 - y["AI"]),
         MASS = as.vector(sum(y[c("AAM","AT", "AF", "AL", "AM")])*RATS),
         CP=as.vector(conc[1]*24450.0/MW)
       ))
})

### Versions that only calculate what is needed for parameter estimation


initparmmx <- function(parms) {
  mx <- matrix(nrow=5, ncol=7)
  mx[1, 6] <- parms["VCH"]
  mx[1, 7] <- parms["MW"]
  mx[4, 6] <- parms["VL"]*parms["PL"]
  mx[5, 6] <- parms["VMAX"]
  mx[5, 7] <- parms["KM"]

  mxx <- matrix(parms[c("QP", "QM", "QT", "QF", "QL")], nrow=5, ncol=5, byrow=TRUE)
  mxx <- sweep(mxx, 2, parms[c("VCH", "VM", "VT", "VF", "VL")], "/")
  mxx <- sweep(mxx, 2, c(1, parms[c("PM", "PT", "PF", "PL")]), "/")

  mxx <- mxx/(parms["QC"] + parms["QP"]/parms["PB"])

  mxx <- sweep(mxx, 1, c(parms["RATS"]*parms["QP"]/parms["PB"],
                         parms[c("QM", "QT", "QF", "QL")]), "*")
  dg <- diag(c(parms["RATS"]*parms["QP"]/parms["VCH"] + parms["KL"],
               parms[c("QM", "QT", "QF", "QL")]/
               (parms[c("PM", "PT", "PF", "PL")]*parms[c("VM", "VT", "VF", "VL")])))
  mxx <- mxx - dg
  mx[1:5, 1:5] <- mxx
  mx

}
  

### Now, include the gradients wrt Vmax, Km, and initial chamber concentration

initstateG <- function(Pm) {
  y <- rep(0, 20)
  names(y) <- c("AI", "AAM", "AT", "AF", "AL",
                "dAIdVm", "dAAMdVm", "dATdVm", "dAFdVm", "dALdVm",
                "dAIdK", "dAAMdK", "dATdK", "dAFdK", "dALdK",
                "dAIdy0", "dAAMdy0", "dATdy0", "dAFdy0", "dALdy0"
                )
  y["AI"] <- Pm["AI0"]
  y["dAIdy0"] <- Pm["VCH"] * Pm["MW"]/24450.0
  y
}

ccl4modelG <- function(t, y, parms) {
  list(c(parms[,1:5] %*% y[1:5] - c(0, 0, 0, 0, parms[5, 6]*y[5] /
                                     ((Kms <- parms[5, 7]*parms[4, 6]) + y[5])),
         parms[, 1:5] %*% y[6:10] - c(0, 0, 0, 0,
                                     y[5]/(Kms + y[5]) +
                                     parms[5, 6]*Kms*y[10]/
                                     (Kms + y[5])^2),
         parms[, 1:5] %*% y[11:15] - c(0, 0, 0, 0,
                                      parms[5, 6]*(y[15]*Kms -
                                                  parms[4, 6]*y[5])/
                                      (Kms + y[5])^2),
         parms[,1:5] %*% y[16:20] - c(0, 0, 0, 0,
                                      parms[5, 6]*Kms*y[20]/(Kms + y[5])^2)
         ),
       c(CP = as.vector(y[1]*(zz <- 24450.0/parms[1, 6]/parms[1, 7])),
         dCPdVm = as.vector(y[6]*zz),
         dCPdK = as.vector(y[11]*zz),
         dCPdy0 = as.vector(y[16]*zz)
         )
       )
}

### Function to use in gnls.  This is more complicated than usual for such
### functions, because each value for each animal depends on the previous
### value for that animal.  Normal vectorization doesn't work.  Work with
### log(Vmax) and log(Km)

ccl4gnls <- function(time, initconc, lVmax, lKm, lconc) {
  Vmax <- if(length(lVmax) == 1) rep(exp(lVmax), length(time)) else exp(lVmax)
  Km <- if (length(lKm) == 1) rep(exp(lKm), length(time)) else exp(lKm)
  conc <- if (length(lconc) == 1) rep(exp(lconc), length(time)) else exp(lconc)
  Concs <- levels(initconc)
  CP <- numeric(length(time))
  .grad <- matrix(nrow=length(time), ncol=3,
                  dimnames=list(NULL, c("lVmax", "lKm", "lconc")))
  ### Run the model once for each unique initial concentration
  for (Conc in Concs) {
    sel <- initconc == Conc

    parms <- initparms(CONC=conc[sel][1], VMAX=Vmax[sel][1], KM=Km[sel][1])
    parmmx <- initparmmx(parms)

    y <- initstateG(parms)

    TTime <- sort(unique(time[sel]))
    if (! 0 %in% TTime) TTime <- c(0, TTime)

    out <- lsoda(y, TTime, ccl4modelG, parmmx, rtol=1e-12, atol=1e-12)
    CP[sel] <- out[match(time[sel], out[,"time"]),"CP"]
    .grad[sel, "lVmax"] <- out[match(time[sel], out[, "time"]), "dCPdVm"]
    .grad[sel, "lKm"] <- out[match(time[sel], out[, "time"]), "dCPdK"]
    .grad[sel, "lconc"] <- out[match(time[sel], out[, "time"]), "dCPdy0"]
  }
  
  .grad <- .grad * cbind(Vmax, Km, conc)
  attr(CP, "gradient") <- .grad
  CP
}

if (require(nlme, quietly=TRUE)) {

  start <- log(c(lVmax = 0.11, lKm=1.3, 25, 100, 250, 1000))

### Data are from:
### Evans, et al. (1994) Applications of sensitivity analysis to a
### physiologically
### based pharmacokinetic model for carbon tetrachloride in rats.
### Toxicology and Applied Pharmacology 128: 36--44.

  data(ccl4data)
  ccl4data.avg<-aggregate(ccl4data$ChamberConc,
                        by=ccl4data[c("time", "initconc")], mean)
  names(ccl4data.avg)[3]<-"ChamberConc"

### Estimate log(Vmax), log(Km), and the logs of the initial
### concentrations with gnls
  cat("\nThis may take a little while ... \n")
  ccl4.gnls <- gnls(ChamberConc ~ ccl4gnls(time, factor(initconc),
                                           lVmax, lKm, lconc),
                    params = list(lVmax + lKm ~ 1, lconc ~ factor(initconc)-1),
                    data=ccl4data.avg, start=start,
                    weights=varPower(fixed=1),
                    verbose=TRUE)

  start <- coef(ccl4.gnls)
  ccl4.gnls2 <- gnls(ChamberConc ~ ccl4gnls(time, factor(initconc),
                                            lVmax, lKm, lconc),
                     params = list(lVmax + lKm ~ 1,
                       lconc ~ factor(initconc)-1),
                     data=ccl4data, start=start,
                     weights=varPower(fixed=1),
                     verbose=TRUE)

  print(summary(ccl4.gnls2))

  ### Now fit a separate initial concentration for each animal
  start <- c(coef(ccl4.gnls))
  cat("\nApprox. 95% Confidence Intervals for Metabolic Parameters:\n")
  tmp <- exp(intervals(ccl4.gnls2)[[1]][1:2,])
  row.names(tmp) <- c("Vmax", "Km")
  print(tmp)
  cat("\nOf course, the statistical model is inappropriate, since\nthe concentrations within animal are pretty highly autocorrelated:\nsee the graph.\n")

  opar <- par(ask=TRUE, no.readonly=TRUE)
  plot(ChamberConc ~ time, data=ccl4data, xlab="Time (hours)",
       xlim=range(c(0, ccl4data$time)),
       ylab="Chamber Concentration (ppm)", log="y")

  out <- predict(ccl4.gnls2, newdata=ccl4data.avg)

  concentrations <- sort(unique(ccl4data$initconc))
  for (conc in concentrations) {
    times <- ccl4data.avg$time[sel <- ccl4data.avg$initconc == conc]
    CP <- out[sel]
    lines(CP ~ times)
  }
  par(opar)
} else {
  cat("This example requires the package nlme\n")
}

