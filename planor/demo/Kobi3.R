library("planor")
#---------------------------------------------------------------------------
# EXAMPLES FROM THE PLANOR MANUAL
#---------------------------------------------------------------------------
# Exemple 3 page 16
#---------------------------------------------------------------------------
cat("\n")
cat("***************** EXEMPLE 3 PAGE 16 *****************\n")
cat("\n")
cat("Five treatment factors at 2 levels \n")
cat("robotFac \n")
cat("Block structure: 4 plates / (2 rows x 4 columns) \n")
cat("Model: nsoil + qsoil + Rug + conc + Tact \n")
cat("Estimate: idem \n")
cat("Hierarchy 1: conc constant in plate \n")
cat("Hierarchy 2: Tact constant in plate \n")
cat("Hierarchy 3: nsoil constant in plate x column \n")
cat("Hierarchy 4: qsoil constant in plate x column \n")
cat("N=32 units\n")
cat("\n")
# 
cat("*** RUN ***\n")
robotFac <- planor.factors( factors=list(
                        conc=c(1,3),
                        Tact=c(15,30),
                        nsoil=c("curd","Saint-Paulin"),
                        qsoil=c("0.01g","0.10g"),
                        Rug=c(0.25,0.73),
                        plate=1:4,
                        row=1:2,
                        col=1:4),
                     hierarchy=list(~conc/plate,
                                    ~Tact/plate,
                                    ~nsoil/(plate*col),
                                    ~qsoil/(plate*col)))
# Obligation for the moment of a complicated model formulation:
# - all factors in the design must appear in the model-estimate formulae
#   (it is not sufficient to declare them as basic factors)
robotMod <- planor.model( model=~nsoil * qsoil * Rug * conc * Tact,
                    listofmodels=list(c(~plate+row+col+Rug,~Rug)) )
robotKey <- planor.designkey(factors=robotFac, model=robotMod,
                             nunits=32,
                       base=~plate+row+col, max.sol=1)
robotDes <- planor.design(robotKey)
print(robotDes@design)
summary(robotKey, maxprint=31)
# --------------------------------------------------


cat("\n ROBOT1B \n")
cat("Block structure: 2 plates / (2 rows x 4 columns) \n")
cat("Model: nsoil + qsoil + Rug + conc \n")
cat("Estimate: idem \n")
cat("Hierarchy 1: conc constant in plate \n")
cat("Hierarchy 2: Tact constant in plate \n")
cat("Hierarchy 3: nsoil constant in plate x column \n")
cat("Hierarchy 4: qsoil constant in plate x column \n")
cat("N=16 units\n")
cat("\n")
# 
cat("*** RUN ***\n")
F3B <- planor.factors( factors=list(
                        conc=c(1,3),
                        Tact=c(15,30),
                        nsoil=c("curd","Saint-Paulin"),
                        qsoil=c("0.01g","0.10g"),
                        Rug=c(0.25,0.73),
                        plate=1:4,
                        row=1:2,
                        col=1:4),
                     hierarchy=list(~conc/plate,
                                    ~Tact/plate,
                                    ~nsoil/(plate*col),
                                    ~qsoil/(plate*col)))
# Obligation for the moment of a complicated model formulation:
# - all factors in the design must appear in the model-estimate formulae
#   (it is not sufficient to declare them as basic factors)
M3B <- planor.model( model=~nsoil * qsoil * Rug * conc * Tact,
                    listofmodels=list(c(~plate+row+col,~Rug)) )
K3B <- planor.designkey(factors=F3B, model=M3B, nunits=32,
                       base=~plate+row+col, max.sol=1)
P3B <- planor.design(key=K3B, select=c(1,1))
resum3B <-summary(K3B[1,1], maxprint=31)
