make.parvec <-
function(emul, fix.betas) {

# Check whether betas are to be fixed #!+
if (fix.betas) emul$beta.vec=NULL
  
# Parameter vector #!+
parvec <- c(emul$rho, emul$kappa, emul$zeta, emul$beta.vec, emul$phi.vec)

# Names for the parameter vector #!+
namevec   <- c("rho", "kappa", "zeta") 
if (!is.null(emul$beta.vec)) { 
   for (i in 1:length(emul$beta.vec)) namevec <- c(namevec, "beta") 
}
for (i in 1:length(emul$phi.vec)) namevec  <- c(namevec, "phi") 
names(parvec) <- namevec 

# Output #!+
parvec
}
