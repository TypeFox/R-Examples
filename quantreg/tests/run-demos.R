library(quantreg)

(dDIR <- system.file("demo", package = "quantreg"))
setwd(dDIR)
set.seed(1) # since some demos randomly generate

cat("Running demos from package 'quantreg' : \n\n")
for(ff in list.files(dDIR, pattern="\\.R$", full.names = TRUE)) {
   f <- basename(ff)
   cat("\n", f," :\n", paste(rep.int("-", nchar(f)), collapse=''),
       "\n", sep='')

   source(ff, echo = TRUE)
}

