library("RUnit")
library("prc")

test.prcsp <- function() {


RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6


# semiparametric fit
t1=Sys.time()
fit=prcsp (mtct.eg$V3_BioV3B_2500[1:50], 2500, mtct.eg$V3_BioV3B_500[1:50], 500, try.additiona.support.sets=FALSE, verbose=2, max.iter=2, grid.density=50, stop.when.dropping=FALSE)
t2=Sys.time()
print(t2-t1)
if (!is.null(fit)) {
    checkEqualsNumeric(coef(fit), c(19.0012791, 28334.5818187,    -2.0655433,     0.4956502 ), tolerance=tolerance)
    checkEqualsNumeric(mixlik(fit), 85.83998, tolerance=tolerance)
}

# structural model fit
t1=Sys.time()
fit=prcstruct (mtct.eg$V3_BioV3B_2500[1:50], 2500, mtct.eg$V3_BioV3B_500[1:50], 500, verbose=2, max.iter=1, grid.density=100)
t2=Sys.time()
print(t2-t1)
if (!is.null(fit)) {
    checkEqualsNumeric(coef(fit), c(19.0012791, 28420.1881214,    -2.0440495,     0.5005353), tolerance=tolerance)
}


}
