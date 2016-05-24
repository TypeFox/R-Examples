# spearman.r
# Time-stamp: c:/x/rpack/corrgram/tests/spearman.r

dat <- data.frame(
a=c(4.427371, 4.652635, 5.117894, 4.648367),
b=c(3.065796, 2.652594, 2.702682, 2.664237),
c=c(1.820023, 1.248203, 1.583516, 0.943176),
d=c(1.0981531, 1.3303935, 1.4942022, 0.9183384),
e=c(2.445383, 2.230344, 2.796065, 2.108716),
f=c(-3.058478, -3.287118, -3.063749, -2.809525),
g=c(3.512448, 3.196442, 3.450073, 3.190861),
h=c(2.16616, 2.601121, 2.175354, 2.87997),
i=c(3.223686, 2.84551, 2.959626, 2.845965),
j=c(2.298442, 1.562345, 1.897684, 1.541668),
k=c(-1.481442, -2.59257, -1.751966, -1.667766),
l=c(-3.321619, -3.422294, -3.501768, -3.175093),
m=c(1.898439, 1.331196, 2.141605, 1.928749),
n=c(2.0098, 1.713265, 2.110812, 2.506027),
o=c(3.180492, 3.101148, 3.488329, 3.329727),
p=c(2.397239, 2.091732, 2.521272, 2.510106),
q=c(1.3390112, 0.8297595, 1.1517592, 1.0869653)
  )

cmat <- cor(dat, method='spearman')
# The cmat matrix has a boundary case in which one of the correlations
# is exactly -1 - .Machine$double.eps
print(min(cmat), digits=17) # -1.0000000000000002
is.matrix(cmat) # TRUE
isSymmetric(cmat) # TRUE
min(cmat, na.rm = TRUE) > -1 - .Machine$double.eps # FALSE
min(cmat, na.rm = TRUE) >= -1 - .Machine$double.eps # TRUE
max(cmat, na.rm = TRUE) < 1 + .Machine$double.eps # TRUE

corrgram(cmat, type="cor", order=FALSE,
         lower.panel=panel.shade, upper.panel=panel.pts,
         text.panel=panel.txt, diag.panel=panel.minmax, main="Data")

# ----------------------------------------------------------------------------
