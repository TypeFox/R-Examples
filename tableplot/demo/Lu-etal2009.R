## Lu-etal2009   Summarizing results from a large simulation study

description <- "
Results from a simulation study by Lu et al. (2009) comparing five methods for
modeling latent variable relationships. The manipulated conditions of the
simulation are: method of modeling (M); number of indicators or items (I) per
latent variable; loading size (L) of the items; size of parameters or betas (B)
in the latent variable relationships; and sample size (S). 

The simulation results are examined in terms of five criteria: bias in parameter
estimates; coverage of confidence intervals for parameters (i.e., proportion of
confidence intervals that capture the true parameter values); mean square error
from the analysis; power for significance testing on parameter values; and
standard error of estimates. 

The results are summarized in terms of the effect sizes (partial eta-squares)
for all terms in the five univariate linear models for the simulations.
"

partial.eta2 <- read.table(yy <- textConnection("
Bias Cover MSE Power SE
(M)ethod     0.977    0.826 0.849 0.299 0.850
(I)tem       0.935    0.481 0.764 0.617 0.570
(L)oading    0.737    0.441 0.865 0.619 0.746
(B)eta       0.870    0.881 0.414 0.992 0.601
(S)ampleSize 0.594    0.716 0.596 0.925 0.922
MxI          0.936    0.531 0.772 0.055 0.518
MxL          0.871    0.438 0.633 0.166 0.547
MxB          0.956    0.887 0.599 0.156 0.820
MxS          0.782    0.742 0.337 0.136 0.475
IxL          0.557    0.012 0.579 0.480 0.323
IxB          0.511    0.547 0.588 0.656 0.048
IxS          0.040    0.035 0.131 0.083 0.214
LxB          0.269    0.463 0.112 0.633 0.033
LxS          0.115    0.004 0.206 0.096 0.385
BxS          0.480    0.611 0.407 0.907 0.268
MxIxL        0.658    0.038 0.470 0.072 0.355
MxIxB        0.780    0.560 0.853 0.121 0.156
MxIxS        0.407    0.035 0.606 0.024 0.239
MxLxB        0.354    0.535 0.239 0.123 0.283
MxLxS        0.262    0.053 0.218 0.088 0.203
MxBxS        0.704    0.665 0.632 0.078 0.520
IxLxB        0.163    0.492 0.250 0.601 0.005
IxLxS        0.054    0.005 0.042 0.243 0.121
IxBxS        0.050    0.463 0.290 0.720 0.015
LxBxS        0.021    0.452 0.125 0.699 0.003
MxIxLxB      0.513    0.497 0.594 0.106 0.059
MxIxLxS      0.232    0.005 0.332 0.038 0.163
MxIxBxS      0.485    0.497 0.622 0.080 0.074
MxLxBxS      0.138    0.501 0.276 0.086 0.064
IxLxBxS      0.082    0.278 0.077 0.472 0.001
MxIxLxBxS    0.245    0.318 0.288 0.065 0.028
Overall      0.999    0.997 0.996 0.999 0.987"), header=T)
close(yy)


specs <- make.specs(scale.max = 100, back.fill="grey80", label=1, label.size=0.5, shape.col="grey50")

pdf(file="Lu2009.pdf", paper="letter")

tableplot(
	values = round(as.matrix(partial.eta2)*100),
	cell.specs = specs,
	assign.sets = matrix(1,32,5),
	h.parts = c(5,10,10,5,1,1), 
	v.parts=c(1,1,1,1,1),
	label.size = 0.6, 
	gap=1.5,
	left.space=22)

dev.off()