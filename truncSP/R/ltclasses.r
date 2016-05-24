setClass("lt", representation(call="call",coefficients="matrix",startcoef="matrix",cvalues="data.frame",value="numeric",counts="integer",convergence="integer",message="character",residuals="matrix",fitted.values="matrix",df.residual="integer",covariance="matrix",bootrepl="matrix"))

setClass("summary.lt", representation(level="numeric",confint="matrix", bootconfint="matrix"), contains="lt")
