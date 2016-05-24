setClass("stls", representation(call="call",coefficients="matrix",startcoef="matrix",value="numeric",counts="integer",convergence="integer",message="character",residuals="matrix",fitted.values="matrix",df.residual="integer",covariance="matrix",bootrepl="matrix"))

setClass("summary.stls", representation(level="numeric",confint="matrix", bootconfint="matrix"), contains="stls")
