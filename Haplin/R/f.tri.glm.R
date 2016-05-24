"f.tri.glm" <- function(observed, design.matrix, maternal = F, info, ...)
{
#
# THE PROGRAM ESTIMATES EFFECTS OF SEVERAL ALLELES IN A CASE-TRIAD,
# CASE-CONTROL-TRIAD OR CASE-CONTROL DESIGN
# ASSUMING HARDY-WEINBERG EQUILIBRIUM (AND RARE DISEASE, IF NECESSARY)
# IMPORTANT: OBSERVED FREQUENCIES (observed) MUST COME FROM A COMPLETE
# GRID IN APPROPRIATE ORDER! 
#
# THE DESIGN MATRIX IS COMPUTED BY f.make.design
#
# THE ... PASSES ON OTHER ARGUMENTS, LIKE start, TO THE GLM
#
## NOTE: ARGUMENTS SUCH AS maternal ETC SHOULD BE EXTRACTED FROM info 
## IN THE FUTURE. BUT: NBNB! DIFFERENT RUNS ARE OFTEN DONE WITH DIFFERENT
## SETTINGS FOR response (AND POSSIBLY ALSO maternal FOR LIKELIHOOD RATIO TEST)
## THIS IS NOT REFLECTED IN info
#
# OBSERVED FREQUENCY, ARRANGED ACCORDING TO GRID:
.o <- observed	#
.n.haplo <- sum(info$haplos$selected.haplotypes) ## BRUKES BARE I OUTPUT
#
## COMPUTE THE CORRECT DESIGN MATRIX (COULD HAVE BEEN DONE ONCE AND FOR ALL, BUT REMEMBER THAT response, AT LEAST, CAN BE DIFFERENT FROM WHAT'S IN info)
.design.matrix <- design.matrix
.d1 <- nrow(.design.matrix)
#
## IN FIRST STEP, INITIALIZE. -1 SIGNALS THIS
if(identical(.o, -1)){
	.o <- rep(1, .d1)
}
#
## JUST AN EXTRA CHECK THAT DESIGN MATCHES OBSERVED FREQUENCY DATA:
if(length(.o) != .d1) stop("Problem with design matrix!")
#
## NUMBER OF TRIADS (OR INDIVIDUALS, FOR CASE-CONTROL)
.ntri <- sum(.o)
#
## ADD OBSERVED FREQUENCY DATA TO DESIGN MATRIX
.design.matrix <- cbind(.o, .design.matrix)
#
## CONSTRUCT GENERAL FORMULA
.formula <- formula(paste(c(".o ~ -1 ", names(.design.matrix)[-1]), collapse = "+"))
#
## ESTIMATION: ##
#
## ACTUAL ESTIMATION (SUPPRESSES WARNINGS SINCE FREQUENCIES MAY BE NON-INTEGER)
.res <- suppressWarnings(glm(.formula, family = poisson, data = .design.matrix, ..., maxit = 20))
###.res1 <- suppressWarnings(glm.fit(x = as.matrix(.design.matrix[,-1]), y = .o, family = poisson(), intercept = F))
# 
## FREQUENCY PREDICTION (NOTE: COULD ALSO HAVE USED FITTED VALUES IN OBJECT): 
.pred <- predict(.res, type = "response")
if(abs(sum(.pred) - .ntri) > 0.001 * .ntri) stop("Potential problem in prediction!")
#
## ADDING INFORMATION TO OUTPUT:
.out <- list(result = .res, pred = .pred, nall = .n.haplo, ntri = .ntri, ref.cat = info$haplos$ref.cat, maternal = maternal, design = info$model$design, orig.call = sys.call(), date = date())
class(.out) <- "tri.glm"
return(invisible(.out))
}
