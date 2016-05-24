f.var.covar <- function(pred, X, data, info){
##
## COMPUTES EXPLICITLY THE VARIANCE-COVARIANCE MATRIX FOR THE LIKELIHOOD WITH MISSING DATA
##
## X IS THE DESIGN MATRIX (COMPLETE GRID), pred ARE PREDICTED FREQUENCIES FROM LAST RUN OF THE EM,
## data ARE THE "ORIGINAL" DATA USED BY EM
#
##
#
## STANDARD POISSON CONTRIBUTION:
.pX <- pred * X
.d2Poisson <- -t(X) %*% (.pX) # THIS IS STANDARD POISSON WITH VALUES PREDICTED ACCORDING TO EM
#
## MIDLERTIDIG SJEKK (SER UT TIL AA VAERE UNODVENDIG): #-#
.sjekk <- tapply(data$orig.lines, data$ind, unique)
.sjekk0 <- table(.sjekk)
.sjekk1 <- table(tapply(data$ind, data$orig.lines, unique))
if(any(.sjekk0 != 1) | any(.sjekk1 != 1)) stop() # SJEKKER AT ind OG orig.lines ER EN-TIL-EN, BURDE IKKE VAERE NODVENDIG
#
## MATCH EVERYTHING TO data
.pos <- f.pos.match(data = data, info = info)
.X <- X[.pos,]
.m <- dim(.X)[1]
.k <- dim(.X)[2]
.l <- pred[.pos]
.colnames <- colnames(.X)
#	
## NORMALIZE PREDICTIONS OVER AMBIGUITY GROUPS:
.orig.lines <- data$orig.lines
.lstar <- .l/f.groupsum(.l, INDICES = .orig.lines) # NORMALIZED
.Xl <- .lstar * .X
#	
## COMPUTE INDIVIDUAL SCORE PARTS:	
.lTX <- tapply(as.numeric(.Xl), list(orig.lines = rep(.orig.lines, .k), col = rep(1:.k, each = .m)), sum) ## THIS IS (gamma_i)^T X, WITH ONE ROW FOR EACH i AND ONE COLUMN FOR EACH COLUMN OF X. EACH ROW IS THE INDIVIDUAL PART OF THE SCORE VECTORS
colnames(.lTX) <- .colnames
#
## COMPUTE FIRST AND SECOND PART OF AMBIGUITY-PART OF SECOND DERIV. MULTINOMIAL LOGLIKE
.d2phi.1 <- t(.lTX) %*% .lTX
.d2phi.2 <- t(.X) %*% (.Xl)
.d2phi.del1 <- - .d2phi.1 + .d2phi.2
#	
## ADD STANDARD POISSON AND PART DUE TO AMBIGUITIES:
.d2.loglike.Poisson <- .d2phi.del1 + .d2Poisson
#
## SCORE COMPUTATION (COMPUTED FROM MULTINOMIAL FORMULA t(gamma_i - gamma) %*% X
## COMMON ELEMENT FOR ALL SCORE VECTORS:
.gammaTX <- (t(.pX/sum(pred)) %*% rep(1, dim(.pX)[1]))[,1] ## LAST SUBSETTING REDUCES MATRIX TO NAMED VECTOR
.score <- t(t(.lTX) - .gammaTX) # SUBTRACT .gammaTX COLUMNWISE
#
## INVERT TO OBTAIN VAR-COVAR:	
.var.covar.Poisson <- -solve(.d2.loglike.Poisson)
#
##
return(list(var.covar = .var.covar.Poisson, score = .score))
}
