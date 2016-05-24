f.pos.match <- function(data, info){
##
## FOR data, COMPUTES THE POSITION OF EACH LINE IN A COMPLETE
## GRID AS USED IN THE GLM
#
## EXTRACT CHARACTERISTICS FROM THE DESIGN MATRIX
.char <- f.make.design(info = info, ret.characteristics = T)
#
##
## MATCH PREDICTED PROBABILITIES TO ORIGINAL DATA:
.pos <- f.pos.in.grid(A = .char, comb = as.matrix(data[, names(.char)]))
#
##
return(.pos)
}
