StepSize <- function(W, psi, dX, stepDirection){

# binary search: to = 2^-M for smallest positive M s.t. H(to) >= 0
#
# OUTPUT
# t : optimal step length in (0,1] in direction "stepDirection"
#     via Hermite interpolation for convex functions
# 
# function 'HelpFunk' from Matlab implementation directly included here

# search for optimal to
to <- 1
Hto <- LikFunk(W, psi + to * stepDirection, dX) - LikFunk(W, psi, dX)
while (Hto < 0){
    to <- to / 2
    Hto <- LikFunk(W, psi + to * stepDirection, dX) - LikFunk(W, psi, dX)
    }

# step length
dH0 <- t(GradientL(W, psi, dX)) %*% stepDirection
dH0to <- dH0 * to
if (Hto >= dH0to / 2){t <- to} else {t <- to * dH0to / (2 * (dH0to + Hto))}

return(as.numeric(t))
}
