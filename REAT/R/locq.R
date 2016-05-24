locq <-
function (e_ij, e_j, e_i, e) {

if (e_ij > e_j) { return (NA) }
if (e_i > e) { return (NA) }
if (e_j > e) { return (NA) }

s_ij <- e_ij/e_i
s_j <- e_j/e  
LQ <- s_ij/s_j

return(LQ)
}
