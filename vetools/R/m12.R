# Verified 1.3.18
m12 <-
function(x) {
        r = x %% 12
        if (r == 0) { return(12) }
        return(r)
}
