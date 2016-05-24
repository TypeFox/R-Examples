"powerF" <-
function(PV, df2, df1=1, alpha=.05)
 {
ncp <- df2 * ( PV / (1-PV))
power <- 1- pf(qf(1-alpha, df1, df2), df1, df2, ncp)
return(power)
}

