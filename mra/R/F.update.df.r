F.update.df <- function( fit, df=NA ){
#
#   Function to update the df of a model
#

if( is.na(df) ){
    df <- df$df.estimated   # use the rank estimated by MRAWIN, stored in fit after version 1.12 of MRA
} else if( df <= 0 ){
    df <- fit$aux$nx + fit$aux$ny  # assume full rank
} # else use the values supplied by user (unchanged from input)
fit$df <- df

# recompute fit statistics
fit$aic <- -2*fit$loglik + 2*df
fit$qaic <- -2*(fit$loglik/fit$vif) + 2*df
fit$aicc <- fit$aic + (2*df*(df+1))/(fit$aux$nan - df - 1)
fit$qaicc <- fit$qaic + (2*df*(df+1))/(fit$aux$nan - df - 1)

fit

}
