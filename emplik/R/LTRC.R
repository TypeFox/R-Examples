LTRC <- function(x,d,w=rep(1, length(d)),y=rep(-Inf, length(x)) ) {

#### This works for Left Truncated and Right Censored data.
#### this is actually Lynden-Bell or WJT estimator.
#### Or the Kaplan-Meier/Nelson Aalen est. for right censor only data.

temp <- Wdataclean2(x,d,w)
dd <- temp$dd
ww <- temp$weight
dd[length(dd)] <- 1
xx <- temp$value

######why not use DnR?

temp <- DnR(xx,dd,ww,y=y)

NelAal <- temp$n.event/temp$n.risk
survP <- cumprod( 1 - NelAal )
NelAal <- cumsum(NelAal)
jumps <- -diff( c(1, survP) )

list(times=xx[dd==1], survjump=jumps, surv=survP, CumHaz=NelAal)
}
