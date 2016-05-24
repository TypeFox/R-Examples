quantile.manual.thresh.scalewise<-function(x, which.levels=c(1,2,3),
                                            hard = TRUE,quantile=.9, ...)
{
wc.shrink <- x
if (hard) {
for (i in names(x)[which.levels]) {
wci <- x[[i]]
unithresh <- quantile(abs(wci),quantile)
wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
}
}
else {
for (i in names(x)[which.levels]) {
wci <- x[[i]]
unithresh <- quantile(abs(wci),quantile)
wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) *
(abs(wci) > unithresh)
}
}
wc.shrink
}
