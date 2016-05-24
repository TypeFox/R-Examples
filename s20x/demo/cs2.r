data (body.df)
ethbod<-crosstabs(~ethnicity+factor(bodyim,levels=c("slight.uw","right","slight.ow","mod.ow","very.ow")),data=body.df)
ethbod$exp
fisher.test(body.df$ethnicity,body.df$bodyim,simulate.p.value=T,B=1e5,alternative="greater")
rowdistr(ethbod)
rowdistr(ethbod,comp="between")

