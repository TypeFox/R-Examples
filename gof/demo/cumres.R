example(cumres.glm)
par(mfrow=c(2,2));
plot(g1)
plot(g1, ci="sim", legend="type2")
plot(g1, col="darkblue", legend=NULL)
plot(g1, col=NULL, ci="sim", col.ci="darkred")

opt <- par(ask=TRUE)
par(mfrow=c(1,1)); plot(g, idx=3, col=NULL, ci=TRUE, legend="type2")
par(opt)
