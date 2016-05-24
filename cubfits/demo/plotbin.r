suppressMessages(library(cubfits, quietly = TRUE))

reu13.df <- ex.train$reu13.df
phi.Obs <- ex.train$phi.Obs
b.Init <- b.Init$roc
aa.names <- names(reu13.df)

# summary data.
ret.bin <- prop.bin.roc(reu13.df, phi.Obs)
phi.Obs.lim <- range(phi.Obs)
ret.model <- prop.model.roc(b.Init, phi.Obs.lim = phi.Obs.lim)

# plot.
par(mfrow = c(1, 3))
for(i.aa in 1:length(aa.names)){
  plotbin(ret.bin[[i.aa]], ret.model = ret.model[[i.aa]], main = aa.names[i.aa])
}
