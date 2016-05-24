data(bmw)
qplot(bmw)
meplot(abs(bmw[bmw > 0]))
meplot(abs(bmw[bmw < 0]))
out <- gpd(-bmw, ne = 200)
plot(out)
tp <- tailplot(out)
gpd.q(tp, 0.999)
gpd.sfall(tp, 0.999)

exindex(bmw, 30)
exindex(-bmw, 30)
out <- pot(-bmw, ne = 200, run = 30)
tp <- plot(out)
if(!is.null(tp)) {
    gpd.q(tp, 0.999)
    gpd.sfall(tp, 0.999)
}
