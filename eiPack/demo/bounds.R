
data(senc)
out <- bounds(cbind(dem, rep, non) ~ cbind(black, white, natam), data
              = senc, rows = c("black", "white"), column = "dem",
              excluded = "non", threshold = 0.9)
par(ask = TRUE)
plot(out, row = "black", column = "dem")
plot(out, row = "white", column = "dem")
par(ask = FALSE)
