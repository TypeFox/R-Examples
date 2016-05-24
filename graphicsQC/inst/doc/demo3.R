plot(pressure^.125 ~ temperature, data=pressure)
plot(lm(pressure^.125 ~ temperature, data=pressure))
