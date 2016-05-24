data(carcass)

car <- as.gmData(carcass)

m1 <- ggm(~.^.,gmData=car)

m1 <- fit(m1)

#dynamic.Graph(m1)

