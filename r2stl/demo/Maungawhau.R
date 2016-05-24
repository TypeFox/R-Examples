# Now let's look at R's Volcano data

z <- volcano

x <- 1:dim(volcano)[1]

y <- 1:dim(volcano)[2]

r2stl(x, y, z, filename="volcano.stl", show.persp=TRUE)
