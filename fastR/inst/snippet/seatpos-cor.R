require(car) # for scatterplot.matrix
scatterplot.matrix(~Age+Arm+hipcenter+Ht+HtShoes+Leg+Seated+Thigh+Weight, 
    reg.line=lm, smooth=TRUE, span=0.5, 
    diagonal = 'density', data=seatpos) 
round(cor(seatpos),2)
