data(coagulation,package="faraway")
summary(coag~diet,data=coagulation,fun=favstats)
coag.xyplot <- xyplot(coag~diet,coagulation)
coag.bwplot <- bwplot(coag~diet,coagulation)
