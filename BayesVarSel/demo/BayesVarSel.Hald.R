#read Hald data
data(Hald)
#run the main function: (in this small example we keep all models)
hald.Bvs<- Bvs(formula="y~x1+x2+x3+x4", data=Hald, n.keep=16)

#print the result
hald.Bvs

#obtain a summary of the result
summary(hald.Bvs)

#plots for the result
#image plot of the conditional probabilities
plotBvs(hald.Bvs, option="conditional")


