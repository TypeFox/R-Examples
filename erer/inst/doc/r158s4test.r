# A. Test the new function for linear model: lm2
library(erer); data(daIns); source("C:/aErer/r157s4.r")
ny <- lm2(data = daIns, name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
slotNames(ny); getSlots("lm2")
ny@sigma; slot(object = ny, name = "sigma")

show(ny); ny; plot(ny); class(ny)
showMethods(classes = "lm2") 
getMethod(f = "show", signature = "lm2") 
getMethod(f = "plot", signature = "lm2")

# Two error trials
lm2(data = daIns, name.y = "Y", name.x = c("Injury", "gender")) # Gender
lm2(data = daIns, name.y = 45,  name.x = c("Injury", "Gender")) # 45

# B. Comparison with 'lm' as implemented in the base R by S3
gg <- lm(formula = Y ~ 1 + Injury + HuntYrs + Nonres + Lspman + Lnong +
  Gender + Age + Race + Marital + Edu + Inc + TownPop, data = daIns)
names(gg); summary(gg); plot(residuals(gg))
summary(gg)$r.squared; summary(gg)$sigma