library(erer); data(daIns); names(daIns)

# Manual creation of a formula object
f1 <- Y ~ Injury + HuntYrs + Nonres
class(f1); typeof(f1); str(f1); f1

# Converting a character string into a formula
st <- paste("Y", 
  paste(c("Injury", "HuntYrs", "Nonres"), collapse = " + "), sep = " ~ ")
f2 <- as.formula(st)
f3 <- formula(st)
st; class(st); f3

# Working on many variable names
k1 <- Y ~ .
k2 <- formula(daIns[, -14])
k3 <- bsFormu(name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
k1; k2; k3

# Extracting elements from a formula object
f1[1]; f1[2]; f1[3]
tf <- terms(f1); tf
attributes(tf)
attr(x = tf, which = "term.labels") 