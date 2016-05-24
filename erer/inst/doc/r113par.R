# Parameters in par()
?par; # example(par)
p72 <- par(no.readonly = FALSE)             # 72 total parameters
p66 <- par(no.readonly = TRUE )             # 66 settable parameters
p06 <- p72[!(names(p72) %in% names(p66))]   # 6 for queries and read-only
str(p72)
p72[1:3]
sort(noquote(names(p72)))
sort(noquote(names(p06)))

# Three uses of par()
# a. querying graphics state
par(c("pch", "col")) 

# b. setting parameters
par(pch = 3, col = "red")

# c. Saving, setting, and retoring
oldPar <- par(mfrow = c(3, 2), col = "green")  # save and set 
oldPar                      # show old parameter values
par(c("mfrow", "col"))      # show new parameter values
for (i in 1:6) {
  plot(rnorm(100), pch = i) # draw graphs with new parameter values
}
par(oldPar)                 # restore old parameter values
par(c("mfrow", "col"))      # show old parameter values again