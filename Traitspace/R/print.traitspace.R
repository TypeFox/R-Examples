print.traitspace <-
function(x,...){

print(x[8])
print(x[9])

cat("Predicted relative abundances:", "\n")
print(data.frame(x$predicted.p$P_level_1_level_3))
cat("Predicted species distributions:", "\n")
print(data.frame(x$predicted.p$P_level_1_level_3_dist))
cat("\n")

#cat("Observed relative abundances:", "\n")
#print(data.frame(x$true.p$p))
#cat("Observed species distributions:", "\n")
#print(data.frame(x$true.p$p_dist))
#cat("** warning:", "the observed relative abundances/species distributions are calculated from the trait data. This may not be the correct value.", "\n")
}
