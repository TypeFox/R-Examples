"dhhrss" <-
function (dhhrobj) 
{
v.hft <- dhhrobj$v.hft
s.hft <- dhhrobj$s.hft
k.hft <- dhhrobj$k.hft
cat("Skewness calculations\n")
cat("_____________________\n")
cat("Summary s.hft:\n")
print(summary(s.hft))
cat("Kurtosis calculations\n")
cat("_____________________\n")
cat("Summary k.hft:\n")
print(summary(k.hft))
cat("Variance calculations\n")
cat("_____________________\n")
mv <- mean(v.hft)
cat("Variance of scale corrected HFT is", var(v.hft/mv), "\n")
}

