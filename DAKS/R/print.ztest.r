print.ztest<-function(x, ...){
if(length(x$diff_value) == 1){cat("\n \t One sample Z-test\n")}
if(length(x$diff_value) == 2){cat("\n \t Two sample Z-test\n")}
cat("\nz = ", round(x$Z.value, digits = 4), " ")
cat("p-value = ", round(x$p.value, digits = 4), " ")

if(x$alternative == "two.sided"){
if(is.null(x$imp_alt)){
cat("\nalternative hypothesis: true mean is not equal", x$mu, " ")
}else{
cat("\nalternative hypothesis: true difference in means is not equal", x$mu, " ")	
}
}
if(x$alternative == "greater"){
if(is.null(x$imp_alt)){
cat("\nalternative hypothesis: true mean is greater", x$mu, " ")
}else{
cat("\nalternative hypothesis: true difference in means is greater", x$mu, " ")	
}
}
if(x$alternative == "less"){
if(is.null(x$imp_alt)){
cat("\nalternative hypothesis: true mean is less", x$mu, " ")
}else{
cat("\nalternative hypothesis: true difference in means is less", x$mu, " ")	
}
}

cat("\n")
cat(x$conf.level*100, "percent confidence interval:\n", " ")
write(x$conf, "")

cat("sample estimates:\n")	
if(length(x$diff_value) == 1){
estimate <- round(x$diff_value, digits = 5)
names(estimate) <- "mean in imp"
print(estimate)
}
if(length(x$diff_value) == 2){
estimate <- round(x$diff_value, digits = 5)
names(estimate) <- c("mean in imp", "mean in imp_alt")
print(estimate)
}
cat("\n")
}
