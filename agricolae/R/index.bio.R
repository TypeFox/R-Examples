`index.bio` <-
function (data, method = c("Margalef","Simpson.Dom", "Simpson.Div", "Berger.Parker",
"McIntosh","Shannon"),level = 95, nboot = 100, console=TRUE)
{
    method <- match.arg(method)
    x <- data
    if (length(x) > 1) {
    if (method == "Margalef") 
        formula1 <- expression((length(x) - 1)/log(sum(x)))
    if (method == "Simpson.Dom")
        formula1 <- expression(sum((x/sum(x))^2))
    if (method == "Simpson.Div") 
        formula1 <- expression(1 - sum((x/sum(x))^2))
    if (method == "Berger.Parker") 
        formula1 <- expression(max(x)/sum(x))
    if (method == "McIntosh") 
        formula1 <- expression((sum(x) - sqrt(sum(x^2)))/(sum(x) - 
            sqrt(sum(x))))
    if (method == "Shannon") 
        formula1 <- expression(-sum((x/sum(x)) * log((x/sum(x))^2)))
    index <- eval(formula1)
    n = length(x)
    estimador <- rep(0, n)
    for (i in 1:nboot) {
        x <- sample(data, n, replace = TRUE)
        estimador[i] <- eval(formula1)
    }
    lic = as.numeric(quantile(estimador, (1 - 0.01 * level)/2, type = 6))
    lsc = as.numeric(quantile(estimador, (1 + 0.01 * level)/2, type = 6))
    if (console) {
    cat("\nMethod:", method, "\n")
    cat("\nThe index:", index, "\n\n")
    cat(level, "percent confidence interval:\n", lic, ";", lsc, 
        "\n\n")
        }
    }
if ( length(x) > 1) return(data.frame(row.names = NULL,method=method,index=index,confidence=level,lic=lic,lsc=lsc,nboot=nboot))
else return("Error data < 2")
}

