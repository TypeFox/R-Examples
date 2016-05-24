ci.boot <- function(x, method = "all", sigma.t = NULL, conf = 0.95){

indices <- c("all", "norm", "basic", "perc", "BCa", "student")
method <- match.arg(method, indices)

mes <- NULL
B <- x$dist 
se <- x$res[4]
est <- x$res[1]
alpha <- 1 - conf
R <- length(B)
uZ <- round((1 - (alpha/2))*R, 0)
lZ <- round((alpha/2)*R, 0)
sB <- sort(B)


if(method == "all" | method == "norm")
{z <- qnorm(1 - (alpha/2))
nciL <- est - z * se
nciU <- est + z * se  
nci <- c(nciL, nciU)
}

else nci <- c(NA, NA)

if(method == "all" | method == "basic")
{
bciL <- 2*est-sB[uZ]
bciU <- 2*est-sB[lZ]
bci <- c(bciL, bciU)
}

else bci <- c(NA, NA)

if(method == "all" | method == "perc")
{
pciL <- sB[lZ]
pciU <- sB[uZ]
pci <- c(pciL, pciU)
}

else pci <- c(NA, NA)

if(method == "all" | method == "BCa")
{
pv <- pseudo.v(x$data, statistic = x$statistic)
jk <- mean(pv[,1])
a <- sum((jk - pv[,1])^3)/(6*((sum((jk - pv[,1])^2))^1.5))# Eq. 3.11 Manly 
p <- length(sB[sB > est])/R
z0 <- qnorm(1 - p)
zU <- (z0 - qnorm(alpha/2))/(1 - a*(z0 - qnorm(alpha/2)))+z0
zL <- (z0 + qnorm(alpha/2))/(1 - a*(z0 + qnorm(alpha/2)))+z0
pL <- pnorm(zL)
pU <- pnorm(zU)
uZ <- round(pL*R, 0)
lZ <- round(pU*R, 0); uZ <- ifelse(uZ == 0, 1, uZ) 
bcciU <- sB[lZ]
bcciL <- sB[uZ]
bcci <- c(bcciL, bcciU)
}

else bcci <- c(NA, NA)

if(method == "all" | method == "student")
{
if(is.null(sigma.t)){
mes <- "Bootstrap SEs req'd for studentized intervals"
sci <- c(NA, NA)
}
else{
t <- (B - est)/sigma.t
sciL <- est - qt(1-(alpha/2))*se
sciU <- est + qt(1-(alpha/2))*se
sci <- c(sciL, sciU)
}
}

else sci = c(NA, NA)

head <- paste(conf*100,"%", " Bootstrap confidence interval(s)", sep = "")
ends <-  c(paste(as.character(c((1 - conf)/2, 1 - ((1 - conf)/2)) * 100), "%", sep = ""))
res <- matrix(nrow = 5, data=rbind(nci, bci, pci, bcci, sci), dimnames = list(c("Normal","Basic","Percentile", "BCa","Studentized"), ends))
out <- list(head = head, res = res, mes = mes, a = ifelse(method == "BCa"|method == "all", a, 0))
class(out) <- "ciboot"
out
} 

print.ciboot <- function(x, digits = max(3, getOption("digits")), ...) 
{
    cat("\n")
    cat(x$head, "\n\n")
    print(x$res, digits = digits, justify = "center")
    cat("\n")
    if(!is.null(x$mes))cat(x$mes,"\n")
    invisible(x)
}
 


 



