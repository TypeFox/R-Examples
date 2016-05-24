library(drfit)
data(pyrithione)
rpyr <- drfit(pyrithione,linlogit=TRUE,linlogitWrong=c("MSPT","MSPHI","PyS"))
print(rpyr,digits=2)
