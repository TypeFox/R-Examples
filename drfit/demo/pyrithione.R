require(drfit)
data(pyrithione)

# MSPT and MSPHI give unreasonable fits with the linlogit model
rpyr <- drfit(pyrithione,linlogit=TRUE,linlogitWrong=c("MSPT","MSPHI"))
rpyr

# Consult the above result list to sort out the colors
drplot(rpyr,pyrithione,dtype="none",overlay=TRUE,bw=FALSE,
    colors=rainbow(14),xlim=c("auto",8))

drplot(rpyr,pyrithione)
