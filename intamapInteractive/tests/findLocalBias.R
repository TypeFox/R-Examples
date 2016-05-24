# Assuming that the soil type is the source of biases
library(intamapInteractive)
data(meuse)
coordinates(meuse) = ~x+y

lb = findLocalBias(meuse,gid = "soil",formulaString=as.formula(zinc~1))
lb$single$bias


meuseUnbias = removeLocalBias(meuse,localBias = lb, gid = "soil",
    formulaString = zinc~1) 
summary(meuseUnbias)