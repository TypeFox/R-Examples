# LVdataGen.R

# Generate the simulated data sets from the Lotka-Volterra model
#  to be used as reference data sets for testing and comparing
#  inference algorithms

require(smfsb)

# fix the RNG seed for reproducability
set.seed(1)

# first simulate some "perfect" data
LVperfect=simTs(c(x1=50,x2=100),0,30,2,stepLVc)
message("LVPerfect")
plot(LVperfect,plot.type="single")
print(LVperfect)
LVprey=LVperfect[,1]
message("LVprey")
print(LVprey)

# some noisy data
noiseSD=10
LVnoise10=LVperfect+rnorm(prod(dim(LVperfect)),0,noiseSD)
message("LVnoise10")
print(LVnoise10)
LVnoise10Scale10=LVnoise10/10
message("LVnoise10Scale10")
print(LVnoise10Scale10)
LVpreyNoise10=LVnoise10[,1]
message("LVpreyNoise10")
print(LVpreyNoise10)
LVpreyNoise10Scale10=LVpreyNoise10/10
message("LVpreyNoise10Scale10")
print(LVpreyNoise10Scale10)
noiseSD=30
LVnoise30=LVperfect+rnorm(prod(dim(LVperfect)),0,noiseSD)
message("LVnoise30")
print(LVnoise30)
LVnoise3010=LVnoise30
LVnoise3010[,2]=LVnoise10[,2]
message("LVnoise3010")
print(LVnoise3010)

LVtimed=as.timedData(LVperfect)
LVnoise10timed=as.timedData(LVnoise10)
rows=c(1,2,3,6,11,16)
LVirregular=LVtimed[rows,]
message("LVirregular")
print(LVirregular)
LVirregularNoise10=LVnoise10timed[rows,]
message("LVirregularNoise10")
print(LVirregularNoise10)

dump(c("LVperfect","LVprey","LVnoise10Scale10","LVpreyNoise10","LVpreyNoise10Scale10","LVnoise10","LVnoise30","LVnoise3010","LVirregular","LVirregularNoise10"),"LVdata.R")



# eof

