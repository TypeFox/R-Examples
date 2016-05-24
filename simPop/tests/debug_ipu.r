## goal: calculate household weights that should also fulfil person-type constraints
rm(list=ls())
library(simPop)
# load sample and population data
data(eusilcS)
data(eusilcP)

# variable generation and preparation
eusilcS$hsize <- factor(eusilcS$hsize)

# make sure, factor levels in sample and population match
eusilcP$region <- factor(eusilcP$region, levels = levels(eusilcS$db040))
eusilcP$gender <- factor(eusilcP$gender, levels = levels(eusilcS$rb090))
eusilcP$hsize  <- factor(eusilcP$hsize , levels = levels(eusilcS$hsize))

# generate input matrix
# we want to adjust to variable "db040" (region) as household variables and
# variable "rb090" (gender) as individual information
samp <- data.table(eusilcS)
pop <-  data.table(eusilcP)
setkeyv(samp, "db030")
hh <- samp[!duplicated(samp$db030),]
hhpop <- pop[!duplicated(pop$hid),]

# reg contains for each region the number of households
reg <- data.table(model.matrix(~db040 +0, data=hh))
# hsize contains for each household size the number of households
hsize <- data.table(model.matrix(~factor(hsize) +0, data=hh))

# aggregate persons-level characteristics per household
# gender contains for each household the number of males and females
gender <- data.table(model.matrix(~db030+rb090 +0, data=samp))
setkeyv(gender, "db030")
gender <- gender[, lapply(.SD, sum), by = key(gender)]

# bind together and use it as input
inp <- cbind(reg,
    hsize,
    gender)
# the totals we want to calibrate to
con <- c(
  as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$region)),
  as.list(xtabs(rep(1, nrow(hhpop)) ~ hhpop$hsize)),
  as.list(xtabs(rep(1, nrow(eusilcP)) ~ eusilcP$gender))
)
# we need to have the same names as in 'inp'
names(con) <- setdiff(names(inp), "db030")

# run ipu und check results
res <- ipu(inp=inp, hid="db030", con=con, verbose=TRUE)

is <- sapply(2:(ncol(res)-1), function(x) { 
  sum(res[,x]*res$weights)
}) 
data.frame(required=unlist(con), is=is)




######Basic not converting example
# basic example
require(simPop)
inp <- as.data.frame(matrix(0, nrow=8, ncol=7))
colnames(inp) <- c("hhid","hh1","hh2","p1","p2","p3","p4")
inp$hhid <- 1:8
inp$hh1[1:3] <- 1
inp$hh2[4:8] <- 1
inp$p1 <- c(1,1,2,1,0,1,2,1)
inp$p2 <- c(1,0,1,0,2,1,1,1)
inp$p3 <- c(1,1,0,2,1,0,2,0)
inp$p4 <- c(0,0,0,0,0,0,0,0)
con <- list(hh1=35, hh2=65, p1=91, p2=65, p3=104,p4=5)
res <- ipu(inp=inp, hid="hhid", con=con, verbose=TRUE)