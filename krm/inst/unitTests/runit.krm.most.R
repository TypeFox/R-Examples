library("RUnit")
library("krm")


test.krm.mos.test <- function() {


tolerance=1e-3
# more stringent tolerance for one system to ensure algorithm accuracy
if(file.exists("D:/gDrive/3software/_checkReproducibility")) {
    tolerance=1e-6
} 
RNGkind("Mersenne-Twister", "Inversion")
dat.file.name=paste(system.file(package="krm")[1],'/misc/y1.txt', sep="") # needed for testing kernel sequences
seq.file.name=paste(system.file(package="krm")[1],'/misc/sim1.fasta', sep="") # needed for testing kernel sequences


##########################################################################################
# Euclidean kernel covariates

#### logistic regression

# parametric bootstrap
data=sim.liu.2008 (n=100, a=.1, seed=1) 
test = krm.most(y~x, data, regression.type="logistic", formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf", n.rho=2, n.mc = 100, range.rho=.99, verbose=TRUE)
checkEqualsNumeric(test$p.values, c(0.91,   0.90,   0.93,   0.91), tolerance = tolerance)

# perturbation
data=sim.liu.2008 (n=100, a=.1, seed=1) 
test = krm.most(y~x, data, regression.type="logistic", formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf", n.rho=2, n.mc = 100, inference.method="perturbation", verbose=TRUE)
# mvrnorm behaves differently between 32 bit and 64 bit
if (R.Version()$system %in% c("x86_64, mingw32")) {
    checkEqualsNumeric(test$p.values, c(0.87,     NA,   0.89,     NA), tolerance = tolerance)
} 

# LGL2008
data=sim.liu.2008 (n=50, a=.1, seed=1) 
test = krm.most(y~x, data, regression.type="logistic", formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf", n.mc = 100, range.rho=.99, inference.method="Davies", verbose=TRUE)
checkEqualsNumeric(test$p.values, 0.1223421, tolerance = tolerance)

# todo: add linear kernel


##########################################################################################
# sequence kernel covariates

#### logistic regression

# parametric bootstrap
dat=read.table(dat.file.name); names(dat)="y"
dat=cbind(dat, seq=unlist(readFastaFile(seq.file.name))); dat$seq=as.character(dat$seq)

test = krm.most (y~1, dat, regression.type="logistic", seq.file.name=seq.file.name, kern.type="mi", n.rho=2, n.mc = 5e1, inference.method="parametric.bootstrap", verbose=TRUE)
checkEqualsNumeric(test$p.values, c(0.68,   0.60,   0.66,   0.60), tolerance = tolerance)

test.2 = krm.most (y~1, dat, regression.type="logistic", formula.kern=~seq, kern.type="mi", n.rho=2, n.mc = 5e1, inference.method="parametric.bootstrap", verbose=TRUE)
checkEqualsNumeric(test.2$p.values, c(0.68,   0.60,   0.66,   0.60), tolerance = tolerance)

test.3 = krm.most (y~1, dat, regression.type="logistic", formula.kern=~seq, kern.type="mi", n.rho=2, n.mc = 5e1, inference.method="parametric.bootstrap", seq.start=1, seq.end=10, verbose=TRUE)
checkEqualsNumeric(test.3$p.values, c(0.62,   0.48,   0.54,   0.46), tolerance = tolerance)


## performance 

#data=sim.liu.2008 (n=50, a=.1, seed=1) 
#system.time({
#    test = krm.most(y~x, data, regression.type="logistic", formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf", inference.method="LGL2008", verbose=FALSE)
#})
## On gizmo
## parametric.bootstrap
## n=50: 117s; n=100: 240s; n=200: s
## perturbation
## n=50: 21s; n=100: 83s; n=200: 447s
## LGL2008
## n=50: .3s; n=100: 1s; n=200: 8s
## to find out the hours needed to run sim_most_batch: 
## one_test_time_in_sec *(100+20*2)/3600  # 140 is the number of tests to do for one core
## here we assume each core needs to run 140 simulations (100 for size and 40 for power)


#dat=read.table(dat.file.name); names(dat)="y"
#system.time({
#    test = krm.most (y~1, dat, regression.type="logistic", seq.file.name=seq.file.name, kern.type="mi", n.rho=2, inference.method="parametric.bootstrap", verbose=TRUE)
#})



}
