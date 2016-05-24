library("RUnit")
library("chngpt")

test.chngpt.test <- function() {


    RNGkind("Mersenne-Twister", "Inversion")    
    data=sim.chngpt("sigmoid4", type="step", n=250, seed=1, beta=0, x.distr="norm", e.=3.4, b.=-Inf)
    tolerance=1e-1
    if(file.exists("D:/gDrive/3software/_checkReproducibility")) tolerance=1e-6

    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="hinge", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="lr"); test
    checkEqualsNumeric(test$p.value, c(.3475), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="step", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="lr"); test
    checkEqualsNumeric(test$p.value, c(.6229), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="segmented", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="lr"); test
    checkEqualsNumeric(test$p.value, c(0.7709), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="stegmented", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="lr"); test
    checkEqualsNumeric(test$p.value, c(0.7077), tolerance=tolerance) 
    
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="hinge", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="score"); test
    checkEqualsNumeric(test$p.value, c(.3334), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="step", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="score"); test
    checkEqualsNumeric(test$p.value, c(.625), tolerance=tolerance) 
    test = chngpt.test (formula.null=y~z, formula.chngpt=~x, data, type="segmented", mc.n=1e4, verbose=0, chngpts.cnt=10, main.method="score"); test
    checkEqualsNumeric(test$p.value, c(0.777), tolerance=tolerance) 


#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x, data, type="step", mc.n=1e3, interaction.method="lr.mc", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.618), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x, data, type="hinge", mc.n=1e3, interaction.method="lr.mc", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.36), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e4, interaction.method="lr.mc", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.6196), tolerance=tolerance) 
#    test = chngpt.test (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e4, verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.6422), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="lr.pastor", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(0.8552631), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="weighted.two.sided", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.646), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="weighted.one.sided", verbose=0, chngpts.cnt=10); test
#    checkEqualsNumeric(test$p.value, c(.557), tolerance=tolerance) 
#    
#    test = chngpt.test.2 (formula.null=y~z, formula.chngpt=~x*z, data, type="step", mc.n=1e3, interaction.method="main.itxn", verbose=0, chngpts.cnt=50); test
#    checkEqualsNumeric(test$p.value, c(.679), tolerance=tolerance) 




}
