### R code from vignette source 'sybil.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sybil.Rnw:433-434
###################################################
library(sybil)


###################################################
### code chunk number 2: sybil.Rnw:443-444
###################################################
library(help = "sybil")


###################################################
### code chunk number 3: sybil.Rnw:448-449
###################################################
help(doubleGeneDel)


###################################################
### code chunk number 4: sybil.Rnw:453-454
###################################################
help.search("flux variability analysis")


###################################################
### code chunk number 5: sybil.Rnw:457-458 (eval = FALSE)
###################################################
## vignette("sybil")


###################################################
### code chunk number 6: sybil.Rnw:479-480
###################################################
mp  <- system.file(package = "sybil", "extdata")


###################################################
### code chunk number 7: sybil.Rnw:483-484
###################################################
mod <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")


###################################################
### code chunk number 8: sybil.Rnw:504-505 (eval = FALSE)
###################################################
## modelorg2tsv(mod, prefix = "Ec_core")


###################################################
### code chunk number 9: sybil.Rnw:510-511
###################################################
data(Ec_core)


###################################################
### code chunk number 10: sybil.Rnw:524-525
###################################################
ex <- findExchReact(Ec_core)


###################################################
### code chunk number 11: sybil.Rnw:528-530
###################################################
upt <- uptReact(ex)
ex[upt]


###################################################
### code chunk number 12: sybil.Rnw:536-538
###################################################
mod <- changeBounds(Ec_core, ex[c("EX_glc(e)", "EX_lac_D(e)")], lb = c(0, -10))
findExchReact(mod)


###################################################
### code chunk number 13: sybil.Rnw:557-558
###################################################
optL <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = FALSE)


###################################################
### code chunk number 14: sybil.Rnw:563-564
###################################################
opt <- optimizeProb(Ec_core, algorithm = "fba", retOptSol = TRUE)


###################################################
### code chunk number 15: sybil.Rnw:570-571
###################################################
lp_obj(opt)


###################################################
### code chunk number 16: sybil.Rnw:575-576
###################################################
checkOptSol(opt)


###################################################
### code chunk number 17: sybil.Rnw:590-591
###################################################
fba <- optimizeProb(Ec_core, algorithm = "fba")


###################################################
### code chunk number 18: sybil.Rnw:594-595
###################################################
mod_obj(fba)


###################################################
### code chunk number 19: sybil.Rnw:599-600
###################################################
mtf <- optimizeProb(Ec_core, algorithm = "mtf", wtobj = mod_obj(fba))


###################################################
### code chunk number 20: sybil.Rnw:603-604
###################################################
lp_obj(mtf)


###################################################
### code chunk number 21: sybil.Rnw:609-611
###################################################
nvar(fluxdist(fba))
nvar(fluxdist(mtf))


###################################################
### code chunk number 22: sybil.Rnw:616-618
###################################################
help("sysBiolAlg_fba-class")
help("sysBiolAlg_mtf-class")


###################################################
### code chunk number 23: sybil.Rnw:621-623
###################################################
?fba
?mtf


###################################################
### code chunk number 24: sybil.Rnw:628-630
###################################################
fl <- getFluxDist(mtf)
length(fl)


###################################################
### code chunk number 25: sybil.Rnw:635-637
###################################################
fd <- getFluxDist(mtf, ex)
getNetFlux(fd)


###################################################
### code chunk number 26: sybil.Rnw:642-643
###################################################
mod_obj(mtf)


###################################################
### code chunk number 27: sybil.Rnw:668-683
###################################################
data(Ec_core)
optimizeProb(Ec_core)
Ec_core_C <- Ec_core # we copy the model to compare later.

# define carbon sources...
CS_List = c('EX_fru(e)' , 'EX_glc(e)', 'EX_fum(e)')
CS_CA = c(6,6,4); # ...and the number of carbon atoms per molekule

cntCAtoms = 6 * abs(lowbnd(Ec_core_C)[react_id(Ec_core_C)=='EX_glc(e)'])
# prohibit excretion of carbon sources
uppbnd(Ec_core_C)[react_id(Ec_core_C) %in% CS_List] = 0
# co2 mustn't act as source, too!
lowbnd(Ec_core_C)[react_id(Ec_core_C) == "EX_co2(e)"] <- 0
lowbnd(Ec_core_C)[react_id(Ec_core_C) %in% CS_List] = -1000
findExchReact(Ec_core_C)


###################################################
### code chunk number 28: sybil.Rnw:687-695
###################################################
help("fbaEasyConstraint")

ec <- list(
	react=list(match(CS_List,react_id(Ec_core_C))), # affected reactions
	x=list(-1*CS_CA), # coefficient
	ub = cntCAtoms, # atoms count is the upper bound
	rtype="U" # type of constraint U => upper bound
)


###################################################
### code chunk number 29: sybil.Rnw:698-701
###################################################
opt <- optimizeProb(Ec_core_C, algorithm=("fbaEasyConstraint"), easyConstraint=ec)
opt2 <- optimizeProb(Ec_core_C, algorithm=("mtfEasyConstraint"), easyConstraint=ec)
mtf_glc <- optimizeProb(Ec_core, algorithm="mtf") # normal opt to compare to.


###################################################
### code chunk number 30: sybil.Rnw:704-719
###################################################
# check fluxes in FBA result:
print(data.frame(Csrc=CS_List, 
		CS_CA, flx=fluxes(opt)[match(CS_List,react_id(Ec_core_C))]))
# check fluxes in MTF result:
print(data.frame(Csrc=CS_List, 
		CS_CA, flx=fluxes(opt2)[match(CS_List,react_id(Ec_core_C))]))
# look at biomass values:
print(data.frame(fbaCA=mod_obj(opt), 
		mtfCA=mod_obj(opt2), mtf_glc=mod_obj(mtf_glc)))
# compare the absolute sum over fluxes
print(data.frame(mtfCA=lp_obj(opt2), mtf_glc=lp_obj(mtf_glc)))
# write problem to file (optional)
prob <- sysBiolAlg(Ec_core_C, algorithm = "fbaEasyConstraint",
			easyConstraint=ec, useNames=TRUE)
writeProb(problem(prob), fname='test_easyCons.lp')


###################################################
### code chunk number 31: sybil.Rnw:729-730
###################################################
ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0)


###################################################
### code chunk number 32: sybil.Rnw:751-753
###################################################
ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0,
                   algorithm = "lmoma", wtflux = getFluxDist(mtf))


###################################################
### code chunk number 33: sybil.Rnw:769-772 (eval = FALSE)
###################################################
## ko <- optimizeProb(Ec_core, gene = "b2276", lb = 0, ub = 0,
##                    algorithm = "room", wtflux = getFluxDist(mtf),
##                    solverParm = list(PRESOLVE = GLP_ON))


###################################################
### code chunk number 34: sybil.Rnw:804-805
###################################################
opt <- oneGeneDel(Ec_core)


###################################################
### code chunk number 35: sybil.Rnw:819-820
###################################################
checkOptSol(opt)


###################################################
### code chunk number 36: sybil.Rnw:824-825
###################################################
plot(opt, nint = 20)


###################################################
### code chunk number 37: sybil.Rnw:831-834
###################################################
opt <- oneGeneDel(Ec_core, algorithm = "lmoma", wtflux = getFluxDist(mtf))
checkOptSol(opt)
plot(opt, nint = 20)


###################################################
### code chunk number 38: sybil.Rnw:840-841
###################################################
opt <- geneDeletion(Ec_core)


###################################################
### code chunk number 39: sybil.Rnw:844-846 (eval = FALSE)
###################################################
## opt2 <- geneDeletion(Ec_core, combinations = 2)
## opt3 <- geneDeletion(Ec_core, combinations = 3)


###################################################
### code chunk number 40: sybil.Rnw:862-864
###################################################
opt <- fluxVar(Ec_core, percentage = 80, verboseMode = 0)
plot(opt)


###################################################
### code chunk number 41: sybil.Rnw:892-894
###################################################
opt <- robAna(Ec_core, ctrlreact = "EX_o2(e)", verboseMode = 0)
plot(opt)


###################################################
### code chunk number 42: sybil.Rnw:911-918
###################################################
Ec_core_wo_glc <- changeUptake(Ec_core, off = "glc_D[e]")
opt <- phpp(Ec_core_wo_glc,
            ctrlreact = c("EX_succ(e)", "EX_o2(e)"),
            redCosts = TRUE,
            numP = 25,
            verboseMode = 0)
plot(opt)


###################################################
### code chunk number 43: phpp_rf
###################################################
plot(opt, "EX_succ(e)")


###################################################
### code chunk number 44: phpp_rs
###################################################
plot(opt, "EX_o2(e)")


###################################################
### code chunk number 45: sybil.Rnw:952-953
###################################################
opt <- oneGeneDel(Ec_core, algorithm = "fba", fld = "all")


###################################################
### code chunk number 46: sybil.Rnw:956-957
###################################################
sum <- summaryOptsol(opt, Ec_core)


###################################################
### code chunk number 47: sybil.Rnw:971-972
###################################################
printExchange(sum, j = c(1:50), dense = TRUE)


###################################################
### code chunk number 48: sybil.Rnw:987-991
###################################################
ref    <- optimizeProb(Ec_core)
opt    <- oneGeneDel(Ec_core)
let    <- lethal(opt, wt = mod_obj(ref))
nletid <- c(1:length(allGenes(Ec_core)))[! let] 


###################################################
### code chunk number 49: sybil.Rnw:1000-1001 (eval = FALSE)
###################################################
## gmat <- combn(nletid, 3)


###################################################
### code chunk number 50: sybil.Rnw:1006-1007 (eval = FALSE)
###################################################
## opt <- multiDel(Ec_core, nProc = 4, todo = "geneDeletion", del1 = gmat)


###################################################
### code chunk number 51: sybil.Rnw:1018-1019 (eval = FALSE)
###################################################
## mapply(checkOptSol, opt)


###################################################
### code chunk number 52: sybil.Rnw:1030-1032
###################################################
opt <- optimizeProb(Ec_core, poCmd = list("getRedCosts"))
postProc(opt)


###################################################
### code chunk number 53: sybil.Rnw:1068-1069 (eval = FALSE)
###################################################
## optimizeProb(Ec_core, method = "exact")


###################################################
### code chunk number 54: sybil.Rnw:1072-1073 (eval = FALSE)
###################################################
## optimizeProb(Ec_core, solver = "cplexAPI", method = "dualopt")


###################################################
### code chunk number 55: sybil.Rnw:1096-1099 (eval = FALSE)
###################################################
## opt <- oneGeneDel(Ec_core,
##                   solverParm = list(TM_LIM = 1000,
##                                     PRESOLVE = GLP_ON))


###################################################
### code chunk number 56: sybil.Rnw:1112-1116 (eval = FALSE)
###################################################
## opt <- optimizeProb(Ec_core,
##                     solverParm = list(CPX_PARAM_SCRIND = CPX_ON,
##                                       CPX_PARAM_EPRHS = 1E-09),
##                     solver = "cplexAPI")


###################################################
### code chunk number 57: sybil.Rnw:1135-1139 (eval = FALSE)
###################################################
## opt <- optimizeProb(Ec_core,
##                     solverParm = list(verbose = "full",
##                                       timeout = 10),
##                     solver = "lpSolveAPI")


###################################################
### code chunk number 58: sybil.Rnw:1153-1154
###################################################
help(SYBIL_SETTINGS)


###################################################
### code chunk number 59: sybil.Rnw:1176-1177 (eval = FALSE)
###################################################
## SYBIL_SETTINGS("parameter name", value)


###################################################
### code chunk number 60: sybil.Rnw:1182-1183 (eval = FALSE)
###################################################
## SYBIL_SETTINGS("parameter name")


###################################################
### code chunk number 61: sybil.Rnw:1188-1189 (eval = FALSE)
###################################################
## SYBIL_SETTINGS()


###################################################
### code chunk number 62: sybil.Rnw:1204-1205
###################################################
SYBIL_SETTINGS("SOLVER", "cplexAPI", loadPackage = FALSE)


###################################################
### code chunk number 63: sybil.Rnw:1211-1212
###################################################
SYBIL_SETTINGS("METHOD")


###################################################
### code chunk number 64: sybil.Rnw:1215-1216
###################################################
SYBIL_SETTINGS("SOLVER", "glpkAPI")


###################################################
### code chunk number 65: sybil.Rnw:1219-1220
###################################################
SYBIL_SETTINGS("METHOD")


###################################################
### code chunk number 66: sybil.Rnw:1255-1257
###################################################
data(Ec_core)
Ec_core


###################################################
### code chunk number 67: sybil.Rnw:1261-1262
###################################################
help("modelorg")


###################################################
### code chunk number 68: sybil.Rnw:1269-1270
###################################################
react_num(Ec_core)


###################################################
### code chunk number 69: sybil.Rnw:1273-1274
###################################################
id <- react_id(Ec_core)


###################################################
### code chunk number 70: sybil.Rnw:1277-1278
###################################################
react_id(Ec_core)[13] <- "biomass"


###################################################
### code chunk number 71: sybil.Rnw:1282-1284
###################################################
cg <- gray(0:8/8)
image(S(Ec_core), col.regions = c(cg, rev(cg)))


###################################################
### code chunk number 72: sybil.Rnw:1297-1298 (eval = FALSE)
###################################################
## mod <- readTSVmod(reactList = "reactionList.txt")


###################################################
### code chunk number 73: sybil.Rnw:1322-1323
###################################################
help("optsol")


###################################################
### code chunk number 74: sybil.Rnw:1327-1329
###################################################
os <- optimizeProb(Ec_core)
is(os)


###################################################
### code chunk number 75: sybil.Rnw:1332-1333
###################################################
lp_obj(os)


###################################################
### code chunk number 76: sybil.Rnw:1336-1337
###################################################
getFluxDist(os)


###################################################
### code chunk number 77: sybil.Rnw:1382-1383
###################################################
lp <- optObj(solver = "glpkAPI", method = "exact")


###################################################
### code chunk number 78: sybil.Rnw:1407-1408
###################################################
lp <- initProb(lp)


###################################################
### code chunk number 79: sybil.Rnw:1412-1417
###################################################
cm <- Matrix(c(0.5, 2, 1, 1), nrow = 2)
loadLPprob(lp, nCols = 2, nRows = 2, mat = cm,
           lb = c(0, 0), ub = rep(1000, 2), obj = c(1, 1),
           rlb = c(0, 0), rub = c(4.5, 9), rtype = c("U", "U"),
           lpdir = "max")


###################################################
### code chunk number 80: sybil.Rnw:1422-1423
###################################################
lp


###################################################
### code chunk number 81: sybil.Rnw:1426-1427
###################################################
status <- solveLp(lp)


###################################################
### code chunk number 82: sybil.Rnw:1430-1431
###################################################
getMeanReturn(code = status, solver = solver(lp))


###################################################
### code chunk number 83: sybil.Rnw:1434-1436
###################################################
status <- getSolStat(lp)
getMeanStatus(code = status, solver = solver(lp))


###################################################
### code chunk number 84: sybil.Rnw:1440-1442
###################################################
getObjVal(lp)
getFluxDist(lp)


###################################################
### code chunk number 85: sybil.Rnw:1445-1446
###################################################
getRedCosts(lp)


###################################################
### code chunk number 86: sybil.Rnw:1449-1451
###################################################
delProb(lp)
lp


###################################################
### code chunk number 87: sybil.Rnw:1482-1484
###################################################
ec <- sysBiolAlg(Ec_core, algorithm = "fba")
is(ec)


###################################################
### code chunk number 88: sybil.Rnw:1489-1490
###################################################
opt <- optimizeProb(ec)


###################################################
### code chunk number 89: sybil.Rnw:1497-1500
###################################################
ecr <- sysBiolAlg(Ec_core, algorithm = "room", wtflux = opt$fluxes)
is(ecr)
ecr


###################################################
### code chunk number 90: sybil.Rnw:1549-1550 (eval = FALSE)
###################################################
## promptSysBiolAlg(algorithm = "foo")


