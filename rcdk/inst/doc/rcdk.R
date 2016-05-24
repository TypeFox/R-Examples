### R code from vignette source 'rcdk.Rnw'

###################################################
### code chunk number 1: rcdk.Rnw:23-25
###################################################
options(width=74)
library(xtable)


###################################################
### code chunk number 2: <
###################################################
library(rcdk)


###################################################
### code chunk number 3: rcdk.Rnw:66-67 (eval = FALSE)
###################################################
## library(help=rcdk)


###################################################
### code chunk number 4: rcdk.Rnw:80-81 (eval = FALSE)
###################################################
## mols <- load.molecules( c('data1.sdf', '/some/path/data2.sdf') )


###################################################
### code chunk number 5: rcdk.Rnw:99-104 (eval = FALSE)
###################################################
## iter <- iload.molecules('verybig.sdf', type='sdf')
## while(hasNext(iter)) {
##  mol <- nextElem(iter)
##  print(get.property(mol, "cdk:Title"))
## }


###################################################
### code chunk number 6: rcdk.Rnw:109-111
###################################################
smile <- 'c1ccccc1CC(=O)C(N)CC1CCCCOC1'
mol <- parse.smiles(smile)[[1]]


###################################################
### code chunk number 7: rcdk.Rnw:115-117
###################################################
smiles <- c('CCC', 'c1ccccc1', 'CCCC(C)(C)CC(=O)NC')
mols <- parse.smiles(smiles)


###################################################
### code chunk number 8: rcdk.Rnw:124-134 (eval = FALSE)
###################################################
## options("java.parameters"=c("-Xmx4000m"))
## library(rcdk)
## 
## for (smile in smiles) {
##   m <- parse.smiles(smile)
##   ## perform operations on this molecule
##   
##   jcall("java/lang/System","V","gc")
##   gc()
## }


###################################################
### code chunk number 9: rcdk.Rnw:139-140 (eval = FALSE)
###################################################
## write.molecules(mols, filename='mymols.sdf')


###################################################
### code chunk number 10: rcdk.Rnw:147-148 (eval = FALSE)
###################################################
## write.molecules(mols, filename='mymols.sdf', together=FALSE)


###################################################
### code chunk number 11: rcdk.Rnw:152-154
###################################################
get.smiles(mols[[1]])
unlist(lapply(mols, get.smiles))


###################################################
### code chunk number 12: rcdk.Rnw:173-179 (eval = FALSE)
###################################################
## smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
##             'c1ccccc1Cc1ccccc1',
##             'C1CCC1CC(CN(C)(C))CC(=O)CC')
## mols <- parse.smiles(smiles)
## view.molecule.2d(mols[[1]])
## view.molecule.2d(mols)


###################################################
### code chunk number 13: rcdk.Rnw:198-202 (eval = FALSE)
###################################################
## dframe <- data.frame(x = runif(4),
##                      toxicity = factor(c('Toxic', 'Toxic', 'Nontoxic', 'Nontoxic')),
##                      solubility = c('yes', 'yes', 'no', 'yes'))
## view.table(mols[1:4], dframe)


###################################################
### code chunk number 14: rcdk.Rnw:231-235
###################################################
mol <- parse.smiles('c1ccccc1')[[1]]
set.property(mol, "title", "Molecule 1")
set.property(mol, "hvyAtomCount", 6)
get.property(mol, "title")


###################################################
### code chunk number 15: rcdk.Rnw:240-241
###################################################
get.properties(mol)


###################################################
### code chunk number 16: rcdk.Rnw:246-247 (eval = FALSE)
###################################################
## write.molecules(mol, 'tagged.sdf', write.props=TRUE)


###################################################
### code chunk number 17: rcdk.Rnw:259-264
###################################################
mol <- parse.smiles('c1ccccc1C(Cl)(Br)c1ccccc1')[[1]]
atoms <- get.atoms(mol)
bonds <- get.bonds(mol)
cat('No. of atoms =', length(atoms), '\n')
cat('No. of bonds =', length(bonds), '\n')


###################################################
### code chunk number 18: rcdk.Rnw:272-273
###################################################
unlist(lapply(atoms, get.symbol))


###################################################
### code chunk number 19: rcdk.Rnw:278-279
###################################################
coords <- get.point3d(atoms[[1]])


###################################################
### code chunk number 20: rcdk.Rnw:284-285 (eval = FALSE)
###################################################
## coords <- do.call('rbind', lapply(atoms, get.point3d))


###################################################
### code chunk number 21: rcdk.Rnw:290-293 (eval = FALSE)
###################################################
## if ( any(apply(coords, 2, function(x) length(unique(x))) == 1) ) {
##   print("molecule is flat")
## }


###################################################
### code chunk number 22: rcdk.Rnw:311-315
###################################################
mols <- parse.smiles(c('CC(C)(C)C','c1ccc(Cl)cc1C(=O)O', 'CCC(N)(N)CC'))
query <- '[#6D2]'
hits <- match(query, mols)
print(hits)


###################################################
### code chunk number 23: rcdk.Rnw:330-332
###################################################
dc <- get.desc.categories()
dc


###################################################
### code chunk number 24: rcdk.Rnw:339-340
###################################################
dn <- get.desc.names(dc[1])


###################################################
### code chunk number 25: rcdk.Rnw:348-350
###################################################
aDesc <- eval.desc(mol, dn[1])
allDescs <- eval.desc(mol, dn)


###################################################
### code chunk number 26: rcdk.Rnw:361-362 (eval = FALSE)
###################################################
## descNames <- unique(unlist(sapply(get.desc.categories(), get.desc.names)))           


###################################################
### code chunk number 27: rcdk.Rnw:366-375
###################################################
data(bpdata)
mols <- parse.smiles(bpdata[,1])
descNames <- c(
     'org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor',
     'org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor',
     'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor')          
descs <- eval.desc(mols, descNames) 
class(descs)
dim(descs)


###################################################
### code chunk number 28: rcdk.Rnw:392-397
###################################################
mol <- parse.smiles('CC(=O)CC(=O)NCN')[[1]]
convert.implicit.to.explicit(mol)
get.tpsa(mol)
get.xlogp(mol)
get.alogp(mol)


###################################################
### code chunk number 29: rcdk.Rnw:406-411 (eval = FALSE)
###################################################
## descs <- descs[, !apply(descs, 2, function(x) any(is.na(x)) )]
## descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )]
## r2 <- which(cor(descs)^2 > .6, arr.ind=TRUE)
## r2 <- r2[ r2[,1] > r2[,2] , ]
## descs <- descs[, -unique(r2[,2])]


###################################################
### code chunk number 30: rcdk.Rnw:421-423
###################################################
model <- lm(BP ~ khs.sCH3 + khs.sF + apol + nHBDon, data.frame(bpdata, descs))
summary(model)


###################################################
### code chunk number 31: rcdk.Rnw:430-435
###################################################
par(mar=c(4.3,4.3,1,1),cex.lab=1.3, pty='s')
plot(bpdata$BP, model$fitted, pch=19, 
     ylab='Predicted BP', xlab='Observed BP',
     xlim=range(bpdata$BP), ylim=range(bpdata$BP))
abline(0,1, col='red')


###################################################
### code chunk number 32: rcdk.Rnw:463-468
###################################################
smiles <- c('CCC', 'CCN', 'CCN(C)(C)', 
            'c1ccccc1Cc1ccccc1',
            'C1CCC1CC(CN(C)(C))CC(=O)CC')
mols <- parse.smiles(smiles)
fp <- get.fingerprint(mols[[1]], type='maccs')


###################################################
### code chunk number 33: rcdk.Rnw:476-480
###################################################
mols <- parse.smiles(bpdata[1:50,1]) 
fps <- lapply(mols, get.fingerprint, type='extended') 
fp.sim <- fp.sim.matrix(fps, method='tanimoto') 
fp.dist <- 1 - fp.sim 


###################################################
### code chunk number 34: rcdk.Rnw:488-490
###################################################
clustering <- hclust(as.dist(fp.dist))
plot(clustering, main='A Clustering of the BP dataset')


###################################################
### code chunk number 35: rcdk.Rnw:510-518
###################################################
query.mol <- parse.smiles('CC(=O)')[[1]]
target.mols <- parse.smiles(bpdata[,1])
query.fp <- get.fingerprint(query.mol, type='maccs')
target.fps <- lapply(target.mols, get.fingerprint, type='maccs')
sims <- unlist(lapply(target.fps, 
                      distance, 
                      fp2=query.fp, method='tanimoto'))
hits <- which(sims > 0.3)


###################################################
### code chunk number 36: rcdk.Rnw:522-527
###################################################
d <- data.frame(SMILES=bpdata[hits,1], Similarity=sims[hits])
row.names(d) <- NULL
d <- d[sort.list(d[,2], dec=TRUE),]
xtable(d, label="tab:sims", 
       caption="Summary of molecules in the BP dataset that are greater than 0.3 similar to acetaldehyde")


###################################################
### code chunk number 37: rcdk.Rnw:553-557
###################################################
sp <- get.smiles.parser()
molecule <- parse.smiles('N')[[1]]
convert.implicit.to.explicit(molecule)
formula <- get.mol2formula(molecule,charge=0)


###################################################
### code chunk number 38: rcdk.Rnw:566-567
###################################################
formula@mass


###################################################
### code chunk number 39: rcdk.Rnw:571-572
###################################################
formula@charge


###################################################
### code chunk number 40: rcdk.Rnw:580-581
###################################################
formula@isotopes


###################################################
### code chunk number 41: rcdk.Rnw:586-587
###################################################
formula@objectJ


###################################################
### code chunk number 42: rcdk.Rnw:591-592
###################################################
formula@string


###################################################
### code chunk number 43: rcdk.Rnw:598-599 (eval = FALSE)
###################################################
## formula <- set.charge.formula(formula, charge=1)


###################################################
### code chunk number 44: rcdk.Rnw:609-611
###################################################
formula <- get.formula('NH4', charge = 1);
formula


###################################################
### code chunk number 45: rcdk.Rnw:627-633
###################################################
mfSet <- generate.formula(18.03383, window=1, 
		elements=list(c("C",0,50),c("H",0,50),c("N",0,50)), 
		validation=FALSE);
for (i in 1:length(mfSet)) {
  print(mfSet[i])
}


###################################################
### code chunk number 46: rcdk.Rnw:640-642
###################################################
formula <- get.formula('NH4', charge = 0);
isvalid.formula(formula,rule=c("nitrogen","RDBE"))


###################################################
### code chunk number 47: rcdk.Rnw:658-661
###################################################
formula <- get.formula('CHCl3', charge = 0) 
isotopes <- get.isotopes.pattern(formula,minAbund=0.1)
isotopes


###################################################
### code chunk number 48: rcdk.Rnw:671-672
###################################################
plot(isotopes, type="h", xlab="m/z", ylab="Intensity")


