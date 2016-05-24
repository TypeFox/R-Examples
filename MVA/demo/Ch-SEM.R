### R code from vignette source 'Ch-SEM.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library("MVA")
set.seed(280875)
library("lattice")
lattice.options(default.theme =
    function()
        standard.theme("pdf", color = FALSE))

if (file.exists("deparse.R")) {
    if (!file.exists("figs")) dir.create("figs")
    source("deparse.R")
    options(prompt = "R> ", continue = "+  ", width = 64,
        digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

    options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           figtwo =   function() {par(mfrow = c(2,1))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           figthree = function() {par(mfrow = c(3,1))},
                           fourfig =  function() {par(mfrow = c(2,2))},
                           sixfig =   function() {par(mfrow = c(3,2))},
                           nomar = function() par("mai" = c(0, 0, 0, 0))))
}


###################################################
### code chunk number 2: setup
###################################################
library("sem")
ab <- c(0.73,
        0.70, 0.68,
        0.58, 0.61, 0.57,      
        0.46, 0.43, 0.40, 0.37,      
        0.56, 0.52, 0.48, 0.41, 0.72)      
ability <- diag(6) / 2
ability[upper.tri(ability)] <- ab
ability <- ability + t(ability)  
rownames(ability) <- colnames(ability) <- 
    c("SCA","PPE","PTE","PFE","EA","CP")


alienation <- matrix(c(11.839, 6.947, 6.819, 4.783, -3.834, -21.899,
                        6.947, 9.364, 5.091, 5.028, -3.889, -18.831,
                        6.819, 5.091,12.532, 7.495, -3.841, -21.748,
                        4.783, 5.028, 7.495, 9.986, -3.625, -18.775,
                       -3.834,-3.889,-3.841,-3.625,  9.60,   35.522,
                      -21.899,-18.831,-21.748,-18.775,35.522,450.283),
                      ncol = 6, byrow = TRUE)
rownames(alienation) <- colnames(alienation) <- c("Anomia67",
    "Powles67", "Anomia71", "Powles71", "Educ", "SEI")

d <-
c(0.447,       
  0.422, 0.619,       
  0.435, 0.604, 0.583,         
  0.114, 0.068, 0.053, 0.115,  
  0.203, 0.146, 0.139, 0.258, 0.349,          
  0.091, 0.103, 0.110, 0.122, 0.209, 0.221,
  0.082, 0.063, 0.066, 0.097, 0.321, 0.355, 0.201,
  0.513, 0.445, 0.365, 0.482, 0.186, 0.315, 0.150, 0.154,
  0.304, 0.318, 0.240, 0.368, 0.303, 0.377, 0.163, 0.219, 0.534,
  0.245, 0.203, 0.183, 0.255, 0.272, 0.323, 0.310, 0.288, 0.301, 0.302,
  0.101, 0.088, 0.074, 0.139, 0.279, 0.367, 0.232, 0.320, 0.204, 0.368, 0.340,
  0.245, 0.199, 0.184, 0.293, 0.278, 0.545, 0.232, 0.314, 0.394, 0.467, 0.392, 0.511)

druguse <- diag(13) / 2
druguse[upper.tri(druguse)] <- d
druguse <- druguse + t(druguse)

rownames(druguse) <- colnames(druguse) <- c("cigarettes", "beer", "wine", "liquor", "cocaine", 
        "tranquillizers", "drug store medication", "heroin", 
        "marijuana", "hashish", "inhalants", "haluucinogenics", "amphetamine")



###################################################
### code chunk number 3: ch:SEM:ability:plot
###################################################
ord <- order.dendrogram(as.dendrogram(hclust(dist(ability))))
panel.corrgram <-    
    function(x, y, z, subscripts, at,  
             level = 0.9, label = FALSE, ...)
{
    require("ellipse", quietly = TRUE)
    x <- as.numeric(x)[subscripts]   
    y <- as.numeric(y)[subscripts]   
    z <- as.numeric(z)[subscripts]   
    zcol <- level.colors(z, at = at, col.regions = grey.colors, ...)
    for (i in seq(along = z)) {
        ell <- ellipse(z[i], level = level, npoints = 50,   
                       scale = c(.2, .2), centre = c(x[i], y[i]))
        panel.polygon(ell, col = zcol[i], border = zcol[i], ...)
    }
    if (label)  
        panel.text(x = x, y = y, lab = 100 * round(z, 2), cex = 0.8,
                   col = ifelse(z < 0, "white", "black"))
}

print(levelplot(ability[ord, ord], at = do.breaks(c(-1.01, 1.01), 20),
          xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
          scales = list(x = list(rot = 90)),
          panel = panel.corrgram, label = TRUE))


###################################################
### code chunk number 4: ch:SEM:ability-model
###################################################
mod <- c("Ability     -> SCA, lambda1, NA",
         "Ability     -> PPE, lambda2, NA",
         "Ability     -> PTE, lambda3, NA",
         "Ability     -> PFE, lambda4, NA",
         "Aspiration  -> EA, lambda5, NA",
         "Aspiration  -> CP, lambda6, NA",
         "Ability    <-> Aspiration, rho, NA",
         "SCA        <-> SCA, theta1, NA",
         "PPE        <-> PPE, theta2, NA",
         "PTE        <-> PTE, theta3, NA",
         "PFE        <-> PFE, theta4, NA",
         "EA         <-> EA, theta5, NA",
         "CP         <-> CP, theta6, NA",
         "Ability    <-> Ability, NA, 1",
         "Aspiration <-> Aspiration, NA, 1")
writeLines(mod, con = "ability_model.txt")


###################################################
### code chunk number 5: ch:SEM:ability-sem
###################################################
ability_model <- specify.model(file = "ability_model.txt")
ability_sem <- sem(ability_model, ability, 556)


###################################################
### code chunk number 6: ch:SEM:ability-sem
###################################################
writeLines(readLines("ability_model.txt"))
r <- file.remove("ability_model.txt")


###################################################
### code chunk number 7: ch:SEM:ability-summary
###################################################
summary(ability_sem)


###################################################
### code chunk number 8: ch:SEM:ability-summary
###################################################
su <- summary(ability_sem)


###################################################
### code chunk number 9: ch:SEM:ability-path
###################################################
path.diagram(ability_sem, file = "ability_sem", 
             ignore.double = FALSE, edge.labels = "both")


###################################################
### code chunk number 10: ch:SEM:ability-files
###################################################
file.remove("ability_sem.dot")
file.copy("ability_sem.pdf", "figs", overwrite = TRUE)
file.remove("ability_sem.pdf")


###################################################
### code chunk number 11: ch:SEM:druguse-model
###################################################
mod <- c("Alcohol   -> Cigs, lambda1, NA",
         "Alcohol   -> Beer, lambda3, NA",
         "Alcohol   -> Wine, lambda4, NA",
         "Alcohol   -> Liqr, lambda6, NA",
         "Cannabis  -> Cigs, lambda2, NA",
         "Cannabis  -> Wine, lambda5, NA",
         "Cannabis  -> Marj, lambda12, NA",
         "Cannabis  -> Hash, lambda13, NA",
         "Hard      -> Liqr, lambda7, NA",
         "Hard      -> Cocn, lambda8, NA",
         "Hard      -> Tran, lambda9, NA",
         "Hard      -> Drug, lambda10, NA",
         "Hard      -> Hern, lambda11, NA",
         "Hard      -> Hash, lambda14, NA",
         "Hard      -> Inhl, lambda15, NA",
         "Hard      -> Hall, lambda16, NA",
         "Hard      -> Amph, lambda17, NA",
         "Cigs     <-> Cigs, theta1, NA",
         "Beer     <-> Beer, theta2, NA",
         "Wine     <-> Wine, theta3, NA",
         "Liqr     <-> Liqr, theta4, NA",
         "Cocn     <-> Cocn, theta5, NA",
         "Tran     <-> Tran, theta6, NA",
         "Drug     <-> Drug, theta7, NA",
         "Hern     <-> Hern, theta8, NA",
         "Marj     <-> Marj, theta9, NA",
         "Hash     <-> Hash, theta10, NA",
         "Inhl     <-> Inhl, theta11, NA",
         "Hall     <-> Hall, theta12, NA",
         "Amph     <-> Amph, theta13, NA",
         "Alcohol  <-> Alcohol, NA, 1",
         "Cannabis <-> Cannabis, NA, 1",
         "Hard     <-> Hard, NA, 1",
         "Alcohol  <-> Cannabis, rho1, NA",
         "Alcohol  <-> Hard, rho2, NA",
         "Cannabis <-> Hard, rho3, NA")
writeLines(mod, con = "druguse_model.txt")


###################################################
### code chunk number 12: ch:SEM:druguse-names
###################################################
rownames(druguse) <- colnames(druguse) <- c("Cigs", 
    "Beer", "Wine", "Liqr", "Cocn", "Tran", "Drug", 
    "Hern", "Marj", "Hash", "Inhl", "Hall", "Amph")


###################################################
### code chunk number 13: ch:SEM:druguse-sem
###################################################
druguse_model <- specify.model(file = "druguse_model.txt")
druguse_sem <- sem(druguse_model, druguse, 1634)


###################################################
### code chunk number 14: ch:SEM:druguse-sem
###################################################
writeLines(readLines("druguse_model.txt"))
r <- file.remove("druguse_model.txt")


###################################################
### code chunk number 15: ch:SEM:druguse-path
###################################################
path.diagram(druguse_sem, file = "druguse_sem", 
             ignore.double = FALSE, edge.labels = "both") 


###################################################
### code chunk number 16: ch:SEM:druguse-files
###################################################
file.remove("druguse_sem.dot")
file.copy("druguse_sem.pdf", "figs", overwrite = TRUE)
file.remove("druguse_sem.pdf")


###################################################
### code chunk number 17: ch:SEM:druguse-summary
###################################################
summary(druguse_sem)


###################################################
### code chunk number 18: ch:SEM:druguse-summary
###################################################
su <- summary(druguse_sem)


###################################################
### code chunk number 19: ch:SEM:druguse-cov
###################################################
round(druguse_sem$S - druguse_sem$C, 3)


###################################################
### code chunk number 20: ch:SEM:alienation:plot
###################################################
a <- cov2cor(alienation)
ord <- order.dendrogram(as.dendrogram(hclust(dist(a))))
panel.corrgram <-    
    function(x, y, z, subscripts, at,  
             level = 0.9, label = FALSE, ...)
{
    require("ellipse", quietly = TRUE)
    x <- as.numeric(x)[subscripts]   
    y <- as.numeric(y)[subscripts]   
    z <- as.numeric(z)[subscripts]   
    zcol <- level.colors(z, at = at, col.regions = grey.colors, ...)
    for (i in seq(along = z)) {
        ell <- ellipse(z[i], level = level, npoints = 50,   
                       scale = c(.2, .2), centre = c(x[i], y[i]))
        panel.polygon(ell, col = zcol[i], border = zcol[i], ...)
    }
    if (label)  
        panel.text(x = x, y = y, lab = 100 * round(z, 2), cex = 0.8,
                   col = ifelse(z < 0, "white", "black"))
}

print(levelplot(a[ord, ord], at = do.breaks(c(-1.01, 1.01), 20),
          xlab = NULL, ylab = NULL, colorkey = list(space = "top"),
          scales = list(x = list(rot = 90)),
          panel = panel.corrgram, label = TRUE))


###################################################
### code chunk number 21: ch:SEM:alienation-model
###################################################
mod <- c("SES           -> Educ, NA, 1",
    "SES           -> SEI, lambda1, NA",
    "Alienation67  -> Anomia67, NA, 1",
    "Alienation67  -> Powles67, lambda2, NA",
    "Alienation71  -> Anomia71, NA, 1",
    "Alienation71  -> Powles71, lambda3, NA",
    "SES           -> Alienation67, beta1, NA",
    "SES           -> Alienation71, beta2, NA",
    "Alienation67  -> Alienation71, beta3, NA",
    "Educ         <-> Educ, theta1, NA",
    "SEI          <-> SEI, theta2, NA",
    "SES          <-> SES, delta0, NA",
    "Anomia67     <-> Anomia67, theta3, NA",
    "Powles67     <-> Powles67, theta4, NA",
    "Anomia71     <-> Anomia71, theta5, NA",
    "Powles71     <-> Powles71, theta6, NA",
    "Alienation67 <-> Alienation67, delta1, NA",
    "Alienation71 <-> Alienation71, delta2, NA")
writeLines(mod, con = "alienation_model.txt")
mod2 <- c(mod, "Anomia67 <-> Anomia71,psi,NA")
writeLines(mod2, con = "alienation_model2.txt")


###################################################
### code chunk number 22: ch:SEM:alienation-sem
###################################################
alienation_model <- specify.model(
    file = "alienation_model.txt")
alienation_sem <- sem(alienation_model, alienation, 932)


###################################################
### code chunk number 23: ch:SEM:alienation-sem
###################################################
writeLines(readLines("alienation_model.txt"))
r <- file.remove("alienation_model.txt")


###################################################
### code chunk number 24: ch:SEM:alienation-path
###################################################
path.diagram(alienation_sem, file = "alienation_sem", 
             ignore.double = FALSE, edge.labels = "both")


###################################################
### code chunk number 25: ch:SEM:alienation-files
###################################################
file.remove("alienation_sem.dot")
file.copy("alienation_sem.pdf", "figs", overwrite = TRUE)
file.remove("alienation_sem.pdf")


###################################################
### code chunk number 26: ch:SEM:alienation-summary
###################################################
summary(alienation_sem)


###################################################
### code chunk number 27: ch:SEM:alienation-summary
###################################################
su <- summary(alienation_sem)


###################################################
### code chunk number 28: ch:SEM:alienation-sem2
###################################################
alienation_model2 <- specify.model(file = "alienation_model2.txt")
alienation_sem2 <- sem(alienation_model2, alienation, 932)
su <- summary(alienation_sem2)
r <- file.remove("alienation_model2.txt")


