### R code from vignette source 'Ch-EFA.Rnw'

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
### code chunk number 2: ch:EFA:data
###################################################
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
        "marijuana", "hashish", "inhalants", "hallucinogenics", "amphetamine")



###################################################
### code chunk number 3: ch:EFA:life:tab
###################################################
"life" <-
structure(.Data = list(c(63., 34., 38., 59., 56., 62., 50., 65., 56., 69., 65., 64., 56., 60., 61., 49., 59., 63., 59., 65., 65., 64.,
        64., 67., 61., 68., 67., 65., 59., 58., 57.)
, c(51., 29., 30., 42., 38., 44., 39., 44., 46., 47., 48., 50., 44., 44., 45., 40., 42., 44., 44., 48., 48., 63.,
        43., 45., 40., 46., 45., 46., 43., 44., 46.)
, c(30., 13., 17., 20., 18., 24., 20., 22., 24., 24., 26., 28., 25., 22., 22., 22., 22., 23., 24., 28., 26., 21.,
        21., 23., 21., 23., 23., 24., 23., 24., 28.)
, c(13., 5., 7., 6., 7., 7., 7., 7., 11., 8., 9., 11., 10., 6., 8., 9., 6., 8., 8., 14., 9., 7., 6., 8., 10., 8.,
        8., 9., 10., 9., 9.)
, c(67., 38., 38., 64., 62., 69., 55., 72., 63., 75., 68., 66., 61., 65., 65., 51., 61., 67., 63., 68., 67., 68.,
        68., 74., 67., 75., 74., 71., 66., 62., 60.)
, c(54., 32., 34., 46., 46., 50., 43., 50., 54., 53., 50., 51., 48., 45., 49., 41., 43., 48., 46., 51., 49., 47.,
        47., 51., 46., 52., 51., 51., 49., 47., 49.)
, c(34., 17., 20., 25., 25., 28., 23., 27., 33., 29., 27., 29., 27., 25., 27., 23., 22., 26., 25., 29., 27., 25.,
        24., 28., 25., 29., 28., 28., 27., 25., 28.)
, c(15., 6., 7., 8., 10., 14., 8., 9., 19., 10., 10., 11., 12., 9., 10., 8., 7., 9., 8., 13., 10., 9., 8., 10., 11.,
        10., 10., 10., 12., 10., 11.)
)
, class = "data.frame" 
, names = c("m0", "m25", "m50", "m75", "w0", "w25", "w50", "w75")
, row.names = c("Algeria", "Cameroon", "Madagascar", "Mauritius", "Reunion", "Seychelles", "South Africa (C)", "South Africa (W)",
        "Tunisia", "Canada", "Costa Rica", "Dominican Rep.", "El Salvador", "Greenland", "Grenada", "Guatemala",
        "Honduras", "Jamaica", "Mexico", "Nicaragua", "Panama", "Trinidad (62)", "Trinidad (67)",
        "United States (66)", "United States (NW66)", "United States (W66)", "United States (67)", "Argentina",
        "Chile", "Colombia", "Ecuador")
)

toLatex(HSAURtable(life), pcol = 1, rownames = TRUE,
    caption = "Life expectancies for different countries by age and gender.",
    label = "ch:EFA:life:tab")


###################################################
### code chunk number 4: ch:EFA:life
###################################################
sapply(1:3, function(f) 
    factanal(life, factors = f, method ="mle")$PVAL)


###################################################
### code chunk number 5: ch:EFA:life3
###################################################
factanal(life, factors = 3, method ="mle")


###################################################
### code chunk number 6: ch:EFA:life:scores
###################################################
(scores <- factanal(life, factors = 3, method = "mle",
                   scores = "regression")$scores)


###################################################
### code chunk number 7: ch:EFA:life:3d
###################################################
cex <- 0.8
plot(scores[,1], scores[,2], type = "n", xlab = "Factor 1", ylab = "Factor 2")
text(scores[,1], scores[,2], abbreviate(rownames(life), 5), cex = cex)
plot(scores[,1], scores[,3], type = "n", xlab = "Factor 1", ylab = "Factor 3")
text(scores[,1], scores[,3], abbreviate(rownames(life), 5), cex = cex)
plot(scores[,2], scores[,3], type = "n", xlab = "Factor 2", ylab = "Factor 3")
text(scores[,2], scores[,3], abbreviate(rownames(life), 5), cex = cex)


###################################################
### code chunk number 8: ch:EFA:druguse:plot
###################################################
ord <- order.dendrogram(as.dendrogram(hclust(dist(druguse))))  
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

print(levelplot(druguse[ord, ord], at = do.breaks(c(-1.01, 1.01), 20),
          xlab = NULL, ylab = NULL, colorkey = list(space = "top"), 
          scales = list(x = list(rot = 90)),
          panel = panel.corrgram, label = TRUE))


###################################################
### code chunk number 9: ch:EFA:drugs
###################################################
sapply(1:6, function(nf)
    factanal(covmat = druguse, factors = nf, 
             method = "mle", n.obs = 1634)$PVAL)


###################################################
### code chunk number 10: ch:EFA:drugs
###################################################
(factanal(covmat = druguse, factors = 6, 
          method = "mle", n.obs = 1634))


###################################################
### code chunk number 11: ch:EFA:drugdiff
###################################################
pfun <- function(nf) {
    fa <- factanal(covmat = druguse, factors = nf, 
                   method = "mle", n.obs = 1634)
    est <- tcrossprod(fa$loadings) + diag(fa$uniquenesses)
    ret <- round(druguse - est, 3)
    colnames(ret) <- rownames(ret) <- 
        abbreviate(rownames(ret), 3)
    ret
}
pfun(6)


###################################################
### code chunk number 12: ch:opt
###################################################
op <- options(width = 150, prompt = "              R> ")
pfun2 <- pfun
pfun <- function(...) {
     x <- pfun2(...)
     rownames(x) <- paste("             ", rownames(x))
     x
}


###################################################
### code chunk number 13: ch:EFA:drufdiff34
###################################################
pfun(3)
pfun(4)


###################################################
### code chunk number 14: ch:opt
###################################################
options(op)


