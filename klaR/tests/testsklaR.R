###############
# preliminaries
library("klaR")
library("MASS")
data(B3)
postscript("testklaR.ps")

############
# classifier

# Naive Bayes
set.seed(123)
print(NB <- NaiveBayes(PHASEN ~ ., data = B3))
predict(NB)

# S-KNN
# numerical too instable to check posteriors
SK <- sknn(PHASEN ~ ., data = B3)
predict(SK, B3)$class
SK <- sknn(PHASEN ~ ., data = B3, gamma = 10, kn = 10)
predict(SK, B3)$class


## SVMlight
## this works on Windows only, hence omitted for the meantime:
#if(class(try(system("svm_learn -?", intern = TRUE))) == "try-error"){
#    cat("SVMlight seems not to be installed, hence svmlight cannot be used", "\n",
#        "and these differences are expected")
#}else{
#    print(SVM <- svmlight(PHASEN ~ ., data = B3))
#    predict(SVM, B3)
#}

# RDA
set.seed(123)
rda(PHASEN ~ ., data = B3)
rB3 <- rda(PHASEN ~ ., gamma = 0.05, lambda = 0.2, data = B3)
print(pB3 <- predict(rB3))
cB3 <- pB3$class
pB3 <- pB3$posterior

# stepclass
set.seed(123)
print(SC <- stepclass(PHASEN ~ ., data = B3, method = "lda", criterion = "AS", output=FALSE))
print(lda(SC$formula, data = B3))

#########
# scaling

## beta scaling, e.scal
pbB3 <- b.scal(pB3, B3$PHASEN, dis = TRUE)
#betascale(pbB3)
#e.scal(pB3)
#e.scal(pB3, tc = B3$PHASEN)

# ucpm
ucpm(pB3, B3$PHASEN)
ucpm(pbB3$member, B3$PHASEN)


##############
# greedy.wilks
data(B3)
gw_obj <- greedy.wilks(PHASEN ~ ., data = B3, niveau = 0.1)
print(gw_obj, digits=4)
print(lda(gw_obj$formula, data = B3))
gw_obj2 <- greedy.wilks(B3[,-1], B3$PHASEN, niveau = 0.1)
identical(all.equal(gw_obj$results, gw_obj2$results), TRUE)


####
# nm
data(B3)
x <- nm(PHASEN ~ ., data = B3)
x$learn
x <- nm(PHASEN ~ ., data = B3, gamma = 0.1)
predict(x)$post


##########
# meclight
data(iris)
set.seed(123)
meclight.obj <- meclight(Species ~ ., data = iris)
print(meclight.obj)
set.seed(123)
meclight.obj2 <- meclight(iris[,1:4], iris[,5])
identical(all.equal(meclight.obj$Proj.matrix, meclight.obj2$Proj.matrix), TRUE)

######
# misc

# calc.trans, hmm.sop, errormatrix
print(trans.matrix <- calc.trans(B3$PHASEN))
errormatrix(B3$PHASEN, apply(pB3, 1, which.max))
print(prior.prob <- hmm.sop("2", trans.matrix, pB3))
errormatrix(B3$PHASEN, apply(prior.prob, 1, which.max))

# friedman.data
set.seed(123)
friedman.data(1, 6, 40)

# dkernel
set.seed(123)
kern <- density(rnorm(50))
x <- seq(-3, 3, len = 100)
print(dkernel(x, kern))



########
## plots

# Naive Bayes, stepclass, RDA
plot(NB)
plot(SC)
plot(rB3)

classscatter(PHASEN ~ BSP91JW + EWAJW + LSTKJW, 
    data = B3, method = "lda")

plineplot(PHASEN ~ ., data = B3, method = "lda", 
    x = "EWAJW", xlab = "EWAJW")


# quadplot
quadtrafo(pB3)
s3d <- quadplot(pB3, col = rainbow(4)[B3$PHASEN], 
        labelpch = 22:25, labelcex = 0.8,
        pch = (22:25)[apply(pB3, 1, which.max)],
        main = "RDA posterior assignments")
quadlines(centerlines(4), sp = s3d, lty = "dashed")
par("mar")


# triplot
triplot(grid = TRUE, frame = FALSE)
some.triangle <- rbind(c(0, 0.65, 0.35), c(0.53, 0.47, 0), 
                       c(0.72, 0, 0.28))[c(1:3, 1), ]
trilines(some.triangle, col = "green", pch = 16, type = "b")
triframe(label = c("left", "top", "right"), col = "blue", 
         label.col = "green3")
triperplines(1/6, 1/3, 1/2)

pred <- predict(lda(Species ~ ., data = iris), iris)
plotchar <- rep(1, 150)
plotchar[pred$class != iris$Species] <- 19
triplot(pred$posterior, label = colnames(pred$posterior), 
        main = "LDA posterior assignments", center = TRUE, 
        pch = plotchar, col = rep(c("blue", "green3", "red"), rep(50, 3)), 
        grid = TRUE)
legend(x = -0.6, y = 0.7, col = c("blue", "green3", "red"), 
    pch = 15, legend = colnames(pred$posterior))
par("mar")

# partimat
partimat(Species ~ ., data = iris, method = "lda", 
    plot.matrix = TRUE, imageplot = FALSE)

dev.off()
