## ----paquete,eval=FALSE--------------------------------------------------
#  library("Anthropometry")

## ----trimowa1,eval=FALSE,tidy=FALSE--------------------------------------
#  dataTrimowa <- sampleSpanishSurvey
#  numVar <- dim(dataTrimowa)[2]
#  bust <- dataTrimowa$bust
#  bustSizes <- bustSizesStandard(seq(74, 102, 4), seq(107, 131, 6))

## ----trimowa2,eval=FALSE,tidy=FALSE--------------------------------------
#  orness <- 0.7
#  weightsTrimowa <- weightsMixtureUB(orness, numVar)

## ----trimowa3,eval=FALSE,tidy=FALSE--------------------------------------
#  numClust <- 3 ; alpha <- 0.01 ; niter <- 10 ; algSteps <- 7
#  ah <- c(23, 28, 20, 25, 25)
#  
#  set.seed(2014)
#  numSizes <- bustSizes$nsizes - 1
#  res_trimowa <- computSizesTrimowa(dataTrimowa, bust, bustSizes$bustCirc,
#                                    numSizes, weightsTrimowa, numClust,
#                                    alpha, niter, algSteps, ah, FALSE)

## ----trimowa4,eval=FALSE,tidy=FALSE--------------------------------------
#  prototypes <- anthrCases(res_trimowa, numSizes)

## ----trimowa5,eval=FALSE,tidy=FALSE--------------------------------------
#  bustVariable <- "bust"
#  xlim <- c(72, 132)
#  color <- c("black", "red", "green", "blue", "cyan", "brown", "gray",
#             "deeppink3", "orange", "springgreen4", "khaki3", "steelblue1")
#  variable <- "necktoground"
#  ylim <- c(116, 156)
#  title <- "Prototypes \n bust vs neck to ground"
#  plotPrototypes(dataTrimowa, prototypes, numSizes, bustVariable,
#                 variable, color, xlim, ylim, title, FALSE)
#  plotPrototypes(dataTrimowa, prototypes, numSizes, bustVariable,
#                 variable, color, xlim, ylim, title, TRUE)

## ----TDDclust,eval=FALSE,tidy=FALSE--------------------------------------
#  dataTDDcl <- sampleSpanishSurvey[1 : 25, c(2, 3, 5)]
#  dataTDDcl_aux <- sampleSpanishSurvey[1 : 25, c(2, 3, 5)]

## ----TDDclust2,eval=FALSE,tidy=FALSE-------------------------------------
#  numClust <- 3 ; alpha <- 0.01 ; lambda <- 0.5 ; niter <- 5
#  Th <- 0 ; T0 <- 0 ; simAnn <- 0.9
#  
#  set.seed(2014)
#  res_TDDcl <- TDDclust(dataTDDcl, numClust, lambda, Th, niter, T0, simAnn,
#                        alpha, dataTDDcl_aux, verbose = FALSE)

## ----TDDclust3,eval=FALSE,tidy=FALSE-------------------------------------
#  table(res_TDDcl$NN[1,])
#  #1  2  3
#  #5 10  9
#  res_TDDcl$Cost
#  #[1] 0.3717631
#  res_TDDcl$klBest
#  #[1] 3

## ----TDDclust4,eval=FALSE,tidy=FALSE-------------------------------------
#  prototypes <- anthrCases(res_TDDcl)
#  trimmed <- trimmOutl(res_TDDcl)

## ----hipam,eval=FALSE,tidy=FALSE-----------------------------------------
#  dataHipam <- sampleSpanishSurvey
#  bust <- dataHipam$bust
#  bustSizes <- bustSizesStandard(seq(74, 102, 4), seq(107, 131, 6))

## ----hipam2,eval=FALSE,tidy=FALSE----------------------------------------
#  type <- "IMO"
#  maxsplit <- 5 ; orness <- 0.7
#  ah <- c(23, 28, 20, 25, 25)
#  
#  set.seed(2013)
#  numSizes <- bustSizes$nsizes - 1
#  res_hipam <- computSizesHipamAnthropom(dataHipam, bust, bustSizes$bustCirc,
#                                         numSizes, maxsplit, orness, type,
#                                         ah, FALSE)

## ----hipam3,eval=FALSE,tidy=FALSE----------------------------------------
#  fitmodels <- anthrCases(res_hipam, numSizes)
#  outliers <- trimmOutl(res_hipam, numSizes)

## ----hipam4,eval=FALSE,tidy=FALSE----------------------------------------
#  bustVariable <- "bust"
#  xlim <- c(72, 132)
#  color <- c("black", "red", "green", "blue", "cyan", "brown", "gray",
#             "deeppink3", "orange", "springgreen4", "khaki3", "steelblue1")
#  variable <- "hip"
#  ylim <- c(83, 153)
#  title <- "Fit models HIPAM_IMO \n bust vs hip"
#  title_outl <- "Outlier women HIPAM_IMO \n bust vs hip"
#  plotPrototypes(dataHipam, fitmodels, numSizes, bustVariable,
#                 variable, color, xlim, ylim, title, FALSE)
#  plotTrimmOutl(dataHipam, outliers, numSizes, bustVariable,
#                variable, color, xlim, ylim, title_outl)

## ----ssa,eval=FALSE,tidy=FALSE-------------------------------------------
#  landmarksNoNa <- na.exclude(landmarksSampleSpaSurv)
#  numLandmarks <- (dim(landmarksNoNa)[2]) / 3
#  landmarksNoNa_First50 <- landmarksNoNa[1 : 50, ]
#  numIndiv <- dim(landmarksNoNa_First50)[1]

## ----ssa1,eval=FALSE,tidy=FALSE------------------------------------------
#  array3D <- array3Dlandm(numLandmarks, numIndiv, landmarksNoNa_First50)

## ----ssa2,eval=FALSE,tidy=FALSE------------------------------------------
#  numClust <- 3 ; alpha <- 0.01 ; algSteps <- 5
#  niter <- 5 ; stopCr <- 0.0001

## ----ssa22,eval=FALSE,tidy=FALSE-----------------------------------------
#  set.seed(2013)
#  res_kmProc <- trimmedLloydShapes(array3D, numIndiv, alpha, numClust,
#                                       algSteps, niter, stopCr,
#                                       verbose = FALSE)

## ----ssa3,eval=FALSE,tidy=FALSE------------------------------------------
#  clust_kmProc <- res_kmProc$asig
#  table(clust_kmProc)
#  #1  2  3
#  #19 18 12

## ----ssa4,eval=FALSE,tidy=FALSE------------------------------------------
#  prototypes <- anthrCases(res_kmProc)
#  trimmed <- trimmOutl(res_kmProc)

## ----ssa5,eval=FALSE,tidy=FALSE------------------------------------------
#  data_First50 <- sampleSpanishSurvey[1 : 50, ]
#  data_First50_notrimm <- data_First50[-trimmed, ]
#  boxplot(data_First50_notrimm$necktoground ~ as.factor(clust_kmProc),
#          main = "Neck to ground")

## ----ssa6,eval=FALSE,tidy=FALSE------------------------------------------
#  projShapes(1, array3D, clust_kmProc, prototypes)
#  legend("topleft", c("Registrated data", "Mean shape"),
#                  pch = 1, col = 1:2, text.col = 1:2)
#  title("Procrustes registrated data for cluster 1 \n
#                  with its mean shape superimposed", sub = "Plane xy")

## ----AA,eval=FALSE,tidy=FALSE--------------------------------------------
#  USAFSurvey_First50 <- USAFSurvey[1 : 50, ]
#  variabl_sel <- c(48, 40, 39, 33, 34, 36)
#  USAFSurvey_First50_inch <- USAFSurvey_First50[,variabl_sel] / (10 * 2.54)
#  USAFSurvey_preproc <- preprocessing(data = USAFSurvey_First50_inch,
#                                      stand = TRUE, percAccomm = 0.95,
#                                      mahal= TRUE)

## ----AA3,eval=FALSE,tidy=FALSE-------------------------------------------
#  set.seed(2010)
#  numArch <- 10 ; numRep <- 20
#  lass <- stepArchetypesRawData(data = USAFSurvey_preproc$data,
#                                numArch=1:numArch, numRep = numRep,
#                                verbose = FALSE)
#  screeplot(lass)

## ----AA4,eval=FALSE,tidy=FALSE-------------------------------------------
#  numArchoid <- 3
#  res_archoids_ns <- archetypoids(numArchoid, USAFSurvey_preproc$data,
#                                  huge = 200, step = FALSE, ArchObj = lass,
#                                  nearest = "cand_ns" , sequ = TRUE)
#  res_archoids_alpha <- archetypoids(numArchoid, USAFSurvey_preproc$data,
#                                     huge = 200, step = FALSE, ArchObj = lass,
#                                     nearest = "cand_alpha", sequ = TRUE)
#  res_archoids_beta <- archetypoids(numArchoid, USAFSurvey_preproc$data,
#                                    huge = 200, step = FALSE, ArchObj = lass,
#                                    nearest = "cand_beta", sequ = TRUE)
#  
#  boundaries_ns <- anthrCases(res_archoids_ns)
#  boundaries_alpha <- anthrCases(res_archoids_alpha)
#  boundaries_beta <- anthrCases(res_archoids_beta)

## ----AA5,eval=FALSE,tidy=FALSE-------------------------------------------
#  matPer <- matPercs(boundaries_ns, USAFSurvey_preproc$data)

## ----AA6,eval=FALSE,tidy=FALSE-------------------------------------------
#  barplot(matPer, beside = TRUE, main = paste(numArchoid,
#                                              " archetypoids", sep = ""),
#          ylim = c(0, 100), ylab = "Percentile",
#          xlab = "Each bar is related to each anthropometric
#                 variable selected")

