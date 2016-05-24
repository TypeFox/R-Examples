## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">", dev = 'pdf')

## ------------------------------------------------------------------------
library(RClone)
data(posidonia)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(posidonia[1:10,1:8], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  data(posidonia)
#  
#  sort_all(posidonia)

## ------------------------------------------------------------------------
#Let's create your example table:
test <- matrix("232/231", ncol = 2, nrow = 2)
colnames(test) <- paste("locus", 1:2, sep = "_")


#Use :
data1 <- convert_GC(as.data.frame(test), 3, "/")

## ---- eval = FALSE-------------------------------------------------------
#  data1

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(data1, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  #library(adegenet)
#  #with data1, a genind object from Adegenet:
#  
#  test <- genind2df(data1)
#  data2 <- convert_GC(test, 3, "/")
#  #only if yours alleles are of length "3"

## ---- eval = FALSE-------------------------------------------------------
#  data(infile)
#  
#  #This is nearly a GenClone file, type:
#  write.table(infile, "infile.csv", col.names = FALSE, row.names = FALSE, sep = ";")
#  
#  #Now you have a formatted GenClone file:
#  res <- transcript_GC("infile.csv", ";", 2, 7, 3)
#  posidonia <- res$data_genet
#  coord_posidonia <- res$data_coord

## ---- eval = FALSE-------------------------------------------------------
#  data(posidonia)
#  
#  list_all_tab(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  list_all_tab(haplodata, haploid = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  list_all_tab(posidonia)

## ---- echo = FALSE-------------------------------------------------------
data(posidonia)
knitr::kable(list_all_tab(posidonia), align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(haplodata)

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(posidonia)

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(MLG_tab(posidonia)[1:5,], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(haplodata, haploid = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(posidonia) #on ramets
#  freq_RR(posidonia, genet = TRUE) #on genets
#  freq_RR(posidonia, RR = TRUE) #Round-Robin methods

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(posidonia)

## ---- echo = FALSE-------------------------------------------------------
res <- cbind(freq_RR(posidonia), freq_RR(posidonia, genet = TRUE)[,3], freq_RR(posidonia, RR = TRUE)[,3])[1:7,]
colnames(res)[3:5] <- c("freq_ramet", "freq_genet", "freq_RR")
knitr::kable(res, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(posidonia, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(haplodata, haploid = TRUE, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(posidonia, nbrepeat = 1000, He = TRUE) #with He results
#  sample_loci(posidonia, nbrepeat = 1000, graph = TRUE) #graph displayed
#  sample_loci(posidonia, nbrepeat = 1000, bar = TRUE) #progression bar
#  													#could be time consuming
#  sample_loci(posidonia, nbrepeat = 1000, export = TRUE) #graph export in .eps

## ---- eval = FALSE-------------------------------------------------------
#  res <- sample_loci(posidonia, nbrepeat = 1000, He = TRUE) #time consuming
#  names(res)

## ---- echo = FALSE-------------------------------------------------------
data(resvigncont)
names(resvigncont$resvigncont$res_SU1)

## ---- eval = FALSE-------------------------------------------------------
#  #Results: MLG
#  res$res_MLG

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$res_SU1$res_MLG, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  #Results: alleles
#  res$res_alleles

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$res_SU1$res_alleles, align = "c")

## ------------------------------------------------------------------------
#Results: raw data
#res$raw_He
#res$raw_MLG
#res$raw_all

## ---- eval = FALSE-------------------------------------------------------
#  boxplot(res$raw_MLG, main = "Genotype accumulation curve",
#  	xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes")

## ---- echo = FALSE-------------------------------------------------------
boxplot(resvigncont$res_SU1$raw_MLG, main = "Genotype accumulation curve", xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes") 

## ---- eval = FALSE-------------------------------------------------------
#  sample_units(posidonia, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_units(haplodata, haploid = TRUE, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  pgen(posidonia)
#  data(factoR) #for psex
#  psex(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  pgen(haplodata, haploid = TRUE)
#  data(factoR) #for psex
#  psex(haplodata, haploid = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  #allelic frequencies computation:
#  psex(posidonia) #psex on ramets
#  psex(posidonia, genet = TRUE) #psex on genets
#  psex(posidonia, RR = TRUE) #psex with Round-Robin method
#  #psex computation
#  psex(posidonia) #psex with one psex per replica
#  psex(posidonia, MLGsim = TRUE) #psex MLGsim method
#  #pvalues:
#  psex(posidonia, nbrepeat = 100) #with p-values
#  psex(posidonia, nbrepeat = 1000, bar = TRUE) #with p-values and a progression bar

## ---- eval = FALSE-------------------------------------------------------
#  data(factoR)
#  res <- psex(posidonia, RR = TRUE, nbrepeat = 1000)
#  res[[1]] #if nbrepeat != 0, res contains a table of psex values
#  									#and a vector of sim-psex values

## ---- echo = FALSE-------------------------------------------------------
data(factoR)

knitr::kable(resvigncont$res_PS2, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  res[[2]] #sim psex values

## ---- echo = FALSE-------------------------------------------------------
resvigncont$res_PS1[[2]]

## ---- eval = FALSE-------------------------------------------------------
#  Fis(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  Fis(posidonia) #Fis on ramets
#  Fis(posidonia, genet = TRUE) #Fis on genets
#  Fis(posidonia, RR = TRUE) #Fis with Round-Robin methods
#  #RR = TRUE contains two results : a table with allelic frequencies
#  							 #and a table with Fis results

## ---- eval = FALSE-------------------------------------------------------
#  Fis(posidonia, RR = TRUE)[[2]]

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(Fis(posidonia, RR = TRUE)[[2]], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  pgen_Fis(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  #allelic frequencies:
#  psex_Fis(posidonia) #psex Fis on ramets
#  psex_Fis(posidonia, genet = TRUE) #psex Fis on genets
#  psex_Fis(posidonia, RR = TRUE) #psex Fis with Round-Robin method
#  #psex computation
#  psex_Fis(posidonia) #psex Fis, one for each replica
#  psex_Fis(posidonia, MLGsim = TRUE) #psex Fis with MLGsim method
#  #pvalues
#  psex_Fis(posidonia, nbrepeat = 100) #with p-values
#  psex_Fis(posidonia, nbrepeat = 1000, bar = TRUE) #with p-values and a progression bar

## ---- eval = FALSE-------------------------------------------------------
#  data(factoR)
#  res <- psex_Fis(posidonia, RR = TRUE, nbrepeat = 1000)
#  res[[1]]
#  #if nbrepeat != 0, res contains a table of psex values
#  						   #and a vector of sim-psex Fis values

## ---- echo = FALSE-------------------------------------------------------

knitr::kable(resvigncont$res_PS4, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  res[[2]] #sim psex Fis values

## ---- echo = FALSE-------------------------------------------------------
resvigncont$res_PS3[[2]]

## ---- eval = FALSE-------------------------------------------------------
#  data(popsim)
#  
#  #genetic distances computation, distance on allele differences:
#  respop <- genet_dist(popsim)
#  ressim <- genet_dist_sim(popsim, nbrepeat = 1000) #theoretical distribution:
#  												  #sexual reproduction
#  ressimWS <- genet_dist_sim(popsim, genet = TRUE, nbrepeat = 1000) #idem, without selfing

## ---- echo = FALSE-------------------------------------------------------
data(popsim)
respop <- resvigncont$respop
ressim <- resvigncont$ressim
ressimWS <- resvigncont$ressimWS

## ---- fig.width = 10, fig.height = 8-------------------------------------
#graph prep.:
p1 <- hist(respop$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), main = "popsim", 
			xlab = "Genetic distances", breaks = seq(0, max(respop$distance_matrix)+1, 1))
p2 <- hist(ressim$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), main = "popSR", 
			xlab = "Genetic distances", breaks = seq(0, max(ressim$distance_matrix)+1, 1))
p3 <- hist(ressimWS$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3), 
			main = "popSRWS", xlab = "Genetic distances", 
			breaks = seq(0, max(ressimWS$distance_matrix)+1, 1))
limx <- max(max(respop$distance_matrix), max(ressim$distance_matrix), 
			max(ressimWS$distance_matrix))

#graph superposition: 
plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), main = "", 
		xlab = "Genetic distances")
plot(p2, col = rgb(0.7,0.9,1,0.5), freq = FALSE, add = TRUE)
plot(p3, col = rgb(0.9,0.5,1,0.3), freq = FALSE, add = TRUE)

#adding a legend:
leg.txt <- c("original data","simulated data", "without selfing")
col <- c(rgb(0,0.4,1,1), rgb(0.7,0.9,1,0.5), rgb(0.9,0.5,1,0.3))
legend("top", fill = col, leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, 
bg = "white")

## ------------------------------------------------------------------------
#determining alpha2
table(respop$distance_matrix)
#alpha2 = 4

## ------------------------------------------------------------------------
#creating MLL list:
MLLlist <- MLL_generator(popsim, alpha2 = 4)
#or
res <- genet_dist(popsim, alpha2 = 4)
MLLlist <- MLL_generator2(res$potential_clones, MLG_list(popsim))

## ---- eval = FALSE-------------------------------------------------------
#  respop <- genet_dist(haplodata, haploid = TRUE)
#  ressim <- genet_dist_sim(haplodata, haploid = TRUE, nbrepeat = 1000)
#  MLLlist <- MLL_generator(haplodata, haploid = TRUE, alpha2 = 4)
#  #or
#  res <- genet_dist(haplodata, haploid = TRUE, alpha2 = 4)
#  MLLlist <- MLL_generator2(res$potential_clones, haploid = TRUE, MLG_list(haplodata))

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(popsim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(haplodata)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(posidonia)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(resvigncont$rescl, align = "c")
data(coord_posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(popsim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(haplodata)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(posidonia, graph = TRUE) #classic graphic
#  Pareto_index(posidonia, legends = 2, export = TRUE) #export option
#  Pareto_index(posidonia, full = TRUE) #all results

## ------------------------------------------------------------------------
res <- Pareto_index(posidonia, full = TRUE, graph = TRUE, legends = 2)
names(res)
res$Pareto
res$c_Pareto
#res$regression_results
#res$coords_Pareto #points coordinates

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(popsim, coords = coord_sim, Loiselle = TRUE, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(haplodata, haploid = TRUE, coords = coord_haplo, Loiselle = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  data(posidonia)
#  data(coord_posidonia)
#  
#  #kinship distances:
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE)
#  autocorrelation(posidonia, coords = coord_posidonia, Ritland = TRUE)
#  
#  #ramets/genets methods:
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE) #ramets
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  					genet = TRUE, central_coords = TRUE)
#  											#genets, central coordinates of each MLG
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  				genet = TRUE, random_unit = TRUE) #genets, one random unit per MLG
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  				genet = TRUE, weighted = TRUE) #genets, with weighted matrix on kinships
#  
#  #distance classes construction:
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE)
#  													#10 equidistant classes
#  distvec <- c(0,10,15,20,30,50,70,76.0411074)
#  						#with 0, min distance and 76.0411074, max distance
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  					vecdist = distvec) #custom distance vector
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  					class1 = TRUE, d = 7) #7 equidistant classes
#  autocorrelation(posidonia, coords = coord_posidonia, Loiselle = TRUE,
#  					class2 = TRUE, d = 7)
#  					#7 distance classes with the same number of units in each
#  
#  #graph options:
#  autocorrelation(posidonia, coords = coord_posidonia, Ritland = TRUE, graph = TRUE)
#  																	#displays graph
#  autocorrelation(posidonia, coords = coord_posidonia, Ritland = TRUE, export = TRUE)
#  																	#export graph
#  
#  #pvalues computation
#  autocorrelation(posidonia, coords = coord_posidonia, Ritland = TRUE, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  res <- autocorrelation(posidonia, coords = coord_posidonia, Ritland = TRUE,
#  						nbrepeat = 1000, graph = TRUE)

## ---- echo = FALSE-------------------------------------------------------
plot(resvigncont$resauto$Main_results[,3], resvigncont$resauto$Main_results[,6], main = "Spatial aucorrelation analysis",
ylim = c(-0.2,0.2), type = "l", xlab = "Spatial distance", ylab = "Coancestry (Fij)")
points(resvigncont$resauto$Main_results[,3], resvigncont$resauto$Main_results[,6], pch = 20)
abline(h = 0, lty = 3)

## ---- eval = FALSE-------------------------------------------------------
#  names(res)

## ---- echo = FALSE-------------------------------------------------------
names(resvigncont$resauto)

## ---- eval = FALSE-------------------------------------------------------
#  res$Main_results #enables graph reproduction

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$resauto$Main_results, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  apply(res$Main_results, 2, mean)[6] #mean Fij

## ---- echo = FALSE-------------------------------------------------------
apply(resvigncont$resauto$Main_results, 2, mean)[6] #mean Fij

## ---- eval = FALSE-------------------------------------------------------
#  res$Slope_and_Sp_index #gives b and Sp indices

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$resauto$Slope_and_Sp_index, align = "c")

## ------------------------------------------------------------------------
#raw data:
#res$Slope_resample
#res$Kinship_resample
#res$Matrix_kinship_results
#res$Class_kinship_results 
#res$Class_distance_results

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(posidonia, coords = coord_posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(haplodata, haploid = TRUE, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(posidonia, coords = coord_posidonia) #basic, with 10 equidistant classes
#  distvec <- c(0,10,15,20,30,50,70,76.0411074)
#  						#with 0, min distance and 76.0411074, max distance
#  clonal_sub(posidonia, coords = coord_posidonia, vecdist = distvec)
#  												#custom distance classes
#  clonal_sub(posidonia, coords = coord_posidonia, class1 = TRUE, d = 7)
#  												#7 equidistant classes
#  clonal_sub(posidonia, coords = coord_posidonia, class1 = TRUE, d = 7)
#  				#7 distance classes with the same number of units in each

## ---- eval = FALSE-------------------------------------------------------
#  res <- clonal_sub(posidonia, coords = coord_posidonia)
#  res[[1]] #Global clonal subrange

## ---- echo = FALSE-------------------------------------------------------
resvigncont$rescs[[1]]

## ---- eval = FALSE-------------------------------------------------------
#  res$clonal_sub_tab  #details per class

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$rescs$clonal_sub_tab, align ="c")

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(posidonia, coords = coord_posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(haplodata, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(posidonia, coords = coord_posidonia, nbrepeat = 100) #pvalue computation
#  agg_index(posidonia, coords = coord_posidonia, nbrepeat = 1000, bar = TRUE)
#  															#could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  res <- agg_index(posidonia, coords = coord_posidonia, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  res$results #Aggregation index

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$resagg$results, align = "c")

## ------------------------------------------------------------------------
#res$simulation #vector of sim aggregation index

## ---- eval = FALSE-------------------------------------------------------
#  #for posidonia, center of quadra is at 40,10
#  edge_effect(posidonia, coords = coord_posidonia, center = c(40,10))

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(popsim, coords = coord_sim, center = c(40,10), listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(haplodata, coords = coord_haplo, center = c(40,10))

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(posidonia, coords = coord_posidonia, center = c(40,10), nbrepeat = 100)
#  																	#pvalue computation
#  edge_effect(posidonia, coords = coord_posidonia, center = c(40,10), nbrepeat = 1000,
#  													bar = TRUE) #could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  res <- edge_effect(posidonia, coords = coord_posidonia, center = c(40,10), nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  res$results #Aggregation index

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont$resee$results, align = "c")

## ------------------------------------------------------------------------
#res$simulation #vector of sim aggregation index

## ---- eval = FALSE-------------------------------------------------------
#  genclone(posidonia, coords = coord_posidonia)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(haplodata, haploid = TRUE, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(posidonia, coords = coord_posidonia, nbrepeat = 100) #pvalues
#  genclone(posidonia, coords = coord_posidonia, nbrepeat = 1000, bar = TRUE)
#  													#could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  genclone(posidonia, coords = coord_posidonia)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(resvigncont$resgen[,1:10], longtable = TRUE, align = "c")
knitr::kable(resvigncont$resgen[,11:17], longtable = TRUE, align = "c")
knitr::kable(resvigncont$resgen[,18:24], longtable = TRUE, align = "c")

