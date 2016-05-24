## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">", dev = 'pdf')

## ------------------------------------------------------------------------
library(RClone)
data(posidonia)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(posidonia[1:10,1:8], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  data(zostera)
#  head(zostera)

## ---- echo = FALSE-------------------------------------------------------
data(zostera)
knitr::kable(head(zostera), align = "c")

## ------------------------------------------------------------------------
popvec <- zostera[,1] #futur vecpop
coord_zostera <- zostera[,2:3] #futur coordinates
zostera <- zostera[,4:ncol(zostera)] #dataset

zostera <- convert_GC(zostera, 3) #We used "3" because this is the length of each allele.

## ---- eval = FALSE-------------------------------------------------------
#  head(zostera)

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(zostera[1:6,1:7], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  sort_all(zostera)

## ---- eval = FALSE-------------------------------------------------------
#  #library(adegenet)
#  #with data1, a genind object from Adegenet:
#  
#  test <- genind2df(data1)
#  data2 <- convert_GC(test, 3, "/")
#  #only if yours alleles are of length "3"

## ---- eval = FALSE-------------------------------------------------------
#  list_all_tab(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  list_all_tab(haplodata, haploid = TRUE, vecpop = haplovec)

## ---- eval = FALSE-------------------------------------------------------
#  list_all_tab(zostera, vecpop = popvec)

## ------------------------------------------------------------------------
#SaintMalo

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(list_all_tab(zostera, vecpop = popvec)[[1]], align = "c")

## ------------------------------------------------------------------------
#Arcouest

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(list_all_tab(zostera, vecpop = popvec)[[2]], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(haplodata, vecpop = haplovec)

## ---- eval = FALSE-------------------------------------------------------
#  MLG_tab(zostera, vecpop = popvec)[[1]]
#  #SaintMalo

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(MLG_tab(zostera, vecpop = popvec)[[1]][1:5,], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(haplodata, haploid = TRUE, vecpop = haplovec)

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(zostera, vecpop = popvec) #on ramets
#  freq_RR(zostera, vecpop = popvec, genet = TRUE) #on genets
#  freq_RR(zostera, vecpop = popvec, RR = TRUE) #Round-Robin methods

## ---- eval = FALSE-------------------------------------------------------
#  freq_RR(zostera, vecpop = popvec)[[1]]
#  #SaintMalo

## ---- echo = FALSE-------------------------------------------------------
res <- cbind(freq_RR(zostera, vecpop = popvec)[[1]], freq_RR(zostera, vecpop = popvec, genet = TRUE)[[1]][,3], freq_RR(zostera, vecpop = popvec, RR = TRUE)[[1]][,3])[1:7,]
colnames(res)[3:5] <- c("freq_ramet", "freq_genet", "freq_RR")
knitr::kable(res, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(zostera, vecpop = popvec, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(haplodata, haploid = TRUE, vecpop = haplovec, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_loci(zostera, vecpop = popvec, nbrepeat = 1000, He = TRUE) #with He results
#  sample_loci(zostera, vecpop = popvec, nbrepeat = 1000, graph = TRUE) #graph displayed
#  sample_loci(zostera, vecpop = popvec, nbrepeat = 1000, bar = TRUE)
#  												#progression bar, could be time consuming
#  sample_loci(zostera, vecpop = popvec, nbrepeat = 1000, export = TRUE)
#  												#graph export in .eps

## ---- eval = FALSE-------------------------------------------------------
#  res <- sample_loci(zostera, vecpop = popvec, nbrepeat = 1000, He = TRUE)
#  names(res)

## ---- echo = FALSE-------------------------------------------------------
data(resvigncont2)
names(resvigncont2$res2_SU1)

## ---- eval = FALSE-------------------------------------------------------
#  names(res$SaintMalo)

## ---- eval = FALSE-------------------------------------------------------
#  names(resvigncont2$res2_SU1$SaintMalo)

## ---- eval = FALSE-------------------------------------------------------
#  #Results: MLG
#  res$Arcouest$res_MLG

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2_SU1$Arcouest$res_MLG, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  #Results: alleles
#  res$Arcouest$res_alleles

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(res$Arcouest$res_alleles, align = "c")

## ------------------------------------------------------------------------
#Results: raw data
#res$Arcouest$raw_He
#res$Arcouest$raw_MLG
#res$Arcouest$raw_all

## ---- eval = FALSE-------------------------------------------------------
#  boxplot(res$SaintMalo$raw_MLG, main = "Genotype accumulation curve",
#  	xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes")

## ---- echo = FALSE-------------------------------------------------------
boxplot(resvigncont2$res2_SU1$SaintMalo$raw_MLG, main = "Genotype accumulation curve", xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes") 

## ---- eval = FALSE-------------------------------------------------------
#  boxplot(res$Arcouest$raw_MLG, main = "Genotype accumulation curve",
#  	xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes")

## ---- echo = FALSE-------------------------------------------------------
boxplot(resvigncont2$res2_SU1$Arcouest$raw_MLG, main = "Genotype accumulation curve", xlab = "Number of loci sampled", ylab = "Number of multilocus genotypes") 

## ---- eval = FALSE-------------------------------------------------------
#  sample_units(zostera, vecpop = popvec, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  sample_units(haplodata, haploid = TRUE, vecpop = haplovec, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  pgen(zostera, vecpop = popvec)
#  data(factoR) #for psex
#  psex(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  pgen(haplodata, haploid = TRUE, vecpop = haplovec)
#  data(factoR) #for psex
#  psex(haplodata, haploid = TRUE, vecpop = haplovec)

## ---- eval = FALSE-------------------------------------------------------
#  #allelic frequencies computation:
#  psex(zostera, vecpop = popvec) #psex on ramets
#  psex(zostera, vecpop = popvec, genet = TRUE) #psex on genets
#  psex(zostera, vecpop = popvec, RR = TRUE) #psex with Round-Robin method
#  #psex computation
#  psex(zostera, vecpop = popvec) #psex with one psex per replica
#  psex(zostera, vecpop = popvec, MLGsim = TRUE) #psex MLGsim method
#  #pvalues:
#  psex(zostera, vecpop = popvec, nbrepeat = 100) #with p-values
#  psex(zostera, vecpop = popvec, nbrepeat = 1000, bar = TRUE)
#  										#with p-values and a progression bar

## ---- eval = FALSE-------------------------------------------------------
#  data(factoR)
#  res <- psex(zostera, vecpop = popvec, RR = TRUE, nbrepeat = 1000)
#  res$Arcouest[[1]]
#  #if nbrepeat != 0, res contains a table of psex values and a vector of sim-psex values

## ---- echo = FALSE-------------------------------------------------------
data(factoR)
knitr::kable(resvigncont2$res2_PS2, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  res$Arcouest[[2]] #a part of sim-psex values

## ---- echo = FALSE-------------------------------------------------------
resvigncont2$res2_PS1$Arcouest[[2]][1:10]

## ---- eval = FALSE-------------------------------------------------------
#  Fis(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  Fis(zostera, vecpop = popvec) #Fis on ramets
#  Fis(zostera, vecpop = popvec, genet = TRUE) #Fis on genets
#  Fis(zostera, vecpop = popvec, RR = TRUE) #Fis with Round-Robin methods
#  #RR = TRUE contains two results : a table with allelic frequencies
#  							 #and a table with Fis results

## ---- eval = FALSE-------------------------------------------------------
#  Fis(zostera, vecpop = popvec, RR = TRUE)$Arcouest[[2]]

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(Fis(zostera, vecpop = popvec, RR = TRUE)$Arcouest[[2]], align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  pgen_Fis(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  #allelic frequencies:
#  psex_Fis(zostera, vecpop = popvec) #psex Fis on ramets
#  psex_Fis(zostera, vecpop = popvec, genet = TRUE) #psex Fis on genets
#  psex_Fis(zostera, vecpop = popvec, RR = TRUE) #psex Fis with Round-Robin method
#  #psex computation
#  psex_Fis(zostera, vecpop = popvec) #psex Fis, one for each replica
#  psex_Fis(zostera, vecpop = popvec, MLGsim = TRUE) #psex Fis with MLGsim method
#  #pvalues
#  psex_Fis(zostera, vecpop = popvec, nbrepeat = 100) #with p-values
#  psex_Fis(zostera, vecpop = popvec, nbrepeat = 1000, bar = TRUE)
#  											#with p-values and a progression bar

## ---- eval = FALSE-------------------------------------------------------
#  data(factoR)
#  res <- psex_Fis(zostera, vecpop = popvec, RR = TRUE, nbrepeat = 1000)
#  res$Arcouest[[1]]
#  #if nbrepeat != 0, res contains a table of psex values and a vector of sim-psex Fis values

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2_PS4, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  res$Arcouest[[2]] #a part of sim psex Fis values

## ---- echo = FALSE-------------------------------------------------------
resvigncont2$res2_PS3$Arcouest[[2]][1:10]

## ------------------------------------------------------------------------
data(popsim)
vecsim <- c(rep(1,50), rep(2,50))

## ---- eval = FALSE-------------------------------------------------------
#  #genetic distances computation, distance on allele differences:
#  respop <- genet_dist(popsim, vecpop = vecsim)
#  ressim <- genet_dist_sim(popsim, vecpop = vecsim , nbrepeat = 1000)
#  								#theoretical distribution: sexual reproduction
#  ressimWS <- genet_dist_sim(popsim, vecpop = vecsim , genet = TRUE, nbrepeat = 1000)
#  															#idem, without selfing

## ---- echo = FALSE-------------------------------------------------------
respop <- resvigncont2$respop
ressim <- resvigncont2$ressim
ressimWS <- resvigncont2$ressimWS

## ---- fig.width = 10, fig.height = 8-------------------------------------
#graph prep.:
#first pop: 
p1 <- hist(respop[[1]]$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), 
			breaks = seq(0, max(respop[[1]]$distance_matrix)+1, 1), 
			main = "pop_1_sim", xlab = "")
p2 <- hist(ressim[[1]]$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), 
			breaks = seq(0, max(ressim[[1]]$distance_matrix)+1, 1), 
			main = "pop_1_SR", xlab = "")
p3 <- hist(ressimWS[[1]]$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3), 
			breaks = seq(0, max(ressimWS[[1]]$distance_matrix)+1, 1), 
			main = "pop_1_SRWS", xlab = "")
limx <- max(max(respop[[1]]$distance_matrix), max(ressim[[1]]$distance_matrix), 
		max(ressimWS[[1]]$distance_matrix))

#graph superposition: 
plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), 
		main = paste("pop", unique(vecsim)[[1]], sep = "_"), 
		xlab = "Genetic distances")
plot(p2, col = rgb(0.7,0.9,1,0.5), freq = FALSE, add = TRUE)
plot(p3, col = rgb(0.9,0.5,1,0.3), freq = FALSE, add = TRUE)

#adding a legend:
leg.txt <- c("original data","simulated data", "without selfing")
col <- c(rgb(0,0.4,1,1), rgb(0.7,0.9,1,0.5), rgb(0.9,0.5,1,0.3))
legend("top", fill = col, leg.txt, plot = TRUE, bty = "o", box.lwd = 1.5, 
bg = "white")

#second pop:
p <- 2 #useful if several populations: just change *p* and run lines

p1 <- hist(respop[[p]]$distance_matrix, freq = FALSE, col = rgb(0,0.4,1,1), 
			breaks = seq(0, max(respop[[p]]$distance_matrix)+1, 1), 
			main = paste("pop", p, sep = "_"), xlab = "")
p2 <- hist(ressim[[p]]$distance_matrix, freq = FALSE, col = rgb(0.7,0.9,1,0.5), 
			breaks = seq(0, max(ressim[[p]]$distance_matrix)+1, 1), 
			main = paste("pop", p, sep = "_"), xlab = "")
p3 <- hist(ressimWS[[p]]$distance_matrix, freq = FALSE, col = rgb(0.9,0.5,1,0.3), 
			breaks = seq(0, max(ressimWS[[p]]$distance_matrix)+1, 1), 
			main = paste("pop", p, sep = "_"), xlab = "")
limx <- max(max(respop[[p]]$distance_matrix), max(ressim[[p]]$distance_matrix), 
			max(ressimWS[[p]]$distance_matrix))

#graph superposition: 
plot(p1, col = rgb(0,0.4,1,1), freq = FALSE, xlim = c(0,limx), 
		main = paste("pop", unique(vecsim)[[p]], sep = "_"), 
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
table(respop[[1]]$distance_matrix)
#alpha2 = 3

## ------------------------------------------------------------------------
#creating MLL list:
MLLlist <- MLL_generator(popsim, vecpop = vecsim, alpha2 = c(3,0))
##This will create a list of MLL (alpha2 = 3) and MLG (alpha2 = 0) !

#or
res <- genet_dist(popsim, vecpop = vecsim, alpha2 = c(3,0))
MLLlist <- MLL_generator2(list(res[[1]]$potential_clones, 
	res[[2]]$potential_clones), MLG_list(popsim, vecpop = vecsim), vecpop = vecsim)

## ---- eval = FALSE-------------------------------------------------------
#  respop <- genet_dist(haplodata, haploid = TRUE, vecpop = vechaplo)
#  ressim <- genet_dist_sim(haplodata, haploid = TRUE, vecpop = vechaplo,
#  							nbrepeat = 1000)
#  MLLlist <- MLL_generator(haplodata, haploid = TRUE, vecpop = vechaplo,
#  							alpha2 = c(3,0))
#  #or
#  res <- genet_dist(haplodata, haploid = TRUE, vecpop = vechaplo, alpha2 = c(3,0))
#  MLLlist <- MLL_generator2(list(res[[1]]$potential_clones, res[[2]]$potential_clones),
#  			haploid = TRUE, MLG_list(haplodata, vecpop = vechaplo), vecpop = vechaplo)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(popsim, vecpop = vecsim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(haplodata, vecpop = vechaplo)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_index(zostera, vecpop = popvec)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(clonal_index(zostera, vecpop = popvec), align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(popsim, vecpop = vecsim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(haplodata, vecpop = vechaplo)

## ---- eval = FALSE-------------------------------------------------------
#  Pareto_index(zostera, vecpop = popvec, graph = TRUE) #classic graphic
#  Pareto_index(zostera, vecpop = popvec, legends = 2, export = TRUE)
#  														#export option
#  Pareto_index(zostera, vecpop = popvec, full = TRUE) #all results

## ---- eval = FALSE-------------------------------------------------------
#  res <- Pareto_index(zostera, vecpop = popvec, full = TRUE, graph = TRUE, legends = 2)

## ---- echo = FALSE-------------------------------------------------------
require(RClone)
resz <- split(zostera, popvec)
res1 <- Pareto_index(resz[[1]],  full = TRUE, graph = TRUE, legends = 2)
res2 <- Pareto_index(resz[[2]],  full = TRUE, graph = TRUE, legends = 2)

## ---- eval = FALSE-------------------------------------------------------
#  names(res$SaintMalo)

## ---- echo = FALSE-------------------------------------------------------
names(res1)

## ---- eval = FALSE-------------------------------------------------------
#  res$SaintMalo$Pareto

## ---- echo = FALSE-------------------------------------------------------
res1$Pareto

## ---- eval = FALSE-------------------------------------------------------
#  res$SaintMalo$c_Pareto

## ---- echo = FALSE-------------------------------------------------------
res1$c_Pareto

## ---- eval = FALSE-------------------------------------------------------
#  #res$SaintMalo$regression_results
#  #res$SaintMalo$coords_Pareto #points coordinates

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(popsim, coords = coord_sim, Loiselle = TRUE, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  autocorrelation(haplodata, haploid = TRUE, coords = coord_haplo, Loiselle = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  #kinship distances:
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE)
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Ritland = TRUE)
#  
#  #ramets/genets methods:
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE)
#  																			#ramets
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec,, Loiselle = TRUE,
#  					genet = TRUE, central_coords = TRUE)
#  											#genets, central coordinates of each MLG
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE,
#  					genet = TRUE, random_unit = TRUE)
#  													#genets, one random unit per MLG
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE,
#  					genet = TRUE, weighted = TRUE)
#  											#genets, with weighted matrix on kinships
#  
#  #distance classes construction:
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE)
#  																#10 equidistant classes
#  distvec <- c(0,10,15,20,30,50,70,76.0411074)
#  									#with 0, min distance and 76.0411074, max distance
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE,
#  					vecdist = distvec) #custom distance vector
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE,
#  					class1 = TRUE, d = 7) #7 equidistant classes
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Loiselle = TRUE,
#  					class2 = TRUE, d = 7)
#  							#7 distance classes with the same number of units in each
#  
#  #graph options:
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Ritland = TRUE,
#  														graph = TRUE) #displays graph
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Ritland = TRUE,
#  														export = TRUE) #export graph
#  
#  #pvalues computation
#  autocorrelation(zostera, coords = coord_zostera, vecpop = popvec, Ritland = TRUE,
#  																	nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  res <- autocorrelation(zostera, coords = coord_zostera, vecpop = popvec,
#  										Ritland = TRUE, nbrepeat = 1000, graph = TRUE)

## ---- echo = FALSE-------------------------------------------------------
plot(resvigncont2$res2auto1$Main_results[,3], resvigncont2$res2auto1$Main_results[,6], main = "Spatial aucorrelation analysis",
ylim = c(-0.2,0.2), type = "l", xlab = "Spatial distance", ylab = "Coancestry (Fij)")
points(resvigncont2$res2auto1$Main_results[,3], resvigncont2$res2auto1$Main_results[,6], pch = 20)
abline(h = 0, lty = 3)

plot(resvigncont2$res2auto2$Main_results[,3], resvigncont2$res2auto2$Main_results[,6], main = "Spatial aucorrelation analysis",
ylim = c(-0.2,0.2), type = "l", xlab = "Spatial distance", ylab = "Coancestry (Fij)")
points(resvigncont2$res2auto2$Main_results[,3], resvigncont2$res2auto2$Main_results[,6], pch = 20)
abline(h = 0, lty = 3)

## ---- eval = FALSE-------------------------------------------------------
#  names(res$Arcouest)

## ---- echo = FALSE-------------------------------------------------------
names(resvigncont2$res2auto2)

## ---- eval = FALSE-------------------------------------------------------
#  res$Arcouest$Main_results #enables graph reproduction

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2auto2$Main_results, align = "c")

## ---- eval = FALSE-------------------------------------------------------
#  apply(res$Arcouest$Main_results, 2, mean)[6] #mean Fij

## ---- echo = FALSE-------------------------------------------------------
apply(resvigncont2$res2auto2$Main_results, 2, mean)[6]

## ---- eval = FALSE-------------------------------------------------------
#  res$Arcouest$Slope_and_Sp_index #gives b and Sp indices

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2auto2$Slope_and_Sp_index, align = "c")

## ------------------------------------------------------------------------
#raw data:
#res$Arcouest$Slope_resample
#res$Arcouest$Kinship_resample
#res$Arcouest$Matrix_kinship_results
#res$Arcouest$Class_kinship_results 
#res$Arcouest$Class_distance_results

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(zostera, coords = coord_zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(haplodata, haploid = TRUE, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  clonal_sub(posidonia, coords = coord_posidonia) #basic, with 10 equidistant classes
#  distvec <- c(0,10,15,20,30,50,70,76.0411074)
#  								#with 0, min distance and 76.0411074, max distance
#  clonal_sub(zostera, coords = coord_zostera, vecpop = popvec, vecdist = distvec)
#  															#custom distance classes
#  clonal_sub(zostera, coords = coord_zostera, vecpop = popvec, class1 = TRUE, d = 7)
#  																#7 equidistant classes
#  clonal_sub(zostera, coords = coord_zostera, vecpop = popvec, class1 = TRUE, d = 7)
#  							#7 distance classes with the same number of units in each

## ------------------------------------------------------------------------
res <- clonal_sub(zostera, coords = coord_zostera, vecpop = popvec)
res$Arcouest[[1]] #Global clonal subrange

## ---- eval = FALSE-------------------------------------------------------
#  res$Arcouest$clonal_sub_tab  #details per class

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(res$Arcouest$clonal_sub_tab, align ="c")

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(zostera, coords = coord_zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(haplodata, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  agg_index(zostera, coords = coord_zostera, vecpop = popvec, nbrepeat = 100)
#  															#pvalue computation
#  agg_index(zostera, coords = coord_zostera, vecpop = popvec, nbrepeat = 1000,
#  											bar = TRUE) #could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  res <- agg_index(zostera, coords = coord_zostera, vecpop = popvec, nbrepeat = 1000)

## ---- eval = FALSE-------------------------------------------------------
#  res$SaintMalo$results #Aggregation index

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2_agg$SaintMalo$results, align = "c")

## ------------------------------------------------------------------------
#res$SaintMalo$simulation #vector of sim aggregation index

## ---- eval = FALSE-------------------------------------------------------
#  #for zostera, centers of quadra is at 15,10
#  edge_effect(zostera, coords = coord_zostera, vecpop = popvec,
#  				center = rep(c(15,10),2))

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(popsim, coords = coord_sim, center = rep(c(15,10),2), listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(haplodata, coords = coord_haplo, center = rep(c(15,10),2))

## ---- eval = FALSE-------------------------------------------------------
#  edge_effect(zostera, coords = coord_zostera, vecpop = popvec, center = rep(c(15,10),2),
#  													nbrepeat = 100) #pvalue computation
#  edge_effect(zostera, coords = coord_zostera, vecpop = popvec, center = rep(c(15,10),2),
#  									nbrepeat = 1000, bar = TRUE) #could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  res <- edge_effect(zostera, coords = coord_zostera, vecpop = popvec,
#  		center = rep(c(15,10),2), nbrepeat = 100) #better put 1000 nbrepeat at least

## ---- eval = FALSE-------------------------------------------------------
#  res$SaintMalo$results #Aggregation index

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(resvigncont2$res2_ee$SaintMalo$results, align = "c")

## ------------------------------------------------------------------------
#res$SaintMalo$simulation #vector of sim aggregation index

## ---- eval = FALSE-------------------------------------------------------
#  genclone(zostera, coords = coord_zostera, vecpop = popvec)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(popsim, coords = coord_sim, listMLL = MLLlist)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(haplodata, haploid = TRUE, coords = coord_haplo)

## ---- eval = FALSE-------------------------------------------------------
#  genclone(zostera, coords = coord_zostera, vecpop = popvec, nbrepeat = 100) #pvalues
#  genclone(zostera, coords = coord_zostera, vecpop = popvec, nbrepeat = 1000, bar = TRUE)
#  																#could be time consuming

## ---- eval = FALSE-------------------------------------------------------
#  genclone(zostera, coords = coord_zostera, vecpop = popvec)

## ---- echo = FALSE, results = 'asis'-------------------------------------
knitr::kable(resvigncont2$res2_gen[,1:10], longtable = TRUE, align = "c")
knitr::kable(resvigncont2$res2_gen[,11:17], longtable = TRUE, align = "c")
knitr::kable(resvigncont2$res2_gen[,18:24], longtable = TRUE, align = "c")

