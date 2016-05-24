### R code from vignette source 'diveRsity.Rnw'

###################################################
### code chunk number 1: diveRsity.Rnw:66-67 (eval = FALSE)
###################################################
## divOnline()


###################################################
### code chunk number 2: diveRsity.Rnw:83-84 (eval = FALSE)
###################################################
## citation("diveRsity")


###################################################
### code chunk number 3: diveRsity.Rnw:127-128 (eval = FALSE)
###################################################
## install.packages("diveRsity")


###################################################
### code chunk number 4: diveRsity.Rnw:133-141 (eval = FALSE)
###################################################
## # install and load devtools
## install.packages("devtools")
## 
## library("devtools")
## 
## # download and install diveRsity
## 
## install_github(name = "diveRsity", subdir = "kkeenan02")


###################################################
### code chunk number 5: diveRsity.Rnw:160-161 (eval = FALSE)
###################################################
## install.packages("package_name")


###################################################
### code chunk number 6: diveRsity.Rnw:169-170 (eval = FALSE)
###################################################
## library("diveRsity")


###################################################
### code chunk number 7: diveRsity.Rnw:174-189 (eval = FALSE)
###################################################
## ?divPart
## ?inCalc
## ?readGenepop
## ?corPlot
## ?difPlot
## ?chiCalc
## ?divOnline
## ?divBasic
## ?fstOnly
## ?divRatio
## ?microPlexer
## ?arp2gen
## ?divMigrate
## ?haploDiv
## ?fastDivPart


###################################################
### code chunk number 8: diveRsity.Rnw:433-436 (eval = FALSE)
###################################################
## divPart(infile = NULL, outfile = NULL, gp = 3, pairwise = FALSE,
##         WC_Fst = FALSE, bs_locus = FALSE, bs_pairwise = FALSE, 
##         bootstraps = 0, plot = FALSE, parallel = FALSE)


###################################################
### code chunk number 9: diveRsity.Rnw:492-498
###################################################
#data(Test_data,package='diveRsity')
#library(diveRsity)
#x<-capture.output(res<-divPart(Test_data,"outtt",3,T,T,T,T,3,F,T))

load("./div_res.RData")
res <- div_results


###################################################
### code chunk number 10: diveRsity.Rnw:508-510
###################################################
options(width=50)
res$standard[c(1:10,nrow(res$standard)),] 


###################################################
### code chunk number 11: diveRsity.Rnw:540-542
###################################################
options(width=50)
res$estimate[c(1:10,nrow(res$estimate)),]


###################################################
### code chunk number 12: diveRsity.Rnw:579-595
###################################################
noquote(names(res$pairwise)[1])
noquote(res$pairwise[[1]][1:4,1:4])
noquote(names(res$pairwise)[2])
noquote(res$pairwise[[2]][1:4,1:4])
noquote(names(res$pairwise)[3])
noquote(res$pairwise[[3]][1:4,1:4])
noquote(names(res$pairwise)[4])
noquote(res$pairwise[[4]][1:4,1:4])
noquote(names(res$pairwise)[5])
noquote(res$pairwise[[5]][1:4,1:4])
noquote(names(res$pairwise)[6])
noquote(res$pairwise[[6]][1:4,1:4])
noquote(names(res$pairwise)[7])
noquote(res$pairwise[[7]][1:4,1:4])
noquote(names(res$pairwise)[8])
noquote(res$pairwise[[8]][1:4,1:4])


###################################################
### code chunk number 13: diveRsity.Rnw:608-625
###################################################
noquote(names(res$bs_locus)[1])
res$bs_locus[[1]][c(1:3,nrow(res$bs_locus[[1]])),]
noquote(names(res$bs_locus)[2])
res$bs_locus[[2]][c(1:3,nrow(res$bs_locus[[2]])),]
noquote(names(res$bs_locus)[3])
res$bs_locus[[3]][c(1:3,nrow(res$bs_locus[[3]])),]
noquote(names(res$bs_locus)[4])
res$bs_locus[[4]][c(1:3,nrow(res$bs_locus[[4]])),]
noquote(names(res$bs_locus)[5])
res$bs_locus[[5]][c(1:3,nrow(res$bs_locus[[5]])),]
noquote(names(res$bs_locus)[6])
res$bs_locus[[6]][c(1:3,nrow(res$bs_locus[[6]])),]
noquote(names(res$bs_locus)[7])
res$bs_locus[[7]][c(1:3,nrow(res$bs_locus[[7]])),]
noquote(names(res$bs_locus)[8])
res$bs_locus[[8]][c(1:3,nrow(res$bs_locus[[8]])),]
#$


###################################################
### code chunk number 14: diveRsity.Rnw:640-657
###################################################
noquote(names(res$bs_pairwise)[1])
res$bs_pairwise[[1]][c(1:3,nrow(res$bs_pairwise[[1]])),]
noquote(names(res$bs_pairwise)[2])
res$bs_pairwise[[2]][c(1:3,nrow(res$bs_pairwise[[2]])),]
noquote(names(res$bs_pairwise)[3])
res$bs_pairwise[[3]][c(1:3,nrow(res$bs_pairwise[[3]])),]
noquote(names(res$bs_pairwise)[4])
res$bs_pairwise[[4]][c(1:3,nrow(res$bs_pairwise[[4]])),]
noquote(names(res$bs_pairwise)[5])
res$bs_pairwise[[5]][c(1:3,nrow(res$bs_pairwise[[5]])),]
noquote(names(res$bs_pairwise)[6])
res$bs_pairwise[[6]][c(1:3,nrow(res$bs_pairwise[[6]])),]
noquote(names(res$bs_pairwise)[7])
res$bs_pairwise[[7]][c(1:3,nrow(res$bs_pairwise[[7]])),]
noquote(names(res$bs_pairwise)[8])
res$bs_pairwise[[8]][c(1:3,nrow(res$bs_pairwise[[8]])),]
#$


###################################################
### code chunk number 15: diveRsity.Rnw:668-671 (eval = FALSE)
###################################################
## inCalc(infile, outfile = NULL, gp = 3, bs_locus = FALSE,
##        bs_pairwise = FALSE, bootstraps = 0, plot = FALSE,
##        parallel = FALSE)


###################################################
### code chunk number 16: diveRsity.Rnw:720-723
###################################################

load("./in_res.RData")
res_in <- in_results


###################################################
### code chunk number 17: diveRsity.Rnw:735-736
###################################################
noquote(res_in$Allele_In[1:10,c(1:5,ncol(res_in$Allele_In))])


###################################################
### code chunk number 18: diveRsity.Rnw:754-756
###################################################
res_in$l_bootstrap[1:10,]
#$


###################################################
### code chunk number 19: diveRsity.Rnw:775-785
###################################################
noquote(names(res_in$PW_bootstrap)[1])
res_in$PW_bootstrap[[1]][1:5,]
noquote(names(res_in$PW_bootstrap)[2])
res_in$PW_bootstrap[[2]][1:5,]
noquote(names(res_in$PW_bootstrap)[3])
res_in$PW_bootstrap[[3]][1:5,]
noquote(names(res_in$PW_bootstrap)[4])
res_in$PW_bootstrap[[4]][1:5,]
noquote(names(res_in$PW_bootstrap)[5])
res_in$PW_bootstrap[[5]][1:5,]


###################################################
### code chunk number 20: diveRsity.Rnw:792-793 (eval = FALSE)
###################################################
## readGenepop(infile = NULL, gp = 3, bootstrap = FALSE)


###################################################
### code chunk number 21: diveRsity.Rnw:853-854 (eval = FALSE)
###################################################
## corPlot(x,y)


###################################################
### code chunk number 22: diveRsity.Rnw:888-889 (eval = FALSE)
###################################################
## difPlot(x, outfile = NULL, interactive = FALSE)


###################################################
### code chunk number 23: diveRsity.Rnw:933-934 (eval = FALSE)
###################################################
## chiCalc(infile = NULL, outfile = NULL, gp = 3, minFreq = NULL)


###################################################
### code chunk number 24: diveRsity.Rnw:961-962 (eval = FALSE)
###################################################
## divOnline()


###################################################
### code chunk number 25: diveRsity.Rnw:969-971 (eval = FALSE)
###################################################
## fstOnly(infile = NULL, outfile = NULL, gp = 3, bs_locus = FALSE,
##         bs_pairwise = FALSE, bootstraps = 0, parallel = FALSE)


###################################################
### code chunk number 26: diveRsity.Rnw:1005-1007 (eval = FALSE)
###################################################
## divRatio(infile = NULL, outfile = NULL, gp = 3, pop_stats = NULL,
##          refPos = NULL, bootstraps = 1000, parallel = FALSE)


###################################################
### code chunk number 27: diveRsity.Rnw:1060-1062 (eval = FALSE)
###################################################
## bigDivPart(infile = NULL, outfile = NULL, WC_Fst = FALSE,
##            format = NULL)


###################################################
### code chunk number 28: diveRsity.Rnw:1089-1090 (eval = FALSE)
###################################################
## microPlexer()


###################################################
### code chunk number 29: diveRsity.Rnw:1110-1111 (eval = FALSE)
###################################################
## arp2gen(infile)


###################################################
### code chunk number 30: diveRsity.Rnw:1119-1120 (eval = FALSE)
###################################################
## divMigrate(infile = NULL, stat = c("gst", "d_jost"))


###################################################
### code chunk number 31: diveRsity.Rnw:1128-1130 (eval = FALSE)
###################################################
## haploDiv(infile = NULL, outfile = NULL, pairwise = FALSE,
##          bootstraps = 0)


###################################################
### code chunk number 32: diveRsity.Rnw:1143-1144
###################################################
library("diveRsity")


###################################################
### code chunk number 33: diveRsity.Rnw:1154-1155 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 34: diveRsity.Rnw:1162-1163
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 35: diveRsity.Rnw:1172-1177 (eval = FALSE)
###################################################
## div_results <- divPart(infile = Test_data, outfile = "Test", 
##                          gp = 3, pairwise = TRUE, 
##                          WC_Fst = TRUE, bs_locus = TRUE, 
##                          bs_pairwise = TRUE, bootstraps = 100, 
##                          plot = FALSE, parallel = TRUE)


###################################################
### code chunk number 36: diveRsity.Rnw:1180-1182
###################################################

load("./div_res.RData")


###################################################
### code chunk number 37: diveRsity.Rnw:1192-1193
###################################################
names(div_results)


###################################################
### code chunk number 38: diveRsity.Rnw:1198-1199
###################################################
typeof(div_results$bs_locus)


###################################################
### code chunk number 39: diveRsity.Rnw:1202-1203
###################################################
#$


###################################################
### code chunk number 40: diveRsity.Rnw:1207-1208
###################################################
names(div_results$bs_locus)


###################################################
### code chunk number 41: diveRsity.Rnw:1210-1211
###################################################
#$


###################################################
### code chunk number 42: diveRsity.Rnw:1218-1219 (eval = FALSE)
###################################################
## mymatrix[5, 1]


###################################################
### code chunk number 43: diveRsity.Rnw:1225-1226
###################################################
div_results$bs_locus$Gst[1:10, ]


###################################################
### code chunk number 44: diveRsity.Rnw:1230-1231 (eval = FALSE)
###################################################
## div_results$bs_locus$Gst[ ,1]


###################################################
### code chunk number 45: diveRsity.Rnw:1243-1244 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 46: diveRsity.Rnw:1251-1252
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 47: diveRsity.Rnw:1261-1265 (eval = FALSE)
###################################################
## in_results <- inCalc (infile = Test_data, outfile = "Test", 
##                          gp = 3, bs_locus = TRUE, 
##                          bs_pairwise = TRUE, bootstraps = 100, 
##                          plot = FALSE, parallel = TRUE)


###################################################
### code chunk number 48: diveRsity.Rnw:1268-1270
###################################################
load("./in_res.RData")



###################################################
### code chunk number 49: diveRsity.Rnw:1283-1284
###################################################
names(in_results)


###################################################
### code chunk number 50: diveRsity.Rnw:1289-1290
###################################################
typeof(in_results$PW_bootstrap)


###################################################
### code chunk number 51: diveRsity.Rnw:1293-1294
###################################################
#$


###################################################
### code chunk number 52: diveRsity.Rnw:1299-1300
###################################################
names(in_results$PW_bootstrap)


###################################################
### code chunk number 53: diveRsity.Rnw:1302-1303
###################################################
#$


###################################################
### code chunk number 54: diveRsity.Rnw:1308-1309 (eval = FALSE)
###################################################
## mymatrix[5, 1]


###################################################
### code chunk number 55: diveRsity.Rnw:1315-1316
###################################################
in_results$PW_bootstrap[["pop1, vs. pop2,"]][1:3, ]


###################################################
### code chunk number 56: diveRsity.Rnw:1318-1319
###################################################
#$


###################################################
### code chunk number 57: diveRsity.Rnw:1323-1324 (eval = FALSE)
###################################################
## in_results$PW_bootstrap[["pop1, vs. pop2,"]][ ,1]


###################################################
### code chunk number 58: diveRsity.Rnw:1335-1336 (eval = FALSE)
###################################################
## setwd("mypath")


###################################################
### code chunk number 59: diveRsity.Rnw:1343-1344
###################################################
data(Test_data, package = "diveRsity")


###################################################
### code chunk number 60: diveRsity.Rnw:1352-1354
###################################################
gp_res <- readGenepop(infile = Test_data, gp = 3,
                      bootstrap = FALSE)


###################################################
### code chunk number 61: diveRsity.Rnw:1362-1363
###################################################
names(gp_res)


###################################################
### code chunk number 62: diveRsity.Rnw:1378-1385
###################################################
locus18_pop1 <- c(gp_res$pop_alleles[[1]][[1]][,18], 
                  gp_res$pop_alleles[[2]][[1]][,18])
# sort alleles by size
allele_sort <- order(locus18_pop1, decreasing = FALSE)
#plot
plot(locus18_pop1[allele_sort], ylab = "allele size", col="blue",
     pch = 16)


###################################################
### code chunk number 63: diveRsity.Rnw:1401-1436
###################################################
# Define a results matrix with 37 columns (loci) and
# 1000 rows (bootstraps)to record allele number per locus

num_all <- matrix(rep(0,(37*10)), ncol = 37)

# Now using readGenepop we can fill the matrix
bs<-10
for(i in 1:bs){
    # first produce a bootstrap file
    
    x <- readGenepop(infile = Test_data, gp = 3,
                     bootstrap = TRUE)
                     
    # Now record the number of alleles at each locus
    
    num_all[i, ] <- x$nalleles                      
}

# Now we can use this data to calculate the mean
# number of alleles per locus as well at their
# 95% confidence intervals

mean_num <- colMeans(num_all)
lower<-vector()
upper<-vector()
for(i in 1:ncol(num_all)){
    lower[i] <- mean_num[i] - (1.96 * sd(num_all[,i]))  
    upper[i] <- mean_num[i] + (1.96 * sd(num_all[,i]))
}

# Now we can create a data frame of these results

bs_res <- data.frame(mean_num, lower, upper)

bs_res[1:10,]


###################################################
### code chunk number 64: diveRsity.Rnw:1461-1519 (eval = FALSE)
###################################################
##  # load the diveRsity package 
## library("diveRsity")
## 
## # We can specify the names of our simulation folders in two ways 
## 
## # manually 
## fold_names <- paste("sim", 1:10, sep = "")
## 
## or
## 
## # automatically (when there is only a single level below the
## # top directory)
## fold_names <- list.dirs(full.names = TRUE, recursive = FALSE) 
## 
## # Now we can determine the names of all genepop files in each folder
## 
## file_names <- lapply(fold_names, function(x){
##   files <- dir(path = paste(x, "/", sep = ""), pattern = "*.gen",
##                full.names = TRUE)
##   return(files) })
## 
## # file_names will be a list of length 10. Each element will contain
## # the names of all .gen files within the respective simulation folder
## 
## # Before we are ready to run the main analyses, we should set up 
## # the parallel environment
## 
## # load the doParallel package
## library("doParallel")
## # set up a cluster of 10 CPUs (one for each batch of files)
## cl <- makeCluster(10)
## # Export the 'divPart' function to the cluster cores 
## clusterExport(cl, "divPart", envir = environment())
## 
## # Now we can run the main analyses
## 
## results <- parLapply(cl, file_names, function(x){
##   sim_res <- sapply(x, function(y){
##     out <- divPart(infile = y, gp = 3, WC_Fst = TRUE)
##     return(out$estimate[nrow(out$estimate), 4:7])
##   })
##   return(t(sim_res))  # transpose sim_res
## })
## 
## # This will generate a list (of length 10), with each element 
## # containing a matrix of 1000 rows (1 per file) and 4 columns 
## # (1 for each diversity statistic) 
## 
## # Example of output for simulation 1
## # G_st_est        G_hed_st_est    D_Jost_est      Fst_WC
## # 0.3905          0.8938          0.8256          0.4010
## # 0.5519          0.8719          0.6986          0.6031
## # 0.5924          0.8880          0.7092          0.6096
## # ...             ...             ...             ...
## # ...             ...             ...             ...
## 
## # these results could then be piped to further analyses or 
## # visualisation tools


###################################################
### code chunk number 65: diveRsity.Rnw:1534-1535
###################################################
print(sessionInfo(), locale = FALSE)


