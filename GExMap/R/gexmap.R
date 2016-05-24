gexmap<-function (genome = "homosapiens", scale = "", source = "", res = "", isGO = FALSE, isMAP = TRUE, lim_chi = 5, global_test_choice = 4, pcorrd = 1, pcorrv = 1)
{
################################################################################
#
#     FUNCTIONS
#
################################################################################

gexload.data <-
function (source, genome)
{
   ###############################################################################
   #                                                                             #
   #                                LOAD.DATA                                    #
   #                                                                             #
   #      ACTION = Load the data.Rdata file                                      #
   #      OUTPUT = Matrix "data" of all ENSEMBL annotations                      #
   #                                                                             #
   ###############################################################################

    write("Loading localization data", file = "")
    data.Rdata.dir = file.path(paste(source, genome, ".Rdata", sep = ""))
    if (!file.exists(data.Rdata.dir))
    {
        write(sprintf("ERROR: file ", data.Rdata.dir, " does not exist in the default directory", source), file = "")
        write(paste("Load manually the ", genome, ".Rdata file", sep = ""), file = "")
        load(file = file.choose())
    } else
    {
        load(file = data.Rdata.dir)
        write(paste(data.Rdata.dir, " file Loaded", sep = ""), file = "")
    }
    return(data)
}

gextest<-function (genome.test, nhazard, results.graph, lim_chi, global_test_choice, pcorrd, pcorrv)
{
    #########################################################################################################
    #	The test take into acount all the chr and compares reel effectives of our list to the presence probability issued from the ENSEMBL genome.
    #	Probabilities are computed chr by chr
    #	Sum of presence probabilities of one chr is 1
    #	The program tests each chr at progressives scaling: from 1 to 10 Mb by unit. 
    #
    #########################################################################################################


######## Global tests to select chromosomes
   
    #round of genome chr by chr

    lev = levels(factor(genome.test[, "chr"]))
    max = 20
    g_test = matrix(0, ncol = max, nrow = length(lev))   #matrix Nchr*11 will stock the results of the tests for each chr and for each scale
    rownames(g_test) = lev
    colnames(g_test) = seq(1, max, 1)
    
    #g_test_graph: only pval to be reported graphically
    g_test_graph = g_test
    
    global_tests = matrix(0, ncol = 5, nrow = length(lev))
    colnames(global_tests) = c("chr", "pval", "scale", "min", "W")
    global_tests[, 1] = as.matrix(lev)
    
    lev = as.data.frame(genome.test)
    
    local_test_results = matrix(0, ncol = 12, nrow = 0)
    colnames(local_test_results) = c("chr", "unit", "total", "list", "up", "down", "hazard", "theo", "regions", "binome.test", "Optmized", "Wilcoxon")
    
    write(paste("Global tests", sep = ""), file = "")
    for (I in levels(lev[, "chr"]))
    {
      write(paste("Chr ", I, sep = ""), file = "")
      gentest = genome.test[genome.test[, "chr"] == I, ]  # extract data of crh
      if(nrow(gentest)!=0)
      {      
          #Wilcoxon global test comparaing real distribution (gentest[,"list"]) to expected distribution (gentest[,"hazard"])
          p = 1
          p = try(as.numeric(wilcox.test(as.matrix(as.numeric(gentest[, "list"])), as.matrix(as.numeric(gentest[, "hazard"])), paired = TRUE, exact = FALSE, correct = FALSE)[3]), silent=TRUE)
          if (p <= 0.05)
          {
              global_tests[global_tests[, "chr"] == I, 5] = p
          }
          
          nhazard = sum(as.numeric(gentest[, "total"]))    #N = toltal ENSEMBL genes of the chr
          theo = matrix(0, ncol = 1, nrow = nrow(gentest))     #matrix of probabilities of presence
          colnames(theo) = "theo"
          for (K in 1:nrow(gentest)) #cration of the matrix of presence probabilities
          {
              theo[K, ] = as.numeric(gentest[K, "total"])/nhazard
          }
          
          gentest = cbind(gentest, theo)
          for (N in 1:max)
          {
              # concatenates units
              if (N >= 2)
              {
                  if (floor(nrow(gentest)/N) * N == nrow(gentest))
                  {
                    x = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo_tmp = as.matrix(gentest[, "theo"])
                    x_tmp = as.matrix(gentest[, "list"])
                  } else     #si le compte n'est pas juste
                  {
                    x = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo = matrix(0, ncol = 1, nrow = (floor(nrow(gentest)/N)))
                    theo_tmp = as.matrix(gentest[, "theo"])
                    theo_tmp[(floor(nrow(gentest)/N) * N), ] = sum(as.numeric(gentest[((floor(nrow(gentest)/N)) * N):nrow(gentest), "theo"]))
                    theo_tmp = as.matrix(theo_tmp[1:((floor(nrow(gentest)/N)) * N), ])
                    x_tmp = as.matrix(gentest[, "list"])
                    x_tmp[(floor(nrow(gentest)/N) * N), ] = sum(as.numeric(gentest[((floor(nrow(gentest)/N)) * N):nrow(gentest), "list"]))
                    x_tmp = as.matrix(x_tmp[1:((floor(nrow(gentest)/N)) * N), ])
                  }
                  for (K in 1:nrow(x)) #concatenates units
                  {
                    x[K, ] = sum(as.numeric(x_tmp[(K * N - N + 1):(K * N), ]))
                    theo[K, ] = sum(as.numeric(theo_tmp[(K * N - N + 1):(K * N), ]))
                  }
              } else
              {
                  x = as.matrix(as.numeric(gentest[, "list"]))
                  theo = as.matrix(as.numeric(gentest[, "theo"]))
              }
              #test
              p = suppressWarnings(try(as.numeric(chisq.test(x, p = theo)[3]), silent=TRUE))
              if (is.nan(p) | (p> 0.05))
              {
                  g_test[I, N] = ""
                  g_test_graph[I, N] = ""
              } else
              {
                  g_test_graph[I, N] = 1
                  g_test_graph[I, N] = suppressWarnings(try(as.numeric(chisq.test(x, p = theo)[3]), silent=TRUE))
                  g_test[I, N] = paste(g_test_graph[I, N], "(", min(x), ")", sep = "")
                  #test of the minimal number of genes by unit
                  if ((min(x) >= lim_chi) & (global_tests[global_tests[, 1] == I, 2] == 0))
                  {
                    global_tests[global_tests[, 1] == I, 2] = g_test[I, N]
                    global_tests[global_tests[, 1] == I, 3] = N
                    global_tests[global_tests[, 1] == I, 4] = min(x)
                  }
              }
          }
          
          #interesting region detection by a local test if the chr is interesting
          local_test = matrix(0, ncol = (ncol(genome.test) + 3), nrow = 0)
          colnames(local_test) = c(colnames(genome.test), "theo", "regions", "binome.test")
          
          # user choice of CHI/Wilcoxon/CHI-Wilcoxon
          pass_to_local = FALSE
          # 1- At least CHI is OK
          if ((global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 1))
          {
              pass_to_local = TRUE
              write(paste("1- CHI is OK", sep = ""), file = "")
          }
          # 2- At least Wilcoxon is OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) & (global_test_choice == 2))
          {
              pass_to_local = TRUE
              write(paste("2- At least Wilcoxon is OK", sep = ""), file = "")
          }
           # 3- CHI & Wilcoxon are OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) & (global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 3))
          {
              pass_to_local = TRUE
              write(paste("3- CHI & Wilcoxon are OK", sep = ""), file = "")
          }
          # 4- CHI OR Wilcoxon is OK
          if ((global_tests[global_tests[, 1] == I, 5] != 0) | (global_tests[global_tests[, 1] == I, 2] != 0) & (global_test_choice == 4)) {
              pass_to_local = TRUE
              write(paste("4- CHI OR Wilcoxon are OK", sep = ""), file = "")
          }
          
          if (pass_to_local)
          {
              #matrix of detected regions of interest
              regions = matrix(0, ncol = 1, nrow = nrow(gentest))
              colnames(regions) = "regions"
              
              #matrix of pval form test of each unit of regions of interest
              test = matrix("", ncol = 1, nrow = nrow(gentest))
              colnames(test) = "test"
              
              totlist = sum(as.numeric(genome.test[, "list"]))   # total number of genes of the tested list
              for (U in 1:(nrow(gentest)))
              {
                  if ((as.numeric(gentest[U, "list"]) > as.numeric(gentest[U, "hazard"])))
                  {
                    regions[U, 1] = 1
                    test[U, ] = 1
                    test[U, ] = try(as.numeric(binom.test(as.numeric(gentest[U, "list"]), totlist, as.numeric(gentest[U, "total"])/nhazard)[3]), silent=TRUE)
                  }
                  else
                  {
                    regions[U, 1] = 0
                    test[U, ] = 1
                  }
              }
              
              #pval correction for binomial test	
              write(paste("global pval correction", sep = ""), file = "")
              library(multtest)
              procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
              # correct results
              test = as.matrix(as.numeric(test))
              pval.adjusted <- mt.rawp2adjp(test, procs)
              pval.adjusted <- pval.adjusted$adjp[order(pval.adjusted$index), ]
              test = pval.adjusted[, pcorrd]
              
              gentest = cbind(gentest, regions, test)
              local_test = rbind(local_test, gentest)
  
        			#########################################################################################################
        			#
       			  #        Validation & optimization
        			#
       	 			#########################################################################################################
              
              write(paste("Validation & optimization", sep = ""), file = "")
              GExRegions = as.matrix(local_test[, "regions"])
              #region extension
              for (U in 1:nrow(GExRegions))
              {
                  if (local_test[U, "hazard"] >= local_test[U, "list"])
                  {
                    GExRegions[U, ] = 0
                  }
              }
              for (U in 2:(nrow(GExRegions) - 1))
              {
                  if ((as.numeric(local_test[U, "binome.test"]) > 0.05))
                  {
                    GExRegions[U, ] = 0
                  }
              }
              GExRegionst = GExRegions
              
              for (U in 2:(nrow(GExRegions) - 1))
              {
                  if (GExRegions[U, ] == 1)
                  {
                    GExRegionst[U - 1, ] = 1
                    GExRegionst[U + 1, ] = 1
                  }
              }
              GExRegions = GExRegionst
              if (GExRegions[1, ] == 1)
              {
                  GExRegions[2, ] = 1
              }
              if (GExRegions[nrow(GExRegions), ] == 1)
              {
                  GExRegions[(nrow(GExRegions) - 1), ] = 1
              }
              
              genome.test.tmp = cbind(local_test[, c("chr", "unit", "list", "hazard")], GExRegions)
              colnames(genome.test.tmp) = c("chr", "unit", "list", "hazard", "regions")
              #matrix for pval of all tested genome
              test = matrix(1, ncol = 1, nrow = 0)
              colnames(test) = "test"
              
              # listing data localisations of interesting regions before binding in gexROI matrix
              for (C in levels(as.data.frame(genome.test.tmp)[, "chr"]))
              {
  				      #########################################################################################################
  				      #
  				      #        censing
  				      #
  				      #########################################################################################################
  				      
                  write(paste("Finding interesting regions for chr ", C, sep = ""), file = "")
                  gexChr = genome.test.tmp[genome.test.tmp[, "chr"] == C, ]
                  start = matrix(0, ncol = 3, nrow = (nrow(gexChr) + 1))
                  start[1:nrow(start) - 1, 1] = as.numeric(gexChr[, "regions"])
                  start[2:nrow(start), 2] = as.numeric(gexChr[, "regions"])
                  
                  start[, 3] = start[, 1] - start[, 2]
                  start = start[1:nrow(start) - 1, ]
                  
                  #correspondences matrix of statrs and units
                  regions = cbind(gexChr[, "unit"], start[, 3])
                  
                  #end shifted matrix (cf to annexes)
                  end = matrix(0, ncol = 3, nrow = (nrow(gexChr) + 1))
                  end[2:nrow(end), 1] = as.numeric(gexChr[, "regions"])
                  end[1:nrow(end) - 1, 2] = as.numeric(gexChr[, "regions"])
                  
                  end[, 3] = end[, 1] - end[, 2]
                  end = end[2:nrow(end), ]
                  
                  regions = cbind(regions, end[, 3])
                  
                  gexROI = matrix(0, ncol = 3, nrow = length(regions[regions[, 2] == 1, 1]))
                  gexROI[, 1] = as.numeric(regions[as.numeric(regions[, 2]) == 1, 1])
                  gexROI[, 2] = as.numeric(regions[regions[, 3] == 1, 1])
                  
                  gexROI[, 3] = gexROI[, 2] - gexROI[, 1]
                  
                  ########################################################################################################
    				      #
    				      #        TEST
    				      #
    				      #########################################################################################################
    				      
    				      #Matrix of test results
    				      #two columns: the first for the chromosome units, the second for the pvalues
                  pval = matrix("", ncol = 2, nrow = nrow(gexChr))
                  pval[, 1] = gexChr[, "unit"]
                  
                  if(nrow(gexROI)!=0)
                  {
                    for (I in 1:nrow(gexROI))
                    {
                      #Test of entire region
                      list = as.vector(as.numeric(gexChr[gexROI[I, 1]:gexROI[I, 2], "list"]))
                      hazard = as.vector(as.numeric(gexChr[gexROI[I, 1]:gexROI[I, 2], "hazard"]))
                      p = wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)
                      p = as.numeric(p[3])
                      #saving results
                      pval[gexROI[I, 1]:gexROI[I, 2], 2] = p
                      
                      				         
    				      #########################################################################################################
    				      #
    				      #        PROGRESSIVE OPTIMISATION 
    				      #
    				      #########################################################################################################
    				         
    				         #localization information of the region
                      l = gexROI[I, 2] - gexROI[I, 1]
                      start = gexROI[I, 1]
                      end = gexROI[I, 2]
                      
                      #While region width>3 and pval is non significative
                      while ((p > 0.05) && (l >= 3))
                      {
    				            #upstream reduction 
                        list = as.vector(as.numeric(gexChr[(start + 1):end, "list"]))
                        hazard = as.vector(as.numeric(gexChr[(start + 1):end, "hazard"]))
                        p = wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)
                        p1 = as.numeric(p[3])
                        
                        #downstream reduction 
                        list = as.vector(as.numeric(gexChr[start:(end - 1), "list"]))
                        hazard = as.vector(as.numeric(gexChr[start:(end - 1), "hazard"]))
                        p2 = 1
                        p2 = try(as.numeric(wilcox.test(list, hazard, paired = TRUE, exact = FALSE, correct = FALSE)[3]), silent=TRUE)
                        
                        #Choice of best pval
                        if (p1 > p2)
                        {
                          p = p2
                          pval[start:(end - 1), 2] = p
                          end = end - 1
                        }
                        if (p2 >= p1)
                        {
                          p = p1
                          pval[(start + 1):end, 2] = p
                          start = start + 1
                        }
                        l = l - 1
                      } #End of optimisation	
                    } #End of tests for one chromosome
                    #pval correction for local wolcoxon tests
                    write(paste("local pval correction", sep = ""), file = "")
                 }else{
                  write(paste("*****No Region of Interest found\n", sep=""), file="")
                 }
                 test = rbind(test, as.matrix(pval[, 2]))
                  test = as.matrix(as.numeric(test))
                  procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY")
                  # correct results
                  pval.adjusted <- mt.rawp2adjp(test, procs)
                  pval.adjusted <- pval.adjusted$adjp[order(pval.adjusted$index), ]
                  test = as.matrix(pval.adjusted[, pcorrv])
              }
              colnames(GExRegions) = "Optmized"
              colnames(test) = "Wilcoxon"
              local_test = cbind(local_test, GExRegions, test)
              
              local_test[local_test[, "regions"] == "0", "regions"] = ""
              local_test[local_test[, "Optmized"] == "0", "Optmized"] = ""
              
              local_test_results = rbind(local_test_results, local_test)
          }
      }
    }
    pdf(file = file.path(paste(results.graph, "test_global.pdf", sep = "")))
    matplot(t(g_test_graph), pch = c(1:ncol(t(g_test_graph))), type = "p", col = c(1:ncol(t(g_test_graph))))
    legend(1, 0.045, pch = c(1:ncol(t(g_test_graph))), colnames(t(g_test_graph)), col = c(1:ncol(t(g_test_graph))), cex = 0.7)
    dev.off()
    write.table(g_test, row.names = TRUE, sep = "\t", file = file.path(paste(results.graph, "Chi_results.txt", sep = "")))
    write.table(global_tests, row.names = TRUE, sep = "\t", file = file.path(paste(results.graph, "global_tests_results.txt", sep = "")))
    return(local_test_results)
}

gex.mapping <-
function (data, list.chr, scale, I, list.ENS)
{
   ###############################################################################
   #                                                                             #
   #                                 MAPPING                                     #
   #                                                                             #
   #                           List, Genome, Hazard                              #
   #                                                                             #
   ###############################################################################
  
    ENS = cbind(as.matrix(as.numeric(data[data[, "chr"] == I, "chr"])), as.matrix(as.numeric(data[data[, "chr"] == I, "m"]))) # extracts data for the current chr
    colnames(ENS) = c("chr", "m")
    ENS = as.matrix(ENS[order(ENS[, "m"]), ])       # sort m by asc
    size = max(ENS[, "m"])                          # max valure of m
    size = floor(size/scale) + 1                    # scaling, heigth = number of unit for the chr
    
    genome.temp = matrix(0, ncol = 3, nrow = size)   # matrix of the number of ENSEMBL genes by unit)
    genome.temp[, 1] = seq(1, size)
    chr = matrix(I, ncol = 1, nrow = size)
    genome.temp = cbind(chr, genome.temp)
    colnames(genome.temp) = c("chr", "scale", "ensembl", "list")
    
    # treatment of the matrix list.chr
    for (J in 1:nrow(ENS))                          # for each gene of the current chr
    {
        genome.temp[floor(as.integer(ENS[J, "m"])/scale) + 1, "ensembl"] = as.numeric(genome.temp[floor(as.integer(ENS[J, "m"])/scale) + 1, "ensembl"]) + 1
    }
    
    for (J in 1:nrow(list.chr))
    {
        genome.temp[floor(as.integer(list.chr[J, "m"])/scale) + 1, "list"] = as.numeric(genome.temp[floor(as.integer(list.chr[J, "m"])/scale) + 1, "list"]) + 1
    }
    up = matrix(0, ncol = 1, nrow = nrow(genome.temp)) #upregulated genes
    colnames(up) = "up"
    down = matrix(0, ncol = 1, nrow = nrow(genome.temp)) #downregulated genes
    colnames(down) = "down"
    for (J in 1:nrow(list.chr))                     # round of current chromosome
    {
        if (as.numeric(list.chr[J, "expression"] > 0))
        {
            up[floor(as.integer(list.chr[J, "m"])/scale) + 1, "up"] = up[floor(as.integer(list.chr[J, "m"])/scale) + 1, "up"] - 1
        }
        if (as.numeric(list.chr[J, "expression"] < 0))
        {
            down[floor(as.integer(list.chr[J, "m"])/scale) + 1, "down"] = down[floor(as.integer(list.chr[J, "m"])/scale) + 1, "down"] - 1
        }
    }
    genome.temp = cbind(genome.temp, up, down)
    return(genome.temp)
}

gexgraph <-
function (genome.test, scale, results.graph, W, cytobands)
{
   ###############################################################################
   #                                                                             #
   #                              GExGraph                                       #
   #                                                                             #
   #        INPUT: Genome.map                                                    #
   #        OUTPUT: graphics in pdf format                                       #
   #                                                                             #
   ###############################################################################

    write("Graph", file = "")
    lev = as.data.frame(genome.test)
    genome.test[as.numeric(genome.test[, "binome.test"]) > 0.05, "binome.test"] = ""
    genome.test[as.numeric(genome.test[, "Wilcoxon"]) > 0.05, "Wilcoxon"] = ""
    for (C in levels(lev[, "chr"]))
    {
        #temporary data matrix of each chr
        genome.temp = genome.test[genome.test[, "chr"] == C, ]
        
        #cytobands matrix
        cytobands_tmp = cytobands[cytobands[,"chr"]==C, c("cytoband", "debut", "fin")]       
        temp = cytobands_tmp
        cytobands_tmp = cbind(as.matrix(as.numeric(temp[,"debut"])), as.matrix(as.numeric(temp[,"fin"])))
        rownames(cytobands_tmp) = temp[,"cytoband"]
        colnames(cytobands_tmp) = c("debut", "fin")
        
        #maximum value of x = maximum value of total genes =>width of the graphic
        max = max(as.numeric(genome.temp[, "total"]))
        #maximum value of y = number of unit for the current chr = height of the graphic
        may = max(as.numeric(genome.temp[, "unit"]))
        
        write(paste("graph ", genome.temp[1, "chr"], sep = ""), file = "")
        
        colours = character()
        COL.ENS = "green"
        COL.list = "orange"
        COL.Haz = "black"
        COL.up = "red"
        COL.down = "blue"
        COL.regions = "red"
        COL.binome.test = "cyan"
        COL.Optmized = "orange"
        COL.Wilcoxon = "green"
        colours = c(COL.ENS, COL.list, COL.Haz, COL.up, COL.down, COL.regions, COL.binome.test, COL.Optmized, COL.Wilcoxon)
        
        y = as.matrix(genome.temp[, "unit"])
        
        #min and max of y axis
        min_up = min(as.numeric(genome.temp[, "up"]))
        min_down = min(as.numeric(genome.temp[, "down"]))
        min = min(min_down, min_up) - 1
        
        #replace the "1" which marks the regions of interest by the min value to draw the line of interesting region under the curves
        genome.temp[genome.temp[, "regions"] == 1, "regions"] = min
        
        #pvals selection <=0.05   
        genome.temp[genome.temp[, "binome.test"] != "", "binome.test"] = min - 0.5
        genome.temp[genome.temp[, "Optmized"] != "", "Optmized"] = min - 1
        genome.temp[genome.temp[, "Wilcoxon"] != "", "Wilcoxon"] = min - 1.5
        
        min = min - 1
        
        x = as.matrix(genome.temp[, c("total", "list", "hazard", "up", "down", "regions", "binome.test", "Optmized", "Wilcoxon")])
        
        pdf(file = paste(results.graph, "Chr", C, ".pdf", sep = ""), width = (may*0.8), height = (abs(min-2)+max)/2.5, family = "Helvetica", title = paste("chr", C, sep = ""))
        #pdf(file = paste(results.graph, "Chr", C, ".pdf", sep = ""), width = (may*1.5), height = (abs(min-2)+max)/2.5, family = "Helvetica", title = paste("chr", C, sep = ""))
        
        matplot(y, x, col = colours, type = "o", pch = 20, axes = FALSE, ylim = c(min, max), ylab = "", xlab = "")
        
        title(main = paste("Chr", C), cex.main = 3, font.main = 2, col.main = "black", font.lab = 2, col.lab = "black", xlab = "", ylab = "")
        title(xlab = paste(scale, " pb", sep = ""), line = 3.5, cex.lab = 3)
        title(ylab = "Gene frequency", line = 1.5, cex.lab = 3)
        axis(side = 2, at = min:max, tick = TRUE, lab = c(rep(min:max, 1)), cex.axis = 1)
        axis(side = 1, at = seq(1, (may + 1), by = 5), tick = TRUE, lab = c(seq(0, may, by = 5)), cex.axis = 1)
        
        #cytobands
        pair=0
        for(K in 1:nrow(cytobands_tmp))
        {
          if(pair == 0)
          {
            rect((cytobands_tmp[K, "debut"]/scale), (min-0.75), (cytobands_tmp[K, "fin"]/scale), (min - 1.25), col="black") 
            text(((((cytobands_tmp[K, "fin"]/scale) - (cytobands_tmp[K, "debut"]/scale))/2) + (cytobands_tmp[K, "debut"]/scale)) , (min-1), rownames(cytobands_tmp)[K], cex = 1, col="white")
            pair = 1
          } else
          {
            rect((cytobands_tmp[K, "debut"]/scale), (min-0.75), (cytobands_tmp[K, "fin"]/scale), (min - 1.25), col="white") 
            text(((((cytobands_tmp[K, "fin"]/scale) - (cytobands_tmp[K, "debut"]/scale))/2) + (cytobands_tmp[K, "debut"]/scale)) , (min-1), rownames(cytobands_tmp)[K], cex = 1, col="black")
           pair = 0
          }
        }
        
        #  Legend of curves
        leg.txt = character()
        leg.txt <- c(leg.txt, "Curves")
        leg.txt <- c(leg.txt, "")
        leg.txt <- c(leg.txt, "Ensembl genome")
        leg.txt <- c(leg.txt, "Tested list")
        leg.txt <- c(leg.txt, "Hazard estimation")
        leg.txt <- c(leg.txt, "Upregulated genes")
        leg.txt <- c(leg.txt, "Downregulated genes")
        COL.blank = "white"
        colours = c(COL.Haz, COL.blank, COL.ENS, COL.list, COL.Haz, COL.up, COL.down)
        legend(1, max, leg.txt, pch = "  -----", col = colours, cex = 3)
        
        #  Legend of statistical results
        leg.txt = character()
        leg.txt <- c(leg.txt, "Statistucal results")
        leg.txt <- c(leg.txt, "")
        leg.txt <- c(leg.txt, "Units where Tested list > Hazard estimation")
        leg.txt <- c(leg.txt, "Statistically validated")
        leg.txt <- c(leg.txt, "Optimized regions")
        leg.txt <- c(leg.txt, "Significant")
        colours = c(COL.Haz, COL.blank, COL.regions, COL.binome.test, COL.Optmized, COL.Wilcoxon)
        legend(15, max, leg.txt, pch = "  ----", col = colours, cex = 3)
        dev.off()
    }
}

gexgo <-
function (list.ENS, source, res)
{
    write(paste("Analyzing GO identifiers", sep = ""), file = "")
    data.Rdata.dir = file.path(paste(source, "GO.Rdata", sep = ""))
    if (!file.exists(data.Rdata.dir))
    {
        write(sprintf("ERROR: file ", data.Rdata.dir, " does not exist in the default directory", source), file = "")
        write(paste("Load manually the GO.Rdata file", sep = ""), file = "")
        load(file = file.choose())
    } else
    {
        load(file = data.Rdata.dir)
        write(paste(data.Rdata.dir, " file Loaded", sep = ""), file = "")
    }
    
    #bind two columns to the GO data matrix
    go = cbind(go, matrix(0, ncol = 1, nrow = nrow(go))) #one column of 0 for the scores of each GO id
    go = cbind(go, matrix("", ncol = 1, nrow = nrow(go))) #one column for the list of genes relatives to the GO id
    colnames(go) = c("id", "ontology", "type", "score", "genes")
    gengo = list.ENS[list.ENS[, "GO"] != "0", ]   #filters only genes with GO id
    write(paste(nrow(gengo), " genes have GO identifiers", sep = ""), file = "")
    
    
    #comput GO id scrores and create corresponding genes lists
    
    tmp = gengo[,c("GO", "name")]
    
    allgo = as.matrix(unlist(list(apply(as.matrix(gengo[, 6]), 1, function(x) unlist(strsplit(x, ","))))))
    allnames = as.matrix(unlist(apply(tmp, 1, function(x) rep(x[2], length(unlist(strsplit(x[1], ",")))))))
    
    #scores
    t=as.matrix(table(factor(allgo)))
    for(I in 1:nrow(t))
    {
       go[go[,"id"]==rownames(t)[I], "score"]=t[I,]
    }
    
    concat<-function(X)
    {
      if(length(X)>=2) return(paste(X,collapse=",",sep=","))
      if(length(X)==1) return(paste(X))
      if(length(X)==0) return("")
    }
      
    #names
    gonames = cbind(allgo, allnames)
    gonames = gonames[order(gonames[,1]),]
    for (I in levels(factor(gonames[, 1])))
    {
       go[go[,"id"]==I, "genes"]= concat(gonames[gonames[,1]==I,2])
    }

    resume <- file(paste(res, "resume.txt", sep = ""), "w+")
    cat(paste("Size of input list: ", nrow(list.ENS), sep = ""), file = resume)
    cat(paste("\n", nrow(gengo), " genes have GO identifiers", sep = ""), file = resume)
    
    #filters GO id with score=0
    resgo = go[go[, 4] != "0", ]
    resgoP = resgo[resgo[, 3] == "P", ]
    resgoF = resgo[resgo[, 3] == "F", ]
    resgoC = resgo[resgo[, 3] == "C", ]
    
    orderP = as.matrix(order(as.matrix(as.numeric(resgoP[, 4])), decreasing = TRUE))  #sort the Biological Processes by score
    tP = as.matrix(unlist(strsplit(resgoP[orderP[1, ], 5], ",")))  #matrix of the longuest list of genes for a GO identifier 
    #tP = as.matrix(tP[2:nrow(tP), ]) #remove the first row which is empty
    detailP = matrix("", ncol = nrow(orderP), nrow = nrow(tP))  #matrix of all genes for each GO id
    sumP = sum(as.numeric(resgoP[, 4]))   #Sum of genes for this GO type
    resP = matrix(0, ncol = 4, nrow = 0)
    prop = matrix(0, ncol = 2, nrow = nrow(orderP))   #part of each GO id in the GO type
    for (I in 1:nrow(orderP))
    {
        resP = rbind(resP, resgoP[orderP[I, ], 1:4])
        prop[I, 1] = round((as.numeric(resP[resP[, 1] == resgoP[orderP[I, ], 1], 4])/sumP) * 100, 2)
        prop[I, 2] = paste(resP[I, 2], " (", prop[I, 1], "%,", resP[I, 4], ")", sep = "")
    }
    resP = cbind(resP, prop)

    orderF = as.matrix(order(as.matrix(as.numeric(resgoF[, 4])), decreasing = TRUE))  #sort the Molecular Function by score
    tF = as.matrix(unlist(strsplit(resgoF[orderF[1, ], 5], ",")))  #matrix of the longuest list of genes for a GO identifier 
    tF = as.matrix(tF[2:nrow(tF), ]) #remove the first row which is empty
    detailF = matrix("", ncol = nrow(orderF), nrow = nrow(tF))  #matrix of all genes for each GO id
    sumF = sum(as.numeric(resgoF[, 4]))   #Sum of genes for this GO type
    resF = matrix(0, ncol = 4, nrow = 0)
    prop = matrix(0, ncol = 2, nrow = nrow(orderF))   #part of each GO id in the GO type
    for (I in 1:nrow(orderF)) 
    {
        resF = rbind(resF, resgoF[orderF[I, ], 1:4])
        prop[I, 1] = round((as.numeric(resF[resF[, 1] == resgoF[orderF[I, ], 1], 4])/sumF) * 100, 2)
        prop[I, 2] = paste(resF[I, 2], " (", prop[I, 1], "%, ", resF[I, 4], ")", sep = "")
    }
    resF = cbind(resF, prop)
    
    orderC = as.matrix(order(as.matrix(as.numeric(resgoC[, 4])), decreasing = TRUE))  #sort the Molecular Function by score
    tC = as.matrix(unlist(strsplit(resgoC[orderC[1, ], 5], ",")))  #matrix of the longuest list of genes for a GO identifier 
    tC = as.matrix(tC[2:nrow(tC), ]) #remove the first row which is empty
    detailC = matrix("", ncol = nrow(orderC), nrow = nrow(tC))  #matrix of all genes for each GO id
    sumC = sum(as.numeric(resgoC[, 4]))   #Sum of genes for this GO type
    resC = matrix(0, ncol = 4, nrow = 0)
    prop = matrix(0, ncol = 2, nrow = nrow(orderC))   #part of each GO id in the GO type
    for (I in 1:nrow(orderC))
    {
        resC = rbind(resC, resgoC[orderC[I, ], 1:4])
        prop[I, 1] = round((as.numeric(resC[resC[, 1] == resgoC[orderC[I, ], 1], 4])/sumC) * 100, 2)
        prop[I, 2] = paste(resC[I, 2], " (", prop[I, 1], "%,", resC[I, 4], ")", sep = "")
        #t = as.matrix(unlist(strsplit(resgoC[orderC[I, ], 5], ",")))
        #temp = matrix("", ncol = 1, nrow = (nrow(detailC) - nrow(t) + 1))
        #t = rbind(as.matrix(t[2:nrow(t), ]), temp)
        #detailC[, I] = t
    }
    resC = cbind(resC, prop)
    #colnames(detailC) = as.vector(resC[, 1])
    
    res = paste(res, "GO.results", sep = "")
    write(paste("Go results folder : ", res, sep = ""), file = "")
    if (!dir.create(res))
    {
        dir.create(res)
    }
    
    resP = cbind(resP[, 1:2], resP[, 4:6])
    colnames(resP) = c("ID GO", "Biological Process", "Score", "Proportion", "labels")
    resF = cbind(resF[, 1:2], resF[, 4:6])
    colnames(resF) = c("ID GO", "Molecular Functions", "Score", "Proportion", "labels")
    resC = cbind(resC[, 1:2], resC[, 4:6])
    colnames(resC) = c("ID GO", "Cellular Component", "Score", "Proportion", "labels")
    
    write.table(resP, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/GO.BP.txt", sep = "")))
    #write.table(detailP, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/detail.BP.txt", sep = "")))
    write.table(resC, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/GO.CC.txt", sep = "")))
    #write.table(detailC, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/detail.CC.txt", sep = "")))
    write.table(resF, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/GO.MF.txt", sep = "")))
    #write.table(detailF, row.names = FALSE, sep = "\t", file = file.path(paste(res, "/detail.MF.txt", sep = "")))
    
    lim = 1   #minimum number genes linked with a GO to be mapped
    x = as.vector(as.numeric(resP[as.numeric(resP[, 4]) >= lim, 4]))
    labels = as.vector(resP[as.numeric(resP[, 4]) >= lim, "labels"])
    
    pdf(file = file.path(paste(res, "/GO_BP.pdf", sep = "")), width = 150, height = 80, family = "Helvetica", title = "Biological Processes")
      pie(x, labels, col = palette(rainbow(100)), cex = 8)
      text(0, 1, paste("Biological Processes (", sum(x), "%)", sep = ""), cex = 15)
    dev.off()
    
    write(paste(nrow(resgoP), " Biological Processes have have been identified in the list", sep = ""), file = "")
    cat("\n", paste(nrow(resgoP), " Biological Processes have have been identified in the list", sep = ""), file = resume)
    cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep = ""), file = resume)
    x = as.vector(as.numeric(resF[as.numeric(resF[, 4]) >= lim, 4]))
    labels = as.vector(resF[as.numeric(resF[, 4]) >= lim, "labels"])
    
    pdf(file = file.path(paste(res, "/GO_MF.pdf", sep = "")), width = 150, height = 80, family = "Helvetica", title = "Molecular Functions")
      pie(x, labels, col = palette(rainbow(100)), cex = 8)
      text(0, 1, paste("Molecular Functions, (", sum(x), "%)", sep = ""), cex = 15)
    dev.off()
    
    write(paste(nrow(resgoF), " Molecular Functions have have been identified in the list", sep = ""), file = "")
    cat("\n", paste(nrow(resgoF), " Molecular Functions have have been identified in the list", sep = ""), file = resume)
    cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep = ""), file = resume)
    x = as.vector(as.numeric(resC[as.numeric(resC[, 4]) >= lim, 4]))
    labels = as.vector(resC[as.numeric(resC[, 4]) >= lim, "labels"])
    
    pdf(file = file.path(paste(res, "/GO_CC.pdf", sep = "")), width = 150, height = 80, family = "Helvetica", title = "Cellular Component")
      pie(x, labels, col = palette(rainbow(100)), cex = 8)
      text(0, 1, paste("Cellular Component (", sum(x), "%)", sep = ""), cex = 15)
    dev.off()
    
    write(paste(nrow(resgoC), " Cellular Component have have been identified in the list", sep = ""), file = "")
    cat("\n", paste(nrow(resgoC), " Cellular Component have have been identified in the list", sep = ""), file = resume)
    cat("\n", paste(length(x), " GO id represent ", sum(x), "%", sep = ""), file = resume)
    close(resume)
    write(paste("GO analyze OK", sep = ""), file = "")
}


################################################################################
################################################################################
################################################################################


  write("###############################################################################", file = "")
  write("#                                                                             #", file = "")
  write("#                              GExMap 1.1.3                                   #", file = "")
  write("#                                                                             #", file = "")
  write("###############################################################################", file = "")
  
  options(warn=-1)

	###############################################################################
	#                                                                             #
	#                              INPUT PARAMETERS                               #
	#                                                                             #
	#    scale = width of units (Mbp)                                             #
	#    path = Source folder                                                     #
	#    res = main results folder                                                #
	#                                                                             #
	###############################################################################
  
  library("multtest")
  if (scale == "")
  {
    scale = 1e+06
    write("Default scale: 1 unit = 1 000 000 bp", file = "")
  } else
  {
    write(paste("Scale unite =", scale, sep = ""), file = "")
  }
  #Default source folder
  if (source == "")
  {
    source = paste(getwd(), "/library/GExMap/data/", sep = "")
    write(paste("Default source folder : ", source, sep = ""), file = "")
  } else
  {
    write(paste("Custom source folder :", source, sep = ""), file = "")
  }
  #Results folder
  if (res == "")
  {
    res = paste(getwd(), "/GExMap", format(Sys.time(), "(%H-%M-%S) %a%d %b%Y"), sep = "")
    write(paste("Default results folder : ", res, sep = ""), file = "")
  } else
  {
    write(paste("Custom results folder :", res, sep = ""), file = "")
  }
  dir.create(res)
  res = paste(res, "/", sep = "")
  
  ###############################################################################
  #                                                                             #
  #                               UPLOAD LISTE                                  #
  #                                                                             #
  ###############################################################################
  
  write("Upload of list to analyze", file = "")
  flush.console()
  if (source == "test")
  {
    data(list)
  } else
  {
    list = read.table(file = file.choose(), sep = "\t", header = FALSE, row.names = 1)
  }
  names = c(rownames(list)[1], as.vector(list[1, ]))
  col = as.matrix(rownames(list))
  col = as.matrix(col[2:nrow(col), ])
  list = cbind(col, as.matrix(list[2:nrow(list), ]))
  colnames(list) = names
  
  write(paste("Loading the ", genome, " genome Rdata file", sep = ""), file = "")
  flush.console()
  if (source == "test")
  {
    data(data)
  } else
  {
    data = gexload.data(source, genome)
  }
  
  #
  ######################################
  # upload of correspondances file
  ######################################
  
  mtype = unlist(strsplit(colnames(list)[1], ","))[2]
  ttype = unlist(strsplit(colnames(list)[1], ","))[1]
  if (length(agrep(",", colnames(list)[1], ignore.case = TRUE)) != 0)
  {
    write(paste("Loading the ", ttype, " ", mtype, " file", sep = ""), file = "")
    Rdata.dir = file.path(paste(source, ttype, "-", mtype, ".Rdata", sep = ""))
  } else
  {
    write(paste("Loading the ", ttype, ".Rdata file", sep = ""), file = "")
    Rdata.dir = file.path(paste(source, ttype, ".Rdata", sep = ""))
  }
  
  if (source == "test")
  {
      corr = data(corr)
  } else
  {
    if (!file.exists(Rdata.dir))
    {
      write(sprintf("ERROR: file %s does not exist in the default directory", source), file = "")
      write(paste("Load manually the ", ttype, "-", mtype,".Rdata file", sep = ""), file = "")
      flush.console()
      load(file = file.choose())
    } else
    {
      load(file = Rdata.dir)
      write("Rdata file Loaded")
    }
  }
  
  ######################################
  # Matrix of correspondences
  ######################################
  
  absents_corr = matrix(0, ncol = 1, nrow = 0)  #matrix of the identifiers not found in the data.Rdata matrix
  write("Looking for correspondences", file = "")
  #deb = format(Sys.time(), "%H:%M:%S")
  id = as.matrix(unlist(apply(list, 1, function(x) as.matrix(corr[corr[, "probes"] == x[1], "ensembl"]))))
  #write(paste("debut:", deb, ", fin:", format(Sys.time(), "%H:%M:%S"), sep=""), file="")
  id = cbind(id, as.matrix(rep(1, nrow(id))))
  
  write("Correspondences OK", file = "")
  if (nrow(absents_corr) != 0)
  {
    write.table(absents_corr, row.names = FALSE, sep = "\t", file = paste(res, "absents_corr.txt", sep = ""))
    write(paste(nrow(absents_corr), " identifiers unknown with no ENSEMBL correspondences have been placed in the absents_corr.txt file", sep = ""), file = "")
  }
  
  # Suppression of the ENSG00000000000 identifiers which corresponds to unidentified probes
  id = as.matrix(id[id[, 1] != "ENSG00000000000", ])
  write(paste("Nbre of identified probes=", nrow(id), sep = ""), file = "")
  
  # looking for localisation information
  write("Collecting localisation data", file = "")
  tbl_scale = matrix(0, ncol = 1, nrow = nrow(id))
  absents_data = matrix(0, ncol = 1, nrow = 0)
  id_map = matrix(0, ncol=ncol(data), nrow=0)
  
  for(K in 1:nrow(id))
  {
    if(is.null(nrow(data[data[,"ensembl"]==id[K,1],])))  id_map = rbind(id_map, t(data[data[,"ensembl"]==id[K,1],]))
  }
  
  #id_map = rbind(id_map, apply(id, 1, function(x) t(as.matrix(data[data[,"ensembl"]==x[1],]))))
  id_map =  cbind(id_map, as.matrix(rep(1, nrow(id_map))))
  colnames(id_map) = c("chr", "cytoband", "ensembl", "name", "m", "GO", "expression")        
  list.ENS = id_map[,c("chr", "m", "cytoband", "name", "ensembl", "GO", "expression")]
  
  write("Localisation data OK", file = "")
  write("Saving data", file = "")
  write.table(id_map, row.names = FALSE, sep = "\t", file = file.path(paste(res, "report.txt", sep = "")))
  
  ######################################################################################################
  #                                                                                                    #
  #                             FORMATTING GENOME TABLE                                                #
  #                                                                                                    #
  #   chr = chromosome                                                                                 #
  #   unit = scale unit (million of bp...)                                                             #
  #   total = Nbre of genes by unit for the ENSEMBL genome                                             #
  #   list = Nbr of genes by unit for the tested list                                                  #
  #   hasard = Nbr of genes by unit for the tested gene list expected by chance                        #
  #                                                                                                    #                                                                           
  #   genome = matrix of all data for graphics                                                         #
  #   ENS = ENSEMBL gross data matrix for the current chromosome                                       #
  #   genome.temp = matrix of the number of ENSEMBL genes of each unit for the current chromosome      #   
  #   liste.chr = matrix of the number of genes from the list of each unit for the current chromosome  #
  #   list.ENS = all data about the gene list                                                          #
  #                                                                                                    #
  ######################################################################################################
  
  if (isMAP)
  {
    genome = matrix(ncol = 6, nrow = 0, dimnames = NULL)
    colnames(genome) = c("chr", "unit", "total", "list", "up", "down")
    
    # matrix of number of ENSEMBL genes by unit
    for (I in levels(factor(data[, "chr"])))   # chromosome par chromosome
    {
      write(paste("Treatment of Chr ", I, sep = ""), file = "")
      list.chr = list.ENS[list.ENS[, "chr"] == I, ]      #extracts list data for the current chromosome
      
      #If the list contains only one gene the defaul class is "character"
      if (!is.matrix(list.chr))
      {
        list.chr = t(as.matrix(list.chr))
      } else   # matrix of only one gene
      {
        list.chr = as.matrix(list.chr[order(list.chr[, "m"]), ])      #sorts conlumn by m
      }
      if (nrow(list.chr) > 5)   #If there is at least five gene on the current chromosome
      {
        genome.temp = gex.mapping(data, list.chr, scale, I, list.ENS)
      } else 
      {
        genome.temp = matrix(ncol = ncol(genome), nrow = 0)
      }
      genome = rbind(genome, genome.temp)
    }
    
    #cytobands matrix creation
    cytobands = matrix(0, ncol=4, nrow=0)
    colnames(cytobands) = c("chr", "cytoband", "debut", "fin")
    for (I in levels(as.data.frame(data)[,"chr"]))
    {
      chr_tmp = as.data.frame(data[data[,"chr"]==I,])
      for (J in levels(factor(chr_tmp[, "cytoband"])))
      {
        cytobands = rbind(cytobands, matrix(0, ncol=4, nrow=1))
        cytobands[nrow(cytobands), "chr"] = I
        cytobands[nrow(cytobands), "cytoband"] = J
        cytobands[nrow(cytobands), "debut"] = min(as.matrix(chr_tmp[chr_tmp[,"cytoband"]==J, "m"]))
        cytobands[nrow(cytobands), "fin"] = max(as.matrix(chr_tmp[chr_tmp[,"cytoband"]==J, "m"]))
      }
    }
    
    #TESTS using two computation of the hazard: global then specific for each one of the two test packages
    #Resulting files are placed in separated folders
   
    #Estimation of hazard global then specific
    
    chance.compute = matrix(0, ncol = 2, nrow = 2)
    colnames(chance.compute) = c("nhazard", "genome.graph")
    rownames(chance.compute) = c("specific method", "global method")
    chance.compute[1, 1] = sum(as.numeric(genome[genome[, "list"] != 0, "total"]))   #number of ENSEMBL genes presents in units where there is at least one gene frome the gene list
    chance.compute[2, 1] = sum(as.numeric(genome[, "total"]))
    chance.compute[1, 2] = paste(res, "specific method/", sep = "")
    chance.compute[2, 2] = paste(res, "global method/", sep = "")
    for (W in 1:2)
    {
      write(paste("Hazard computation ", rownames(chance.compute)[W], sep = ""), file = "")
      hazard = matrix(0, ncol = 1, nrow = nrow(genome))
      colnames(hazard) = "hazard"
      nhazard = chance.compute[W, 1]
      results.graph = chance.compute[W, 2]
      for (H in 1:nrow(genome))
      {
        hazard[H, "hazard"] = (sum(as.numeric(genome[, "list"]))/as.numeric(nhazard)) * as.numeric(genome[H, "total"])  #computes number of genes expected by chance
      }
      
      if (!dir.create(results.graph))
      {
        dir.create(results.graph)
      }
      
      genome.test = cbind(genome, hazard)
      
      #statistic method
      genome.test = gextest(genome.test, nhazard, results.graph, lim_chi, global_test_choice, pcorrd, pcorrv)
      write.table(genome.test, row.names = FALSE, sep = "\t", file = file.path(paste(results.graph, "genome.txt", sep = "")))
      gexgraph(genome.test, scale, results.graph, W, cytobands)
      write(paste("Results have been saved in the folder: ", results.graph, sep = ""), file = "")
    }
  }
  #GO
  if (isGO)
  {
    gexgo(list.ENS, source, res)
  }
}

#gexmap(genome="homosapiens", scale="", source="", res="", isGO=F, isMAP=TRUE, lim_chi=5, global_test_choice=4, pcorrd=1, pcorrv=1)
#genome="homosapiens"; scale=""; source=""; res=""; isGO=F; isMAP=TRUE; lim_chi=5; global_test_choice=4; pcorrd=1; pcorrv=1
