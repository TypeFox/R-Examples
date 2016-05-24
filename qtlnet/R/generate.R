#####################################################################
##
## $Id: generate.R,v 2007/11/28 byandell Exp $
##
##     Copyright (C) 2007 Elias Chaibub Neto and Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
## Routines: generate.qtl.markers, generate.qtl.pheno
##############################################################################
## Generating data for some DAGs examples
## This code is specific for particular graphs, and is not meant
## to provide general tools for generating directed graphs.

#########################################
generate.qtl.markers <- function(cross, n.phe, nqtl = 3)
{
  ## randomly selects 2 or 3 markers (per phenotype) 
  nqtl <- array(nqtl, n.phe)
  allqtls <- list()
  markers <- list()
  for(i in 1:n.phe){
    if(nqtl[i]==2){
      chrm <- sample(c(1:20), 2, replace = TRUE)
      position <- sample(c(1:10), 2, replace = FALSE)
      position <- c(cross$geno[[chrm[1]]]$map[position[1]],
                    cross$geno[[chrm[2]]]$map[position[2]])
    }
    else{
      chrm <- sample(c(1:20), 3, replace = TRUE)
      position <- sample(c(1:10), 3, replace = FALSE)
      position <- c(cross$geno[[chrm[1]]]$map[position[1]],
                    cross$geno[[chrm[2]]]$map[position[2]],
                    cross$geno[[chrm[3]]]$map[position[3]])	
    }
    allqtls[[i]] <- makeqtl(cross, chr = chrm, pos = position)
    markers[[i]] <- find.marker(cross, chr = chrm, pos = position)
  }
  names(allqtls) <- paste("y", 1:n.phe, sep = "")
  names(markers) <- paste("y", 1:n.phe, sep = "")
  list(allqtl = allqtls, markers = markers)
}

##################################################################
generate.qtl.pheno <- function(name = c("acyclic","acyc2or3","cyclica","cyclicb","cyclicc"),
                               cross,
                               bp, bq, stdev, allqtl,
                               burnin = 2000, geno)
{
  name <- match.arg(name)
  switch(name,
         acyclic = { generate.data(cross, bp, bq, stdev, allqtl) },
         acyc2or3 = { generate.data.2or3(cross, bp, bq, stdev, allqtl) },
         cyclica = { generate.data.graph.a(cross, burnin, bq, bp, stdev, geno) },
         cyclicb = { generate.data.graph.b(cross, burnin, bq, bp, stdev, geno) },
         cyclicc = { generate.data.graph.c(cross, burnin, bq, bp, stdev, geno) })
}

##################################################################
geno.effect <- function(nodes, bq, bp, stdev, allqtl, n, y)
{
  i <- nodes[1]
  out <- bq[i, allqtl[[i]]$geno[,1,]] + bq[i, allqtl[[i]]$geno[,2,]] +
    bq[i, allqtl[[i]]$geno[,3,]] + rnorm(n, 0, stdev[i])
  len <- length(nodes)
  if(len > 1) {
    for(j in seq(2, len))
      out <- out + bp[i] * y[,j]
  }
  out
}
##################################################################
## Acyclic example (100 phenotypes network)
##
generate.data <- function(cross,bp,bq,stdev,allqtl)
{
  n <- length(cross$pheno[,1])
  y <- matrix(0,n,100)

  ## Ordering is important in terms of the dependencies in the graph.
  ## This is set up for a specific graph used in the paper.

  node.parents <- list(c(6), ## founder nodes
                       c(10),
                       c(21),
                       c(22),
                       c(24),
                       c(27),
                       c(38),
                       c(48),
                       c(89),
                       c(98),
                       c(100),
                       c(1),
                       c(16),
                       c(53),
                       c(43),
                       c(11),
                       c(8),
                       c(18),
                       c(26),
                       c(14),
                       c(80),
                       c(3),
                       c(2),
                       c(4),
                       c(85),
                       c(50),
                       c(9),
                       c(44),
                       c(57),
                       c(7),
                       c(49),
                       c(64),
                       c(5),
                       c(32),
                       c(47),
                       c(59),
                       c(61,6), ## node, parents
                       c(13,6),
                       c(34,27),
                       c(87,1),
                       c(31,1),
                       c(29,16),
                       c(23,16),
                       c(70,1,53,43),
                       c(19,11),
                       c(55,11,13),
                       c(60,34),
                       c(17,8),
                       c(88,8),
                       c(81,8,31),
                       c(40,29),
                       c(28,23),
                       c(75,16,70,11),
                       c(65,29,18),
                       c(25,18,19),
                       c(46,19),
                       c(33,26),
                       c(37,26),
                       c(56,1,40),
                       c(30,17),
                       c(91,14,80,65),
                       c(39,29),
                       c(71,4),
                       c(15,3),
                       c(63,40,2,4),
                       c(94,65,4,85,19),
                       c(54,50),
                       c(74,50),
                       c(51,18,46),
                       c(68,46),
                       c(82,4,9,44),
                       c(45,44),
                       c(42,33),
                       c(41,37),
                       c(90,57,56),
                       c(36,30),
                       c(96,7),
                       c(35,7),
                       c(83,30,7,39,71,49,54),
                       c(20,3),
                       c(78,15,43,51),
                       c(73,63),
                       c(86,63,64),
                       c(62,51),
                       c(52,49),
                       c(12,7,5),
                       c(72,5),
                       c(67,54,19),
                       c(77,45),
                       c(66,36),
                       c(97,41,40,96,83),
                       c(84,40,20,73),
                       c(93,39,78),
                       c(79,32,62),
                       c(76,47,52),
                       c(58,12),
                       c(99,4,79),
                       c(95,76,50),
                       c(69,58,50,43,59),
                       c(92,59))
  for(i in seq(length(node.parents)))
    y[,i] <- geno.effect(node.parents[[i]], bq, bp, stdev, allqtl, n, y)

  y <- data.frame(y)
  names(y) <- paste("y",1:100,sep="")
  cross$pheno <- y
  return(cross)
}

###########################################################
## Actual example used in paper with two or three QTL per node.
##
generate.data.2or3 <- function(cross, bp, bq, stdev, allqtl)
{
  n <- length(cross$pheno[,1])
  y <- matrix(0,n,100)

  ## Ordering is important in terms of the dependencies in the graph.
  ## This is set up for a specific graph used in the paper.

  node.parents <- list(c(6), ## founder nodes
                       c(10),
                       c(21),
                       c(22),
                       c(24),
                       c(27),
                       c(38),
                       c(48),
                       c(89),
                       c(98),
                       c(100),
                       c(1),
                       c(16),
                       c(53),
                       c(43),
                       c(11),
                       c(8),
                       c(18),
                       c(26),
                       c(14),
                       c(80),
                       c(3),
                       c(2),
                       c(4),
                       c(85),
                       c(50),
                       c(9),
                       c(44),
                       c(57),
                       c(7),
                       c(49),
                       c(64),
                       c(5),
                       c(32),
                       c(47),
                       c(59),
                       c(61,6), ## node, parents
                       c(13,6),
                       c(34,27),
                       c(87,1),
                       c(31,1),
                       c(29,16),
                       c(23,16),
                       c(70,1,53,43),
                       c(19,11),
                       c(55,11,13),
                       c(60,34),
                       c(17,8),
                       c(88,8),
                       c(81,8,31),
                       c(40,29),
                       c(28,23),
                       c(75,16,70,11),
                       c(65,29,18),
                       c(25,18,19),
                       c(46,19),
                       c(33,26),
                       c(37,26),
                       c(56,1,40),
                       c(30,17),
                       c(91,14,80,65),
                       c(39,29),
                       c(71,4),
                       c(15,3),
                       c(63,40,2,4),
                       c(94,65,4,85,19),
                       c(54,50),
                       c(74,50),
                       c(51,18,46),
                       c(68,46),
                       c(82,4,9,44),
                       c(45,44),
                       c(42,33),
                       c(41,37),
                       c(90,57,56),
                       c(36,30),
                       c(96,7),
                       c(35,7),
                       c(83,30,7,39,71,49,54),
                       c(20,3),
                       c(78,15,43,51),
                       c(73,63),
                       c(86,63,64),
                       c(62,51),
                       c(52,49),
                       c(12,7,5),
                       c(72,5),
                       c(67,54,19),
                       c(77,45),
                       c(66,36),
                       c(97,41,40,96,83),
                       c(84,40,20,73),
                       c(93,39,78),
                       c(79,32,62),
                       c(76,47,52),
                       c(58,12),
                       c(99,4,79),
                       c(95,76,50),
                       c(69,58,50,43,59),
                       c(92,59))
  
  for(i in seq(length(node.parents)))
    y[,i] <- geno.effect(node.parents[[i]], bq, bp, stdev, allqtl)

  y <- data.frame(y)
  names(y) <- paste("y",1:100,sep="")
  cross$pheno <- y
  cross
}

###############
###############
## cyclic graphs 
###############
###############

###################################################
## compute the means (that depend on the genotypes) 
## of each individual. These means are used in the
## generating of the phenotype data 
##
compute.mu <- function(ind.geno,bq){
  mu <- rep(0, 6)
  for(i in 1:6)
    mu[i] <- sum(bq[ind.geno[(3 * (i - 1)) + (1:3)]])
  mu
}

############################################################
## generate the phenotype data. Each data point is generated 
## by a separate Markov chain
##
generate.data.graph.a <- function(cross,burnin,bq,bp,stdev,geno)
{
  ## Gibbs sampler for graph (a) 
  gibbs.graph.a <- function(n, burnin, bp, stdev, mu){
    phi6 <- stdev[6]
    aux2 <- stdev[4] + stdev[2] * bp[4,2]^2
    phi2 <- stdev[2] * stdev[4] / aux2
    aux4 <- stdev[5] + stdev[4] * bp[5,4]^2
    phi4 <- stdev[4] * stdev[5] / aux4
    aux5 <- stdev[2] * stdev[6] +
            stdev[5] * stdev[6] * (bp[2,5]^2) +
            stdev[2] * stdev[5] * (bp[6,5]^2)
    phi5 <- stdev[2] * stdev[5] * stdev[6] / aux5
    aux1 <- stdev[2] + stdev[1] * bp[2,1]^2
    phi1 <- stdev[1] * stdev[2] / aux1
    aux3 <- stdev[4] + stdev[3] * bp[4,3]^2
    phi3 <- stdev[3] * stdev[4] / aux3
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0, n + burnin)
    for(i in 2:(n + burnin)){
      m1 <- (stdev[2] * mu[1] +
             stdev[1] * bp[2,1] * (y2[i-1] - mu[2] - bp[2,5] * y5[i-1])) / aux1
      y1[i] <- rnorm(1,m1,sqrt(phi1))

      m2 <- (stdev[4] * (mu[2] + bp[2,1] * y1[i] + bp[2,5] * y5[i-1]) +
             stdev[2] * bp[4,2] * (y4[i-1] - mu[4] - bp[4,3] * y3[i-1])) / aux2
      y2[i] <- rnorm(1,m2,sqrt(phi2))

      m3 <- (stdev[4] * mu[3] +
             stdev[3] * bp[4,2] * (y4[i-1] - mu[4] - bp[4,2] * y2[i])) / aux3
      y3[i] <- rnorm(1,m3,sqrt(phi3))

      m4 <- (stdev[5] * (mu[4] + bp[4,2] * y2[i] + bp[4,3] * y3[i]) +
             stdev[4] * bp[5,4] * (y5[i-1] - mu[5])) / aux4
      y4[i] <- rnorm(1,m4,sqrt(phi4))

      m5 <- (stdev[2] * stdev[6] * (mu[5] + bp[5,4] * y4[i]) +
             bp[2,5] * stdev[5] * stdev[6] * (y2[i] - mu[2] - bp[2,1] * mu[1]) +
             bp[6,5] * stdev[2] * stdev[5] * (y6[i-1] - mu[6])) / aux5
      y5[i] <- rnorm(1,m5,sqrt(phi5))

      m6 <- mu[6]+bp[6,5]*y5[i]
      y6[i] <- rnorm(1,m6,sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    y[i,] <- as.vector(gibbs.graph.a(n=1,burnin=burnin,bp,stdev,mu=mu))
  }
  y <- data.frame(y)
  names(y) <- c("y1","y2","y3","y4","y5","y6")
  cross$pheno <- y
  return(cross)
} 


###############################
## Cyclic toy example (graph b)
generate.data.graph.b <- function(cross,burnin,bq,bp,stdev,geno)
{
  gibbs.graph.b <- function(n,burnin,bp,stdev,mu){
    aux1 <- stdev[2]*stdev[3]+stdev[1]*stdev[3]*(bp[2,1]^1)+stdev[1]*stdev[2]*(bp[3,1]^2)
    phi1 <- stdev[1]*stdev[2]*stdev[3]/aux1
    aux2 <- stdev[4]+stdev[2]*bp[4,2]^2
    phi2 <- stdev[2]*stdev[4]/aux2
    aux3 <- stdev[6]+stdev[3]*bp[6,3]^2
    phi3 <- stdev[3]*stdev[6]/aux3
    aux4 <- stdev[5]+stdev[4]*bp[5,4]^2
    phi4 <- stdev[4]*stdev[5]/aux4
    aux5 <- stdev[1]+stdev[5]*bp[1,5]^2
    phi5 <- stdev[1]*stdev[5]/aux5
    aux6 <- stdev[5]+stdev[6]*bp[5,6]^2
    phi6 <- stdev[5]*stdev[6]/aux6
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0,n+burnin)
    for(i in 2:(n+burnin)){
      m1 <- (stdev[2]*stdev[3]*(mu[1]+bp[1,5]*y5[i-1])+stdev[1]*stdev[3]*bp[2,1]*(y2[i-1]-mu[2])+stdev[1]*stdev[2]*bp[3,1]*(y3[i-1]-mu[3]))/aux1
      y1[i] <- rnorm(1,m1, sqrt(phi1))

      m2 <- (stdev[4]*(mu[2]+bp[2,1]*y1[i])+stdev[2]*bp[4,2]*(y4[i-1]-mu[4]))/aux2
      y2[i] <- rnorm(1, m2, sqrt(phi2))

      m3 <- (stdev[6]*(mu[3]+bp[3,1]*y1[i])+stdev[3]*bp[6,3]*(y6[i-1]-mu[6]))/aux3
      y3[i] <- rnorm(1, m3, sqrt(phi3))

      m4 <- (stdev[5]*(mu[4]+bp[4,2]*y2[i])+stdev[4]*bp[5,4]*(y5[i-1]-mu[5]-bp[5,6]*y6[i-1]))/aux4
      y4[i] <- rnorm(1, m4, sqrt(phi4))

      m5 <- (stdev[1]*(mu[5]+bp[5,4]*y4[i]+bp[5,6]*y6[i-1])+stdev[5]*bp[1,5]*(y1[i]-mu[1]))/aux5
      y5[i] <- rnorm(1, m5, sqrt(phi5))

      m6 <- (stdev[5]*(mu[6]+bp[6,3]*y3[i])+stdev[6]*bp[5,6]*(y5[i]-mu[5]-bp[5,4]*y4[i]))/aux6
      y6[i] <- rnorm(1, m6, sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    y[i,] <- as.vector(gibbs.graph.b(n=1,burnin=burnin,bp,stdev,mu=mu))
  }
  y <- data.frame(y)
  names(y) <- c("y1","y2","y3","y4","y5","y6")
  cross$pheno <- y
  return(cross)
} 


###############################
## Cyclic toy example (graph c)
generate.data.graph.c <- function(cross,burnin,bq,bp,stdev,geno)
{
  gibbs.graph.c <- function(n,burnin,bp,stdev,mu){
    aux1 <- stdev[2]+stdev[1]*bp[2,1]^2
    phi1 <- stdev[1]*stdev[2]/aux1
    aux2 <- stdev[3]*stdev[5]+stdev[2]*stdev[5]*(bp[3,2]^2)+stdev[2]*stdev[3]*(bp[5,2]^2)
    phi2 <- stdev[2]*stdev[3]*stdev[5]/aux2
    phi3 <- stdev[3]
    aux4 <- stdev[5]+stdev[4]*bp[5,4]^2
    phi4 <- stdev[4]*stdev[6]/aux4
    aux5 <- stdev[2]*stdev[6]+stdev[5]*stdev[6]*(bp[2,5]^2)+stdev[2]*stdev[5]*(bp[6,5]^2)
    phi5 <- stdev[2]*stdev[5]*stdev[6]/aux5
    phi6 <- stdev[6]
    y1 <- y2 <- y3 <- y4 <- y5 <- y6 <- rep(0,n+burnin)
    for(i in 2:(n+burnin)){
      m1 <- (stdev[2]*mu[1]+stdev[1]*bp[2,1]*(y2[i-1]-mu[2]-bp[2,5]*y5[i-1]))/aux1
      y1[i] <- rnorm(1, m1, sqrt(phi1))
      
      m2 <- (stdev[3]*stdev[5]*(mu[2]+bp[2,1]*y1[i]+bp[2,5]*y5[i-1])+stdev[2]*stdev[5]*bp[3,2]*(y3[i-1]-mu[3])+stdev[2]*stdev[3]*bp[5,2]*(y5[i-1]-mu[5]-bp[5,4]*y4[i-1]))/aux2
      y2[i] <- rnorm(1, m2, sqrt(phi2))

      m3 <- mu[3]+bp[3,2]*y2[i]
      y3[i] <- rnorm(1, m3, sqrt(phi3))

      m4 <- (stdev[5]*mu[4]+stdev[4]*bp[5,4]*(y5[i-1]-mu[5]-bp[5,2]*y2[i]))/aux4
      y4[i] <- rnorm(1, m4, sqrt(phi4))

      m5 <- (stdev[2]*stdev[6]*(mu[5]+bp[5,4]*y4[i]+bp[5,2]*y2[i])+stdev[5]*stdev[6]*bp[2,5]*(y2[i]-mu[2]-bp[2,1]*y1[i])+stdev[2]*stdev[5]*bp[6,5]*(y6[i-1]-mu[6]))/aux5
      y5[i] <- rnorm(1, m5, sqrt(phi5))

      m6 <- mu[6]+bp[6,5]*y5[i]
      y6[i] <- rnorm(1, m6, sqrt(phi6))
    }	
    y <- cbind(y1,y2,y3,y4,y5,y6)
    return(y[-c(1:burnin),])
  }

  n <- length(cross$pheno[,1])
  y <- matrix(0,n,6)
  for(i in 1:n){
    mu <- compute.mu(ind.geno=geno[i,],bq=bq)
    tmp <- gibbs.graph.c(n=1, burnin, bp, stdev, mu)
    y[i,] <- as.vector(gibbs.graph.c(n=1, burnin, bp, stdev, mu))
  }
  y <- data.frame(y)
  names(y) <- paste("y", 1:6, sep = "")
  cross$pheno <- y
  return(cross)
} 
