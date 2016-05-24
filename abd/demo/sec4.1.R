# Figure 4.1-1
histogram(~gene.length, HumanGeneLengths, 
	subset=gene.length < 15000, 
	xlab="Gene length (number of nucleotides)",
	n=30)

# Table 4.1-1
mean(HumanGeneLengths$gene.length)
# this computes the SAMPLE standard deviation
sd(HumanGeneLengths$gene.length)
# For the population variance, SS should divide by n rather than n-1
n <- nrow(HumanGeneLengths)
pop.sd <- sd(HumanGeneLengths$gene.length) * sqrt((n-1)/n)
pop.sd

# Figure 4.1-2
HundredGenes <- sample(HumanGeneLengths,100)
histogram(~gene.length, HundredGenes, 
	subset=gene.length < 15000, 
	xlab="Gene length (number of nucleotides)",
	n=30)

# Table 4.1-2 (won't match exactly because of random sampling)
print(favstats(HundredGenes$gene.length))

# Figure 4.1-3
sampleMeans <- replicate(5000, mean(sample(HumanGeneLengths,100)$gene.length))
histogram(sampleMeans, type="percent", n = 40)

# Figure 4.1-4
sampleMeans20 <- replicate (5000, mean(sample(HumanGeneLengths,20)$gene.length))
sampleMeans500 <- replicate (5000, mean(sample(HumanGeneLengths,500)$gene.length))
means <- data.frame(
	mean = c( sampleMeans, sampleMeans20, sampleMeans500 ),
	size = rep(c(100,20,500), each=5000)
	)
histogram(~mean | paste("n =",size), data=means, 
          type="density", n=100, xlim=c(1200,5200),
          layout=c(1,3))
