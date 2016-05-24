`tileplot.double` <-
function(genesonchip, array1data, array2data, annotationslist, cutoff=-1, cutoff_multiplier=3, outputfile, graphdirectory, outputtable, array1name = "Array 1", array2name="Array 2", loess=TRUE, smoothing_factor = 6)
{
	
#The cluster file is the CD-HIT output piped through a python script to make each line a cluster (a list of identifiers)
#clusters = scan(file=clusterfile, sep="\n", what="raw")

#This for loop performs some further manipulation on that cluster file
#for(i in 1:length(clusters))
#	{
#	clusters[i] = strsplit(as.character(clusters[i]), "\t")
#	}

#The genes file should be a simple list of unique gene identifiers without a header. I use IMG identifiers as they are common for both DNA and protein sequences, making life easier.
genes = read.table(file=genesonchip)

#The array data file is a simple, two-column file containing the probe identifiers in the first column and the median hybridization intensities in the second column.  Probe identifiers should take the form "geneidentifier-probenumber", and be sorted according to probe number so that they land in the correct 5' to 3' order.
array1 = read.table(file=array1data)
array2 = read.table(file=array2data)

if(cutoff==-1)
{
	cutoff=cutoff_multiplier*median(array1[,2])
}
cat("Cutoff calculated as",cutoff,"\n")

#These next lines perform a loess normalization of the loess data (i.e. find the polynomial function that best fits the data, straightens it out to a linear relationship with a slope of 1, then adjusts all data points to that slope)
setwd(graphdirectory)

if(loess==TRUE)
{
array2.loess <- loess(log(y) ~ log(x), span=0.2, degree=2, data.frame(x=array1[,2], y=array2[,2]))
array2.predict <- predict(array2.loess, data.frame(x=array1[,2]))
array2.predict.notlog = 2.71828183^array2.predict
array2adjusted = (array1[,2]/array2.predict.notlog)*array2[,2]
array2[,2] <- array2adjusted
write.table(array2, file="normalized_array2")
}



#The annotations file is a list of annotations containing the gene identifier somewhere in each annotation.
annotations = scan(file=annotationslist, what="list", sep="\n")

#The sequences should each cover a single line and come in the same order as the genes file. 
#sequences = read.table(file=allsequences)

#This is the cutoff hybridization intensity value used for calculating the bright probe fraction.
#cutoff = 150

#The following creates a vector of the number of probes for each gene.
probenums = vector()
for(i in 1:dim(genes)[1])
	{
	probes = grep(as.character(genes[i,]), array1[,1])
	probenums[i] = length(probes)
	}

# In the "probe_matrix" each row is a gene, and each column corresponds to the vertical location of the probes for that gene in the correct order in the array matrix. 
probe_matrix = matrix(nrow=dim(genes)[1],ncol=max(probenums))

for(i in 1:dim(genes)[1])
	{
	probes = grep(as.character(genes[i,]), array1[,1])
	for(j in 1:length(probes))
		{
		probe_matrix[i,j] = probes[j]
		}
	}

#This creates a "chunk matrix" ready to fill with the chunk (bright segment) lengths at each location, but starts out with filling it all with zeros
chunk_matrix1 = matrix(nrow=dim(genes)[1],ncol=max(probenums))

for(i in 1:dim(chunk_matrix1)[1])
	{
	for(j in 1:dim(chunk_matrix1)[2])
		{
		chunk_matrix1[i,j] = 0
		}
	}

chunk_matrix2 = chunk_matrix1

bright_probe_fraction1 = vector()
bright_probe_fraction2 = vector()

bright_gene_means1 = vector()
bright_gene_means2 = vector()

bright_gene_medians1 = vector()
bright_gene_medians2 = vector()

#For each chunk, the starting and stopping point along the sequence is defined in these two matrices...
start_chunk_blast1 = matrix(nrow=dim(genes)[1],ncol=max(probenums))
start_chunk_blast2 = matrix(nrow=dim(genes)[1],ncol=max(probenums))
stop_chunk_blast1 = matrix(nrow=dim(genes)[1],ncol=max(probenums))
stop_chunk_blast2 = matrix(nrow=dim(genes)[1],ncol=max(probenums))

#...so that this sequence_chunks matrix can be filled based on it
#sequence_chunks1 = matrix(nrow=dim(genes)[1],ncol=max(probenums))
#sequence_chunks2 = matrix(nrow=dim(genes)[1],ncol=max(probenums))


#This for loop goes through each gene one by one for array1...
for(i in 1:dim(genes)[1])
	{
	nochunk = 1
	gene = head(array1[probe_matrix[i,],2], n=probenums[i])
	gene_binary =vector(length=length(gene))
#... and then creates a "gene_binary" vector, where each probe is represented by a 1 or a 0 - 1 if it's above the cutoff line, 0 if it's below.
		for(j in 1:length(gene))
			{
				gene_binary[j] = ifelse(gene[j]>cutoff,1,0)
			}
if(sum(gene_binary>0))
{

bright_gene_values = vector(length=sum(gene_binary))
k=1
		for(j in 1:length(gene))
		{
			if(gene[j]>cutoff)
			{
				bright_gene_values[k] = as.numeric(gene[j])
				k=k+1
			}
		}
		bright_gene_means1[i] = mean(bright_gene_values)
		bright_gene_medians1[i] = median(bright_gene_values)
}

else
{
	bright_gene_means1[i] = 0
	bright_gene_medians1[i] = 0	
}
		chunk_coord=1

#Calculating the bright probe fraction is a simple question of summing the gene_binary vector and dividing by the total probe number length of the gene.
		bright_probe_fraction1[i] = sum(gene_binary)/probenums[i]

#This for loop fills out the chunk matrix.  It uses a "chunk_coord" or chunk coordinate to figure out how far along the chunk matrix it is.		
		for(j in 1:length(gene))
			{
#If the probe is bright...
			if(gene_binary[j] == 1)
				{
				chunk_matrix1[i,chunk_coord] = chunk_matrix1[i,chunk_coord] +1
#and a chunk is not currently being recorded, it adds one to the current chunk based on the chunk coordinate.					
					if(nochunk==1)
						{
						start_chunk_blast1[i, chunk_coord] = j*30 - 30
						nochunk = 0
						}
				} else
#If the probe is dark, it moves on in the gene and sets the nochunk variable back to 1 to stop recording a chunk.  In this way every probe below the line is recorded as a 0, but probes above the line are recorded in "chunks" or the length of continuous bright segments.
				{
				chunk_coord = chunk_coord + 1
					if(nochunk==0)
					{
					stop_chunk_blast1[i, chunk_coord-1] = (j-1)*30 + 30
					nochunk = 1
					} 
				}
			}

#If the last probe in the gene is bright, then it wraps up the chunk here for the stop_chunk_blast - otherwise there is not stop value inserted for these cases as the chunk seems to never end.
if(gene_binary[length(gene)] == 1)
{
stop_chunk_blast1[i, chunk_coord] = length(gene)*30
}			

#This fills out the sequence_chunks matrix, a matrix whose first column is gene identifiers, followed by the actual chunks of gene hybridization that have caused brightness.
#sequence_chunks1[i,1] = genes[i,]
#p=2
#for(j in 1:length(gene))
#	{
#	if(is.na(start_chunk_blast1[i,j]))
#	{next} else
#	{
#	sequence_chunks1[i,p] = substr(sequences[grep(genes[i,],sequences[,1]),2],start_chunk_blast1[i,j],stop_chunk_blast1[i,j])
#	p = p+1					
#	}	
#		}			
	}

#This next bit is vital - it takes the length of every chunk, and squares it! This is where the magic happens - chunks that are just one probe long will remain as one, but those that are 2 or greater will incraese in a non-linear fashion.
chunk_matrix_square1 = chunk_matrix1^2

#The squared chunk lengths are then used to make a chunk score
chunk_score1 = vector(length=dim(genes)[1])

for(i in 1:dim(genes)[1])
	{
	chunk_score1[i] <- sum(chunk_matrix_square1[i,])
	}

#The probe matrix is then sorted based on its chunk score - higher chunk scores go up the top.
probe_matrix_chunksort1 = probe_matrix[order(chunk_score1, decreasing=TRUE),]
probenums_chunksort1 = probenums[order(chunk_score1, decreasing=TRUE)]


#The chunk count matrix is used to figure out what number each chunk is, ie so if there are three bright segments in a gene, they can be labeled 1, 2, and 3.
chunk_count_matrix1 = matrix(nrow = dim(chunk_matrix1)[1], ncol = dim(chunk_matrix1)[2])

for(i in 1:dim(chunk_count_matrix1)[1])
	{
		for(j in 1:dim(chunk_count_matrix1)[2])
			{
			if(chunk_matrix1[i,j] == 0)
				{
				chunk_count_matrix1[i,j] = 0
				} else
					{
					chunk_count_matrix1[i,j] = 1
					}
			}
	}

chunk_number1 = sum(chunk_count_matrix1)


#This for loop goes through each gene one by one for array2...
for(i in 1:dim(genes)[1])
	{
	nochunk = 1
	gene = head(array2[probe_matrix[i,],2], n=probenums[i])
	gene_binary =vector(length=length(gene))
#... and then creates a "gene_binary" vector, where each probe is represented by a 1 or a 0 - 1 if it's above the cutoff line, 0 if it's below.
		for(j in 1:length(gene))
			{
				gene_binary[j] = ifelse(gene[j]>cutoff,1,0)
			}
	
			if(sum(gene_binary>0))
			{

			bright_gene_values = vector(length=sum(gene_binary))
			k=1
					for(j in 1:length(gene))
					{
						if(gene[j]>cutoff)
						{
							bright_gene_values[k] = as.numeric(gene[j])
							k=k+1
						}
					}
					bright_gene_means2[i] = mean(bright_gene_values)
					bright_gene_medians2[i] = median(bright_gene_values)
			}

			else
			{
				bright_gene_means2[i] = 0
				bright_gene_medians2[i] = 0	
			}
	
		chunk_coord=1

#Calculating the bright probe fraction is a simple question of summing the gene_binary vector and dividing by the total probe number length of the gene.
		bright_probe_fraction2[i] = sum(gene_binary)/probenums[i]

#This for loop fills out the chunk matrix.  It uses a "chunk_coord" or chunk coordinate to figure out how far along the chunk matrix it is.		
		for(j in 1:length(gene))
			{
#If the probe is bright...
			if(gene_binary[j] == 1)
				{
				chunk_matrix2[i,chunk_coord] = chunk_matrix2[i,chunk_coord] +1
#and a chunk is not currently being recorded, it adds one to the current chunk based on the chunk coordinate.					
					if(nochunk==1)
						{
						start_chunk_blast2[i, chunk_coord] = j*30 - 30
						nochunk = 0
						}
				} else
#If the probe is dark, it moves on in the gene and sets the nochunk variable back to 1 to stop recording a chunk.  In this way every probe below the line is recorded as a 0, but probes above the line are recorded in "chunks" or the length of continuous bright segments.
				{
				chunk_coord = chunk_coord + 1
					if(nochunk==0)
					{
					stop_chunk_blast2[i, chunk_coord-1] = (j-1)*30 + 30
					nochunk = 1
					} 
				}
			}

#If the last probe in the gene is bright, then it wraps up the chunk here for the stop_chunk_blast - otherwise there is not stop value inserted for these cases as the chunk seems to never end.
if(gene_binary[length(gene)] == 1)
{
stop_chunk_blast2[i, chunk_coord] = length(gene)*30
}			

#This fills out the sequence_chunks matrix, a matrix whose first column is gene identifiers, followed by the actual chunks of gene hybridization that have caused brightness.
#sequence_chunks2[i,1] = genes[i,]
#p=2
#for(j in 1:length(gene))
#	{
#	if(is.na(start_chunk_blast2[i,j]))
#	{next} else
#	{
#	sequence_chunks2[i,p] = substr(sequences[grep(genes[i,],sequences[,1]),2],start_chunk_blast2[i,j],stop_chunk_blast2[i,j])
#	p = p+1					
#	}	
#		}			
	}

#This next bit is vital - it takes the length of every chunk, and squares it! This is where the magic happens - chunks that are just one probe long will remain as one, but those that are 2 or greater will incraese in a non-linear fashion.
chunk_matrix_square2 = chunk_matrix2^2

#The squared chunk lengths are then used to make a chunk score
chunk_score2 = vector(length=dim(genes)[1])

for(i in 1:dim(genes)[1])
	{
	chunk_score2[i] <- sum(chunk_matrix_square2[i,])
	}

#The probe matrix is then sorted based on its chunk score for array 1 - higher chunk scores go up the top.
probe_matrix_chunksort2 = probe_matrix[order(chunk_score1, decreasing=TRUE),]
probenums_chunksort2 = probenums[order(chunk_score1, decreasing=TRUE)]


#The chunk count matrix is used to figure out what number each chunk is, ie so if there are three bright segments in a gene, they can be labeled 1, 2, and 3.
chunk_count_matrix2 = matrix(nrow = dim(chunk_matrix2)[1], ncol = dim(chunk_matrix2)[2])

for(i in 1:dim(chunk_count_matrix2)[1])
	{
		for(j in 1:dim(chunk_count_matrix2)[2])
			{
			if(chunk_matrix2[i,j] == 0)
				{
				chunk_count_matrix2[i,j] = 0
				} else
					{
					chunk_count_matrix2[i,j] = 1
					}
			}
	}

chunk_number2 = sum(chunk_count_matrix2)

mean_probe_intensity1 = vector()
mean_probe_intensity2 = vector()

ordered_BPF = bright_probe_fraction1[order(bright_probe_fraction1, decreasing="TRUE")]
BPF_slopes = vector()

BPF_slopes_identifier = vector()

for(i in 1:length(ordered_BPF))
{
	BPF_slopes[i] = ordered_BPF[i] - ordered_BPF[i+smoothing_factor]
	BPF_slopes_identifier[i] = (ordered_BPF[i] + ordered_BPF[i+smoothing_factor])/2
	if(i==length(ordered_BPF)-smoothing_factor)
	{
		break
	}
}

maximum_slope = max(BPF_slopes)

for(i in 1:length(BPF_slopes))
{
	if(BPF_slopes[i]==maximum_slope)
	{
		cat("Recommended BPF threshold is",BPF_slopes_identifier[i],"based on a slope of",BPF_slopes[i]/smoothing_factor)
	}
}


#The rest is just data output, first of all plotting all the hybridization patterns and chunk scores to give an idea of the diversity in the sample.
cat("\\documentclass{article}\n\\usepackage{graphicx}\n\\usepackage{epstopdf}\n\\usepackage{color}\n\\usepackage{fullpage}\n\\begin{document}\n\\ttfamily", file = outputfile)
cat("Gene\tAnnotation\tMean Probe Intensity", array1name, "\tMedian Probe Intensity", array1name, "\tBright Segment Length Dependent Score", array1name, "\tBright Probe Fraction", array1name, "\tMean Bright Probe Intensity", array1name, "\tMedian Bright Probe Intensity", array1name, "\tMean Probe Intensity", array2name, "\tMedian Probe Intensity", array2name, "\tBright Segment Length Dependent Score", array2name, "\tBright Probe Fraction", array2name, "\tMean Bright Probe Intensity", array2name, "\tMedian Bright Probe Intensity", array2name, "\t", array2name, "/", array1name, "Median\t", array2name, "/", array1name, " Median Absolute Deviation\n", file = outputtable)

postscript(file = paste(graphdirectory,"bright_probe_fraction_plot.eps",sep="/"), width=9, height=5)
bpf_matrix = matrix(nrow =length(bright_probe_fraction1), ncol=2)
bpf_matrix[,1] = bright_probe_fraction1[order(bright_probe_fraction1)]
bpf_matrix[,2] = bright_probe_fraction2[order(bright_probe_fraction2)]
matplot(bpf_matrix, type="l", ylab="Bright probe fraction", xlab="Genes on the chip", main="Bright Probe Fraction Distribution")
dev.off()
cat("\n\\includegraphics[angle=-90,width=15cm]{bright_probe_fraction_plot.eps}\\\\", file = outputfile, append=TRUE)


postscript(file = paste(graphdirectory,"chunk_score_plot.eps",sep="/"), width=9, height=5)
cs_matrix = matrix(nrow =length(chunk_score1), ncol=2)
cs_matrix[,1] = chunk_score1[order(chunk_score1)]
cs_matrix[,2] = chunk_score2[order(chunk_score2)]
matplot(cs_matrix, type="l", ylab="Bright probe fraction", xlab="Genes on the chip", main="Bright Segment Length Dependent Score Distribution")
dev.off()
cat("\n\n\\includegraphics[angle=-90,width=15cm]{chunk_score_plot.eps}\\\\\n\\clearpage", file = outputfile, append=TRUE)
cat("\\begin{center}\n", file = outputfile, append=TRUE)



#Next up, the hybridization pattern for each gene is plotted individually...
for(i in 1:dim(genes)[1])
	{	
	straightline=vector(length=probenums_chunksort1[i])
	for(j in 1:probenums_chunksort1[i])
		{
			straightline[j]=cutoff
		}

temp_matrix = matrix(nrow =probenums_chunksort1[i], ncol=3)
temp_matrix[,1] = head(array1[probe_matrix_chunksort1[i,],2], n=probenums_chunksort1[i])
temp_matrix[,2] = head(array2[probe_matrix_chunksort1[i,],2], n=probenums_chunksort1[i])
temp_matrix[,3] = straightline

#The next two lines calculate the difference in probe intensity between each of the two arrays, and then finds the sum of the absolute value for these numbers. Basically, the larger the difference between the two samples, the higher the difference_vector will be.
#difference_vector = abs(temp_matrix[,1] - temp_matrix[,2])
#difference_value = sum(difference_vector)

array_comparison_vector = temp_matrix[,2]/temp_matrix[,1]
median_array_comparison = median(log(array_comparison_vector))
mad_array_comparison = mad(log(array_comparison_vector))

below_zero = 0
above_zero = 0

for(t in 1:length(array_comparison_vector))
{
	if(log(array_comparison_vector)[t]>0)
	{
		above_zero = above_zero + 1
	}
	else
	{
		below_zero = below_zero + 1
	}
}

p_value = binom.test(above_zero,above_zero+below_zero,(0.5))$p.value
estimated_probability = binom.test(above_zero,above_zero+below_zero,(0.5))$estimate
upper_CI = binom.test(above_zero,above_zero+below_zero,(0.5))$conf.int[2]
lower_CI = binom.test(above_zero,above_zero+below_zero,(0.5))$conf.int[1]

postscript(file = paste(graphdirectory,paste(i,".eps",sep=""),sep="/"), width=9, height=5)

par(col="black")
matplot(log(temp_matrix), type="l", lwd="3", ylim = c(2,12),  xlab = "Distance along gene (probes)", ylab = "log hybridization intensity", lty=c(1,1,2), col=c("black","red","green"))
#title(main=annotations[grep(genes[order(chunk_score1, decreasing=TRUE)[i],],annotations)], font.main=1, ps=5)
#par(font.lab=2, ps=14, oma=c(1,1,1,1), col="black")
#text(2,1,labels="Bright Probe Fraction:", pos=4)
#text(probenums_chunksort1[i]/1.5,0.5,labels=round(100*bright_probe_fraction1[order(chunk_score1, decreasing=TRUE)[i]],digits=2), pos=4)
#text(probenums_chunksort1[i]/1.3,0.5,labels="%", pos=4)
#par(col="red")
#text(probenums_chunksort1[i]/1.2,0.5,labels=round(100*bright_probe_fraction2[order(chunk_score1, decreasing=TRUE)[i]],digits=2), pos=4)
#text(probenums_chunksort1[i]/1.1,0.5,labels="%", pos=4)
#par(col="black")
#text(2,2,labels="Bright Segment Length Dependent Score:", pos=4)
#text(probenums_chunksort1[i]/1.5,1.5,labels=chunk_score1[order(chunk_score1, decreasing=TRUE)[i]], pos=4)
#par(col="red")
#text(probenums_chunksort1[i]/1.2,1.5,labels=chunk_score2[order(chunk_score1, decreasing=TRUE)[i]], pos=4)
#par(col="black")
#text(2,3,labels="Mean Probe Intensity:", pos=4)
#text(probenums_chunksort1[i]/1.5,2.5,labels=round(mean(temp_matrix[,1]), digits=2), pos=4)
mean_probe_intensity1[i] = mean(temp_matrix[,1])
#par(col="red")
#text(probenums_chunksort1[i]/1.2,2.5,labels=round(mean(temp_matrix[,2]), digits=2), pos=4)
mean_probe_intensity2[i] = mean(temp_matrix[,2])
#par(col="black")
#text(probenums_chunksort1[i]/1.5,3.5,labels="Array 1", pos=4)
#par(col="red")
#text(probenums_chunksort1[i]/1.2,3.5,labels="Array 2", pos=4)

par(col="black")
par(ps=10)
chunknum = 1
for(j in 1:probenums_chunksort1[i])
	{
	if(is.na(start_chunk_blast1[order(chunk_score1, decreasing=TRUE)[i], j]))
		{next}else
		{
		text((start_chunk_blast1[order(chunk_score1, decreasing=TRUE)[i], j]+30)/30,12, labels=chunknum)
		chunknum = chunknum+1
		}
	}
	
	chunknum = 1
	for(j in 1:probenums_chunksort2[i])
		{
		if(is.na(start_chunk_blast2[order(chunk_score1, decreasing=TRUE)[i], j]))
			{next}else
			{
			text((start_chunk_blast2[order(chunk_score1, decreasing=TRUE)[i], j]+30)/30,11, labels=chunknum, col="red")
			chunknum = chunknum+1
			}
		}
dev.off()
cat(annotations[grep(genes[order(chunk_score1, decreasing=TRUE)[i],],annotations)], file = outputfile, append=TRUE)
cat("\\\\\n", file = outputfile, append=TRUE)
cat(paste("\n\\includegraphics[angle=-90,width=15cm]{",i,".eps}\\\\",sep=""), file = outputfile, append=TRUE)

cat("\\begin{tabular}{| l | l | l |}\n", file = outputfile, append=TRUE)
cat("\\hline\n", file = outputfile, append=TRUE)
cat(paste(genes[order(chunk_score1, decreasing=TRUE)[i],]," & ", array1name, " & \\textcolor{red}{", array2name, "} \\\\ \\hline\n"), file = outputfile, append=TRUE)
cat(paste("Mean Probe Intensity: &",round(mean(temp_matrix[,1]), digits=2)," & \\textcolor{red}{", round(mean(temp_matrix[,2]), digits=2), "}\\\\\n"), file = outputfile, append=TRUE) 
cat(paste("Median Probe Intensity: &",round(median(temp_matrix[,1]), digits=2)," & \\textcolor{red}{", round(median(temp_matrix[,2]), digits=2), "}\\\\\n"), file = outputfile, append=TRUE) 
cat(paste("Bright Segment Length Dependent Score: &",chunk_score1[order(chunk_score1, decreasing=TRUE)[i]]," & \\textcolor{red}{", chunk_score2[order(chunk_score1, decreasing=TRUE)[i]], "}\\\\\n"), file = outputfile, append=TRUE) 
cat(paste("Bright Probe Fraction: &",round(100*bright_probe_fraction1[order(chunk_score1, decreasing=TRUE)[i]],digits=2),"\\% & \\textcolor{red}{", round(100*bright_probe_fraction2[order(chunk_score1, decreasing=TRUE)[i]],digits=2), "\\%}\\\\\n"), file = outputfile, append=TRUE) 
cat(paste("Mean of Bright Probes: &",round(bright_gene_means1[order(chunk_score1, decreasing=TRUE)[i]], digits=2)," & \\textcolor{red}{", round(bright_gene_means2[order(chunk_score1, decreasing=TRUE)[i]], digits=2), "}\\\\\n"), file = outputfile, append=TRUE)
cat(paste("Median of Bright Probes: &",round(bright_gene_medians1[order(chunk_score1, decreasing=TRUE)[i]], digits=2)," & \\textcolor{red}{", round(bright_gene_medians2[order(chunk_score1, decreasing=TRUE)[i]], digits=2), "}\\\\\n\\hline\n"), file = outputfile, append=TRUE)
cat(paste("Median of ", array2name, " / ", array1name, ": &\\multicolumn{2}{c|}{",round(median_array_comparison, digits=2), "}\\\\\n"), file = outputfile, append=TRUE)
cat(paste("MAD of ", array2name, " / ", array1name, ": &\\multicolumn{2}{c|}{",round(mad_array_comparison, digits=2), "}\\\\\n"), file = outputfile, append=TRUE)
cat(paste("P-value Binomial Test of ", array2name, "=", array1name, ": &\\multicolumn{2}{c|}{",p_value, "}\\\\\n"), file = outputfile, append=TRUE)
cat(paste("Estimated Probability of ", array2name, ">", array1name, ": &\\multicolumn{2}{c|}{",round(lower_CI, digits=3),"<",round(estimated_probability, digits=3),"<",round(upper_CI, digits=3), "}\\\\\n"), file = outputfile, append=TRUE)

cat("\\hline\n\\end{tabular}\n", file = outputfile, append=TRUE)
#cat("\\end{center}\n\\small\n\\sffamily\nOther genes in this cluster:\\\\\n", file = outputfile, append=TRUE)

cat(paste(genes[order(chunk_score1, decreasing=TRUE)[i],], annotations[grep(genes[order(chunk_score1, decreasing=TRUE)[i],],annotations)], mean(temp_matrix[,1]), median(temp_matrix[,1]), chunk_score1[order(chunk_score1, decreasing=TRUE)[i]], bright_probe_fraction1[order(chunk_score1, decreasing=TRUE)[i]], bright_gene_means1[order(chunk_score1, decreasing=TRUE)[i]], bright_gene_medians1[order(chunk_score1, decreasing=TRUE)[i]], mean(temp_matrix[,2]), median(temp_matrix[,2]), chunk_score2[order(chunk_score1, decreasing=TRUE)[i]], bright_probe_fraction2[order(chunk_score1, decreasing=TRUE)[i]], bright_gene_means2[order(chunk_score1, decreasing=TRUE)[i]],bright_gene_medians2[order(chunk_score1, decreasing=TRUE)[i]],median_array_comparison,mad_array_comparison, "\n", sep="\t"), file = outputtable, append=TRUE)

#cluster = clusters[[grep(genes[order(chunk_score1, decreasing=TRUE)[i],], clusters)]]
#for(q in 1:length(cluster))
#	{
#	if(cluster[q]!=genes[order(chunk_score1, decreasing=TRUE)[i],])
#		{
#		cat(annotations[grep(cluster[q], annotations)], file = outputfile, append=TRUE)
#		cat("\\\\\n", file = outputfile, append=TRUE)
#		}
#	}

#cat("\\ttfamily\n", file = outputfile, append=TRUE)

#followed by the gene annotation.
#cat(as.character(annotations[grep(genes[order(chunk_score1, decreasing=TRUE)[i],],annotations)][1]), file = outputfile, append=TRUE)
#cat("<br><!--BPF=", file = outputfile, append=TRUE)
#cat(as.character(round(100*bright_probe_fraction1[order(chunk_score1, decreasing=TRUE)[i]],digits=2)), file = outputfile, append=TRUE)
#cat(":BSLDP=", file = outputfile, append=TRUE)
#cat(as.character(chunk_score1[order(chunk_score1, decreasing=TRUE)[i]]), file = outputfile, append=TRUE)
#cat(":MPI=", file = outputfile, append=TRUE)
#cat(as.character(round(mean(temp_matrix[,1]))), file = outputfile, append=TRUE)
#cat(":DV=", file = outputfile, append=TRUE)
#cat(as.character(difference_value), file = outputfile, append=TRUE)
#cat(":>", file = outputfile, append=TRUE)



#	for(j in 1:sum(chunk_count_matrix1[order(chunk_score1, decreasing=TRUE)[i],]))
#		{
#followed by the chunk number		
#		cat(j, file = outputfile, append=TRUE)	
#		cat("\\\\\n", file = outputfile, append=TRUE)
#followed by the sequence of the bright segments from the sequence_chunks matrix.		
#		sequence = sequence_chunks1[order(chunk_score1, decreasing=TRUE)[i], j+1]
#		for(k in 1:floor(nchar(sequence)/80))
#			{
#				sequence = paste(substring(sequence,1,k*80),"\\\\",substring(sequence,1+k*80,nchar(sequence)),sep="")
#			}
#		cat(sequence, file = outputfile, append=TRUE)
#		cat("\\\\\n", file = outputfile, append=TRUE)	
#		}
#cat("\\color{red}\n", file = outputfile, append=TRUE)	
#		for(j in 1:sum(chunk_count_matrix2[order(chunk_score1, decreasing=TRUE)[i],]))
#			{
#followed by the chunk number
#			cat(j, file = outputfile, append=TRUE)	
#			cat("\\\\\n", file = outputfile, append=TRUE)	
	#followed by the sequence of the bright segments from the sequence_chunks matrix.		
#	sequence = sequence_chunks2[order(chunk_score1, decreasing=TRUE)[i], j+1]
#	for(k in 1:floor(nchar(sequence)/80))
#		{
#			sequence = paste(substring(sequence,1,k*80),"\\\\",substring(sequence,1+k*80,nchar(sequence)),sep="")
#		}
#	cat(sequence, file = outputfile, append=TRUE)
#	cat("\\\\\n", file = outputfile, append=TRUE)	
#			}
#cat("\\color{black}", file = outputfile, append=TRUE)	
	
		cat("\n\\clearpage\n", file = outputfile, append=TRUE)	
	}

cat("\n\\end{center}\n", file = outputfile, append=TRUE)

postscript(file = paste(graphdirectory,"probe_intensity_comparison_plot.eps",sep="/"), width=9, height=9)
compareplot(arrayset1=mean_probe_intensity1, arrayset2=mean_probe_intensity2, array1label="Array 1 Log Mean Probe Intensity", array2label="Array 2 Log Mean Probe Intensity", title="Log Mean Probe Intensity")
dev.off()
cat(paste("\n\n\\includegraphics[angle=-90,width=15cm]{probe_intensity_comparison_plot.eps}"), file = outputfile, append=TRUE)

postscript(file = paste(graphdirectory,"BSLDS_comparison_plot.eps",sep="/"), width=9, height=9)
compareplot(arrayset1=chunk_score1, arrayset2=chunk_score2, array1label="Array 1 Log BSLDS", array2label="Array 2 Log BSLDS", title="Log Bright Segment Length Dependent Score")
dev.off()
cat(paste("\n\n\\includegraphics[angle=-90,width=15cm]{BSLDS_comparison_plot.eps}"), file = outputfile, append=TRUE)

postscript(file = paste(graphdirectory,"BPF_comparison_plot.eps",sep="/"), width=9, height=9)
compareplot(arrayset1=bright_probe_fraction1, arrayset2=bright_probe_fraction2, array1label="Array 1 Log Bright Probe Fraction", array2label="Array 2 Log Bright Probe Fraction", title="Log Bright Probe Fraction")
dev.off()
cat(paste("\n\n\\includegraphics[angle=-90,width=15cm]{BPF_comparison_plot.eps}"), file = outputfile, append=TRUE)

postscript(file = paste(graphdirectory,"individual_probe_intensity_comparison_plot.eps",sep="/"), width=9, height=9)
compareplot(arrayset1=array1[,2], arrayset2=array2[,2], array1label="Array 1 Log Probe Intensity", array2label="Array 2 Log Probe Intensity", title="Log Individual Probe Intensity")
dev.off()
cat(paste("\n\n\\includegraphics[angle=-90,width=15cm]{individual_probe_intensity_comparison_plot.eps}"), file = outputfile, append=TRUE)


cat("\n\\end{document}", file = outputfile, append=TRUE)
}

