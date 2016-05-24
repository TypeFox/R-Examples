##########
# This script is developped for the Pasha package. It contains a function which aims to plot different information about the alignment performed.
# According to gff files and a "AlignedData" object given by the user, the "plottingPasha" function performs:
# 1) A barplot of chromosomes composition: For each chromosomes it gives the number of reads aligned, the length of the chromosome and the number of annotations by gff file.
# 2) Pie charts of:
#                   - The percentage of coverage of each type of annotation contained in gff files.
#                   - The percentage of the number of reads by gff files.
# It gives also some reports below the pie charts.
# 3) A bar plot representation of the information of the pie charts.
# Descostes 03/10/2013
##########




##################
# FUNCTIONS
##################




.restrict_label <- function(index_vec, chrom_vec)
{
    if(length(index_vec) > 0)  # if some chromosomes names were not fitting the length criterion
    {
        restricted_list <- lapply(strsplit(chrom_vec,""), function(x){
                    
                    if(length(x) > 5) {    return(c(x[1:5], "..")) } else return(x) 
                })
        
        return(unlist(lapply(restricted_list, function(x){return(paste(x,collapse=""))})))
    }
    else
    {
        return(chrom_vec)
    }
}



.format_to_scientific <- function(index_vec, number_vec){
    
    if(length(index_vec) > 0)
    {
        number_vec[index_vec] <- format(number_vec[index_vec], scientific=TRUE, digits=1)
    }
    
    return(number_vec)
}




plotReadsRepartitionAnnotations <- function(alignedData,  gff_names_vec, expName, pdfFileName, genomeReferenceFile)
{
    
    #Reading gff files
    list_anno <- mclapply(gff_names_vec, read.table)
    
    
    ##### PART 1: Barplot of genomic information
    
    cat("Plotting information about:\n")
    
    #1) Retrieve size of chromosomes
    
    cat("\t size of chromosomes...\n")
    
    size_chromosome <- readLines(con=genomeReferenceFile) 
    chrom_vec_alignedData <- as.character(unique(seqnames(alignedData)))
    length_chrom <- as.numeric(unlist(strsplit(size_chromosome[2], " ")))
    sizeofgenome <- sum(length_chrom)
    
    names(length_chrom) <- unlist(strsplit(size_chromosome[3], " "))
    length_chrom <- length_chrom[match(chrom_vec_alignedData, names(length_chrom))]
    
    
    #2) Number of reads aligned on each chromosome
    
    cat("\t Number of reads aligned by chromosome...\n")
    
    nb_reads_byChrom <- sapply(chrom_vec_alignedData, function(x,alignedData){
                
                return(length(which(seqnames(alignedData) == x)))
                
            }, alignedData)
    
    #3) Number of annotations by gff file containing a particular chromosome
    
    cat("\t Number of annotations per file and chromosome...\n")
    
    nb_anno_byGFF_byChrom <- sapply(chrom_vec_alignedData, function(chrom, list_anno){
                
                number_anno <- mclapply(list_anno, function(current_anno_type, chrom){
                            
                            return(length(which(current_anno_type == chrom)))
                            
                        }, chrom)
                
                return(unlist(number_anno))
                
            },list_anno)
    
    
    gff_names_vec <- unlist(strsplit(unlist(lapply(strsplit(gff_names_vec, "/"), function(x){return(x[length(x)])})),".gff"))
    row.names(nb_anno_byGFF_byChrom) <- gff_names_vec
    
    #4) Number of annotations by gff files
    
    cat("\t Number of annotations by files...\n\n")
    
    nb_anno_by_gff <- unlist(lapply(list_anno, nrow))
    names(nb_anno_by_gff) <- gff_names_vec
    
    
    # 5) Checking for size of labels, trimming if necessary and output warning message that some names were trimmed and that user should shorten some element names
    
    nb_anno_byGFF_byChrom_full <- nb_anno_byGFF_byChrom
    
    # checking the names of chromosomes that should not be longer than 5 characters
    index_nbAnnoByGFFByChrom_col <- which(nchar(colnames(nb_anno_byGFF_byChrom)) > 5)
    index_nbReadsByChrom <- which(nchar(names(nb_reads_byChrom)) > 5)
    index_lengthChrom <- which(nchar(names(length_chrom)) > 5)
    
    # checking the length of the names of each annotation
    index_names_gff <- which(nchar(names(nb_anno_by_gff)) > 5)
    index_nbAnnoByGFFByChrom_row <- which(nchar(rownames(nb_anno_byGFF_byChrom)) > 5)
    
    # Checking the number of digits for the number of annotation to know if it has to be converted to a scientific format
    index_number_anno <- which(nchar(as.numeric(nb_anno_by_gff)) > 7)
    index_nb_byChrom <- which(nchar(as.numeric(nb_reads_byChrom)) > 7)
    index_length_chrom <- which(nchar(as.numeric(length_chrom)) > 7)
    
    # Trimming and formatting to scientific annotations when needed
    colnames(nb_anno_byGFF_byChrom) <- .restrict_label(index_nbAnnoByGFFByChrom_col, colnames(nb_anno_byGFF_byChrom))
    rownames(nb_anno_byGFF_byChrom) <- .restrict_label(index_nbAnnoByGFFByChrom_row, rownames(nb_anno_byGFF_byChrom))
    names(nb_reads_byChrom) <- .restrict_label(index_nbReadsByChrom, names(nb_reads_byChrom))
    names(length_chrom) <- .restrict_label(index_lengthChrom, names(length_chrom))
    names(nb_anno_by_gff) <- .restrict_label(index_names_gff, names(nb_anno_by_gff))
    
    nb_anno_by_gff <- .format_to_scientific(index_number_anno, nb_anno_by_gff)
    nb_reads_byChrom <- .format_to_scientific(index_nb_byChrom, nb_reads_byChrom)
    length_chrom <- .format_to_scientific(index_length_chrom, length_chrom)
    
    
    
    # 6) Plotting of the different information retrieved: They are plotted separately since the differences of values are too big to see all information
    
    pdf(file= pdfFileName, onefile=TRUE, paper="a4", colormodel="cmyk", title="pashaReport", height=20,width=10)
    
    
    layout(matrix(c(1,2,3,4,5,5,6,6),2,4, byrow=TRUE), heights=c(5,1.5), TRUE)
    
    #### Figure 1
    
    # Ordering by chromosome names
    chr_names <- unlist(lapply(strsplit(colnames(nb_anno_byGFF_byChrom),"chr"),"[[",2))
    ordered_names_index <- mixedorder(chr_names)
    colnames(nb_anno_byGFF_byChrom) <- chr_names
    
    # plotting
    color_vec_plot1 <- rainbow(nrow(nb_anno_byGFF_byChrom))
    barplot(nb_anno_byGFF_byChrom[,ordered_names_index], beside=TRUE, horiz=TRUE, col= color_vec_plot1, main="Number of annotations\nby chromosomes", xlab="number of annotations", ylab="chromosome", cex.main=0.99, las=1)
    
    #### Figure 2 on right side of 1
    
    #remove names from vector
    nb_anno_by_gff_nonames <- nb_anno_by_gff
    names(nb_anno_by_gff_nonames) <- NULL
    
    ycoord <- barplot(nb_anno_by_gff_nonames, horiz=TRUE, col=color_vec_plot1, main="Number of annotations\nby file", xlab="number of annotations", ylab="file type", cex.main=0.99)
    text(rep(max(nb_anno_by_gff)/2, length(nb_anno_by_gff)), ycoord,labels= paste(names(nb_anno_by_gff), as.numeric(nb_anno_by_gff),sep=": "))
    
    
    #### Figure 3 on right side of 2
    
    # Ordering by chromosome names
    chr_names <- unlist(strsplit(names(nb_reads_byChrom),"chr"))
    chr_names <- chr_names[which(chr_names != "")]
    ordered_names_index <- mixedorder(chr_names)
    names(nb_reads_byChrom) <- chr_names
    
    color_vec_plot3 <- rainbow(length(nb_reads_byChrom))
    ycoord2 <- barplot(nb_reads_byChrom[ordered_names_index], horiz=TRUE, col= color_vec_plot3, main="Number of reads\nby chromosome", xlab="number of reads", ylab="chromosome", cex.main=0.99, las=1)
    text(rep(max(nb_reads_byChrom)/2, length(nb_reads_byChrom)), ycoord2, labels = nb_reads_byChrom)
    
    #### Figure 4 on right side of 3
    
    # Ordering by chromosome names
    chr_names <- unlist(strsplit(names(length_chrom),"chr"))
    chr_names <- chr_names[which(chr_names != "")]
    ordered_names_index <- mixedorder(chr_names)
    names(length_chrom) <- chr_names
    
    color_vec_plot4 <- rainbow(length(length_chrom))
    ycoord3 <- barplot(as.numeric(length_chrom), horiz=TRUE, col= color_vec_plot4, main="Length of\nchromosome", xlab="length", cex.main=0.99, names.arg = names(length_chrom)[ordered_names_index],las=1)
    text(rep(max(as.numeric(length_chrom))/2, length(length_chrom)), ycoord3, labels= length_chrom)
    
    # Figure 5 below 1
    plot(NULL, xlim=c(0,1),ylim=c(0,1), axes=FALSE, xlab="", ylab="")
    legend("left", legend = rev(row.names(nb_anno_byGFF_byChrom_full)), fill = rev(color_vec_plot1), bty="n")
    
    
    
    ##### PART 2: pie chart of coverage and percentage of aligned reads by gff files
    
    color_vec <- rainbow(length(list_anno)+1) # +1 stands for the "others" categorie
    
    # 1) Computation of the overlap of reads by annotations submitted
    
    cat("Computing overlap of annotations with reads:\n")
    
    # Convert each gff file and aligned reads to ranged data (This enables to perform intersections with IRanges)
    list_rangedData_anno <- mclapply(list_anno, function(x){ vector_start <- x[,4]
                vector_end <- x[,5]
                vector_seqname <- x[,1]
                return(RangedData(IRanges(start = vector_start, end = vector_end), space = vector_seqname))})
    
    alignedData_rangedData <- RangedData(IRanges(start = position(alignedData), end = position(alignedData) + qwidth(alignedData)), space = seqnames(alignedData))
    
    
    hits_result_list <- mclapply(list_rangedData_anno, function(x, reads_aligned){ 
                return(findOverlaps(reads_aligned,x))
            }, alignedData_rangedData)
    
    # 2) Compute the coverage of each annotation regarding the genome
    
    cat("\t Computing the coverage of annotations on genome...\n")
    
    coverage_vec <- unlist(mclapply(list_anno, function(x, sizeofgenome){
                        
                        start_vec <- as.numeric(x[,4])
                        end_vec <- as.numeric(x[,5])
                        length_vec <- end_vec - start_vec
                        return((sum(length_vec)*100)/sizeofgenome)
                        
                    }, sizeofgenome))
    
    if(sum(coverage_vec) < 100)
    {
        coverage_vec <- c(coverage_vec, (100 - sum(coverage_vec)))
        names(coverage_vec) <- paste(paste(c(names(nb_anno_by_gff), "others"), round(coverage_vec, digits=1), sep="-"), "%", sep="")
    }
    else
    {
        names(coverage_vec) <- paste(paste(names(nb_anno_by_gff), round(coverage_vec, digits=1), sep="-"), "%", sep="")
    }
    
    
    
    # 3) Retrieve the reads which overlapped each annotations and plot a pie chart of the percentage.
    
    cat("\t Computing the number of reads overlapping each annotation file...\n\n")
    
    total_number_reads <- length(seqnames(alignedData))
    nb_reads_overlappingPercent_vec <- unlist(mclapply(hits_result_list, function(x, total_number_reads){return((length(unique(queryHits(x)))*100)/total_number_reads)}, total_number_reads))
    
    if(sum(nb_reads_overlappingPercent_vec) < 100)
    {
        other_percent <- 100 - sum(nb_reads_overlappingPercent_vec)
        nb_reads_overlappingPercent_vec <- c(nb_reads_overlappingPercent_vec, other_percent)
        names(nb_reads_overlappingPercent_vec) <- paste(paste(c(names(nb_anno_by_gff), "others"), round(nb_reads_overlappingPercent_vec, digits=1), sep="-"), "%", sep="")
    }
    else
    {
        names(nb_reads_overlappingPercent_vec) <- paste(paste(names(nb_anno_by_gff), round(nb_reads_overlappingPercent_vec, digits=1), sep="-"), "%", sep="")
    }
    
    
    # 4) Plotting the pie charts and legends
    
    
    layout(matrix(c(1,2,3,4),2,2, byrow=TRUE))
    
    index_ordered <- order(coverage_vec, decreasing = TRUE)
    percent_vec <- unlist(lapply(strsplit(names(coverage_vec),"-"),"[",2))
    pie(coverage_vec, main="Coverage of the genome by annotations", labels= percent_vec, col = color_vec, border=NA)
    plot(NULL, xlim=c(0,1),ylim=c(0,1), axes=FALSE, xlab="", ylab="")
    legend("top", names(coverage_vec)[index_ordered], cex = 0.9, col = color_vec[index_ordered], bty="n", pch=15)
    
    index_ordered <- order(nb_reads_overlappingPercent_vec, decreasing = TRUE)
    percent_vec <- unlist(lapply(strsplit(names(nb_reads_overlappingPercent_vec),"-"),"[",2))
    pie(nb_reads_overlappingPercent_vec, main=paste("Number of reads of ", expName, "\noverlapping annotations", sep=""), col = color_vec, labels = percent_vec, border=NA)
    plot(NULL, xlim=c(0,1),ylim=c(0,1), axes=FALSE, xlab="", ylab="")
    legend("top", names(nb_reads_overlappingPercent_vec)[index_ordered], cex = 0.9, col = color_vec[index_ordered], bty="n", pch=15)
    
    
    #perform a chi square test
    
    if(length(coverage_vec) != length(nb_reads_overlappingPercent_vec))
    {
        # A section "others" could have been inserted in the vector "nb_reads_overlappingPercent_vec" in this case the last element is not considered
        chisquare_test <- chisq.test(coverage_vec, nb_reads_overlappingPercent_vec[1:(length(nb_reads_overlappingPercent_vec)-1)], simulate.p.value=TRUE)
    }
    else
    {
        chisquare_test <- chisq.test(coverage_vec, nb_reads_overlappingPercent_vec, simulate.p.value=TRUE)
    }
    cat("\nchisquare test p.value = ", chisquare_test$p.value, "\n")
    
    
    
    ##### PART 3: Representing information of pie charts with barplot
    
    cat("Plotting the details of pie charts...\n")
    
    number_row <- round(length(coverage_vec)/2)
    par(mfrow=c(number_row, length(coverage_vec)-number_row))
    
    for(i in 1:length(coverage_vec))
    {
        add_ylim <- 2 * max(coverage_vec[i], nb_reads_overlappingPercent_vec[i])
        namesarg_percent <- c(round(coverage_vec[i], digits=3), round(nb_reads_overlappingPercent_vec[i], digits=3))
        title_main <- strsplit(names(coverage_vec[i]), "-")[[1]][1]
        
        
        barplot(c(coverage_vec[i], nb_reads_overlappingPercent_vec[i]), main = title_main, names.arg = namesarg_percent, col=c("blue","red"), ylim=c(0, max(coverage_vec[i], nb_reads_overlappingPercent_vec[i]) + add_ylim), ylab="percentage")
        legend("topleft", legend=c("genome", expName),bty="n", pch=15, col=c("blue", "red"))
    }
    
    dev.off()
}
