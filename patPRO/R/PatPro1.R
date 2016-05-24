#PatPro1.R
#Geoffrey D Hannigan & Brendan P Hodkinson
#University of Pennsylvania
#Grice Lab

#******************************************************************************
# Copyright 2015 Elizabeth Grice

# This file is part of the patPRO R package.

# The patPRO R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# The patPRO R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with the patPRO R package.  If not, see <http://www.gnu.org/licenses/>.
#******************************************************************************

#===
#Merge the necessary mapping file columns (subject ID and time point which could be days, etc) into the datasets for alpha diversity, bacterial load, and taxa relative abundance.
mergeMapMetaData <- function(map.file, merging.file, map.sub.id="SubjectID", map.tmpt="Time_point", map.smpl.id="SampleID", merge.smpl.id, sample.id.col="SampleID", subject.id.col="SubjectID", tmpt.id.col="Time_point") {
  trimmed_map <- map.file[, c(map.smpl.id, map.sub.id, map.tmpt)]
  colnames(trimmed_map) <- c(sample.id.col, subject.id.col, tmpt.id.col)
  merged_with_map <- merge(merging.file, trimmed_map, by.x = map.smpl.id, by.y = sample.id.col)
  ordered_result <- merged_with_map[ order(merged_with_map[,subject.id.col], merged_with_map[,tmpt.id.col]), ]
  return(ordered_result)
}

#===
#Transpose a dataframe in a way that converts the column headers to row names.
#Use this when dealing with relative abundance calculations using files straight from the Qiime taxonomy output.
transposeRelAbund <- function(x, sample.id.col="SampleID"){
  func1 <- x[,-1]
  row.names(func1) <- x[,1]
  result <- as.data.frame(t(func1))
  result[,sample.id.col] <- c(row.names(result))
  result[,sample.id.col] <- c(gsub("X","",result[,sample.id.col]))
  return(result)
}

#===
#Prepare for graph of the relative abundances of the defined top taxa using the transposed output (originally from Qiime) using 'transpose.rel.abund'
topRelAbundDataFrame <- function(x, top.taxa.num, subject.id.col="SubjectID", sample.id.col="SampleID", tmpt.id.col="Time_point"){
  #Remove the columns that contain the SubjectID and SampleID column headings
  no_id <- x[, !names(x) %in% c(subject.id.col, sample.id.col, tmpt.id.col)]
  #Calculate means of each bacterial taxa (each column is a taxa)
  col_means <- colMeans(no_id)
  #Order the means with the most abundant taxa at the top of the list
  means_order <- order(col_means, decreasing=TRUE)
  #Generate a list of the most relatively abundant taxa based on the list of means
  top_taxa_names <- no_id[,c(means_order[1:top.taxa.num])]
  #Generate a dataframe containing only the columns of the top taxa and the corresponding subject IDs, time points, and sample IDs
  top_taxa <- x[, c(names(top_taxa_names), subject.id.col, tmpt.id.col, sample.id.col)]
  #Use melt to format the dataframe in a way ggplot can use
  melted_top <- melt(top_taxa)
  #Genearte a data frame containing repeats of the subject IDs that correspond to the order of the melted_top data frame
  #This will be two because it includes the sample IDs as well as the time points
  num_plus_two <- top.taxa.num+2
  subject_id <- data.frame(rep(melted_top[melted_top$variable==subject.id.col, "value"], num_plus_two))
  #Also generate a column for the time points
  time_point <- data.frame(rep(melted_top[melted_top$variable==tmpt.id.col, "value"], num_plus_two))
  #Use cbind to add the sample ID column to the rest of the data frame (must use cbind or else the formatting will be screwed up)
  melted_top_ids <- cbind(melted_top, subject_id, time_point)
  #Change the names of the melted dataframe columns
  colnames(melted_top_ids) <- c(sample.id.col,"Bacteria", "Abundance", subject.id.col, tmpt.id.col)
  #Remove the bottom rows of the data frame which contained the subject IDs and were used for generating the list above, but are no longer useful
  no_sub_id <- melted_top_ids[c(!melted_top_ids$Bacteria==subject.id.col), ]
  result <- no_sub_id[c(!no_sub_id$Bacteria==tmpt.id.col), ]
  return(result)
}

#===
# The results from top.rel.abund.data.frame can also be easily used to calculate the mean relative abundance using the following function
# This should be able to be used in for the other data types as well, including the alpha diversity and bacterial load values.
patientMean <- function(x, sub.range, subject.id.col="SubjectID", tmpt.id.col="Time_point", metric.col="Bacteria", abundance.val="Abundance"){
  #First we need to subsample the dataset by removing the subjects that we do not want to include
  taxaSubsampled <- x[c(x[,subject.id.col] %in% sub.range),]
  myDdply <- function(data, group, col){
    result <- ddply(data, group, .fun = function(i){
      Abundance <- c(mean = mean(i[,col]))
      N <- c(n = length(i[,col]))
      StdDev <- c(StdDev = sd(i[,col]))
      StdErr <- StdDev/sqrt(N)
      cbind(Abundance, N, StdDev, StdErr)
    })
    return(result)
  }
  if(metric.col=="") {
    topTaxaMean <- myDdply(taxaSubsampled, c(tmpt.id.col), abundance.val)
    colnames(topTaxaMean) <- c(tmpt.id.col, abundance.val, "N", "StdDev", "StdErr")
  } else {
    topTaxaMean <- myDdply(taxaSubsampled, c(tmpt.id.col,metric.col), abundance.val)
    colnames(topTaxaMean) <- c(tmpt.id.col, metric.col, abundance.val, "N", "StdDev", "StdErr")
  }
  return(topTaxaMean)
}

#===
#Use the results from the top.rel.abund.data.frame function for graphing the top taxa relative abundance
# Right now the value for average counts being TRUE or FALSE is not working. Stragnely it works when the function is fed directly into the script using it, so the issue is using it within the package
plotTopTaxa <- function(top.taxa.data.frame, pat.id, subject.id.col="SubjectID", tmpt.id.col="Time_point", y.lab="Percent Relative Abundance", x.lab="Time Point", plot.title, mark.events=FALSE, mark.times, mark.text="", color.brewer.set="", color.manual.set="", legend.text.size = 7){
  sub_specific_taxa <- top.taxa.data.frame[c(top.taxa.data.frame[,subject.id.col]==pat.id), ]
  #Rename the time point column of the dataframe so that they are easily called through ggplot
  names(sub_specific_taxa)[names(sub_specific_taxa) == tmpt.id.col] <- 'TP'
  # Get the number of time points to be used, in a sorted vector for specifying x axis
  TimePointCount <- unique(sort(sub_specific_taxa$TP))
  # Also make the abundance values up to 100 so that percent goes up to 100
  sub_specific_taxa$Abundance <- sub_specific_taxa$Abundance * 100
  if(length(TimePointCount)>1 & mark.events==FALSE) {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw() + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100))
  } else if(mark.events==TRUE) {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw() + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100))
    for(counter in mark.times) {
      top_taxa_plot <- top_taxa_plot + geom_segment(aes_string(x=counter, y=5, xend=counter, yend=0), arrow = arrow(length = unit(0.3, "cm"))) + geom_text(x=counter, y=7.5, label=mark.text, size=3)
    }
  } else {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_bar(stat = "identity", aes(width = 0)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100) )
  }
  if(color.brewer.set=="" & color.manual.set=="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set=="" & color.manual.set!="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_manual(values=color.manual.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set!="" & color.manual.set=="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_brewer(palette=color.brewer.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  }
}


#===
# This function is similar to the other plot taxa function but uses informaiton for the mean relative abundance information.
plotTopTaxaMean <- function(top.taxa.data.frame, tmpt.id.col="Time_point", y.lab="Percent Relative Abundance", x.lab="Time Point", plot.title, mark.events=FALSE, mark.times, mark.text="", color.brewer.set="", color.manual.set="", legend.text.size = 7){
  sub_specific_taxa <- top.taxa.data.frame
  #Rename the time point column of the dataframe so that they are easily called through ggplot
  names(sub_specific_taxa)[names(sub_specific_taxa) == tmpt.id.col] <- 'TP'
  # Get the number of time points to be used, in a sorted vector for specifying x axis
  TimePointCount <- unique(sort(sub_specific_taxa$TP))
  # Also make the abundance values up to 100 so that percent goes up to 100
  sub_specific_taxa$Abundance <- sub_specific_taxa$Abundance * 100
  if(length(TimePointCount)>1 & mark.events==FALSE) {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100))
  } else if(mark.events==TRUE) {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100))
    for(counter in mark.times) {
      top_taxa_plot <- top_taxa_plot + geom_segment(aes_string(x=counter, y=5, xend=counter, yend=0), arrow = arrow(length = unit(0.3, "cm"))) + geom_text(x=counter, y=7.5, label=mark.text, size=3)
    }
  } else {
    top_taxa_plot <- ggplot(sub_specific_taxa, aes_string(x='TP',y='Abundance',group='Bacteria',fill='Bacteria')) + geom_bar(stat = "identity", aes(width = 0)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + coord_cartesian(ylim = c(0,100))
  }
  if(color.brewer.set=="" & color.manual.set=="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set=="" & color.manual.set!="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_manual(values=color.manual.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set!="" & color.manual.set=="") {
    return(top_taxa_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_brewer(palette=color.brewer.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  }
}

#===
#Use this function to normalize each of the alpha diversity metric measurements to the first time point of each patient (each first time point of each patient will be 100%, and following time points will be relative to that first time point).
#**NOTE** Have user specify first time point in case it is not 1. Default will be 1.
normalizeAlphaDiv <- function(alpha.div.input, alpha.div.metric, subject.id.range, subject.id.col="SubjectID", tmpt.id.col="Time_point"){
  #This script needs to be able to deal with multiple alpha diversity matrices at a time.
  #Use a "for loop" type application of lapply to do this for each subject number identified for i
  total_metrics <- lapply(alpha.div.metric, function(k) {
    metric.id <- k
    total_numbers <- lapply(subject.id.range, function(i) {
      #Change the variable name i to a more specific name
      subject.id <- i
      #Specify the time point which will be used to calculate the relative percents for the subject specified by the for loop
      subj_time_one <- alpha.div.input[c(alpha.div.input[,subject.id.col]==subject.id & alpha.div.input[,tmpt.id.col]==1), metric.id]
      #Get a data frame subset of the data for only subject 1
      subj_one <- alpha.div.input[c(alpha.div.input[,subject.id.col]==subject.id), ]
      #Use sapply to normalize all of the alpha diversity calculations of that patient to the first time point
      normalized_div <- data.frame(sapply(subj_one[ ,metric.id], function(x) x/subj_time_one))
      #Generate a variable name that is specific for the subject ID specified in the for loop.
      variable_name <- paste("final_df_", subject.id, sep="")
      #Assign the output to the variable name that is specific for this iteration of the for loop
      assign(variable_name, normalized_div)
    })
    #Use rbind to cat the data frames together that were specified by the for loop.
    number_cat <- do.call(rbind, total_numbers)
    colnames(number_cat) <- c(paste(metric.id, "Normalized", sep="-"))
    cat_metric_name <- paste("metric_", metric.id, sep="")
    assign(cat_metric_name, number_cat)
  })
  final_metric_data <- do.call(cbind, total_metrics)
  trimmed_input <- alpha.div.input[c(alpha.div.input[,subject.id.col] %in% subject.id.range), ]
  final_df <- cbind(trimmed_input, final_metric_data)
  
  subject_metrics <- lapply(alpha.div.metric, function(i) {
    variable_metric <- paste(i, "Normalized", sep="-")
    specific_metric <- data.frame(final_df[ ,variable_metric])
    colnames(specific_metric) <- c(paste(i, "Normalized", sep="-"))
    cat_metric_name <- paste("metric_", i, sep="")
    assign(cat_metric_name, specific_metric)
  })
  all_metrics <- data.frame(do.call(cbind, subject_metrics))
  melted_metrics <- melt(all_metrics)
  melted_metrics$tmpt.id.col <- c(final_df[,tmpt.id.col])
  melted_metrics$subject.id.col <- c(final_df[,subject.id.col])
  colnames(melted_metrics) <- c("variable", "value", tmpt.id.col, subject.id.col)
  return(melted_metrics)
}

#===
plotNormalizedAlphaDiv <- function(input.df, alpha.div.metrics, tmpt.id.col="Time_point", y.lab="Normalized Alpha Diversity Value", x.lab="Time Point", plot.title, color.brewer.set="", color.manual.set="", mean.mark=FALSE, legend.text.size = 7){
  # Get the number of time points to be used, in a sorted vector for specifying x axis
  TimePointCount <- unique(sort(input.df[,c(tmpt.id.col)]))
  alpha_div_plot <- ggplot(input.df, aes_string(x='Time_point', y='value', group='variable', colour='variable')) + theme_bw() + geom_line() + geom_point()+ theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7))
  if(mean.mark==FALSE) {
    if(color.brewer.set=="" & color.manual.set=="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
    } else if(color.brewer.set=="" & color.manual.set!="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_manual(values=color.manual.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
    } else if(color.brewer.set!="" & color.manual.set=="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_colour_brewer(palette=color.brewer.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
    }
  } else if(mean.mark==TRUE) {
    if(color.brewer.set=="" & color.manual.set=="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_x_continuous(breaks=TimePointCount) + geom_errorbar(aes_string(ymin="value-StdErr", ymax="value+StdErr"), width=0.2) + theme(legend.text = element_text(size = legend.text.size)))
    } else if(color.brewer.set=="" & color.manual.set!="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_manual(values=color.manual.set) + scale_x_continuous(breaks=TimePointCount) + geom_errorbar(aes_string(ymin="value-StdErr", ymax="value+StdErr"), width=0.2) + theme(legend.text = element_text(size = legend.text.size)))
    } else if(color.brewer.set!="" & color.manual.set=="") {
      return(alpha_div_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_colour_brewer(palette=color.brewer.set) + scale_x_continuous(breaks=TimePointCount) + geom_errorbar(aes_string(ymin="value-StdErr", ymax="value+StdErr"), width=0.2) + theme(legend.text = element_text(size = legend.text.size)))
    }
  }
}

#===
plotBacterialLoad <- function(input.bac.load, subject.id, subject.id.col="SubjectID", bac.load.col="bac_load", tmpt.id.col="Time_point", y.lab="Bacterial Load", x.lab="Time Point", plot.title, mean.mark=FALSE){
  if(subject.id=="") {
    TimePointCount <- unique(sort(input.bac.load[,c(tmpt.id.col)]))
    bac_load_plot <- ggplot(input.bac.load, aes_string(x=tmpt.id.col, y=bac.load.col)) + theme_bw() + geom_line() + geom_point() + theme(legend.position="none")  + scale_x_continuous(breaks=TimePointCount)
  } else {
  num_bacteria <- input.bac.load[,bac.load.col][input.bac.load[,subject.id.col]==subject.id]
  time_point <- input.bac.load[,tmpt.id.col][input.bac.load[,subject.id.col]==subject.id]
  num_bind <- cbind(num_bacteria)
  melt_bind <- melt(num_bind)
  melt_bind$Time_point <- c(input.bac.load[,tmpt.id.col][input.bac.load[,subject.id.col]==subject.id])
  TimePointCount <- unique(sort(melt_bind[,c("Time_point")]))
  bac_load_plot <- ggplot(melt_bind, aes_string(x='Time_point', y='value')) + theme_bw() + geom_line() + geom_point() + theme(legend.position="none") + scale_x_continuous(breaks=TimePointCount)
  }
  if(mean.mark==FALSE) {
    return(bac_load_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title))
  } else if(mean.mark==TRUE) {
    std.error.var <- "StdErr"
    return(bac_load_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + geom_errorbar(aes_string(ymin=paste(c(bac.load.col,"StdErr"), collapse='-'), ymax=paste(c(bac.load.col,"StdErr"), collapse='+')), width=0.2))
  }
}

#===
#This function will graph the relative abundance (top taxa like other function) under the values for bacterial load quantification
#First we need to prepare the dataframe for graphing
topAbsAbundDataFrame <- function(rel.abund, bac.load, bac.load.id="bac_load", abund.id.col="Abundance", subject.id.col="SubjectID", tmpt.id.col="Time_point", bacteria.id.col="Bacteria") {
  #Take in the output from top.rel.abund.data.frame and use it in this function as "rel.abund" and the bacterial load quantification as bac.load
  bac_rel_merge <- as.data.frame(merge(rel.abund, bac.load))
  #Add column for calculation of relative abundance as percent of bacterial load (for estimating absolute abundance)
  bac_rel_merge$absolute_abund <- bac_rel_merge[,bac.load.id] * bac_rel_merge[,abund.id.col]
  #Order the resulting dataframe so that it can be properly graphed with ggplot
  bac_rel_abund_ordered <- bac_rel_merge[order(bac_rel_merge[,subject.id.col], bac_rel_merge[,tmpt.id.col], bac_rel_merge[,bacteria.id.col]), ]
  return(bac_rel_abund_ordered)
}

topAbsAbundPlot <- function(rel.abund.df, patient.id, subject.id.col="SubjectID", tmpt.id.col="Time_point", abs.abund.id.col="absolute_abund", bac.id.col="Bacteria", bac.load.col="bac_load", y.lab="Relative Abundance of Bacterial Load", x.lab="Time Point", plot.title, mark.events=FALSE, mark.times, mark.text="", color.brewer.set="", color.manual.set="", legend.text.size = 7) {
  subject_specific_df <- rel.abund.df[c(rel.abund.df[,subject.id.col]==patient.id), ]
  names(subject_specific_df)[names(subject_specific_df) == tmpt.id.col] <- 'TP'
  names(subject_specific_df)[names(subject_specific_df) == abs.abund.id.col] <- 'ABS_ABUND'
  names(subject_specific_df)[names(subject_specific_df) == bac.id.col] <- 'BACTERIA'
  names(subject_specific_df)[names(subject_specific_df) == bac.load.col] <- 'BAC_LOAD'
  # Get the number of time points to be used, in a sorted vector for specifying x axis
  TimePointCount <- unique(sort(subject_specific_df[,c("TP")]))
  ArrowY <- 0.05 * max(subject_specific_df$ABS_ABUND)
  TextY <- 0.08 * max(subject_specific_df$ABS_ABUND)
  if(mark.events==FALSE) {
    top_abund_plot <- ggplot(subject_specific_df, aes_string(x='TP',y='ABS_ABUND',group='BACTERIA',fill='BACTERIA')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + geom_line(aes_string(x='TP', y='BAC_LOAD'), linetype="dashed")
  } else if(mark.events==TRUE){
    top_abund_plot <- ggplot(subject_specific_df, aes_string(x='TP',y='ABS_ABUND',group='BACTERIA',fill='BACTERIA')) + geom_area(position = "stack", stat = "identity", aes(width = 1)) + theme_bw()  + theme(legend.position="right", legend.direction="vertical", legend.box="horizontal", legend.text = element_text(size = 7)) + geom_line(aes_string(x='TP', y='BAC_LOAD'), linetype="dashed")
    for(counter in mark.times) {
      top_abund_plot <- top_abund_plot + geom_segment(aes_string(x=counter, y=ArrowY, xend=counter, yend=0), arrow = arrow(length = unit(0.3, "cm"))) + geom_text(x=counter, y=TextY, label=mark.text, size=3)
    }
  }
  if(color.brewer.set=="" & color.manual.set=="") {
    return(top_abund_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set=="" & color.manual.set!="") {
    return(top_abund_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_manual(values=color.manual.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  } else if(color.brewer.set!="" & color.manual.set=="") {
    return(top_abund_plot + ylab(y.lab) + xlab(x.lab) + ggtitle(plot.title) + scale_fill_brewer(palette=color.brewer.set) + scale_x_continuous(breaks=TimePointCount) + theme(legend.text = element_text(size = legend.text.size)))
  }
}

patproPlotTwo <- function(alpha.div.plot, norm.top.taxa.plot, patpro.plot.title, legend.one.h=0.70, legend.two.h=0.30) {
  # Extract the legend from the figure of interest using this subroutine
  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  # Use the subroutine to get the legends from the two figures
  legendTT <- g_legend(norm.top.taxa.plot)
  legendAD <- g_legend(alpha.div.plot)
  # Create the viewports, push them, draw, and go up
  grid.newpage()
  vp1 <- viewport(width = 0.6, height = 0.9, x = 0.3, y = .50)
  vpleg1 <- viewport(width = 0.25, height = 0.4, x = 0.7, y = legend.one.h)
  vpleg3 <- viewport(width = 0.25, height = 0.4, x = 0.8, y = legend.two.h)
  
  p1 <- norm.top.taxa.plot + theme(legend.position = "none", plot.title=element_blank(), plot.margin = unit(c(-1,0.5,0.5,0.5), "lines"))
  p2 <- alpha.div.plot + theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), plot.title=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(0.5,0.5,-0.3,0.5), "lines"))
  gp1<- ggplot_gtable(ggplot_build(p1))
  gp2<- ggplot_gtable(ggplot_build(p2))
  maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
  gp1$widths[2:3] <- maxWidth
  gp2$widths[2:3] <- maxWidth
  finalGPlot <- arrangeGrob(gp2, gp1)
  
  pushViewport(vp1)
  grid.draw(finalGPlot)
  upViewport(0)
  pushViewport(vpleg1)
  grid.draw(legendAD)
  upViewport(0)
  pushViewport(vpleg3)
  grid.draw(legendTT)
  upViewport(0)
  grid.text(patpro.plot.title, vp = viewport(x = 0.3, y = 0.95))
}

patproPlotThree <- function(alpha.div.plot, bac.load.plot, top.taxa.plot, patpro.plot.title, legend.one.h=0.77, legend.two.h=0.25) {
  # Extract the legend from the figure of interest using this subroutine
  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  # Use the subroutine to get the legends from the two figures
  legendTT <- g_legend(top.taxa.plot)
  legendAD <- g_legend(alpha.div.plot)
  # Create the viewports, push them, draw, and go up
  grid.newpage()
  vp1 <- viewport(width = 0.6, height = 0.9, x = 0.3, y = .50)
  vpleg1 <- viewport(width = 0.25, height = 0.5, x = 0.7, y = legend.one.h)
  vpleg3 <- viewport(width = 0.25, height = 0.5, x = 0.8, y = legend.two.h)
  
  p1 <- top.taxa.plot + theme(legend.position = "none", plot.title=element_blank(), plot.margin = unit(c(-1,0.5,0.5,0.5), "lines"))
  p2 <- bac.load.plot + theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), plot.title=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(-0.45,0.5,-0.3,0.5), "lines"))
  p3 <- alpha.div.plot + theme(legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_blank(), plot.title=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(0.5,0.5,-0.9,0.5), "lines"))
  gp1<- ggplot_gtable(ggplot_build(p1))
  gp2<- ggplot_gtable(ggplot_build(p2))
  gp3<- ggplot_gtable(ggplot_build(p3))
  maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3], gp3$widths[2:3])
  gp1$widths[2:3] <- maxWidth
  gp2$widths[2:3] <- maxWidth
  gp3$widths[2:3] <- maxWidth
  #grid.arrange(gp2, gp1)
  finalGPlot <- arrangeGrob(gp3, gp2, gp1)
  
  pushViewport(vp1)
  grid.draw(finalGPlot)
  upViewport(0)
  pushViewport(vpleg1)
  grid.draw(legendAD)
  upViewport(0)
  pushViewport(vpleg3)
  grid.draw(legendTT)
  upViewport(0)
  grid.text(patpro.plot.title, vp = viewport(x = 0.3, y = 0.95))
}
