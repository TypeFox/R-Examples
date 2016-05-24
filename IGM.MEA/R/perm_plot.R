###############################################################################
# Purpose:  Functions for computing p-values via mann whit and permut. test   #
#              as well as generating plots                                    #
# Author:   Ryan Dhindsa                                                      #
###############################################################################
.extract.rows<-function(data, type){
  # Extracts rows for particular treatment 
  #
  # Args:
  #   data = a dataframe
  #   type = treatment
  #           
  # Returns:
  #   list of data frames containing spike data
  #x = data[data$treatment == type&!is.na(data$treatment), ]
  x = with(data,data[treatment == type&!is.na(treatment), ])
  x$treatment = NULL
  x$well = NULL
  
  return(x)
}

.mann.whit.perm <- function(df, wt, trt, np){
  # Calculates mann-whit p-value and permutation test for a dataframe
  #     must specify wt and the treatment you are testing
  #
  # Args:
  #   df = a dataframe with featre data
  #   wt = the treatment that will be considred wt
  #   np = number of permutations
  #           
  # Returns:
  #   a df that contains a permutation p-value and a mann whit p-value
  wt.df <- .extract.rows(df, wt)
  if (nrow(wt.df) > 0){
    wt.df[wt.df == "NaN"] = NA
  }
  
  trt.df <- .extract.rows(df, trt)
  if (nrow(trt.df) > 0){
    trt.df[trt.df == "NaN"] = NA  #convert NaN to NA
  }
  
  if (is.na(trt) | is.na(wt) |nrow(wt.df)==0 | nrow(trt.df)==0 | all(is.na(trt.df)) | all(is.na(wt.df))){
    #create data frame with empty results if one of the treatments is empty
    result = c(NA, NA)
    result.df = data.frame(result)
    row.names(result.df) = c("perm.p", "data.p")
    colnames(result.df) = trt
    result.df = t(result.df)
    return(result.df)
  }

   #pool wt and trt data
  pool <- rbind(wt.df, trt.df)
  pool = as.matrix(pool)
  
  n.wt <- nrow(wt.df)
  n <- nrow(pool)
  
  # randomly sample the data, compute p-values and store in outp
  outp = matrix(0,np,1)
  for (i in 1:np) {
    wt.smpl = sample(n, n.wt)
    perm.wt = as.numeric(as.vector(pool[wt.smpl,]))
    perm.trt = as.numeric(as.vector(pool[-wt.smpl,]))
    outp[i] <- wilcox.test(perm.wt, perm.trt)$p.value
  }
  outp = sort(outp)
  
  # calculate actual p-val from data
  data.wt = as.numeric(unlist(wt.df))
  data.trt = as.numeric(unlist(trt.df))
  data.p <- wilcox.test(data.wt, data.trt)$p.value
  data.p.rounded <- round(data.p)
  if (data.p.rounded == 0& !is.na(data.p.rounded)){
    data.p.rounded <- signif(data.p,3)
  }
  
  if (data.p == "NaN"){
    perm.p <- "NaN"
  } else{
    perm.p <- length(which(outp<data.p)) / np
    if(perm.p == 0){
      perm.p = paste("<", 1/np)
    }
  }

  
  #create data frame with results
  result = c(perm.p, data.p)
  result.df = data.frame(result)
  row.names(result.df) = c("perm.p", "data.p")
  
  #result <- list(perm.p = perm.p, data.p = data.p)
  result.df <-t(result.df)
  return(result.df)
}

get.wt <- function(s){
  # Uses tcltk user input to specify which treatment should be considered wt
  #     later on, we calculate p-value for every other treatment vs. this wt specification
  #
  # Args:
  #   s object
  #           
  # Returns:
  #   wt
  choices = as.vector(unique(s[[1]]$treatment[!is.na(s[[1]]$treatment)&s[[1]]$treatment!=""]))
  wt = tk_select.list(choices, preselect = NULL, multiple = FALSE,
                      title = "Choose wildtype/reference for permutation test")
  return(wt)
}

.plot.feature <- function(df, feature, platename){
  # Uses ggplot to plot feature data
  #
  # Args:
  #   df = a datafraem containing feature data
  #   feature = name of feature
  #   platename
  #           
  # Returns:
  #   plot
  
  
  df$well = NULL
  df[,-1] = sapply(df[,-1], as.numeric)
  
  melted.df = melt(df, id.vars = "treatment")
  #melted.df = ddply(melted.df, .(treatment, variable), summarize, 
  #                  mean=mean(value),
  #                  sd = sd(value))
  melted.df = with(melted.df, {ddply(melted.df, c("treatment", "variable"), summarize, 
                    mean=mean(value),
                    sem = sqrt(var(value)/length(value)))})
  
  
  pd <- position_dodge(width = 0.1)
  title <- paste(platename, "_", feature, sep="")
  
  x = with(melted.df, { ggplot(melted.df, aes(x=variable, y=mean, group = treatment))+
    geom_point(aes(color=factor(treatment)), position = pd)+
    geom_line(aes(color=factor(treatment)), position = pd)+
    geom_errorbar(aes(x=variable, y=sem, ymin = mean-sem, ymax=mean+sem, 
                      color = factor(treatment)), width = 0.1, position = pd)+
    ggtitle(paste0("\n", title,"\n"))+
    xlab("")+
    ylab(paste0("\n", feature, "\n"))})
  
  x+labs(color="Treatment")
  
  return(x)
}

.apply.perm.and.plot <- function(wt, df, np, feature, platename){
  # Calls .mann.whit.perm() and .plot.feature() to create a single gtable containing
  #     the plot and table
  #
  # Args:
  #   wt, df, feature, and platename
  #           
  # Returns:
  #   gtable
  
  all.treatments = as.vector(unique(df$treatment))
  max(all.treatments, na.rm=TRUE)  # remove non-existent treatments
  
  #extract all treatments that aren't wt
  test.treatments = all.treatments[!all.treatments==wt] 
  
  #calculate mann whit p values for each combination
  num.of.trts = length(test.treatments)

  all.p.values <- data.frame(Treatment=character(),
                 perm.pval=character(), 
                 original.pval=character(), 
                 stringsAsFactors=FALSE) 

  if (num.of.trts==0){ all.p.values<-c(-1,-1) }
  else{
    for (i in 1:num.of.trts){
      trt = test.treatments[i]
      vals = .mann.whit.perm(df, wt, trt, np)
      all.p.values[i,"Treatment"]=paste(wt, " vs. ", trt)
      all.p.values[i,"perm.pval"]= vals[,"perm.p"]
      all.p.values[i,"original.pval"]= vals[,"data.p"]
    }}
  
  colnames(all.p.values)[1] <- paste("Treatment/Genotype")
  
  
  p.value.table = tableGrob(all.p.values)
  feature.plot = .plot.feature(df, feature, platename)
  
  table.and.plot <- arrangeGrob(feature.plot, p.value.table, nrow=2)
  return(table.and.plot)
}

write.pdf <- function(s, wt, np, features.list, type, output.dir){
  # Calls .apply.perm.and.plot() and writes PDF--each page contains a plot and 
  #       table of p-values 
  #
  # Args:
  #   s object
  #   wt
  #   features.list = list of dataframes containing feature data
  #   type = spikes, ns, or bursts
  #   
  #           
  # Returns:
  #   Writes a PDF
  
  out.folder <- paste0(output.dir,"/",type)
  dir.create(out.folder, showWarnings = FALSE)
  
  platename <- get.project.plate.name(s[[1]]$file)
  fname <- paste(platename, "_", type, "_", "analysis", ".pdf", sep='')
  fpath <- paste(out.folder, "/",fname, sep='')
  
  if (length(unique(na.omit(features.list[[1]]$treatment)))==0){
    # No treatments
    return
  }
  
  x = list()
  for(i in 1:length(features.list)){
    feature = names(features.list[i])
    perm.and.plot = .apply.perm.and.plot(wt, features.list[[i]], np, feature, platename)
    x[[feature]] = perm.and.plot
  }
  
  
  all.plots = marrangeGrob(x, nrow=1, ncol=1)
  
  # create pdf
  ggsave(fpath, all.plots, width = 8.5, height = 11, units = "in")
}
