IGM.plot.distributions <- function(s, minVals=1, xlimit=25, binsInSec=5, 
                               feature="non", filterValuesByMin=0, minValues=0,
                               perWell=0, outputdir=getwd()){
# Plot distributions of selected bursting features and print csv output for stats
#
# Args:
#   s        :       the recorded object
#   minVals  : minimum values number per electrode, electrodes with a smaller number of values than that are discarded
#   xlimit   :  max limit of values, for example: xlimit = 25 for IBI analysis means that IBIs longer than 25 seconds 
#     will not be part of distribution calculations
#   binsInSec:  how many bins to cut each of the segments. For example: IBI analysis has 25 seconds as xlimit, 
#     to analyse in a 0.1 sec resolution binsInSec should be set to 10, for 1 sec resolution set binsInsec to 1   
#   feature  :  what feature to analyze, options are "IBI", "ISI, "nspikesInBurst", "duration", "spikesDensityInBurst"
#   filterValuesByMin:  should analysis disregard values with lower then filterValuesByMin number of values ? (0/1, default is 0)
#     for example, if set to 1 for duration analysis, should analysis consider also bursts shorter than filterValuesByMin ?
#   minValues:  disregards values with lower then filterValuesByMin , only if filterValuesByMin set to 1 
#   perWell  : should distribution analysis be performed by testing treatment differences on well level means (1) or electrode level means(0) 
#   outputdir: output directory
#
# Output:
#   Plots the burst distributions.
#   Writes a burst distribution csv to be used for permutation test and plotting
  
  basename <- get.file.basename(s$file)
  logFile <- paste(outputdir, "/", get.project.plate.name(s$file), 
                   "_distributions_log.txt", sep="")
  
  stopifnot(feature != "non")
  ma5 <- c(rep(1, 5)) / 5
  duration <- s$rec.time[2] - s$rec.time[1]
  
  write(paste("->->-> Analysing ", s$file, sep=""), file=logFile, append = TRUE)
  
  cat(paste("Arguments: minVals=", minVals, "; xlimit=", xlimit, "; binsInSec=",
              binsInSec, ";perWell=", perWell, "; duration=", duration, "; feature='", feature, "'; filterValuesByMin=", filterValuesByMin, "; minValues=", minValues,"\n",sep=""))#
  
  treatments <- unique(s$treatment)
  for (tr in 1:length(treatments)){
    if (nchar(treatments[tr])==0){
      treatments <- treatments[-tr]
      next}}
  
  fVals <- NULL
  if (feature == "IBI") {
    fVals <- .calc.all.ibi(s, s$allb)
  } else if (feature == "ISI") {
    fVals <- .calc.all.isi(s, s$allb)
  } else if (feature == "nspikesInBurst") {
    for (j in 1:length(s$channels)) {
      ## Number of spikes in burst !
      fVals[j] <- burstinfo(s$allb[j], "len")
    }} else if (feature == "duration") {
      for (j in 1:length(s$channels)) {
        ## Duration of burst
        fVals[j] <- burstinfo(s$allb[j], "durn")
      }} else if (feature == "spikesDensityInBurst") {
        for (j in 1:length(s$channels)) {
          if (length(burstinfo(s$allb[j], "durn")[[1]]) > 1) {
            ## Number of spikes in burst divided by duration of burst
            fVals[j] <- mapply("/", burstinfo(s$allb[j], "len"), burstinfo(s$allb[j], "durn"), SIMPLIFY = FALSE)
          } else {
            fVals[j] <- 0
            if (burstinfo(s$allb[j], "durn")==0) {
              fVals[j] <- 0
            } else {
              fVals[j] <- mapply("/", burstinfo(s$allb[j], "len"), burstinfo(s$allb[j], "durn"), SIMPLIFY = FALSE)
            }
          }
        }
      }
  
  jump <- binsInSec
  
  # Calculate active E per well
  CperWell <- NULL
  wells <- unique(s$cw)
  if (length(wells) > 0) {
    for (well in wells) {
      treat <- s$treatment[well][[1]]      
      if (treat == "") {
        write (paste("Skipping well-", well, " that has active electrodes since treatment property is empty.", sep=""), file=logFile, append = TRUE)
      }
      active.electrodes <- which(s$cw == well & as.vector(unlist(lapply(s$isis, length))) > 0)
      CperWell <- rbind(CperWell, data.frame(well, newCount=length(active.electrodes)))
    }} else {
      write("No wells found!", file=logFile, append = TRUE)
      return
    }
  
  old <- ""
  wellData <- NULL
  treat <- NULL
  gmeanNormHist <- NULL
  gmeanNormHistByWell <- NULL
  cNum <- 0
  
  # Making the full table
  #  post=c((1/jump) * (minValues * jump):(duration * jump))
  post <- c((1 / jump) * (minValues * jump):(xlimit * jump))
  posT <- post[1:length(post) - 1] + (post[1] + post[2]) / 2
  num_rows <- length(posT) * length(treatments)
  
  # Creating the electrode / well matrix
  allPos <- rep(posT, times=length(treatments))
  allTrt <- rep(treatments, each = length(posT))
  all <- data.frame(pos=allPos, treat=allTrt)
  firstDataColumn <- 3
  
  for (j in 1:length(s$channels)){
    
    #indicator of current well
    icurrentwell <- (s$channels[j])
    well <- strsplit(icurrentwell, "_")[[1]][1]
    
    #    print(paste("well-", icurrentwell, " num of ele-", (CpserWell$newCount[CperWell$well==well])))
    # skip channel if well has less the 4 active E or treatment is ""
    
    treat <- s$treatment[well][[1]]
    # Skip electrode if it does not have a valid treatment property from csv file
    if (!(treat %in% treatments))
    {next}
    if (treat == ""){ 
      next}
    
    # Check number of aE
    if (CperWell$newCount[CperWell$well==well] <= 3){
      write (paste("Skipped well- ", well, " less than 4 electrodes", sep=""), file=logFile, append = TRUE)
      next}else{
        data <- fVals[[j]]
      }
    
    #     Filter data only if over a min
    if (filterValuesByMin==1)
    {  data <- data[data > minValues]}
    
    # Remove from analysis all values that are beyond xLimit, to make sure that
    # the output distribution will have a total on 1 (all values within range)
    data <- data[data <= xlimit]
    histData <- data.frame()
    
    # If electrode has enough values then calculate normal histogram
    if (length(data) > minVals){
      e <- hist(data, plot=FALSE, breaks = c((1 / jump) * (minValues * jump):(xlimit * jump)))
      
      # normalize to number of values in electrode
      histData <- data.frame(normHist=(e$counts / length(data)), pos=e$mids, treat=treat)
      names(histData)[1] <- icurrentwell
      
      if (perWell == 0)
      {
        # add norm hist data to the overall table
        if (treat %in% treatments){
          all[, names(histData)[1]] <- NA
          all$pos <- as.character(all$pos)
          all[all$treat==treat&all$pos==as.character(histData$pos), names(histData)[1]]=histData[, 1]
        }
      } # END perWell ==0 
      
      # for per well analysis only
      if (perWell ==1){
        # new well
        if (well != old ){
          if (old!=""){
            # save old well if there are electrodes to show
            if (cNum != 0){
              valid <- 1
              if (cNum > 3)
              {
                # calculate average values per electrodes if more than 1
                rmeans <- apply(wellData[, 1:cNum], 1, mean)}else{
                  write (paste("Less than four electrodes passed filters for ", old, ", it will be removed from the analysis."), file=logFile, append = TRUE)
                  valid <- 0
                }
              # Continue if treatment is valid
              if (unique(wellData$treat) %in% treatments && valid==1){
                # make a new data frame from the well data
                meansDF <- data.frame(rmeans, pos=wellData$pos, treat=unique(wellData$treat))
                names(meansDF)[1] <- old #change to name of the well that was just analyzed
                
                # Fill in the matrix with well data
                all[, old] <- NA
                all$pos <- as.character(all$pos)
                all[as.character(all$treat)==as.character(unique(wellData$treat))&all$pos==as.character(meansDF$pos), old]=meansDF[, 1]
              }} # more than 0 electrodes
          }
          write(paste("Analyzing well ", well, ",  channel number is:", j, sep=""), file=logFile, append = TRUE)
          old <- well
          cNum <- 1
          wellData <- histData
        }else{ # well == old
          if (cNum==0){
            wellData <- histData}else{
              wellData <- cbind(histData[,1], wellData)}
          cNum <- cNum + 1}}else{ write(paste("Analysis per electrode, electrode-", icurrentwell, "Analyzed."), file=logFile, append = TRUE)}
    } #    if (length(data) > minVals){
  } ## End for j in 1:s$channles
  
  
  # for per well analysis only
  if (perWell ==1){
    if (cNum > 1)
    {
      #save last well 
      rmeans <- apply(wellData[, 1:cNum], 1, mean)#median)}
      # Continue if treatment is valid
      if (unique(wellData$treat) %in% treatments){
        
        meansDF <- data.frame(rmeans, pos=wellData$pos, treat=unique(wellData$treat))
        names(meansDF)[1]=old
        # Fill in the matrix with last well data
        all[, old] <- NA
        all$pos <- as.character(all$pos)
        all[as.character(all$treat)==as.character(unique(wellData$treat))&all$pos==as.character(meansDF$pos), old] <- meansDF[,1]
      }}}
  
  par(mfrow=c(1,1))  
  colors <- c("red", "blue", "green", "orange", "pink", "brown", "yellow", "black")
  first <- 1#
  
  #sahar -time
  table <- all#NULL
  
  if (perWell ==1){
    write("Collapsing genotype data by well.", file=logFile, append = TRUE)
  }else{
    write("Collapsing genotype data by electrode.", file=logFile, append = TRUE)
  }# perwell ?
  
  if (is.null(table))
  {
    write (paste("Based on current parameters there are no cases to plot,  skipping file ", s$file, sep=""), file=logFile, append = TRUE)
    return(1)
  }
  if ((dim(table)[2] - (firstDataColumn - 1)) <= 1)
  {
    write (paste("Not enough electrodes (<=1) to plot, skipping file ",s$file,sep=""),file=logFile,append = TRUE)
    return(1)
  }
  write(paste("Number of electrodes in table : ", (dim(table)[2] - (firstDataColumn - 1)), sep=""), file=logFile, append = TRUE)
  
  # set the ylimit on the top 2% values of the left most dist values (1st fifth 0 - 0.2 (xlimit / 5))
  lim <- apply(table[, (firstDataColumn:dim(table)[2])], 1, mean, na.rm=TRUE)
  
  # Set ylimit for graph based on first columns of data (small values)
  ylimit <- max(lim, na.rm=T)
  
  #  ylimit=max(table[, 4:dim(table)[2]], na.rm=T)
  gmedsTable <- NULL
  gmeansTable <- NULL
  goodTreatments <- treatments
  # Check for empty genotypes
  for (t in treatments){
    if (nrow(table[table$treat==t, ])==0){
      goodTreatments <- goodTreatments[goodTreatments != t]
      write(paste("Treatment", t, "has no electrodes to contribute data after filters and is excluded from analysis.", sep=" "), file=logFile, append = TRUE)
      next}
    tto_plot <- table[table$treat==t, ]
    # Check how many electrodes remain after all the analysis for the treatment
    cleanTreatment <- tto_plot[, colSums(is.na(tto_plot)) != nrow(tto_plot)]
    if (ncol(cleanTreatment) < 7) # 3 are treat, pos, id + 4 electrodes
    {
      goodTreatments <- goodTreatments[goodTreatments != t]
      write(paste("Treatment", t, "has less than 4 electrodes to contribute data after filters and is excluded from analysis.", sep=" "), file=logFile, append = TRUE)
      next}
  }
  Alldistributions <- NULL
  if (length(goodTreatments) <= 1)
  {
    write (paste("Not enough treatments (<=1) with enough data to plot, skipping file ",s$file,sep=""),file=logFile,append = TRUE)
    return(1)
  }  
  # plot per genotype data
  for (tr in 1:length(goodTreatments)){
    t  <- goodTreatments[tr]
    #  select a treatment to plot and calculate means for all the 
    # columns of each row (all electrodes of each genotype)
    tto_plot <- table[table$treat==t, ]
    # Prepare a table for permutations analysis
    cleanTreatment <- tto_plot[, colSums(is.na(tto_plot)) != nrow(tto_plot)]
    
    transposed <- as.data.frame(t(cleanTreatment[1:ceiling(xlimit * jump), firstDataColumn:length(cleanTreatment)]))
    namesT <- rownames(transposed)
    for (pos in 1:length(namesT)){
      namesT[pos] <- strsplit(namesT[pos],  "_")[[1]][1]}
    transposed <- cbind(genotype=t, well=namesT, transposed)
    
    names(transposed) <- c("genotype", 1:(length(transposed) - 1))
    Alldistributions <- rbind(Alldistributions, transposed)
    gmeans <- (apply(tto_plot[, (firstDataColumn:dim(tto_plot)[2])], 1, mean, na.rm=TRUE))
    # smooth the values for representatoin only
    data <- gmeans
    if (xlimit > 5)
    {    data <- filter(gmeans, ma5)}
    pos <- table$pos
    if (first){
      gmeansTable <- data.frame(treat=as.character(t), data=gmeans)
      bottom <- factorial(length(goodTreatments)) / (2 * (factorial(length(goodTreatments)-2))) + 3
      par(mar=c(bottom, 3, 3, 2))
      plot(data[1:ceiling(xlimit * jump)]~pos[1:ceiling(xlimit * jump)], type = "l", col = colors[tr], ylim=c(0, ylimit), xlim=c(0, xlimit),
           ylab = paste(feature, " normalized histogram over genotypes", sep=""), main= paste( feature, " by treatment", sep=""),
           xlab="", lwd=3)
      points(data[1:ceiling(xlimit * jump)]~pos[1:ceiling(xlimit * jump)], type = "l", col = colors[tr], lwd=4)
      first <- 0
    }else{
      gmeansTable <- rbind(gmeansTable, data.frame(treat=t, data=gmeans))
      points(data[1:ceiling(xlimit * jump)]~pos[1:ceiling(xlimit * jump)], type = "l", col = colors[tr], lwd=4)}}
  # print legend for all genotypes and treatments
  legend(xlimit * 0.6, ylimit * 0.8, legend=c(goodTreatments), lty=1, lwd=5, col=colors, bg="white", cex=0.9, y.intersp=0.7)  
  # print values of all ks tests
  line <- 3
  if (length(goodTreatments) > 1){
    for (t2 in 1:(length(goodTreatments) - 1)){
      if (is.na(sum(gmeansTable[gmeansTable$treat==goodTreatments[t2], "data"][1:(xlimit * jump)], na.rm=T))){
        write(paste("Treatment ", t2, " is empty for this recording.", sep=""), file=logFile, append = TRUE)
        next
      }
      for (t3 in (t2 + 1):length(goodTreatments)){
        if (is.na(sum(gmeansTable[gmeansTable$treat==goodTreatments[t3], "data"][1:(xlimit * jump)], na.rm=T))){
          write(paste("Treatment ", t3, " is empty for this recording.", sep=""), file=logFile, append = TRUE)
          next
        }
        
        w <- ks.test(gmeansTable[gmeansTable$treat==goodTreatments[t2], "data"][1:(xlimit * jump)], gmeansTable[gmeansTable$treat==goodTreatments[t3], "data"][1:(xlimit * jump)], alternative ="two.sided")
        write(paste("Kolmogorov-Smirnov test p.value for treatments ", goodTreatments[t2], " vs. ", goodTreatments[t3], " : ", format(w$p.value, digits = 2), ", for: ", feature, " max ", xlimit, " seconds", sep=""), file=logFile, append = TRUE)
        mtext(side = 1, at = 0, line = line,
              text = paste("K-S test for ", goodTreatments[t2], " vs. ", goodTreatments[t3], " : ", format(w$p.value,  digits = 2), ",  for: ", feature, sep=""), col = "black", cex= 0.9, adj=0)
        line <- line + 1}}}
  # write all distributions for a permutation test
  tablePath <- paste(outputdir, "/", basename, "_", feature, "_distributions.csv", sep="")
  csvwell <- paste(outputdir, "/", get.project.plate.name(s$file), "_", feature, "_distributions.csv", sep="")
  #  write.table(  Alldistributions, file=tablePath, sep=",", append = F, col.names=F, row.names=F )
  write.table( Alldistributions, file=csvwell, sep=",", append = T, col.names=F, row.names=F )
}
