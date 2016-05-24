###########################################################
###########################################################
###
### Function to show frequencies in a manner similar to
### the way SPSS' "FREQUENCIES" command does.
###
### File created by Gjalt-Jorn Peters. Questions? You can
### contact me through http://behaviorchange.eu
###
### This file is licensed under Creative Commons BY-SA 3.0
### (Attribution-ShareAlike, which means that you can
### freely use and distribute this file, and you're
### allowed to alter it as long as you release the edited
### version using the same license (i.e. again freely
### available). This license is used to promote Open
### Science and Full Disclosure. For the complete
### license, see
### http://creativecommons.org/licenses/by-sa/3.0/deed.en_US
### For more information about Full Disclosure, see
### http://sciencerep.org/fulldisclosure
###
###########################################################
###########################################################

freq <- function(vector, digits = 1, nsmall=1, transposed=FALSE, round=1,
                 plot=FALSE) {
  
  if(length(vector)<2) {
    stop("The first argument is not a vector! Did you make a typing error? Remember that R is case sensitive!");
  }
  
  if (length(unique(vector)) < 2) {
    stop("There are less than two unique elements in the vector you supplied: ",
         vecTxt(unique(vector), useQuote='"'), ".");
  }
  
  suppressWarnings(variance <- var(as.numeric(vector), na.rm=TRUE));
  
  if (is.null(variance) || is.na(variance)) {
    res <- paste0("I'm unable to look at variance in the variable you specified ('",
                  deparse(substitute(vector)),
                  "') - I get NA or NULL. Trying to convert to a factor.\n");
    vector <- factor(vector);
  }

  suppressWarnings(variance <- var(as.numeric(vector), na.rm=TRUE));
  
  if (is.null(variance) || is.na(variance)) {
    res <- paste0("Sorry, the conversion didn't solve anything. I'm aborting.\n");
    return(res);
  } #else if (duplicated(vector)[-1L]) {
#     res <- paste0("The variable you specified ('", deparse(substitute(vector)),
#         "') has no variance; all values are the same (they are all ",
#         mean(vector), ").");
#     return(res);
#   }
  
  ### Create object to store results
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  ### Store input data
  res$input$vector <- as.factor(vector);  
  
  ### Store category names
  res$intermediate$categoryNames <- levels(res$input$vector);
  ### Store data without missing values
  res$intermediate$vector.valid <- res$input$vector[!is.na(res$input$vector)];
  ### Store frequencies based on full data
  res$intermediate$frequencies.raw <- summary(res$input$vector);
  ### Store frequencies based on data without missing values
  res$intermediate$frequencies.raw.valid <- summary(res$intermediate$vector.valid);
  ### Store proportions based on full data
  res$intermediate$frequencies.prop <- summary(res$input$vector) /
                                       length(res$input$vector);
  ### Store proportions based on data without missing values
  res$intermediate$frequencies.prop.valid <- summary(res$intermediate$vector.valid) /
                                             length(res$intermediate$vector.valid);
  ### Compute cumulative percentages
  res$intermediate$frequencies.prop.cum <- res$intermediate$frequencies.prop.valid;
  for (currentPropIndex in
         2:length(res$intermediate$frequencies.prop.valid)) {
    res$intermediate$frequencies.prop.cum[currentPropIndex] <-
      res$intermediate$frequencies.prop.cum[currentPropIndex - 1] +
      res$intermediate$frequencies.prop.cum[currentPropIndex];
  }

  ### Now we integrate this in a dataframe to show the users. First
  ### ignoring the missing values.
  
  res$intermediate$frequencies.prop.clipped <- res$intermediate$frequencies.prop;
  if (length(res$intermediate$frequencies.prop) > length(res$intermediate$frequencies.prop.valid)) {
    res$intermediate$frequencies.prop.clipped <- res$intermediate$frequencies.prop[1:length(res$intermediate$frequencies.prop)-1];
  }
    
  res$dat <- data.frame(Frequencies = res$intermediate$frequencies.raw.valid,
                        Perc.Total = 100*res$intermediate$frequencies.prop.clipped,
                        Perc.Valid = 100*res$intermediate$frequencies.prop.valid,
                        Cumulative = 100*res$intermediate$frequencies.prop.cum);
  
  ### We then add a row with the totals.
  res$dat <- rbind(res$dat, c(sum(res$intermediate$frequencies.raw.valid),
                              100*sum(res$intermediate$frequencies.prop.clipped),
                              100*sum(res$intermediate$frequencies.prop.valid),
                              NA));
  rownames(res$dat)[nrow(res$dat)] <- "Total valid"
  
  ### Then, if we have missing values, we add the number and percentage of missing values,
  ### as well as the totals for these two columns.
  if (length(res$intermediate$frequencies.prop) > length(res$intermediate$frequencies.prop.valid)) {
    res$dat <- rbind(res$dat, c(res$intermediate$frequencies.raw[length(res$intermediate$frequencies.raw)],
                                100*res$intermediate$frequencies.prop[length(res$intermediate$frequencies.prop)],
                                NA,
                                NA));
    res$dat <- rbind(res$dat, c(sum(res$intermediate$frequencies.raw),
                                100*sum(res$intermediate$frequencies.prop),
                                NA,
                                NA));
    rownames(res$dat)[nrow(res$dat)-1] <- "NA (missing)"
    rownames(res$dat)[nrow(res$dat)] <- "Total"
  }
  
  if (!is.null(round)) {
    tempRowNames <- rownames(res$dat);
    res$dat <- data.frame(lapply(res$dat, function(x) {return(round(x, digits=round));}));
    rownames(res$dat) <- tempRowNames;
  }
  
  if (plot) {
    res$barChart <- plot(res$intermediate$vector.valid);
  }
  
  ## Set object class;
  class(res) <- c("freq");
  return(res);
}

print.freq <- function(x, digits=x$input$digits, nsmall=x$input$nsmall,
                       transposed=x$input$transposed, ...) {
  if (transposed) {
    ### Transpose dataframe
    x$dat <- data.frame(t(x$dat));
    ### Round frequencies and percentages and convert to character vector
    prettyDat <- format(x$dat, digits=digits, nsmall=nsmall);
    ### Replace missing values with a space
    prettyDat <- data.frame(lapply(prettyDat,
                                   function(x) {return(sub("NA", "  ", x));}));
    ### Replace formatted first row with original first row (without
    ### decimals)
    prettyDat[1, ] <- x$dat[1, ];
    ### Add column and row names again
    rownames(prettyDat) <- rownames(x$dat);
    colnames(prettyDat) <- colnames(x$dat);
    ### Print result
    print(prettyDat, ...);
  }
  else {
    ### Round frequencies and percentages and convert to character vector
    prettyDat <- format(x$dat, digits=digits, nsmall=nsmall);
    ### Replace missing values with a space
    prettyDat <- data.frame(lapply(prettyDat,
                                   function(x) {return(sub("NA", "  ", x));}));
    ### Replace formatted first column with original first column (without
    ### decimals)
    prettyDat$Frequencies <- x$dat$Frequencies;
    ### Add row names again
    rownames(prettyDat) <- rownames(x$dat);
    ### Print result
    print(prettyDat, ...);
  }
  if (x$input$plot) {
    x$barChart;
  }
  invisible();
}
