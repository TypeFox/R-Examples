STEPCAM_ABC <- function(data_abundances, data_species, numParticles, n_traits, plot_number, stopRate){
  
  # empty vector for number of species in each community    
  nbsp <- c() 
  
  # number of plots
  n_plots <- nrow(data_abundances) 

  # species present in the community
  esppres <- which(data_abundances[plot_number, ] > 0) ; 
  
  # species richness
  S <- length(esppres); 

  # total number of species of species pool
  taxa <- nrow(data_species); 
  
  # species that are removed
  species_fallout <- taxa - S 

  # frequencies (number of plots of occurrence) of each species
  data_frequencies <- generateFrequencies(data_abundances); 

  # calculate summary statistics of observed focal plot
  # scale trait values to mean = 0, SD = 1
  scaled_species <- scaleSpeciesvalues(data_species,n_traits) 
  observed_traits <- scaled_species[, ] ; row.names(observed_traits) <- scaled_species[, 1]
  observed_abundances <- data_abundances[plot_number, ]
  present_species <- which(observed_abundances > 0)
  observed_traits <- observed_traits[present_species, ]
  observed_abundances <- observed_abundances[present_species]
  observed_presences <- observed_abundances ; observed_presences[1:length(observed_presences)] <- 1
  observed_traits <- observed_traits[, -1]
  
  # calculate FD values observed communities
  Ord <- ordinationAxes(x = scaled_species[,-1], stand.x = FALSE)
  FD_output <- strippedDbFd(Ord, ifelse(data_abundances > 0, 1, 0)) 
  
  trait_means <- c()
  traitvalues <- c()
  for(i in 1:n_traits){
    for(j in 1:length(observed_abundances)){
      traitvalues[j] <- observed_traits[j, i]
    }
    # calculate CTM value
    trait_means[i] <- mean(traitvalues) 
  }
  # bind FD values and CTM together
  summary_stats <- cbind(FD_output$FRic[plot_number],FD_output$FEve[plot_number],FD_output$FDiv[plot_number], t(trait_means)) # bind FD values and CTM together


  # calculate the SD of FD/CTM values: this is used to asses STEPCAM model fit
  sd_vals <- calcSD(scaled_species, data_abundances, n_plots, n_traits); 

  scaled_species <- as.data.frame(cbind(scaled_species[, c(1:(n_traits + 1))], data_frequencies)) ;
  traitnames <- names(data_species)[-1] 
  names(scaled_species) <- c("sp", traitnames[1:n_traits], "freq")
  row.names(scaled_species) <- c(1:taxa)

  

  output <- ABC_SMC(numParticles, species_fallout, taxa,esppres, n_traits,
  sd_vals, summary_stats, plot_number, scaled_species, data_abundances,
  data_frequencies, stopRate, Ord);
 
  return(output);  
}


plotSTEPCAM <- function(output){
  total <- output$DA[1] + output$HF[1] + output$LS[1];
  par(mfrow=c(1, 3));
  hist(output$DA / total, col="grey", main = "Dispersal Assembly", xlim = c(0, 1), ylab = "", xlab = "");
  hist(output$HF / total, col="grey" , main = "Habitat Filtering", xlim = c(0, 1), ylab = "", xlab = "");
  hist(output$LS / total, col="grey" , main = "Limiting Similarity", xlim = c(0, 1), ylab = "", xlab = "");
}

plotElement <- function(d, xmin, xmax, index, maxTime, parameter){
  topLabels <- c("Dispersal Assembly", "Habitat Filtering", "Limiting Similarity", "Richness", "Evenness",
  "Diversity", "Optimum", "Fit");
  mean_data <- mean(d) 
  stdev_data <- sd(d)
 
  x <- seq(xmin, xmax, length = 200);
  y <- dnorm(x, mean = mean_data, sd = stdev_data);
  medY <- 0.8 * max(y);
  title <- ""
  if(index %% maxTime == 1) title <- topLabels[parameter]

  if(index %% maxTime == 0 ) 
  {    
    hist(d, xlim = c(xmin, xmax),main = title, col = "grey", ylab = "",
    xlab = "", yaxt = "n", cex.main = 1.5);
  }
  else  
  {
     hist(d, xlim = c(xmin, xmax), main = title, col = "grey", ylab = "",
     xlab = "", yaxt = "n", xaxt = "n", cex.main = 1.5);
  } 
  par(new = TRUE) 
  plot(y~x, col = "red", type = "l", axes = FALSE, ylab = "", xlab = "", main = "", yaxt = "n");
}

plotSMC <- function(path){  
  end <- ".txt";
  val <- "particles_t="
  calctotal <- read.table(paste(path, val, 1, end, sep = "", collapse = NULL))
  total <- calctotal$V1[1] + calctotal$V2[1] + calctotal$V3[1];
  
  maxTime <- 1;
    
  # determine how many iterations need to be plotted
  for(k in 1:200){
    file_name <- paste(path, val, 200 - k,end, sep = "", collapse = NULL);
    if(file.exists(file_name)){
      maxTime <- 200 - k;
      break;
    }
  }
  maxTime <- maxTime-1;
  numCols <- 8;
  numRows <- maxTime;

  k <- matrix(nrow = numRows, ncol = numCols); 
  cnt <- 1;
  for(i in 1:numCols)
  {
    for(j in 1:numRows)
    {
       k[j, i] <- cnt
       cnt <- cnt + 1;
    }
  } 
  layout(k)
  par(mar = c(0.5, 0, 0.5, 0));
  
  for(c in 1:8)
  {
    fulldata <- c();

    for(i in 1:(maxTime))
    {
      data_name <- paste(path, val, i, end, sep = "", collapse = NULL)
      if(file.exists(data_name))
      {
        data <- read.table(data_name, header = F);
        fulldata <- cbind(fulldata, data[, c]);
      }
    }

    xmin <- min(fulldata);
    xmax <- max(fulldata);
    if(xmin == xmax){xmin = 0.9 * xmin; xmax = 1.1 * xmax;}
    if(xmin == xmax && xmin == 0){xmin = -1; xmax = 1;}
    if(c < 4)
    { 
      fulldata <- fulldata / total
      xmin <- 0;
      xmax <- 1;
    }
    for(i in 1:length(fulldata[1,]))
    {
      plotElement(fulldata[,i], xmin, xmax, i, maxTime, c)
    }
  }
}


  
  
  
  






