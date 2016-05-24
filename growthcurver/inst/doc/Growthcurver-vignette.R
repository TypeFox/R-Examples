## ---- echo = FALSE, eval = TRUE------------------------------------------
# First, load the package and the dataset. 
library(growthcurver)

# Load the sample growth curve data provided with the package 
# The first column is the time in hours, and there is one column 
# for each well in a 96-well plate.
d <- growthdata
knitr::kable(d[1:10, 1:8])

## ---- eval = FALSE-------------------------------------------------------
#  # Replace the next line with the location and name of your input data file.
#  file_name <- "the/path/to/my/data/myfilename.txt"
#  d <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## ---- eval = FALSE-------------------------------------------------------
#  # Convert the "time" column from hours to minutes
#  d$time <- d$time * 60
#  
#  # Convert the "time" column from minutes to seconds
#  d$time <- d$time * 60
#  
#  # Convert the "time" column from seconds to hours
#  d$time <- d$time / 60 / 60

## ---- eval = TRUE--------------------------------------------------------
# First, load the package. 
library(growthcurver)

# Load the sample growth curve data provided in the Growthcurver package.
# The first column is the time in hours, and there is one column 
# for each well in a 96-well plate.
d <- growthdata

# Now, we'll use Growthcurver to summarize the growth curve data using the 
# simple background correction method (minimum value correction). This is the 
# default method, so we don't need to specify it in the command.
# This returns an object of type "gcfit" that holds information about
# the best parameters, the model fit, and additional metrics summarizing
# the growth curve.
gc_fit <- SummarizeGrowth(d$time, d$A1)

# It is easy to get the most useful metrics from a gcfit object, just type:
gc_fit

# And it is easy to plot the raw data and the best fit logistic curve
plot(gc_fit)

## ---- eval = FALSE-------------------------------------------------------
#  # The gcfit object returned from SummarizeGrowth also contains further metrics
#  # summarizing the growth curve data.
#  gc_fit$vals
#  
#  # look at the structure of the gc_fit object
#  str(gc_fit)

## ---- eval = TRUE--------------------------------------------------------
# To see all the available metrics 
str(gc_fit$vals)

# To access a single metric (for example the growth rate r)
gc_fit$vals$r


## ---- eval = TRUE--------------------------------------------------------
# First, load the package and the sample dataset. 
library(growthcurver)
d <- growthdata

## ---- eval = FALSE-------------------------------------------------------
#  # To analyze your data from Excel, you should read your data into the variable
#  # called d. To do so, replace the next line with the name and location of
#  # your input data file.
#  file_name <- "the/path/to/my/data/myfilename.txt"
#  d <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#  
#  # Make sure that you have a column called "time" (and a column called "blank"
#  # if you are using "blanks" for your background correction). See the
#  # "Input Data" data section of the Vignette if you need help with this.

## ---- eval = TRUE--------------------------------------------------------
# Now, we'll use Growthcurver to summarize the growth curve data for the entire
# plate using the default background correction method ("min").
gc_out <- SummarizeGrowthByPlate(d)

## ---- eval = FALSE-------------------------------------------------------
#  # If you would like to use the "blank" background correction, then call
#  # Growthcurver as follows
#  gc_out <- SummarizeGrowthByPlate(d, bg_correct = "blank")
#  
#  # If you would like to generate plots for all of the growth curves in your
#  # plate, then call Growthcurver as follows. You can change the name of
#  # the output file "gc_plots.pdf" to something that makes sense for you.
#  gc_out <- SummarizeGrowthByPlate(d, plot_fit = TRUE,
#                                   plot_file = "gc_plots.pdf")
#  
#  # The summary information for each well is listed as a row in the output
#  # data frame called gc_out.
#  
#  # We can look at the first few rows in the output using the head command.
#  head(gc_out)

## ---- eval = TRUE, echo = FALSE------------------------------------------
knitr::kable(gc_out[1:5, ])

## ---- eval = FALSE-------------------------------------------------------
#  # Or, you can save the entire data table to a tab-separated file that can be
#  # imported into Excel.
#  output_file_name <- "the/path/to/my/data/myfilename.txt"
#  write.table(gc_out, file = output_file_name,
#              quote = FALSE, sep = "\t", row.names = FALSE)

## ---- message = FALSE, fig.width = 7-------------------------------------
# As in the simple example, load the package and the data. 
library(growthcurver)
d <- growthdata

# Let's create an output data frame to store the results in. 
# We'll create it so that it is the right size (it's faster this way!), 
# but leave it empty.
num_analyses <- length(names(d)) - 1
d_gc <- data.frame(sample = character(num_analyses),
                   k = numeric(num_analyses),
                   n0  = numeric(num_analyses),
                   r = numeric(num_analyses),
                   t_mid = numeric(num_analyses),
                   t_gen = numeric(num_analyses),
                   auc_l = numeric(num_analyses),
                   auc_e = numeric(num_analyses),
                   sigma = numeric(num_analyses),
                   stringsAsFactors = FALSE)

# Truncate or trim the input data to observations occuring in the first 20 hours.
# Remember that the times in these sample data are reported in hours. To use  
# minutes (or to trim at a different time), change the next line of code. 
# For example, if you still would like to trim at 20 hours, but your time data 
# are reported in minutes use: trim_at_time <- 20 * 60
trim_at_time <- 20   

# Now, loop through all of the columns in the data frame. For each column,
# run Growthcurver, save the most useful metrics in the output data frame,
# and make a plot of all the growth curve data and their best fits.

# First, create a plot for each of the wells in the 96-well plate.
# Uncomment the next line to save the plots from your 96-well plate to a 
# pdf file in the working directory.
# pdf("growthcurver.pdf", height = 8.5, width = 11)
par(mfcol = c(8,12))
par(mar = c(0.25,0.25,0.25,0.25))
y_lim_max <- max(d[,setdiff(names(d), "time")]) - min(d[,setdiff(names(d), "time")])

n <- 1    # keeps track of the current row in the output data frame
for (col_name in names(d)) {
  
  # Don't process the column called "time". 
  # It contains time and not absorbance data.
  if (col_name != "time") {

    # Create a temporary data frame that contains just the time and current col
    d_loop <- d[, c("time", col_name)]
    
    # Do the background correction.
    # Background correction option 1: subtract the minimum value in a column
    #                                 from all measurements in that column
        min_value <- min(d_loop[, col_name])
    d_loop[, col_name] <- d_loop[, col_name] - min_value
    # Background correction option 2: subtract the mean value of blank wells
    #                                 over the course the experiment
    #                                 (Replace B2, D8, G11 with the column
    #                                  names of your media-only wells)
    #d$blank <- apply(d[, c("B2", "D8", "G11")], 1, mean)
    #d$A1 <- d$A1 - d$blank
    
    # Now, call Growthcurver to calculate the metrics using SummarizeGrowth
    gc_fit <- SummarizeGrowth(data_t = d_loop[, "time"], 
                              data_n = d_loop[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    # Now, add the metrics from this column to the next row (n) in the 
    # output data frame, and increment the row counter (n)
    d_gc$sample[n] <- col_name
    d_gc[n, 2:9] <- c(gc_fit$vals$k,
                      gc_fit$vals$n0,
                      gc_fit$vals$r,
                      gc_fit$vals$t_mid,
                      gc_fit$vals$t_gen,
                      gc_fit$vals$auc_l,
                      gc_fit$vals$auc_e,
                      gc_fit$vals$sigma)
    n <- n + 1
    
    # Finally, plot the raw data and the fitted curve
    # Here, I'll just print some of the data points to keep the file size smaller
    n_obs <- length(gc_fit$data$t)
    idx_to_plot <- 1:20 / 20 * n_obs
    plot(gc_fit$data$t[idx_to_plot], gc_fit$data$N[idx_to_plot], 
         pch = 20, 
         xlim = c(0, trim_at_time), 
         ylim = c(0, y_lim_max),
         cex = 0.6, xaxt = "n", yaxt = "n")
     text(x = trim_at_time / 4, y = y_lim_max, labels = col_name, pos = 1)
     lines(gc_fit$data$t, predict(gc_fit$model), col = "red")
  }
}
# Uncomment the next line to save the plots from your 96-well plate to a file
# dev.off()

## ---- eval = FALSE-------------------------------------------------------
#  # Look at the first few rows (samples) of data in the output data frame.
#  # (I'm only showing the first 4 rows of results, but you may want to see more.
#  #  You can either look at everything using the command "d_gc", or adjust the
#  #  number of rows displayed by changing the 4 to something else,
#  #  e.g., "d_gc[1:15,]").
#  d_gc[1:4, ]

## ---- eval = TRUE, message = FALSE, echo = FALSE-------------------------
library(dplyr)
d_gc[1:4, ] %>% 
    mutate(k = round(k, digits = 5),
         n0 = round(n0, digits = 5), 
         r = round(r, digits = 5),
         t_mid = round(t_mid, digits = 5),
         t_gen = round(t_gen, digits = 5),
         auc_l = round(auc_l, digits = 5),
         auc_e = round(auc_e, digits = 5), 
         sigma = round(sigma, digits = 5))
  

## ---- eval = FALSE, message = FALSE--------------------------------------
#  # Check if Growthcurver provided any notes in a plate of growthcurves returned
#  # from SummarizeGrowthByPlate
#  gc_out %>% filter(note != "")
#  
#  # Check if Growthcurver provided any notes in a single growthcurve returned
#  # from SummarizeGrowth
#  gc_fit$vals$note

## ---- eval = TRUE, message = FALSE---------------------------------------
# Load dplyr and the sample output data
library(dplyr)
gc_out <- as_data_frame(gc_out)

# Plot a histogram of the sigma values in order to check for outliers
hist(gc_out$sigma, main = "Histogram of sigma values", xlab = "sigma")


## ---- eval = FALSE, message = FALSE--------------------------------------
#  # Show the top 5 samples with the largest sigma value
#  # (with the worst model fit to the growth curve data)
#  gc_out %>% top_n(5, sigma) %>% arrange(desc(sigma))

## ---- eval = TRUE, echo = FALSE, message = FALSE-------------------------
gc_out %>%  
  mutate(k = round(k, digits = 5),
         n0 = round(n0, digits = 5), 
         r = round(r, digits = 5),
         t_mid = round(t_mid, digits = 5),
         t_gen = round(t_gen, digits = 5),
         auc_l = round(auc_l, digits = 5),
         auc_e = round(auc_e, digits = 5), 
         sigma = round(sigma, digits = 5)) %>%
  top_n(5, sigma) %>% arrange(desc(sigma))

## ---- eval = TRUE, message = FALSE---------------------------------------
# Load dplyr, ggplot2, and the sample data
library(dplyr)
library(ggplot2)
pca_gc_out <- as_data_frame(gc_out) 

# Prepare the gc_out data for the PCA
rownames(pca_gc_out) <- pca_gc_out$sample

# Do the PCA
pca.res <- prcomp(pca_gc_out %>% select(k:sigma), center=TRUE, scale=TRUE)

# Plot the results
as_data_frame(list(PC1=pca.res$x[,1],
                   PC2=pca.res$x[,2],
                   samples = rownames(pca.res$x))) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)

