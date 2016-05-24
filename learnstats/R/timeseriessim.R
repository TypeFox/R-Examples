#' Run the Time Series Simulation script.
#'
#' @export
#' @examples
#' \dontrun{timeseriessim()}
timeseressim<-function()
{
  ### time series script for R.
  ### by Omar Chavez
  ### 1.0   - First version as cribbed from prob.R

  ##' Given a prompt and acceptable inputs, query user for input until they enter an acceptable value.
  ##' @author Omar Chavez
  getvalidinput <- function(prompt, error, validinput = NULL, low=0, high=0)
  {
    if(is.null(validinput) && low == high) stop("Must either input validinput or both low and high!")
    while(TRUE)
    {
      input <- readline(prompt)
      if(input %in% validinput)
      {
        return(input)
      }
      if(!is.na(as.numeric(input)))
      {
        if(low < high && low < high && low <= as.numeric(input) && as.numeric(input) <= high)
        {
          return(input)
        }
      }
      cat(error)
    }
  }

  options(warn=-1) # make this 0 for testing, -1 for production
  options(scipen=4) #Don't give numbers like 1e-4, instead do .0001

  cat("Are you ready for some nice graphics to help illustrate stability in time plots?\n")
  cat("At each step, type your response and press \"Enter\"\n\n")

  while(TRUE)
  {
    cat("\nHow many data points would you like to simulate?\n")
    # getvalidinput is a function that does exactly what it sounds like. See documentation above.
    data.size <- getvalidinput("Enter a number between 1 and 10000 or 'q'to Quit: ", "Invalid selection! Choose a number from 1 to 10000!\n",
                               validinput = 'q', low = 1, high = 1e4);
    if(data.size == 'q') { cat("Quitting!\n"); break() }
    data.size <- as.numeric(data.size)
    stable.unstable <- getvalidinput("\nWould you like to generate stable or unstable plots?\n1) stable plots\n2) unstable plots \nPlease enter your selection or q to quit: ", "Invalid selection! enter1, 2 or q \n", validinput = c('1', '2', 'q'))
    if(stable.unstable == 'q') { cat("\nQuitting!\n"); break() }
    stable.unstable <- as.numeric(stable.unstable)

    #if unstable plot selected get info about mean and variance
    if(stable.unstable == 2){

      mean = 0
      variance = 1
      mean.and.variance = getvalidinput("\nWould you like to have \n1) changing mean only? \n2) changing variance only? \n3) both changing mean AND changing variance? \nPlease enter your selection or q to quit: ",   "Invalid selection! enter 1, 2, 3 or q.\n",
                                        validinput = c('1', '2', '3', 'q'))
      if( mean.and.variance == 'q') { cat("\nQuitting!\n"); break() }
      mean.and.variance <- as.numeric(mean.and.variance)

      #changing mean
      if(mean.and.variance == 1 | mean.and.variance == 3){

        mean.selection <- getvalidinput("\nWould you like the mean to \n1) Increase?\n2) Decrease?\n3) Oscillate (fluctuate)? \nPlease enter your selection or q to quit: ", "Invalid selection! enter 1,2,3 or q \n", validinput = c('1', '2', '3', 'q'))
        if(mean.selection == 'q') { cat("\nQuitting!\n"); break() }
        mean.selection <- as.numeric(mean.selection)
        #increasing mean
        if(mean.selection == 1){

          mean = c(1:data.size)* (3/data.size)
          if(mean.and.variance == 1){
            description = "increasing mean, constant variance"
          }
          if(mean.and.variance == 3){
            description = "increasing mean"
          }

        }
        #decreasing mean
        if(mean.selection == 2){

          mean = c(data.size:1)* (3/data.size)
          if(mean.and.variance == 1){
            description = "decreasing mean, constant variance"
          }
          if(mean.and.variance == 3){
            description = "decreasing mean"
          }

        }
        #oscillating mean
        if(mean.selection == 3){

          mean = 2 * cos( c(1:data.size)*4*pi/data.size )
          if(mean.and.variance == 1){
            description = "oscillating (fluctuating) mean, constant variance"
          }
          if(mean.and.variance == 3){
            description = "oscillating (fluctuating) mean"
          }

        }
      }

      #changing variance
      if(mean.and.variance == 2 | mean.and.variance == 3){

        variance.selection <- getvalidinput("\nWould you like the variance to \n1) Increase?\n2) Decrease?\n3) Oscillate (fluctuate)? \nPlease enter your selection or q to quit: ", "Invalid selection! enter 1,2,3 or q \n", validinput = c('1', '2', '3', 'q'))
        if(variance.selection == 'q') { cat("\nQuitting!\n"); break() }
        variance.selection <- as.numeric(variance.selection)

        #increasing variance
        if(variance.selection == 1){

          variance = c(1:data.size)* (3/data.size)
          if(mean.and.variance == 2){
            description = "constant mean, increasing variance"
          }
          if(mean.and.variance == 3){
            description = paste(description, ", increasing variance")
          }

        }
        #decreasing variance
        if(variance.selection == 2){

          variance = c(data.size:1)* (3/data.size)
          if(mean.and.variance == 2){
            description = "constant mean, decreasing variance"
          }
          if(mean.and.variance == 3){
            description = paste(description, ", decreasing variance")
          }

        }
        #oscillating variance
        if(variance.selection == 3){

          variance = 2 * cos( c(1:data.size)*4*pi/data.size )
          if(mean.and.variance == 2){
            description = "constant mean, oscillating (fluctuating) variance"
          }
          if(mean.and.variance == 3){
            description = paste(description, ", oscillating (fluctuating) variance")
          }

        }
      }

    }

    switch(stable.unstable,
           {
             #stable plots
             series1 = rnorm(data.size)
             series2 = rnorm(data.size)
             series3 = rnorm(data.size)
             series4 = rnorm(data.size)
           },
  {
    #unstable plots
    series1 = variance * rnorm(data.size) + mean
    series2 = variance * rnorm(data.size) + mean
    series3 = variance * rnorm(data.size) + mean
    series4 = variance * rnorm(data.size) + mean
  }) # end switch


    #generate
    par(mfrow = c(2,2), omi = c(0,0, .5,0), mar = c(4,4, 9,1))
    if(stable.unstable == 1){
      plot(series1, type = 'l', main = "Simulation 1")
      plot(series2, type = 'l', main = "Simulation 2")
      plot(series3, type = 'l', main = "Simulation 3")
      plot(series4, type = 'l', main = "Simulation 4")
      mtext(paste("Stable Plots"), side = 3, line = -1, cex = 2, outer = TRUE)
      mtext("constant mean, constant variance", side = 3, line = -3, cex = 1, outer = TRUE)
    }
    if(stable.unstable == 2){
      plot(series1, type = 'l', main = "Simulation 1")
      plot(series2, type = 'l', main = "Simulation 2")
      plot(series3, type = 'l', main = "Simulation 3")
      plot(series4, type = 'l', main = "Simulation 4")
      mtext(paste("Unstable Plots"), side = 3, line = -1, cex = 2, outer = TRUE)
      mtext(paste(description), side = 3, line = -3, cex = 1, outer = TRUE)
    }
    ################################################
    #                                              #
    #           option to save plots here          #
    #                                              #
    ################################################
    saveplots <- getvalidinput("Do you want to save your plot to the desktop? y/n: ", "Enter (y)es or (n)o!\n", validinput=c('y','n', 'q'))
    if(saveplots == 'q') { cat("Quitting!\n"); break() }
    if(saveplots == 'y')
    {
      titlename <- readline("Please enter your name for the plot: ")
      while(TRUE)
      {
        if(nchar(titlename) > 0) break()
        titlename <- readline("Must enter a name! Enter your name: ")
      }
      os <- .Platform$OS
      td <- gsub(":", ".", as.character(Sys.time()))
      td <- gsub(" ", "_", td)
      if(os == 'unix') #unix == linux or mac
      {
        name <- paste("~/Desktop/time_series_Rplot_", td, sep="")
      }
      if(os == 'windows')
      {
        name <- strsplit(as.character(Sys.getenv('HOME')), 'Documents')
        name <- paste(name, "Desktop\\time_series_Rplot_", td, sep="")
      }
      name <- paste(name, ".jpg", sep="")
      jpeg(name)

      #generate plots to save
      if(stable.unstable == 1){
        par(mfrow = c(2,2), omi = c( 0,0, .25,0), mar = c(4,4, 9,1))
        plot(series1, type = 'l', main = "Simulation 1",xlab = "index", ylab = "time series value")
        plot(series2, type = 'l', main = "Simulation 2",xlab = "index", ylab = "time series value")
        plot(series3, type = 'l', main = "Simulation 3",xlab = "index", ylab = "time series value")
        plot(series4, type = 'l', main = "Simulation 4",xlab = "index", ylab = "time series value")
        mtext(paste("Stable Plots - ", "constant mean, constant variance"), side = 3, line = -1, cex = 1, outer = TRUE)
        mtext(paste("by ",titlename), side = 3, line = -3, cex = 1, outer = TRUE)
      }
      if(stable.unstable == 2){
        par(mfrow = c(2,2), omi = c(0,0, .25,0), mar = c(4,4, 9,1))
        plot(series1, type = 'l', main = "Simulation 1",xlab = "index", ylab = "time series value")
        plot(series2, type = 'l', main = "Simulation 2",xlab = "index", ylab = "time series value")
        plot(series3, type = 'l', main = "Simulation 3",xlab = "index", ylab = "time series value")
        plot(series4, type = 'l', main = "Simulation 4",xlab = "index", ylab = "time series value")
        mtext(paste("Unstable Plots - ", description), side = 3, line = -1, cex = 1, outer = TRUE)
        mtext(paste("by ", titlename), side = 3, line = -3, cex = 1, outer = TRUE)
      }

      dev.off()

    }

    ################################################
  } # end while loop
}
