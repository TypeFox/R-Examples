#' run the qqplot simulator script
#'
#' @export
#' @examples
#' \dontrun{qqsim()}
qqsim<-function()
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

  cat("Are you ready for some nice graphics to help illustrate qq plots for different kinds of distributions?\n", fill = TRUE)
  cat("At each step, type your response and press \"Enter\"\n")

  while(TRUE)
  {
    cat("\nHow many data points would you like to simulate?\n")
    # getvalidinput is a function that does exactly what it sounds like. See documentation above.
    data.size <- getvalidinput("Enter a number between 1 and 10000 or 'q'to Quit: ", "Invalid selection! Choose a number from 1 to 10000!\n",
                               validinput = 'q', low = 1, high = 1e4);
    if(data.size == 'q') { cat("Quitting!\n"); break() }
    data.size <- as.numeric(data.size)
    cat("\nWould you like to see qq-plots and histograms for data sampled from a\n")
    cat("1) normal distribution?\n")
    cat("2) uniform distribution?\n")
    cat("3) skewed right distribution?\n")
    cat("4) skewed left distribution?\n")
    cat("5) bimodal distribution?\n")
    cat("6) bell shaped distribution with outliers?\n")

    distribution <- getvalidinput("Please enter your selection or q to quit: ", "Invalid selection! enter 1, 2, 3, 4, 5, 6 or q \n", validinput = c('1', '2', '3', '4', '5', '6','q'))
    if(distribution == 'q') { cat("\nQuitting!\n"); break() }
    distribution <- as.numeric(distribution)



    switch(distribution,
           {
             #normal distribution
             data1 = rnorm(data.size)
             data2 = rnorm(data.size)
           },
  {
    #uniform distribution
    data1 = runif(data.size)
    data2 = runif(data.size)
  },
  {
    #skewed right distribution
    data1 = rchisq(data.size, 8, ncp = 0)
    data2 = rchisq(data.size, 8, ncp = 0)
  },
  {
    #skewed left distribution
    data1 = -rchisq(data.size, 8, ncp = 0)
    data2 = -rchisq(data.size, 8, ncp = 0)
  },
  {
    #bimodal distribution
    data1 = c(rnorm(data.size/2), rnorm(data.size/2, 4))
    data2 = c(rnorm(data.size/2), rnorm(data.size/2, 4))
  },
  {
    #outliers/heavy tail distribution
    data1 = c(rnorm(data.size*.8), rt(data.size *.2, 3))
    data2 = c(rnorm(data.size*.8), rnorm(data.size *.2, 0, 3))
  }) # end switch


    #generate plots

    par(mfrow = c(2,2), omi = c(0,0, .5,0))
    qq.title1 = "Normal Q-Q Plot of Data Set 1"
    qq.title2 = "Normal Q-Q Plot of Data Set 2"
    hist.title1 = "Histogram of Data Set 1"
    hist.title2 = "Histogram of Data Set 2"

    if(distribution == 1){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1,, main = hist.title1, freq = FALSE)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, main = hist.title2, freq = FALSE)
      mtext(paste("Normal Distribution"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    if(distribution == 2){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1, main = hist.title1, freq = FALSE)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, main = hist.title2, freq = FALSE)
      mtext(paste("Uniform Distribution"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    if(distribution == 3){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1, freq = FALSE, main = hist.title1)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, freq = FALSE, main = hist.title2)
      mtext(paste("Skewed Right Distribution"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    if(distribution == 4){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1, freq = FALSE, main = hist.title1)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, freq = FALSE, main = hist.title2)
      mtext(paste("Skewed Left Distribution"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    if(distribution == 5){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1, freq = FALSE, main = hist.title1)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, freq = FALSE, main = hist.title2)
      mtext(paste("Bimodal Distibution"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    if(distribution == 6){

      qqnorm(data1, main = qq.title1)
      qqline(data1)
      hist(data1, freq = FALSE, main = hist.title1)
      qqnorm(data2, main = qq.title2)
      qqline(data2)
      hist(data2, freq = FALSE, main = hist.title2)
      mtext(paste("Normal Distribution with Possible Outliers"), side = 3, line = -1, cex = 2, outer = TRUE)

    }
    ################################################
    #                                              #
    #           option to save plots here          #
    #                                              #
    ################################################
    saveplots <- getvalidinput("Do you want to save your plots to the desktop? y/n: ", "Enter (y)es or (n)o!\n", validinput=c('y','n', 'q'))
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
        name <- paste("~/Desktop/qq_plot_Rplot_", td, sep="")
      }
      if(os == 'windows')
      {
        name <- strsplit(as.character(Sys.getenv('HOME')), 'Documents')
        name <- paste(name, "Desktop\\qq_plot_Rplot_", td, sep="")
      }
      name <- paste(name, ".jpg", sep="")
      jpeg(name)

      #generate plots to save

      par(mfrow = c(2,2), omi = c( 0,0, .25,0), mar = c(4,4, 9,1))
      if(distribution == 1){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1,, main = hist.title1, freq = FALSE)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, main = hist.title2, freq = FALSE)
        mtext(paste("Normal Distribution"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      if(distribution == 2){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1, main = hist.title1, freq = FALSE)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, main = hist.title2, freq = FALSE)
        mtext(paste("Uniform Distribution"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      if(distribution == 3){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1, freq = FALSE, main = hist.title1)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, freq = FALSE, main = hist.title2)
        mtext(paste("Skewed Right Distribution"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      if(distribution == 4){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1, freq = FALSE, main = hist.title1)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, freq = FALSE, main = hist.title2)
        mtext(paste("Skewed Left Distribution"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      if(distribution == 5){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1, freq = FALSE, main = hist.title1)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, freq = FALSE, main = hist.title2)
        mtext(paste("Bimodal Distibution"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      if(distribution == 6){

        qqnorm(data1, main = qq.title1)
        qqline(data1)
        hist(data1, freq = FALSE, main = hist.title1)
        qqnorm(data2, main = qq.title2)
        qqline(data2)
        hist(data2, freq = FALSE, main = hist.title2)
        mtext(paste("Normal Distribution with Possible Outliers"), side = 3, line = -1, cex = 1.5, outer = TRUE)
        mtext(paste("by", titlename), side = 3, line = -3, cex = 0.8, outer = TRUE)

      }
      dev.off()

    }
    #
    #################################################
  } # end while loop
}
