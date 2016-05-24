#' Run the Probability Walkthrough script.
#'
#' @export
#' @examples
#' \dontrun{prob()}
prob<-function()
{
  ### Probability script for R.
  ### by Josh Errickson
  ### 1.4   - Changed title to always N(0,1), give prob in legend.
  ### 1.3.4 - Never display probabilies as 0 or 1... just say <.0001 or >.9999
  ### 1.3.3 - Display instructions when script is loaded.
  ### 1.3.2 - Added name to plots.
  ### 1.3.1 - Fixes to user input handling.
  ### 1.3   - Helper function to handle user input - should no longer ever crash! Just keep requesting new input.
  ### 1.2.1 - Added better error checking for first choice.
  ### 1.2   - Auto-saving
  ### 1.1   - Added looping action
  ### 1.0.2 - Minor typo (double lower bound)
  ### 1.0.1 - Clarified some language in prompts.
  ### 1.0   - First version as cribbed from pval.R

  ##' Given a prompt and acceptable inputs, query user for input until they enter an acceptable value.
  ##'
  ##'
  ##' @title
  ##' @param prompt String shown to users to prompt for input.
  ##' @param error Error string shown when user enters an invalid selection.
  ##' @param validinput If non-NULL, only accept input in this list.
  ##' @param low Only accept numeric values above this value.
  ##' @param high Only accept numeric values below this value. If low>=high,
  ##' do not accept any numeric input not explictly in validinput.
  ##' @return Acceptable user input.
  ##' @author Josh Errickson
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

  cat("Are you ready for some nice graphics to help illustrate probabilities?\n")
  cat("At each step, type your response and press \"Enter\"\n")

  while(TRUE)
  {
    cat("What is the form of your probability? Choose from the following:\n1: P(Z < c)\n2: P(Z > c)\n3: P(a < Z < b)\n4: P(Z < a or Z > b)\nYou can press 'q' at any input to quit.\n(Remember: < is computed the same as <= and similarly for > and >=)\n")
    # getvalidinput is a function that does exactly what it sounds like. See documentation above.
    dir <- getvalidinput("Enter your selection: ", "Invalid selection! Choose between 1, 2, 3 or 4!\n", validinput = c('1', '2', '3', '4', 'q'))
    if(dir == 'q') { cat("Quitting!\n"); break() }
    dir <- as.numeric(dir)

    switch(dir,
           {
             # less than
             stat <- getvalidinput("Enter the value for the upper bound c: ", "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8);
             if(stat == 'q') { cat("Quitting!\n"); break() }
             stat <- round(as.numeric(stat), 4);

             pv <- round(pnorm(stat), 4); # p-value

             max <- max(3,abs(stat)+1); # width for x-axis.

             # generate the shaded area.
             xcord <- c(-max,seq(-max,stat,.01),stat);
             ycord <- c(0,dnorm(seq(-max,stat,.01)),0);
             if(pv == 0)
               ti <- substitute(paste("P(Z < ", z, ") < .0001", sep=''), list(z=stat)) # create title.
             else if(pv == 1)
               ti <- substitute(paste("P(Z < ", z, ") > .9999 ", sep=''), list(z=stat)) # create title.
             else
               ti <- substitute(paste("P(Z < ", z, ") = ", p, sep=''), list(z=stat, p=pv)) # create title.
           },
  {
    # greater than
    stat <- getvalidinput("Enter the value for the lower bound c: ", "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8);
    if(stat == 'q') { cat("Quitting!\n"); break("Quitting!\n") }
    stat <- round(as.numeric(stat), 4);

    pv <- round(1 - pnorm(stat), 4);

    max <- max(3,abs(stat)+1);

    xcord <- c(stat,seq(stat,max,.01),max);
    ycord <- c(0,dnorm(seq(stat,max,.01)),0);

    if(pv == 0)
      ti <- substitute(paste("P(Z > ", z, ") < .0001", sep=''), list(z=stat)) # create title.
    else if(pv == 1)
      ti <- substitute(paste("P(Z > ", z, ") > .9999 ", sep=''), list(z=stat)) # create title.
    else
      ti <- substitute(paste("P(Z > ", z, ") = ", p, sep=''), list(z=stat, p=pv)) # create title.
  },
  {
    # between
    stat1 <- getvalidinput("First, enter the value for the lower bound a:" , "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8)
    if(stat1 == 'q') { cat("Quitting!\n"); break() }
    stat1 <- round(as.numeric(stat1), 4);

    cat("Next, enter a value for the upper bound b (must be greater than a).\n");
    while(TRUE)
    {
      stat2 <- getvalidinput("Upper bound: " , "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8)
      if(stat2 == 'q') { cat("Quitting!\n"); break() }
      stat2 <- round(as.numeric(stat2), 4);
      if(stat1 < stat2) { break() } else { cat(paste("Upper bound must be greater than lower bound of ", stat1, "!\n", sep='')) }
    }
    if(stat2 == 'q') { break() } #need second before first only breaks out of inner while

    pv <- round(pnorm(stat2) - pnorm(stat1), 4);

    max <- max(3,abs(stat1)+1,abs(stat2)+1);

    xcord <- c(stat1,seq(stat1,stat2,.01),stat2);
    ycord <- c(0,dnorm(seq(stat1,stat2,.01)),0)

    if(pv == 0)
      ti <- substitute(paste("P(", z1, "< Z < ", z2, ") < .0001", sep=''), list(z1=stat1, z2=stat2))
    else if(pv == 1)
      ti <- substitute(paste("P(", z1, "< Z < ", z2, ") > .9999 ", sep=''), list(z1=stat1, z2=stat2))
    else
      ti <- substitute(paste("P(", z1, "< Z < ", z2, ") = ", p, sep=''), list(z1=stat1, z2=stat2, p=pv))
  },
  {
    # outside
    stat1 <- getvalidinput("First, enter the value for the lower bound a:" , "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8)
    if(stat1 == 'q') { cat("Quitting!\n"); break() }
    stat1 <- round(as.numeric(stat1), 4);

    cat("Next, enter a value for the upper bound b (must be greater than a).\n");
    bounds_ok <- FALSE;
    while(!bounds_ok)
    {
      stat2 <- getvalidinput("Upper bound: " , "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8)
      if(stat2 == 'q') { cat("Quitting!\n"); break() }
      stat2 <- round(as.numeric(stat2), 4);
      if(stat1 < stat2) { break() } else { cat(paste("Upper bound must be greater than lower bound of ", stat1, "!\n", sep='')) }
    }
    if(stat2 == 'q') { break() }

    pv <- round(1 - pnorm(stat2) + pnorm(stat1), 4);

    max <- max(3,abs(stat1)+1,abs(stat2)+1);

    xcord1 <- c(-max,seq(-max,stat1,.01),stat1);
    ycord1 <- c(0,dnorm(seq(-max,stat1,.01)),0)
    xcord2 <- c(stat2,seq(stat2,max,.01),max);
    ycord2 <- c(0,dnorm(seq(stat2,max,.01)),0)

    ti <- substitute(paste("P(Z < ", z1, " or Z >", z2, ") = ", p, sep=''), list(z1=stat1, z2=stat2, p=pv))
    if(pv == 0)
      ti <- substitute(paste("P(Z < ", z1, " or Z >", z2, ") < .0001", sep=''), list(z1=stat1, z2=stat2))
    else if(pv == 1)
      ti <- substitute(paste("P(Z < ", z1, " or Z >", z2, ") > .9999 ", sep=''), list(z1=stat1, z2=stat2))
    else
      ti <- substitute(paste("P(Z < ", z1, " or Z >", z2, ") =  ", p, sep=''), list(z1=stat1, z2=stat2, p=pv))
  })

    bound <- floor(max*10)/10; # bounds for axis, "floored" to nearest .1.
    plot(function(x) dnorm(x),-max,max,
         ylim=c(0,.4),
         axes=F,
         xlab="Z values",
         ylab='Density',
         main="N(0,1) Distribution",
         xaxs='i', #to make the axes start exactly at 0,0
         yaxs='i')
    axis(1,at=-ceiling(bound):ceiling(bound))
    switch(dir,
           {
             axis(1,at=stat,lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
             polygon(xcord,ycord,col='red')
           },
  {
    axis(1,at=stat,lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord,ycord,col='red')
  },
  {
    axis(1,at=c(stat1, stat2),lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord,ycord,col='red')
  },
  {
    axis(1,at=c(stat1, stat2),lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord1,ycord1,col='red')
    polygon(xcord2,ycord2,col='red')
  })
    axis(2,at=c(0,1))
    legend('topright', eval(ti), col='red', pch=15, bty='n', pt.cex=2)

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
        name <- paste("~/Desktop/probability_Rplot_", td, sep="")
      }
      if(os == 'windows')
      {
        name <- strsplit(as.character(Sys.getenv('HOME')), 'Documents')
        name <- paste(name, "Desktop\\probability_Rplot_", td, sep="")
      }
      name <- paste(name, ".jpg", sep="")
      jpeg(name)
      plot(function(x) dnorm(x),-max,max,
           ylim=c(0,.4),
           axes=F,
           xlab="Z values",
           ylab='Density',
           main=paste("N(0,1) Distribution",titlename, sep='\nBy '),
           xaxs='i', #to make the axes start exactly at 0,0
           yaxs='i')
      axis(1,at=-ceiling(bound):ceiling(bound))
      switch(dir,
             {
               axis(1,at=stat,lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
               polygon(xcord,ycord,col='red')
             },
  {
    axis(1,at=stat,lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord,ycord,col='red')
  },
  {
    axis(1,at=c(stat1, stat2),lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord,ycord,col='red')
  },
  {
    axis(1,at=c(stat1, stat2),lwd.ticks=2,col.ticks='red',col.axis='red', font=2,lwd=0, tck=-.04, padj=1.25);
    polygon(xcord1,ycord1,col='red')
    polygon(xcord2,ycord2,col='red')
  })
      axis(2,at=c(0,1))
      legend('topright', eval(ti), col='red', pch=15, bty='n', pt.cex=2)
      dev.off()
    }


  }
}
