#' Run the P-val Walkthrough script.
#'
#' @export
#' @examples
#' \dontrun{pval()}
pval<-function()
{
  ### P-value script for R.
  ### by Josh Errickson
  ### Original version by Tom Brown
  ### Version History:
  ### 2.4.2 - Minor spelling issues.
  ### 2.4.1 - Never display probabilies as 0 or 1... just say <.0001 or >.9999
  ### 2.4.0 - Added ability to save plots. Other misc changes mirroring prob().
  ### 2.3.0 - Helper function to handle user input - should no longer ever crash! Just keep requesting new input.
  ### 2.2.0 - Added looping action
  ### 2.1.3 - I can spell distribution, I swear....
  ### 2.1.2 - Fixed bug with negative test statistic & two-sided t/z tests
  ### 2.1.1 - Fixed rounding in p-val
  ### 2.1   - Minor output cleanups
  ###       - Added errors for invalid inputs
  ### 2.0   - Massive restructuring
  ###       - Added mean, chi-sq
  ### 1.0   - Tom Brown's original flavor

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

  options(warn=-1) #make this 0 for testing, -1 for production
  options(scipen=4) #Don't give numbers like 1e-4, instead do .0001

  cat("Are you ready for some nice graphics to help illustrate p-values?\n")
  cat("At each step, type your response and press \"Enter\"\n")
  while(TRUE)
  {
    dis <- cat("What is the distribution of the test statistic? Choose from the follwing:\n1: Z\n2: t\n3: Chi-square\n4: F\nYou can press 'q' at any input to quit.\n")
    # getvalidinput is a function that does exactly what it sounds like. See documentation above.
    dis <- getvalidinput("Enter your selection: ", "Invalid selection! Choose between 1, 2, 3 or 4!\n", validinput = c('1', '2', '3', '4', 'q'))
    if(dis == 'q') {cat("Quitting!\n"); break()}
    if(dis=='1') dis <- 'z'
    if(dis=='2') dis <- 't'
    if(dis=='3') dis <- 'c'
    if(dis=='4') dis <- 'f'

    #      dis <- as.numeric(dis)

    switch(dis,
           'z'=, 't'= {

             stat <- getvalidinput(paste("Enter the ",dis," test statistic: ", sep=''), "Test statistic must be a number!\n", validinput='q', low=-1e8, high=1e8);
             if(stat == 'q') { cat("Quitting!\n"); break() }
             stat <- round(as.numeric(stat), 4);


             # mm is the maximum (for max, assuming stat is reasonable)
             # pv is the left-tailed p-value
             # dl is the density function
             switch(dis,
                    t = {
                      df <- getvalidinput('Enter the degrees of freedom: ', 'Degrees of Freedom must be 1 or greater!\n', validinput='q', low=1, high=1e8)
                      if(df == 'q') { cat("Quitting!\n"); break() }
                      df <- as.numeric(df)
                      mm <- qt(.99,df);
                      pv <- pt(stat,df);
                      dl <- function(x) { dt(x,df) };
                    },
                    z = {
                      mm <- 3;
                      pv <- pnorm(stat);
                      dl <- function(x) { dnorm(x) };
                    }
             )

             cat("What is the direction of the alternative hypothesis? Choose from the following:\n1: greater than\n2: less than\n3: two-sided\n")
             tail <- getvalidinput("Enter your selection: ", "Invalid selection! Choose between 1, 2, or 3 !\n", validinput = c('1', '2', '3', 'q'))
             if(tail == 'q') {cat("Quitting!\n"); break()}
             tail <- as.numeric(tail)

             max <- max(mm,abs(stat)+1)
             # Rightmost X value. Either give me 3, or more if the stat is very large
             bound <- floor(max*10)/10 # bounds for axis, "floored" to nearest .1.

             switch(tail,
                    { # upper tail
                      pv <- 1-pv;
                      xcord <- c(stat,seq(stat,max,.01),max);
                      ycord <- c(0,dl(seq(stat,max,.01)),0)
                    },
{ # lower tail
  xcord <- c(-max,seq(-max,stat,.01),stat);
  ycord <- c(0,dl(seq(-max,stat,.01)),0)
},
{ # two-tailed
  if(stat < 0) { pv <- 2*pv       # correct calc for neg & pos test statistics
  } else pv <- 2*(1-pv)
  stat <- abs(stat)
  xcord <- c(stat,seq(stat,max,.01),max,-max,seq(-max,-stat,.01),-stat);
  ycord <- c(0,dl(seq(stat,max,.01)),0,0,dl(seq(-max,-stat,.01)),0)
}
             )
             pv <- round(pv,4)
             if(pv == 0)
               pvtxt <- "< .0001"
             else if(pv == 1)
               pvtxt <- "> .9999"
             else
               pvtxt <- paste("=", pv)

             switch(dis,
                    t = ti <- paste('t(', df, ') Distribution', sep=''),
                    z = ti <-'N(0,1) Distribution'
             )

             plot(function(x) dl(x),-max,max,
                  ylim=c(0,.4),
                  axes=F,
                  xlab=substitute(paste(dis,"values"),list(dis=dis)),
                  ylab='Density',
                  main=ti,
                  xaxs='i', #to make the axes start exactly at 0,0
                  yaxs='i')
             polygon(xcord,ycord,col='red')
             axis(1,at=-ceiling(bound):ceiling(bound))
             axis(1,at=round(stat,4),lwd.ticks=2,col.ticks='red',col.axis='red',
                  font=2,lwd=0, tck=-.04, padj=1.25)
             if(tail=='t') axis(1,at=-round(stat,4),lwd.ticks=2,col.ticks='red',
                                col.axis='red',font=2,lwd=0, tck=-.04, padj=1.25)
             axis(2,at=c(0,1))
             segments(0,0,0,dl(0),lty=3)  # plot mean
             switch(dis,
                    t = legend('topright',c(paste('p-value',pvtxt),
                                            expression(paste("E(T) under ",H[0]))),
                               pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2),
                    z = legend('topright',c(paste('p-value',pvtxt),
                                            expression(paste("E(Z) under ",H[0]))),
                               pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2)
             )

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
               plot(function(x) dl(x),-max,max,
                    ylim=c(0,.4),
                    axes=F,
                    xlab=substitute(paste(dis," values"),list(dis=dis)),
                    ylab='Density',
                    main='',
                    xaxs='i', #to make the axes start exactly at 0,0
                    yaxs='i')
               title(ti, line=3)
               title(paste("By", titlename), line=1)
               polygon(xcord,ycord,col='red')
               axis(1,at=-ceiling(bound):ceiling(bound))
               axis(1,at=round(stat,4),lwd.ticks=2,col.ticks='red',col.axis='red',
                    font=2,lwd=0, tck=-.04, padj=1.25)
               if(tail=='t') axis(1,at=-round(stat,4),lwd.ticks=2,col.ticks='red',
                                  col.axis='red',font=2,lwd=0, tck=-.04, padj=1.25)
               axis(2,at=c(0,1))
               segments(0,0,0,dl(0),lty=3)  # plot mean
               switch(dis,
                      t = legend('topright',c(paste('p-value',pvtxt),
                                              expression(paste("E(T) under ",H[0]))),
                                 pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2),
                      z = legend('topright',c(paste('p-value',pvtxt),
                                              expression(paste("E(Z) under ",H[0]))),
                                 pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2)
               )
               dev.off()
             }
           },
'c'=,'f'={
  # mm is the maximum (for max, assuming stat is reasonable)
  # pv is the left-tailed p-value
  # dl is the density function
  # xl & ti are xlabel and title, respectively
  # min is to chop off the lower tail of chi-sq with high df

  switch(dis,
         c = {
           stat <- getvalidinput("Enter the Chi-squared test statistic: ", "Chi-squared statistic  must be a non-negative number!\n", validinput='q', low=0, high=1e8);
           if(stat == 'q') { cat("Quitting!\n"); break() }
           stat <- round(as.numeric(stat), 4);
           df <- getvalidinput("Enter the degrees of freedom: ", "Degrees of freedom must be at least 1!\n", validinput='q', low=1, high=1e8);
           if(df == 'q') { cat("Quitting!\n"); break() }
           df <- round(as.numeric(df), 4);
           mm <- qchisq(.999,df);
           min <- min(floor(qchisq(.001,df)),abs(stat)-1);
           pv <- 1-pchisq(stat,df);
           dl <- function(x) { dchisq(x,df) };
           xl <- expression(paste(chi^2,' values', sep=''));
           ti <- substitute(paste(chi^2,"(",a,") Distribution"), list(a=df))
         },
         f = {
           stat <- getvalidinput("Enter the F test statistic: ", "F statistic  must be a non-negative number!\n", validinput='q', low=0, high=1e8);
           if(stat == 'q') { cat("Quitting!\n"); break() }
           stat <- round(as.numeric(stat), 4);
           df1 <- getvalidinput("Enter the first degree of freedom: ", "Degrees of freedom must be at least 1!\n", validinput='q', low=1, high=1e8);
           if(df1 == 'q') { cat("Quitting!\n"); break() }
           df1 <- round(as.numeric(df1), 4);
           df2 <- getvalidinput("Enter the first degree of freedom: ", "Degrees of freedom must be at least 1!\n", validinput='q', low=1, high=1e8);
           if(df2 == 'q') { cat("Quitting!\n"); break() }
           df2 <- round(as.numeric(df2), 4);
           mm <- qf(.99,df1,df2);
           min <- 0;
           pv <- 1-pf(stat,df1,df2);
           dl <- function(x) { df(x,df1,df2) };
           xl <- 'F values';
           ti <- paste("F(",df1,",",df2,") Distribution", sep='')
         }
  )
  pv <- round(pv,4)
  if(pv == 0)
    pvtxt <- "< .0001"
  else if(pv == 1)
    pvtxt <- "> .9999"
  else
    pvtxt <- paste("=", pv)
  max <- max(mm,stat+1)
  # Rightmost X value. Either give me 99.9% of the plot, or more if the stat is very large

  plot(function(x) dl(x),min,max,
       axes=F,
       xlab=xl,
       ylab='Density',
       main=ti,
       xaxs='i', #to make the axes start exactly at 0,0
       yaxs='i')
  xcord <- c(stat,seq(stat,max,.01),max)
  ycord <- c(0,dl(seq(stat,max,.01)),0)
  polygon(xcord,ycord,col='red')

  axis(1,at=0:ceiling(max))
  axis(1,at=round(stat,4), lwd.ticks=2, col.ticks='red', col.axis='red',
       font=2, lwd=0, tck=-.04, padj=1.25)
  axis(2,at=c(0,1))
  switch(dis, # plot mean
         c = segments(df,0,df,dl(df),lty=3),
         f = { m <- df2/(df2+2); segments(m,0,m,dl(m),lty=3) }
  )
  switch(dis,
         c = legend('topright',c(paste('p-value',substitute(p,list(p=pvtxt))),
                                 expression(paste("E(",chi^2,") under ",H[0]))),
                    pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2),
         f = legend('topright',c(paste('p-value',substitute(p,list(p=pvtxt))),
                                 expression(paste("E(F) under ",H[0]))),
                    pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2)
  )

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
    plot(function(x) dl(x),min,max,
         axes=F,
         xlab=xl,
         ylab='Density',
         main="",
         xaxs='i', #to make the axes start exactly at 0,0
         yaxs='i')
    title(ti, line=3)
    title(paste("By", titlename), line=1)
    xcord <- c(stat,seq(stat,max,.01),max)
    ycord <- c(0,dl(seq(stat,max,.01)),0)
    polygon(xcord,ycord,col='red')

    axis(1,at=0:ceiling(max))
    axis(1,at=round(stat,4), lwd.ticks=2, col.ticks='red', col.axis='red',
         font=2, lwd=0, tck=-.04, padj=1.25)
    axis(2,at=c(0,1))
    switch(dis, # plot mean
           c = segments(df,0,df,dl(df),lty=3),
           f = { m <- df2/(df2+2); segments(m,0,m,dl(m),lty=3) }
    )
    switch(dis,
           c = legend('topright',c(paste('p-value',substitute(p,list(p=pvtxt))),
                                   expression(paste("E(",chi^2,") under ",H[0]))),
                      pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2),
           f = legend('topright',c(paste('p-value',substitute(p,list(p=pvtxt))),
                                   expression(paste("E(F) under ",H[0]))),
                      pch=c(15,-1),lty=c(0,3),col=c('red','black'),bty='n',pt.cex=2)
    )
    dev.off()
  }

}
    )
  }
}
