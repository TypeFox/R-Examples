Help <- 
function(topic=NULL) {


  # set up plot window
  set.up.plot <- function(nlines) {
    par(mar=c(.5,.5,.5,.5), bg=rgb(255,253,250,maxColorValue=255),
                            fg=rgb(20,15,15,maxColorValue=255), cex=.72)
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,100))
    if (!missing(nlines)) {
      ybot <- 100 - ((nlines+1) * 4)
      rect(-1, ybot, 95, 96, col=col.rect, lwd=.75, border=col.line)
    }
  }


  help.more <- function(fname, yline) {
    h1 <- "Complete list of Help topics, enter:  Help()"
    h2 <- paste("For more help on a function, enter ? in front of its name:  ?", fname, sep="")

    if (getOption("colors") != "gray") 
      col.sep <- "lightsteelblue"
    else
      col.sep <- "gray50"
    lines(c(5,90), c(yline,yline), col=col.sep)
    text(0,yline-5, label=h1, adj=0)
    text(0,yline-10, label=h2, adj=0)
  }



  # ------------------------------------------------
  # ------------------------------------------------

  if (missing(topic))
    topic <- NULL
  else {
    topic <- deparse(substitute(topic))  # make a char string if not
    topic <- gsub("\"", "", topic)  # if already a char string, remove quotes
  }

  # convert topic to all lowercase letters
  if (!is.null(topic))
    for (i in 1:nchar(topic)) {
      xc <- substr(topic,i,i)
      if (xc %in% LETTERS) substr(topic,i,i) <- letters[which(LETTERS==xc)] 
    }

  if (getOption("colors") != "gray") {
    col.rect <- getOption("col.fill.pt")
    col.line <- "lightsteelblue"
  }
  else {
    col.rect <- "gray90"
    col.line <- "gray30"
  }

  # set graphic parameters
  if (!is.null(topic)) {
    if (topic != "lessr") {
      old.par <- par("mar", "cex", "bg", "fg")
      on.exit(par(old.par))
    }
  }
  else {  # topic is null, Help()
    old.par <- par("mar", "cex", "bg", "fg")
    on.exit(par(old.par))
  }


  if (is.null(topic)) {

  t0 <- "Help Topics for lessR"

  fsys <- bquote(paste(bold("Help(theme)"), "  System level settings, such as a color theme for graphics"))
  fcsv <- bquote(paste(bold("Help(data)"), "  Create a data file from Excel or similar application"))
  frw <- bquote(paste(bold("Help(Read)"), " and ", bold("Help(Write)"), "  Read or write data to or from a file"))
  flib <- bquote(paste(bold("Help(library)"), "  Access libraries of functions called packages"))
  ftrans <- bquote(paste(bold("Help(edit)"), "  Edit data and create new variables from existing variables"))

  fhist <- bquote(paste(bold("Help(Histogram)"), "  Histogram, box plot, dot plot, density curve"))
  fbar <- bquote(paste(bold("Help(BarChart)"), "  Bar chart, pie chart"))
  fline <- bquote(paste(bold("Help(LineChart)"), "  Line chart, such as a run chart or time series chart"))
  fplot <- bquote(paste(bold("Help(ScatterPlot)"), "  Scatterplot for one or two variables, a function plot"))

  fstat <- bquote(paste(bold("Help(SummaryStats)"), "  Summary statistics for one or two variables"))
  fone <- bquote(paste(bold("Help(one.sample)"), "  Analysis of a single sample of data"))
  fmean <- bquote(paste(bold("Help(ttest)"), "  Compare two groups by their mean difference"))
  faov <- bquote(paste(bold("Help(ANOVA)"), "  Compare mean differences for many groups"))
  fpwr <- bquote(paste(bold("Help(power)"), "  Power analysis for the t-test"))
  fcor <- bquote(paste(bold("Help(Correlation)"), "  Correlation analysis"))
  freg <- bquote(paste(bold("Help(Regression)"), " and ", bold("Help(Logit)"), " Regression analysis, logit analysis"))
  ffac <- bquote(paste(bold("Help(factor.analysis)"), "  Confirmatory and exploratory factor analysis"))

  fprob <- bquote(paste(bold("Help(prob)"), "  Probabilities for normal and t-distributions"))
  frnsm <- bquote(paste(bold("Help(random)"), " and ", bold("Help(sample)"), "  Create random numbers or samples"))

  fpdf <- bquote(paste(bold("Help(help.to.pdf)"), "  Obtain a printable pdf of all of the contents"))
  fpck <- bquote(paste(bold("Help(lessR)"), "  lessR manual and list of updates to current version"))

  set.up.plot() 
  pos1 <- 93; pos2 <- 69; pos3 <- 49; pos4 <- 14; pos5 <- 7
  text(50,100, label=t0, font=4)
  text(0,pos1, label=fsys, adj=0)
  text(0,pos1-4, label=fcsv, adj=0)
  text(0,pos1-8, label=frw, adj=0)
  text(0,pos1-12, label=flib, adj=0)
  text(0,pos1-16, label=ftrans, adj=0)
  lines(c(5,90), c(74,74), col=col.line)
  text(0,pos2, label=fhist, adj=0)
  text(0,pos2-4, label=fbar, adj=0)
  text(0,pos2-8, label=fline, adj=0)
  text(0,pos2-12, label=fplot, adj=0)
  lines(c(5,90), c(53,53), col=col.line)
  text(0,pos3, label=fstat, adj=0)
  text(0,pos3-4, label=fone, adj=0)
  text(0,pos3-8, label=fmean, adj=0)
  text(0,pos3-12, label=faov, adj=0)
  text(0,pos3-16, label=fpwr, adj=0)
  text(0,pos3-20, label=fcor, adj=0)
  text(0,pos3-24, label=freg, adj=0)
  text(0,pos3-28, label=ffac, adj=0)
  lines(c(5,90), c(18,18), col=col.line)
  text(0,pos4, label=fprob, adj=0)
  text(0,pos4-4, label=frnsm, adj=0)
  lines(c(5,90), c(7,7), col=col.line)
  #text(0,pos5, label=fagain, adj=0)
  #text(0,pos5-4, label=fpdf, adj=0)
  text(0,pos5-4, label=fpck, adj=0)

  }


  else if (topic %in% c("data", "file", "csv", "sav", "rda")) {
  t0 <- "Data Files"

  t1 <- "
  A data file organizes the data into a table, with variables in the
  columns and the data for a single person, company, etc. in a row.
  By default, include the name of each variable in the first row.
  After the first row, only data values are included. 

  The lessR function Read can read data files in many file formats,
  including MS Excel. The most generic format is the csv format, for
  \"comma separated values\". A csv file is a text file with commas
  that separate adjacent values in each row. Usually the variable
  names are in the first row and each remaining row contains the data
  for one case, such as one person or one company, etc. Each column
  contains the data for the corresponding variable.

  To create the csv file from Excel, do a Save As and choose the csv
  format. With the free, open source LibreOffice Calc, click the arrow
  in the left margin towards the bottom labeled File type. From the
  available options, choose Text CSV. Then Save button. 

  A fixed width format text data file is where each column of data
  values is assigned a specific width, often with no spaces between
  the data values. R can write a data file [see Help(Write)], what is
  called an native R data file with a default file type of .rda. Data files
  written by the SPSS system have the default file type of .sav. Both
  of these files can be read into R, as well as Excel files directly.
  "

  set.up.plot(0)
  text(50,100, label=t0, font=4)
  text(0,52, label=t1, adj=0)

  help.more("Read", 8)
  }


  else if (topic %in% c("rd", "read")) {
  t0 <- "Read Data into R and Prepare for Analysis"

  f1 <- bquote(paste(bold("Read, rd"), "  Read a data file into an R data frame for analysis"))

  t1 <- "
  Browse for a csv, tab-delimited, Excel, R, SAS, or SPSS data file on 
  your file system and read the information into the specified data table
  with Read, or its abbreviation, rd. To browse, use an empty ().
      > mydata <- Read()
  The  <-, the assignment operator, instructs R to assign what was read
  to the data table, here named mydata. This is the default name that
  the lessR data analysis functions assume will be analyzed. Read
  variable labels with the labels option to specify the file of labels.

  Or, specify the full data path name in quotes, such as from the web.
      > mydata <- Read(\"http://lessrstats.com/data/twogroup.csv\")
  To read a text file with a comma for a decimal point, use Read2().

  If you wish to view more information about the data table, do
      > details()
  If the file is not read into mydata, include the name: details(name).

  To read a text file in which each column of data values assigned a
  specific width, add the widths option that specifies the width of each
  column according to the order of the variables.  Enclose the list
  with the c function, here read two variables with widths of 4 and 1.
      > mydata <- Read(widths=c(4,1), col.names=c(\"ID\", \"Gender\"))
  "

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(89,89), col=col.line)
  text(0,49, label=t1, adj=0)

  help.more("Read", 9)
  }


  else if (topic %in% c("wrt", "write")) {
  t0 <- "Write Contents of Data Frame mydata into a Data File"

  f1 <- bquote(paste(bold("Write, wrt"), "  Write a data file called mydata into an R data frame"))

  t1 <- "
  The name of the entire rectangular table of data, called a data frame in R, can 
  be named mydata within R.  This is the default name of the data table assumed
  by the lessR data analysis functions.

  Here is how to write the contents of mydata to a csv data file with the name of 
  mydata.csv.
      > Write()
  The file is written to the default working directory.  The Write function displays
  this location after writing the file.

  Or, explicitly specify the file name.
      > Write(\"mybestdata\")
  The file type of .csv is automatically appended to the file name.

  To write a data file in native R format, use the type=\"R\" option, or the 
  abbreviation for the function name  wrt.r.
      > wrt.r(\"mybestdata\")

  The lessR Write function relies upon the R function write.table, which is
  is quite general, with many options.  For more information, enter ?write.table."

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(89,89), col=col.line)
  text(0,52, label=t1, adj=0)

  help.more("Write", 11)
  }


  else if (topic %in% c("library", "package", "install", "update"))  {
  t0 <- "Contributed Packages"

  f1 <- bquote(paste(bold("install.packages"), "  Download a contributed package"))
  f2 <- bquote(paste(bold("library"), "  Load an installed package from the library into R for access"))
  f3 <- bquote(paste(bold(update.packages), "  Update contributed packages to current versions"))

  t1 <- "
  R works with functions and each function is contained in a specific
  package. Some packages, such as graphics, are included with the
  default installation of R, and are pre-loaded each time R is run.
  Other packages must be explicitly downloaded from the R servers. 

  The example here is for the contributed package lessR. Install one
  time only for a specific computer, with quotes.
      > install.packages(\"lessR\")

  Each time the R application is started, including after the install,
  load the package from the library, without using quotes.
      > library(lessR)
  To see the description of the package and a list of its functions,
      > library(help=lessR)
  To see a list of all installed packages in the library, 
      > library()

  To access updated versions of all installed packages, 
      > update.packages()"

  set.up.plot(3)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  #lines(c(5,90), c(80,80), col=col.line)
  text(0,48, label=t1, adj=0)

  help.more("install.packages", 10)
  }


  else if (topic %in% c("edit", "transform", "trans", "factor", "recode", "rec", 
                        "subset", "subs"))  {
  t0 <- "Edit Data"

  f1 <- bquote(paste(bold("fix"), "  Use a graphical interface to edit data values, add or delete variables"))
  f2 <- bquote(paste(bold("Transform"), "  Transform the values of a variable with a formula"))
  f3 <- bquote(paste(bold("Recode"), "  Recode the values of a variable by specifying the new values"))
  f4 <- bquote(paste(bold("factor"), "  Explicitly define the values of a categorical variable"))
  f5 <- bquote(paste(bold("Subset"), "  Extract a subset of data, variables (columns) and/or rows"))

  t1 <- "
  R function fix provides a graphical/mouse interface for editing data.
      > fix(mydata)
  lessR function Transform creates a new variable or rewrites over existing.
      > mydata <- Transform(SalaryDiv=Salary/1000)
  lessR function Recode changes individual values. 
      > mydata <- Recode(Scores, old=c(1:4), new=c(10,15,20,25))
  R function factor creates a new variable with non-numeric categories.
  Severity was encoded with a 1 for Mild, 2 for Moderate and 3 for Severe.
      > mydata <- Transform(ordered=TRUE, Severity.f= 
               factor(Severity, levels=c(1,2,3), labels=c(\"Mild\", \"Mod\", \"Severe\")))
  Here the values of the new variable are also ordered, from Mild to Severe. 
  Extract subsets of data from a data frame with the lessR Subset function.
      > mydata <- Subset(rows=Gender==\"M\", columns=c(Years, Salary))
  The data frame, mydata, now consists only of data for Males limited to
  the variables Years and Salary. To display a subset, drop the mydata <-.
  "

  set.up.plot(5)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  text(0,82, label=f4, adj=0)
  text(0,78, label=f5, adj=0)
  #lines(c(5,90), c(78,78), col=col.line)
  text(0,44, label=t1, adj=0)

  help.more("Recode", 11)
  }


  else if (topic %in% c("system", "set", "theme")) {
  t0 <- "Global Settings"

  f1 <- bquote(paste(bold("theme"), "  all lessR system settings such as a color theme"))
  f2 <- bquote(paste(bold("showColors"), "  lessR function to illustrate all color names"))

  t1 <- "
  The lessR function theme provides global settings for lessR functions.
  Set the color theme for the graphics functions. The default color theme is
  \"dodgerblue\", with possibilities of \"gray\", \"green\", \"gold\", \"rose\", \"red\", 
  \"purple\", \"sienna\", \"white\", \"orange.black\" and \"gray.black\".  Example:
      > theme(colors=\"gray\") 
  Set the transparency level of bars and plotted points with the
  trans.fill.bar and trans.fill.pt options. Turn off grid lines with col.grid=\"off\".
  Change or remove background color with col.bg, such as col.bg=\"off\".

  Levels of a categorical variable may be encoded with numerical digits,
  such as 0 for Male and 1 for Female. R is obliged to interpret numerical
  variables as numeric.  One option is to redefine these variables as
  factors [see Help(edit)]. Or set the value of the lessR option n.cat.
      > theme(n.cat=3)
  Here any numerical variable with just 3 unique, equally spaced interval 
  values or less is interpreted as a categorical variable. The default
  value of n.cat is 8, and applies to ScatterPlot and SummaryStats.

  To see all available theme options, enter the following.
      > theme(show=TRUE)
  To see all the R named colors, enter the following.
      > showColors()
  "

  set.up.plot(2)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  #lines(c(5,90), c(80,80), col=col.line)
  text(0,47, label=t1, adj=0)

  help.more("theme", 7)
  }


  else if (topic %in% c("histogram", "hs", "hst", "hist", "boxplot", "box", "bx", 
    "dotplot", "dp", "dot", "density", "dn", "dens", "distribution", "dist",
    "univariate")) {
  t0 <- "Histogram, etc."

  f1 <- bquote(paste(bold("Histogram, hs"), "  Histogram"))
  f2 <- bquote(paste(bold("Density, dn"), "  Density curve over histogram"))
  f3 <- bquote(paste(bold("BoxPlot, bx"), "  Box plot"))
  f4 <- bquote(paste(bold("ScatterPlot, sp"), "  Scatter plot of 1 variable"))

  t1 <- "
  Plot a distribution of data values for a continuous variable, here for
  variable Y with the current color theme. Use Histogram or hs.
      > Histogram(Y)
  Specify the gray scale color theme, a title, and a label for the x axis.
      > theme(colors=\"gray\")
      > Histogram(Y, main=\"My Title\", xlab=\"Y (mm)\")
  Specify bins, begin at 60 with a bin width of 10. Can also specify bin.end.
      > Histogram(Y, bin.start=60, bin.width=10)
  Save the output of the function for later analysis and view separately.
      > h <- hs(Y)
      > h

  Density curve superimposed on the underlying histogram, abbreviated
  dn, a BoxPlot or bx, and a one variable ScatterPlot, or sp.
      > Density(Y)   or   > BoxPlot(Y)   or   > ScatterPlot(Y)

  These functions, except sp, can also replace the variable name such as Y 
  with a list of multiple variables, such as c(Salary, Years) or Salary:Years,
  or an entire mydata data frame by passing no argument.
      > Histogram()
  "

  set.up.plot(4)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  text(0,82, label=f4, adj=0)
  #lines(c(5,90), c(78,78), col=col.line)
  text(0,44, label=t1, adj=0)

  help.more("Histogram",8)
  }


  else if (topic %in% c("barchart", "bc", "piechart", "pc", "pareto")) {
  t0 <- "BarChart, PieChart and Pareto Chart"

  f1 <- bquote(paste(bold("BarChart, bc"), "  Bar chart for one or more categorical variables"))
  f2 <- bquote(paste(bold("PieChart, pc"), "  Pie chart for a categorical variable"))
  f3 <- bquote(paste(bold("pareto.chart"), "  Produce a Pareto chart"))

  t1 <- "
  Default bar chart with lessR function BarChart, or bc, as well as the
  frequency table, for one or two variables, here Y and X.
      > BarChart(Y)
      > BarChart(Y, by=X)
      
  With lessR function PieChart or pc, generate a pie chart.
      > PieChart(Y)
      
  The pareto.chart function is part of the external library called gcc
  (see Help(\"libraries\"). This function is not from lessR, so the name
  of the variable must be preceded by the data frame name and a $. 
      > library(gcc)
      > Ycount <- table(mydata$Y)
      > pareto.chart(Ycount$freq)

  Can replace the variable name such as Y with a list of multiple variables,
  such as c(Salary, Years) or Salary:Years, or an entire data frame. The
  default data frame is mydata. Here do a bar chart of all non-numeric
  variables in mydata.
      > BarChart()
  "

  set.up.plot(3)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  #lines(c(5,90), c(81,81), col=col.line)
  text(0,46, label=t1, adj=0)

  help.more("BarChart", 9)
  }


  else if (topic %in% c("linechart", "lc")) {
  t0 <- "Line Chart"

  f1 <- bquote(paste(bold("LineChart, lc"), "  A line chart, such as a run chart or time series chart"))

  t1 <- "
  The lessR function LineChart, or lc, generates a line chart with values
  ordered along some dimension such as time. If the data do not have a 
  pronounced trend, a centerline is automatically provided.
      > LineChart(Y)
  Also provided is a list of all the runs in the data.

  The line chart becomes a time series chart with times/dates on the
  horizontal axis.  Use the time.start and time.by options.
      > LineChart(Y, time.start=\"2005/09/01\", time.by=\"month\")
  Additional options are explained in the R help files for functions par,
  title, points and lines. 

  Color themes are available with the colors option, which can be invoked
  from a specific call to LineChart or system wide for all graphics output with
  the function theme. Here all subsequent graphics output is in gray scale.
      > theme(colors=\"gray\")
      > LineChart(Y)

  Can replace the variable name with a list of multiple variables, such as
  c(Salary, Years) or Salary:Years, or an entire data frame. The default data
  frame is mydata. Here do a line chart of all numerical variables in mydata.
      > LineChart()
  "

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(90,90), col=col.line)
  text(0,48, label=t1, adj=0)

  help.more("LineChart", 8)
  }


  else if  (topic %in% c("scatterplot", "sp", "plot", "scatter")) {
  t0 <- "Scatterplot"

  f1 <- bquote(paste(bold("ScatterPlot, sp"), "  A scatterplot for one or two variables"))

  t1 <- "
  ScatterPlot, or sp, generates a scatter plot for one or two variables with
  the current color theme,  with an optional 0.95 data ellipse. The points
  have a default transparency, which can be set from transparent to oblique.
      > ScatterPlot(X, Y, ellipse=TRUE)
  If the values of X are sorted, a function plot results so that the points
  are not individually displayed and are connected by line segments. If the
  number of response values is less than 10, a bubble plot is produced.

  Here generate a one-dimensional scatterplot, that is, a dot plot. 
      > ScatterPlot(Y)

  ScatterPlot can also provide for plotting two variables with different
  symbols and/or colors for each level of a third variable.
      > ScatterPlot(X, Y, by=Z)
  to better display multiple replications of the same point.

  The colors option specifies color themes, from a call to ScatterPlot
  or system wide for all graphics output with the function: theme. Here all
  subsequent graphics are with the sienna color theme, no transparency.
      > theme(colors=\"sienna\", trans.fill.pt=0)
      > ScatterPlot(X, Y)"

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(90,90), col=col.line)
  text(0,50, label=t1, adj=0)

  help.more("Plot", 8)
  }


  else if  (topic %in% c("summarystats", "ss", "standard score", "z-score")) {
  t0 <- "Summary Statistics"

  f1 <- bquote(paste(bold("SummaryStats, ss"), "  Summarize the values of a variable"))

  t1 <- "
  Summarize the variable Y with lessR SummaryStats, or just ss.  If numerical, 
  sample size, number of  missing data values, mean, sd, skew, kurtosis,
  minimum, maximum, quartiles and interquartile range are provided. If
  categorical, includes cell counts and proportions, plus the chi-square test.
      > SummaryStats(Y)
  A version for abbreviated output also exists, here saving the output for
  analysis of different pieces, both text output and statistics.
      > s <- ss.brief(Y)
      > s

  Summarize all variables in the data frame mydata.
      > SummaryStats()
      
  For a numerical variable Y, provide an optional grouping variable, X, to
  summarize at each level of the grouping variable. Or, if Y is categorical,
  a cross- tabulation table is generated.
      > SummaryStats(Y, by=X)

  SummaryStats also works with two categorical variables, here X and Y.
      > ss(X, Y)"

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(85,85), col=col.line)
  text(0,50, label=t1, adj=0)
  help.more("SummaryStats", 11)
  }


  else if (topic %in% c("one.sample", "one sample", "proportion", "prop")) {
  t0 <- "Inference for a Single Variable"

  f1 <- bquote(paste(bold("ttest, tt"), "  Inference for a mean"))
  f2 <- bquote(paste(bold("binom.test"), "  Inference for a proportion from exact binomial probability"))
  f3 <- bquote(paste(bold("prop.test"), "  Inference for a proportion from approximate normal probability"))

  t1 <- "
  These inference tests analyze the mean of a numeric variable or the
  proportion of a value of a categorical variable with a hypothesis
  test and confidence interval.

  This example uses the lessR function ttest, or tt, to evaluate a variable
  named Y and a null hypothesis of mu=100. The brief version is tt.brief.
      > ttest(Y, mu0=100)

  This example uses ttest to do the analysis from the sample statistics.
      > ttest(n=20, m=47.2, s=8.5, mu0=50)
      
  Here test for a fair coin after getting 53 out of 100 Heads. The R function
  binom.test is based on the exact binomial distribution.  The R prop.test
  function returns a chi-square value based on the normal approximation
  of the binomial.
      > binom.test(53,100, p=.5)
      > prop.test(53,100, p=.5)

  The prop.test function can be specified with or without the Yate's
  correction for continuity factor. The default is to include the correction."

  set.up.plot(3)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  #lines(c(5,90), c(82,82), col=col.line)
  text(0,48, label=t1, adj=0)

  help.more("ttest", 9)
  }


  else if (topic %in% c("ttest", "t-test", "tt")) {
  t0 <- "Compare Two Group Means"

  f1 <- bquote(paste(bold("ttest, tt"), "  An enhanced version of t.test to compare two group means"))
  f2 <- bquote(paste(bold("Model, model"), "  The t-test if Y is numerical and X has two values"))

  t1 <- "
  When responses to a variable are organized into two or more groups,
  compare the group means with a t-test.  For example, suppose the
  response variable is Salary and the grouping variable is Gender, with
  two values, M and F.

  Here the numerical response variable is named Y and the grouping
  variable, also called a factor, is named X, with exactly two values.
      >  ttest(Y ~ X)
  The tilde, ~, expresses the relationship between two or more variables.
  R refers to this expression as a formula, read as: Y is described by X.

  Sometimes the data for a t-test are arranged so that the responses for 
  each group, Y, already are in separate columns called vectors. Here
  calculate the t-test directly from two vectors called Y1 and Y2.
      > ttest(Y1, Y2)
  Add the paired=TRUE option to specify a dependent groups analysis.

  Or, do the analysis directly from summary statistics, the sample size
  (n), sample mean (m) and sample standard deviation (s). Ynm is the
  name of the response variable.
      > ttest(n=34, m=8.92, s=1.67, Ynm=\"Time\")" 
  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #text(0,90, label=f2, adj=0)
  #lines(c(5,90), c(86,86), col=col.line)
  text(0,52, label=t1, adj=0)

  help.more("ttest", 12)
  }

  else if (topic %in% c("anova", "av")) {
  t0 <- "Compare Means of Two or More Groups"

  f1 <- bquote(paste(bold("ANOVA, av"), "  Analysis of variance to compare two or more group means"))
  t1 <- "
  When responses to a variable are organized into exactly two groups,
  either the t-test function, ttest, or the lessR analysis of variance
  function, ANOVA, or simply av, can compare the group means. With
  more than two groups, ANOVA is required. Here the numerical response
  variable is named Y and the grouping variable, or factor, is X. 
      > ANOVA(Y ~ X)
  This is called one-way ANOVA because there is only a single factor, X.

  If the ANOaVA with more than two levels is significant, then a post-hoc
  examination of the mean differences with a controlled error rate will help
  uncover where the differences occurred. The ANOVA function relies
  upon the Tukey HSD procedure, providing both text and graphics output.

  For a randomized block ANOVA invoke a blocking factor with a + .
      > ANOVA(Y ~ X + Blck)

  For a two way between groups ANOVA, specify two factors with a * .
      > ANOVA(Y ~ X1 * X2)
  "

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(85,85), col=col.line)
  text(0,53, label=t1, adj=0)

  help.more("av", 15)
  }



  else if (topic == "power") {
  t0 <- "Power"

  f1 <- bquote(paste(bold("ttestPower, ttp"), "  Power analysis of the t-test"))

  t1 <- "
  The lessR function, ttestPower, uses the standard R function, power.t.test, to 
  calculate a range of power values and automatically provide a power curve. 

  To obtain a power curve with power.t.test requires setting up the range of
  alternative mean or mean difference values, usually by trial and error, 
  invoking ttestPower, saving the results, and then invoking the plot function,
  including the labeling of each axis. Then to analyze related results such 
  as power at a different sample size, the ttestPower function must be run
   several more times. 

  The enhanced function, ttestPower, does all of this automatically for one 
  or two sample t-tests, and also plots the power curve in color. This example is 
  for the default power curve for a sample size of 20 in each group and 
  a within-group or pooled standard deviation of 5.
      > ttestPower(n=20, s=5)
  Related analysis is also provided."

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5lessR ,90), c(89,89), col=col.line)
  text(0,59, label=t1, adj=0)

  help.more("ttp", 25)
  }


  else if (topic %in% c("correlation", "cr", "cor", "corr")) {
  t0 <- "Correlation and Related Graphics"

  f1 <- bquote(paste(bold("Correlation, cr"), "  Correlations between two or more variables"))
  f3 <- bquote(paste(bold("ScatterPlot, sp"), "  Graphics, generate a scatterplot for two or more variables"))

  t1 <- "
  The lessR function Correlation, or cr, can compute a correlation for two
  variables. Or for a data frame, mydata by default, the correlation matrix
  is computed, with pairwise deletion of missing data by default. A heat map
  and scatter plot matrix can also be generated with  graphics=TRUE. The
  matrix is displayed and also is stored as mycor such as for a subsequent
  factor analysis. Set the method option to \"spearman\" or \"kendall\" to get
  these correlations.

  The lessR function, ScatterPlot, or just sp, displays a scatterplot for two
  variables or a scatterplot matrix for a data frame. The corresponding
  correlation or correlation matrix is also displayed. See Help(ScatterPlot)
  for more information.
      > mycor <- Correlation(X,Y)
  The brief form for the correlation analysis for two variables also exists.
      > mycor <- cr.brief(X,Y)

  Or, analyze many correlations at once, such as for Y, X1, X2 and X3 in 
  the data frame called mydata.
      > mycor <- Correlation(c(Y,X1:X3))

  "

  set.up.plot(2)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f3, adj=0)
  #lines(c(5,90), c(83,83), col=col.line)
  text(0,47, label=t1, adj=0)

  help.more("Correlation", 12)
  }


  else if (topic %in% c("regression", "reg")) {
  t0 <- "Regression Analysis"

  f1 <- bquote(paste(bold("Regression, reg, reg.brief"), "  Regression analysis"))
  f2 <- bquote(paste(bold("Model, model"), "  Regression analysis if the variables are numerical"))

  t1 <- "
  The full text output of Regression, or reg, is comprehensive. Here
  specify a multiple regression of response Y and two predictors,
  X1 and X2.
       > Regression(Y ~ X1 + X2)
  Obtain abbreviated output with brief=TRUE, or use the abbreviation,
       > reg.brief(Y ~ X1 + X2)

  Can save the output for later analysis and viewing, such as with knitr,
  and create a text file of R Markdown instructions with file type .Rmd.
       > r <- reg(Y ~ X1 + X2, Rmd=\"reg_out\")
       > r                       # to see all the output
       > r$out_coefs     # to see this one segment of output
       > names(r)        # to see the names of the output segments
  Many types of output are contained in r, which consists of text output for
  viewing, statistics in numerical format, and also the markdown instructions
  to generate the corresponding html, pdf or Word document from RStudio.
  
  To obtain specified prediction intervals for new data, for example,
       > reg(Y ~ X1 + X2, X1.new=c(10,20), X2.new=c(100:110))
  X1.new, X2.new, etc. always specify the values of the predictor
  variables for the prediction intervals regardless of their names."

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,54, label=t1, adj=0)

  help.more("Regression", 10)

  }


  else if (topic %in% c("logit", "lr")) {
  t0 <- "Logit Regression Analysis"

  f1 <- bquote(paste(bold("Logit, lr"), "  Logit regression analysis"))
  f2 <- bquote(paste(bold("Model, model"), "  Logit analysis if a binary response variable"))

  t1 <- "
  Logit preforms a logit analysis. This example specifies a multiple
  regression model with a response variable named Y that has only two
  values, and two predictor variables, X1 and X2.
      > Logit(Y ~ X1 + X2)
  The standard R formula function specifies the model, which uses the
  tilde, ~, to mean 'depends on', and then the plus sign, +, to
  separate terms.

  The input values are not limited to 0 and 1. The output follows the
  general format of the Regression function, but also includes the
  classification table of correct and incorrect predictions.

  The abbreviated form of the function is lr, such as
       > lr(Y ~ X1 + X2)
  "

  set.up.plot(2)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  #lines(c(5,90), c(87,87), col=col.line)
  text(0,56, label=t1, adj=0)

  help.more("Logit", 30)
  }


  else if (topic %in% c("factor.analysis", "fa", "factors", "corcfa", "cfa", "corefa", "efa", "corScree", "scree")) {
  t0 <- "Confirmatory and Exploratory Factor Analysis"

  f1 <- bquote(paste(bold("corCFA, cfa"), "  Confirmatory factor analysis"))
  f2 <- bquote(paste(bold("corEFA, efa"), "  Exploratory factor analysis"))
  f3 <- bquote(paste(bold("corRead, rd.cor"), "  Read a correlation matrix"))
  f4 <- bquote(paste(bold("corScree, scree"), "  Scree plot of eigenvalues of the correlation matrix"))
  f5 <- bquote(paste(bold("corReorder, reord"), "  Reorder the variables in the correlation matrix"))

  t1 <- "
  Several lessR functions analyze data in the form of a correlation matrix,
  by default called: mycor.  Read mycor with corRead, often with the lessR
  function, to, which names a string of sequential variables (items).
      > corRead(names=to(\"m\",20))
  Here browse for the file that contains the matrix, name the 20 variables
  from m01, m02 to m20. Can also compute mycor with Correlation.

  The function corCFA, or cfa, does a confirmatory factor analysis of a
  multiple indicator measurement model. Use lavaan notation, or specify each
  group of items by listing the factor name, Fn, and the corresponding 
  variable names, with items separated by commas or a colon.
      > cfa(F1=V1:V3, F2=V4:V6)

  Accomplish exploratory factor analysis with CorEFA, or just efa.
      > efa(n.fact=2)
  The output includes a specification of the multiple indicator measurement
  model derived from the analysis, plus the corCFA code to analyze.
  "

  set.up.plot(5)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  text(0,82, label=f4, adj=0)
  text(0,78, label=f5, adj=0)
  #lines(c(5,90), c(75,75), col=col.line)
  text(0,40, label=t1, adj=0)

  help.more("cfa", 8)
  }


  else if (topic %in% c("prob", "norm", "pt", "qnorm", "prob.tcut")) {
  t0 <- "Probabilities for Normal and t-distributions"

  f1 <- bquote(paste(bold("prob.norm"), "  Normal distribution probability over a range of values"))
  f2 <- bquote(paste(bold("pt"), "  Probability for a t-distribution related to a specified t-value"))
  f3 <- bquote(paste(bold("qnorm"), "  Quantile for a normal distribution"))
  f4 <- bquote(paste(bold("prob.tcut"), "  Cutoff values for a t-distribution based on corresponding quantiles"))


  t1 <- "
  By default, the lessR function prob.norm, provides the corresponding probability
  of obtaining a randomly sampled normal value, Y, in a range of specified values, as
  well as a plot of the normal curve. The R function pt provides the corresponding
  probability for the t-distribution.

  Upper tail probability for t=1.627, df=24:  > pt(1.627, df=24, lower.tail=FALSE)
  Two-tailed p-value for  t=1.627, df=24:     > 2*pt(1.627, df=24, lower.tail=FALSE)
  Probability and curve for a value between 80 and 120 for, mu=100, sigma=15: 
      > prob.norm(lo=80, hi=120, mu=100, sigma=15)

  The quantile functions are the inverse of the probability functions. For a given 
  probability or area under the curve, the corresponding quantile is the 
  corresponding value of the distribution, Y or t. Usually for t, get the cutoff value.

  The lessR prob.tcut also provides a graph of the cutoff t-value. Here for  df=24.
      > prob.tcut(df=24)
      
  Value from the standard normal distribution that cuts off the top 2.5% of the 
  distribution.  Without specifying mu and sigma, the respective defaults are 0 and 1.
      > qnorm(0.025, lower.tail=FALSE)"


  set.up.plot(4)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  text(0,86, label=f3, adj=0)
  text(0,82, label=f4, adj=0)
  #lines(c(5,90), c(78,78), col=col.line)
  text(0,45, label=t1, adj=0)

  help.more("prob.norm", 11)
  }


  else if (topic %in% c("random", "rnorm", "rbinom")) {
  t0 <- "Normal and Binomial Random Values"

  f1 <- bquote(paste(bold("rnorm"), "  Generate randomly sampled values from a normal distribution"))
  f2 <- bquote(paste(bold("rbinom"), "  Generate randomly sampled values from a binomial distribution"))

  t1 <- "
  R can generate simulated sampling from many different population
  distributions, including the normal and the binomial.

  This example generates 50 randomly sampled values from the standard normal 
  distribution, with a default mu of 0 and sigma of 1.
      > rnorm(50)
      
  The generated data can be stored for further analysis.  Here, generate 100 
  values from a normal distribution with a mean of 50 and a standard deviation 
  of 10, store in the vector Y, and then display the resulting histogram.
      > Y <- rnorm(100, mean=50, sd=10)
      > Histogram(Y)
      
  The binomial distribution describes the process of a binary outcome over 
  many different trials, such as flipping a coin.  In this example, flip a fair 
  coin 20 times with a probability of a Head at 0.5.  Then repeat this set of 20 
  flips 10 times to get the number of Heads obtained on each set of 20 flips.
      > rbinom(10, 20, .5)"

  set.up.plot(2)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  text(0,90, label=f2, adj=0)
  #lines(c(5,90), c(85,85), col=col.line)
  text(0,52, label=t1, adj=0)

  help.more("rnorm", 15)
  }


  else if (topic == "sample") {
  t0 <- "Generate Random Samples"

  f1 <- bquote(paste(bold("sample"), "  Generate random samples"))

  t1 <- "
  To use the sample function, first specify the population from which to randomly 
  sample. The population can be defined from the values of a specified variable,
  or the values can be directly listed. Next specify the size to specify the
  number of elements to sample. By default, sampling is done without replacement,
  each value in the population can only appear once in the resulting sample. 

  The following randomly samples 5 values of the variable Y without replacement.
      > sample(Y, size=5)
      
  If the size of the resulting list of sample values is larger than the available 
  number of values from which to sample, then sampling must be done with 
  replacement. To allow sampling with replacement, invoke replace=TRUE.
      > Y <- sample(c(\"Head\",\"Tail\"), size=10, replace=TRUE)
  Here 10 coin flips are simulated, yielding 10 values of Head or Tail. The
  values are stored in the vector Y for further analysis, such as BarChart(Y).
   
  The following randomly samples 10 numbers from the first 100 integers, without 
  replacement.
      > sample(1:100, size=10)"

  set.up.plot(1)
  text(50,100, label=t0, font=4)
  text(0,94, label=f1, adj=0)
  #lines(c(5,90), c(90,90), col=col.line)
  text(0,54, label=t1, adj=0)

  help.more("sample", 15)
  }


  else if (topic == "help.to.pdf") {
  pdf("R_help.pdf")

  t1 <- "Contents of the Help Files for R Function Help()"
  t2 <- "from the R Contributed Package:"
  t3 <- "lessR"
  t4 <- "David W. Gerbing"
  t5 <- "School of Business Administration"
  t6 <- "Portland State University"
  #t7 <- "Version 2.4, July 10, 2012"
  #set.up.plot()
    plot.new()
    plot.window(xlim=c(0,100), ylim=c(0,100))
  text(50,84, label=t1)
  text(50,80, label=t2)
  text(50,76, label=t3)
  text(50,58, label=t4)
  text(50,54, label=t5)
  text(50,50, label=t6)
  #text(50,24, label=t7)

  #Help()
  Help("data")
  Help("Read")
  Help("Write")
  Help("library")
  Help("transform")
  Help("system")
  Help("Histogram")
  Help("BarChart")
  Help("Plot")
  Help("LineChart")
  Help("SummaryStats")
  Help("one.sample")
  Help("ttest")
  Help("ANOVA")
  Help("power")
  Help("Correlation")
  Help("Regression")
  Help("factor.analysis")
  Help("prob")
  Help("random")
  Help("sample")
  dev.off()

  if (getwd() =="/")
    workdir <- "top level of your file system"
  else
    workdir <- getwd()
  cat("PDF file of Help contents located at current working directory.\n")
  cat("   R_help.pdf at: ", workdir, "\n")

  }



  else if (topic == "lessr") {

    help(package="lessR")
  }


  else {
  cat("
  Value ", topic," for Help not recognized.\n
  Complete list of Help topics, enter:  Help()\n
  \n")

  }

}
