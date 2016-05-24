##
# SISUS: Stable Isotope Sourcing using Sampling
#   Erik Barry Erhardt, Began: December 2006
#
#   Current: 3/10/2014
#
#   Program description
#     Package: sisus
#     Version: 3.09-012  [release].[spreadsheet_version]-[code_revision]
#
#
##
## to do:
#
#
#
#
#
#   rev v.0.08-011: Proportions of Sources Consumed (remove concentration/ass.eff variations)
#     write.model.settings() finish adding switch settings
#     Timeseries needs fixing in .txt and .csv file -- maybe problem only when no solution for one of the mixtures
#     Make sure additional.linear.constraints() handles the empty "NA" cells within a linear constraint correctly (as 0s).
#     Backsolve only if non-BMM
#     Resample for PSCMM using independent Beta distributions for each proportion in the cube
#     \delta to appear, look at plotmath and expression(), something like "expression(delta^{13}C)"
#
#   web:
#     site counter for visitors
#     Online bibliography for the SISUS software
#
#   R-package:
#     arguments in sisus.run() function: output.dir = "where results should be output to", defaulting to current directory
#     SVN to updating R-forge
#     populating subdirectories: tests, src, po, inst, exec, demo, data
#     write documentation in man/
#   bugs:
            #### 12/2/2007 8:40PM bug here in cbind  -FIXED
#     In unique solutions and write-out samples, then error (no error for BMM, but error for AECDMM)
#     If small M, mcmc.diags fail
#     If 1 source, mcmc.diags fail at r.raftery.diag$resmatrix[, 3]
#     as package, running against Getting Started example:
#       Warning messages:
#       1: In sample(seq(n.vertices + 1, M.re), M.re, replace = TRUE, prob = weights.p.dirichlet) :
#         Walker's alias method used: results are different from R < 2.2.0
#   samples:
#     chi-square test for uniform sample by comparing to true volume slices of polytope vinci http://www.lix.polytechnique.fr/Labo/Andreas.Enge/Vinci.html
#     only allow dirichlet resample if simplex constraint is specified
#   plots:
#     Improve scatterplot histograms on the diagonal
#     Scatterplot has different par(mar and oma) settings, make standard somehow
#     Convex Hull plot distinct colors for isotope mvn samples
#   output:
#     In log file output isotope table with all input values formatted neatly as I present in my sisus_1 manuscript.
#     Suggested Citations for paper and software
#     Commas in output so easily read into spreadsheet program
#     Also format in "vertical" format so source names are verticle, allow columns of results -- useful for pubs (just transpose original results)
#     Provide all numerical output that I mean to
#     HPD intervals
#     Provide all numerical output also in a format that's easy to copy/paste into a document (eg. CI's side-by-side)
#     Output assumptions specified by data, eg. uniform concentrations if all concentrations are the same
#     Output canned interpretation for results
#     write.input() "quote" non numeric fields so that it can be read as a csv (not important, since have xls)
#   functionality:
#     Unique solution -- improve determination of this given linear constraints. (needed?)
#       recoding in sample.from.polytope.R, see comment in function
#     Robust to bad input
#     Check return code of all system requests, such as file.remove() (not important now that use new [hex] dir each time)
#     Refer to spreadsheets by name in spreadsheet, rather than order of sheets (not important, since fixed order)
#     What to do with overconstrained system?  Least squares?
#     Unified errors:  (some info from Christian Gunning)
#       Common list of information/warning messages: warning$[label], don't use word "error"
#       Building a list of lists strikes me as a sensible first shot:
#       #begin example
#       msglist = list()
#       this.msg = "msg1"
#       msglist[[this.msg]] = list(ret="return this to user",
#       doc="self-documenting code rocks!")
#       print(msglist$msg1$ret) # the default method?
#       print(msglist$msg1$doc) # print doc upon request
#       #end example
#       fancy S4 classes, slots, etc.
#   spreadsheet modifications:
#     Include sample sizes for variance components for calculation of standard errors (is this necessary? because I want population sd)
#   statistics:
#     tests of hypotheses (frequentist and bayes)
#   legal:
#     register StatAcumen as service mark, and SISUS as trademark http://www.uspto.gov/main/trademarks.htm http://www.uspto.gov/go/tac/doc/basic/
#
#   manuscript:
#     Add references from used packages and functions to a software references section
#
##
## completed:
#
#  v.0.09-011: R 2.7
#     5/23/2008 12:46PM Incorporate a switch for diagnostics
#     5/23/2008 12:46PM R 2.7 updates:
#     5/23/2008 12:46PM   plot formats: bmp, jpeg, tiff directly (no gif)
#     5/23/2008 12:46PM   use chull function in place of chplot package
#     5/23/2008 12:46PM remove "prior" worksheet for now since that functionality is not developed
#     5/23/2008 12:46PM worksheet: M validated to be between 1 and 30000
#     3/23/2008 12:46PM Convex hull plot to have tick marks on all sides (top and right, too)
#   v.0.08-010 cont: MCMC diagnostics
#     3/19/2008 12:46PM R MCMC page (http://cran.r-project.org/src/contrib/Views/Bayesian.html)
#     3/19/2008 12:46PM Provide a collection of diagnostic output with interpretations and suggestions if failed
#     3/19/2008 12:46PM   R coda package: (http://cran.r-project.org/src/contrib/Descriptions/coda.html)
#     3/19/2008 12:46PM     as.mcmc() to create as an mcmc object to use by the functions below
#     3/19/2008 12:46PM     autocorr() for autocorrelation
#     3/19/2008 12:46PM     autocorr.plot() to plot matrix of autocorrelations
#     3/19/2008 12:46PM       levelplot(line[[2]]) also for crosscorrelation plot
#     3/19/2008 12:46PM     crosscorr() for correlations between variables
#     3/19/2008 12:46PM     geweke.diag() for testing mean of first 0.1 same as last 0.5
#     3/19/2008 12:46PM     heidel.diag() for relative accuracy for estimating the mean
#     3/19/2008 12:46PM     plot() of mcmc object does trace and density plots for variables
#     3/19/2008 12:46PM     raftery.diag() determines sample size required for estimating quantiles with specified precision
#     3/19/2008 12:46PM     traceplot() for chains
#     3/19/2008 12:46PM     HPDinterval() for giving most likely credible intervals
#     3/19/2008 12:46PM     Consider crosscorr() as supplimental information to the scatterplot matrix
#  v.0.08-010: Proportions of Sources Consumed Mixing Model enhancement (remove concentration/ass.eff variations)
##    5/28/2007  9:52PM workbook switch for including density smoother on marginal histogram hist.density.sw
##    5/28/2007  3:49PM Improve look of number of seconds in process_info.txt
##    5/27/2007  3:49PM Larger Process Log window (emailed RWUI for help 5/24/2007 11:03AM)
##    5/26/2007  3:48PM exclude polytope vertices from most numerical summaries calculations
##    5/24/2007  9:55PM Using RCDD to determine no solution or unique solution
##    5/24/2007  0:02AM Output a more human-readable log
##    5/24/2007  0:02AM Workbook: Add 1 Biomass worksheet with 2 columns
##    5/24/2007  0:02AM Workbook: add switches, make turn purple if 1
##    5/24/2007  0:02AM Workbook: "Output Settings", instead of "Plot Settings"
##    5/24/2007  0:02AM Workbook: Switch for including solution polytope vertices as part of numerical summaries
##    5/24/2007  0:02AM Workbook: Switch for creating each output file
#  v.0.07-009: Convex Hull Plotting error bars
##    5/22/2007 10:42PM Convex Hull plot std err bars on all mixtures on group plot
##    5/22/2007  9:37PM make unix and windows plotting match up
#  v.0.07-008: Changed version numbering syntax: [release].[spreadsheet_version]-[code_revision]
#  v.0.01-07: Multiple Mixtures
##    5/19/2007  5:32PM Combine two convex.hull[.matrix].R functions into one
##    5/19/2007  3:27PM Convex hull include StdDev bars for each element if isotope.stddev.include.sw=1
##    5/19/2007  1:18PM Check n.mixtures, n.sources, n.isotopes against number in spreadsheet and give warning if different
##    5/19/2007 12:36PM Scatterplot matrix for Biomass AND for elements
##    5/19/2007 12:06PM Resampling -- do not resample vertices, include at top of resam
##    5/19/2007 11:19AM Tightenup scatterplot 2d density estimates
##    5/18/2007  7:42PM More colorful marginal histograms
##    5/18/2007  6:53PM Date/Time in process_info.txt
##    5/18/2007  5:06PM Parameter for SW$time.series.sw switch, to output means in a txt and csv file
##    5/13/2007  5:33PM Worksheet color cells that are different than "no effect", validation, protection
##    5/12/2007  6:14PM Always include the solution polytope extremities in the samples so that the extremes are represented (package rcdd)
##    5/11/2007  5:54PM Convex Hull - Individual plots named *_[indy mixture]_[element1]_[element2]
##    5/11/2007  5:54PM Convex Hull - one plot all mixtures together named *_[mixture group name]
##    5/11/2007  4:23PM Determine which sheets need to be read based on the switches in the first sheet (reading takes a while with a large workbook)
##    5/11/2007  2:54PM Name marginal histograms with "Biomass" and isotope names instead of 0, 1, 2, ...
##    5/11/2007  2:54PM Include all mixtures together in Convex Hull plot
##    5/11/2007  2:54PM For all filenames, replace space in names with "_" and other bad characters with "-"
##    5/11/2007  2:54PM Perform seperate calculations for each mixture
##    5/11/2007  2:54PM Parameter for group name if number of mixtures is more than 1
##    5/11/2007  2:54PM Parameter for the number of mixtures
##    5/11/2007  2:54PM Additional worksheets to allow input of multiple mixtures, mixtures and sources on seperate sheets
#  v.0.01-06: Categorical Parameters and Variation in Concentrations and Efficiencies
##    5/10/2007  9:49PM Modified workbook to work on OpenOffice.org (set numeric fields to Numeric format -- can't be General format)
##    5/10/2007  3:09PM Write out sample and resample of p
##    5/10/2007  1:41PM Write model assumptions based on input (beginning)
##    5/10/2007  1:15PM For no concentration or efficiency, only output Biomass table
##    5/10/2007  0:00AM Stopping execution if no feasible solutions
##    5/10/2007  0:00AM Incorporating switches from input
##    5/ 9/2007  9:14PM Updated workbook: 10 sheets, new sheets for Variation in Concentrations and Efficiencies
##    5/ 9/2007  9:14PM Updated workbook: updated Parameters sheet with categorical layout
#  v.0.01-05: Simplex switch, unix plots
##    5/ 9/2007  6:24PM For resampling, provide numerical summaries for sample and resample
##    5/ 9/2007  5:46PM Copy process_info.txt to [filename.prefix]_process_info.txt at end
##    5/ 8/2007 10:48AM Count number of solutions that fall outside the polytope because of variation in isotope values
##    5/ 5/2007  7:15PM plot formats under "unix", in s.plot.settings.begin.end.R
##    5/ 5/2007  7:15PM replace density.pr() in s.plot.marginal.histogram.R with another smoothing function -- package reldist
##    5/ 5/2007  7:15PM For unique solutions (M=1) change output to provide "Unique solution [p1 p2 ... pn]"
##    4/21/2007 11:00AM Simplex switch for use with Fraction-used problems (ex: boreal forest fire)
##    4/21/2007 10:16AM Marginal Histograms, more height for the plots and less space between plots
#  v.0.01-04: Named Analysis, 2 source, convex hull lines between all pairs
##    3/17/2007  1:47PM Convex Hull - points and labels on top of lines, remove line from diagonal matrix plot
##    3/17/2007  1:36PM Convex Hull - draw all concentration lines between all pairs of sources, make concentration lines heavier than the original lines (light original lines)
##    3/17/2007 12:39PM Filename prefix for specific analysis filenames
##    3/17/2007 11:53AM Spreadsheet 1 change:  Add a "name of analysis" parameter in first spreadsheet to include on plots and tables
##    3/17/2007 11:53AM Spreadsheet 1 change:  Allow selection of plot types desired -- perhaps list in parameters sheet with a 0 or 1 for no and yes -- pass vector of yeses to plot functions.
##    3/17/2007  1:26AM Automate burnin by running a number of iterations, and if still on the boundary, continue burn in until off boundary.
##    3/17/2007  0:59AM Line up output column headers for short names
##    3/17/2007  0:38AM 2 sources only, allow complex hull to draw by 1-dim only
##    3/17/2007  0:08AM If alpha.dirichlet are all 1's then don't resample or plot second scatterplot
##    3/16/2007 11:47PM Generalize filenames using variables for names -- not done for file process_info.txt
##    3/16/2007 11:47PM Delete output filenames before starting processing
#  v.0.01-03: Concentration Convex Hull
##    1/ 3/2007 11:01PM Quit if no feasible solution
##    1/ 3/2007 10:46PM Convex Hull plot curved lines based on concentrations
##    1/ 3/2007  2:47PM Fix plot text around borders for plots to files
#  v.0.01-02
##    1/ 3/2007  0:25AM Output plots and tables automatically to files that can be read into LaTeX and Word
#  v.0.01-01
##   12/31/2006  9:05PM Output input settings
##   12/31/2006  8:30PM Titles on tables
##   12/31/2006  4:21PM Convex hull plot as a matrix of figures in the same plot
##   12/31/2006  4:07PM Titles and axis labels for all plots
##   12/31/2006  3:09PM Convex Hull plot to include all names -- expand region a little
##   12/31/2006 11:46AM Put all routines into seperate functions
##   12/30/2006  1:22PM Narrower line in marginal histogram density estimation
##   12/30/2006 11:16AM Read input values from spreadsheet
##   12/29/2006 11:30PM MVN is polytope
##   12/29/2006 10:30PM Draw from MVN for variation in isotope measurements -- plot iterates on convex hull plot so little clouds around each isotope value
##   12/29/2006  4:50PM Print what functions its performing as it goes:  "Reading data [filename]", "Drawing samples", etc.
##   12/29/2006  2:31PM improve consistent axes on scatterplots
##   12/29/2006 12:11PM Dirichlet prior on p, using resampling:  by sample() and ddirichlet() in gtools
##   12/29/2006 11:26AM convex hull before and after discrimination
##   12/28/2006 10:48PM density estimates using library("fdrtool")
##   12/28/2006 10:48PM 2D scatterplots of all pairs using gclus
##   12/20-28/2006      Core development
#
#
## study completed:
#       *  5/14/2007  4:07PM Run all isosource examples, and others from literature (ask Blair)
#       *  3/17/2007 12:39PM theoretical upper bound for how dimension effects sample size required
#
## web completed:
#       *  5/14/2007  4:07PM logo for SISUSWiki -- maybe just scaled image of website sideways flowers
#       *  5/13/2007 10:07PM donations page
#       *  5/14/2007  8:22PM bug reporting link to SISUSWiki
##
#
#
##
#   RWUI for creating SISUS application on web
##
#     Enter a title for the application
#       SISUS:  Stable Isotope Sourcing using Sampling
#     Enter an introductory explanation (optional)
#       SISUS is a Bayesian statistical model and software aiming to provide a comprehensive solution to stable isotope sourcing inference and prediction problems.
#       http://statacumen.com/sisus/
#     Choose a variable input item
#       Current insert position (see facsimile below for positions of existing items)
#         1
#       File Upload Box
#       [Enter]
#     File upload box
#     Enter the name the variable has in your R script
#       filename
#     Enter explanatory text for the user (optional)
#       Upload SISUS Excel workbook for analysis
#       [Enter]
#     [Finished composing page]
#     Validation
#       Validation [radio button]
#       [Enter]
#     Enter a name for the application
#     Appears at end of URL, so choose a simple name with no spaces
#       sisus
#     Enter the name of a results file to be displayed
#       [Finished entering filenames] -- don't define any since list of files depends on options selected
#     Layout of results files
#       [Enter]
#     Number of columns in the layout tables
#       [Enter]
#     Display results on the analysis page?
#       [Enter]
#     Process information
#       [tick box]
#     Dimensions of box to display process information (in pixels)
#       800 800
#       [Enter]
#     Upload R script
#       sisus_run.R
#       [Enter]
#     Upload subsidiary R scripts and data sources (optional)
#       [Finished entering sibsidiary files]
#     Add an initial Login page?
#       No login page [radio button]
#       [Enter]
#     Create application
#       [Create web application]
#     Completed
#     The completed web application for you to download (Right-Click, Save As ...) :
#     sisus.tgz   or    sisus.zip
#       [sisus.zip], download
##
#     Installing on StatAcumen.com/sisus website.
#       Extract sisus.war
#       If current sisus is running at StatAcumen.com,
#         cd /home/erik/website/apache-tomcat-5.5.23/webapps
#         rm sisus.war
#           This will delete the sisus directory
#       copy new sisus.war to that directory
#         Tomcat will automatically deploy sisus.war when upload is complete
#       chmod 755 ./sisus/WEB-INF/srcmd.sh  (don't know better set of permissions, but this works)
#       copy all F:\USERS\Erik\StatAcumen\sisus\SISUS\*.R programs to /home/erik/website/apache-tomcat-5.5.23/webapps/sisus/WEB-INF
##
#     Resizing process_info.txt window
#       Edit these two files:
#         /home/erik/website/apache-tomcat-5.5.23/webapps/sisus/EnterData.jsp
#         /home/erik/website/apache-tomcat-5.5.23/webapps/sisus/message_sisus.jsp
#       Replace "400" by "1600"
#         Detail: appears in these locations
#           EnterData.jsp
#             function disableFormAll() {
#             messagewindow = window.open ("message_sisus.jsp", "message_window_sisus", "dependent=yes,menubar=yes,resizable=yes,scrollbars=yes,width=1600,height=1600");
#             messagewindow.moveTo(1600,0);
#           message_sisus.jsp
#             <object type="text/plain" data="displaypinfo.do" width="1600" height="1600"></object>
#
#     Google Analytics for submission page.
#       Add the following 6 lines of java script to the bottom of the EnterData.jsp file
#          <script src="http://www.google-analytics.com/urchin.js" type="text/javascript">
#          </script>
#          <script type="text/javascript">
#          _uacct = "UA-2046347-1";
#          urchinTracker();
#          </script>
#
#
##
#   Excel workbook creation
#
##
#     When making changes to the workbook, need to:
#       1. get.data.R define sheet names
#       2. assign.variables.R define variables
#       3. sisus.R assign variables
#       4. write.input.R define sheet names
##
#     Data validation for input parameter fields
#       Menu: Data/Validation
#         Settings:
#           Allow: Whole Number, Data: between, Minimum: 0, Maximum: 1
#         Input Message: leave blank
#         Error Alert:
#           Check "Show error..."
#           Style: Stop, Title: "Valid input is 0 or 1 only", Error message: "Please enter 0 or 1"
##
#     Data validation for input data fields
#       Menu: Data/Validation
#         Settings:
#           Allow: Decimal , Data: between, Minimum: -999999, Maximum: 999999
#         Input Message: leave blank
#         Error Alert:
#           Check "Show error..."
#           Style: Stop, Title: "Number", Error message: "Please enter number between -999999 and 999999"
##
#     Conditional format color of cells different from "no effect"
#       Menu: Format/Conditional Formatting
#         Cell Value Is, equal to, [0 or 1], Format [color],
#         Add to do other values
##
#     Conditional format color of cells for input (names)
#       Menu: Format/Conditional Formatting
#         Cell Value Is, equal to, ="" (blank), Format [color],
#         Cell Value Is, not equal to, ="" (blank), Format [color],
##
#     Protecting Worksheet Cells
#       Select worksheet by clicking upper-left-most row/column header
#         Menu: Format/Cells, Protection: Locked
#       Select range for user input
#         Menu: Format/Cells, Protection: Not Locked
#       Menu: Tools/Protection/Protect Sheet, Select unlocked cells only, no password
##
#     Protecting Workwork
#       Menu: Tools/Protection/Protect Workbook, Structure only, no password
#
#
#
##

## devel stuff
##detach(package:sisus); # remove package from R memory
# library(sisus)
# rm(list=ls());                                          # remove all variables
# setwd("F:\\USERS\\Erik\\StatAcumen\\sisus\\SISUS\\work");  # set working directory
# #library(sisus)
# #setwd("F:\\USERS\\Erik\\StatAcumen\\sisus\\SISUS\\work\\test");  # set working directory
# #filename = "SISUS_v0_09_template_devel.xls";
# #sisus.run(filename)
# #filename = "SISUS_v0_08_Hobson2003_Geese_3-2ExampleAE_Albumen_CNsimplexPlot_AECDMM_corrected_2.xls";
# setwd("F:\\USERS\\Erik\\StatAcumen\\sisus\\SISUS\\work\\TimeStudy");  # set working directory
# #filename = "SISUS_v0_08_Time_Study.xls";
# #filename = "SISUS_v0_08_Time_Study2.xls";
# setwd("F:\\USERS\\Erik\\StatAcumen\\sisus\\SISUS\\work\\mink_test");  # set working directory
# filename = "SISUS_v0_08_mink.xls";
# setwd("F:\\USERS\\Erik\\StatAcumen\\sisus\\SISUS\\work\\KP2002BearAE");  # set working directory
# filename = "SISUS_v0_08_KP2002BearAE.xls";
# filename = "SISUS_v0_08_KP2002BearAE_thrash_for_plots.xls";
# sisus.run(filename)
#
## 3/4/2014
# library(sisus)
# setwd("C:/Dropbox/StatAcumen/sisus/SISUS/R-package/bugs/20140304_upgrade");
# filename = "SISUS_v0_09_template_S33only.xls";
# sisus.run(filename)
##
## 8/26/2014 8:23PM
# per https://mail.google.com/mail/u/0/#search/cran+sisus/14786840e8d42927
# Prof Brian Ripley via statacumen.com
# Jul 30
#
# to Erik, CRAN
# This misuses the 'exec' directory:
#
# sisus/exec:
# SISUS_v0_09_template.xls
#
# which is intended for *scripts* (see 'Writing R Extensions').
# It looks like that file is not in fact used but downloaded from http://statacumen.com/sw/sisus/examples/SISUS_v0_09_template.xls .
# Please arrange to install from 'inst', e.g. 'inst/extdata', and use system.file() to locate it.
##
# Moved from /exec to /inst/extdata, removed /exec, incremented version, recompiled

sisus.run <-
structure( # structure to add an example attribute at the end
function# driver program for SISUS
### runs the stable isotope analysis using the specified Excel-like workbook input dataset
##references<< Erhardt, Erik Barry. SISUS: Stable Isotope Sourcing using Sampling, Getting Started. May 30, 2007
##references<< \url{http://statacumen.com/sisus/}
(filename
### The Excel workbook name.  A modified version of sisus_*_templateX.xls
### template workbook available in the package sisus/inst/extdata folder and at http://statacumen.com/sisus/
)
# DRIVER FUNCTION -------------------------------------------------------------
{

  sisus.version = "3.9-013 (2014-08-26)";   # current version of SISUS, from completed list above
  R.version = R.Version()$version.string; # version of R running

  # Output filenames (more after variables to include filename.prefix)
  process.filename = "process_info.txt";  # hardcoded in write.out.R and write.progress.R
  #file.remove(process.filename);   # delete old process_info file

  time.start = proc.time()[3];    # start timer
        p.o = paste("SISUS: Stable Isotope Sourcing using Sampling", "\n"); write.progress(p.o, time.start);
        p.o = paste("      http://StatAcumen.com/sisus, by Erik Barry Erhardt", "\n"); write.out(p.o);
        p.o = paste("      SISUS version ", sisus.version, "running on", R.version, "\n"); write.out(p.o);
        p.o = paste("      Starting ", Sys.time(), "\n"); write.out(p.o);

  ########################################
        p.o = paste("Reading workbook: ", filename, "\n"); write.progress(p.o, time.start);
  DATA = get.data(filename);   # read data in

  ########################################
        p.o = paste("Process input data and assign variables", "\n"); write.progress(p.o, time.start);
  VARIABLES = assign.variables(DATA); # assign variables from input data

  # Analysis Labels
  analysis.name                          = VARIABLES$analysis.name;
  filename.prefix                        = VARIABLES$filename.prefix;
  mixtures.group.name                    = VARIABLES$mixtures.group.name;

  # Data Dimensions
  n.mixtures                             = VARIABLES$n.mixtures;
  n.sources                              = VARIABLES$n.sources;
  n.isotopes                             = VARIABLES$n.isotopes;

  # Sampling
  M                                      = VARIABLES$M;
  n.samples.isotope.mvn                  = VARIABLES$n.samples.isotope.mvn;
  skip                                   = VARIABLES$skip;
  burnin                                 = VARIABLES$burnin;
  seed                                   = VARIABLES$seed;

  # Mixing Model
  SW                                     = VARIABLES$SW;

  # Plot Settings
  hist.bins                              = VARIABLES$hist.bins;

  # Plot Formats
  plot.format.list                       = VARIABLES$plot.format.list;

  # names of source, isotopes, and sources
  names.isotopes                         = VARIABLES$names.isotopes;
  names.mixtures                         = VARIABLES$names.mixtures;
  names.sources                          = VARIABLES$names.sources;

  # discrimination corrected isotope signatures
  isotopes.mixtures                      = VARIABLES$isotopes.mixtures;
  isotopes.sources                       = VARIABLES$isotopes.sources;

  # concentration values
  concentration.sources                  = VARIABLES$concentration.sources;

  # assimilation efficiency values
  efficiency.sources                     = VARIABLES$efficiency.sources;

  # define additional linear constraints
  lc                                     = VARIABLES$lc;

  # priors
  priors.cols                            = VARIABLES$priors.cols;
  priors.type                            = VARIABLES$priors.type;
  priors.precision                       = VARIABLES$priors.precision;
  priors.sources                         = VARIABLES$priors.sources;

  # construct mean vector and covariance matrix for isotopes
  isotope.mean                           = VARIABLES$isotope.mean;
  isotope.sigma                          = VARIABLES$isotope.sigma;

  # construct mean vector and covariance matrix for concentration
  concentration.mean                     = VARIABLES$concentration.mean;
  concentration.sigma                    = VARIABLES$concentration.sigma;

  # construct mean vector and covariance matrix for efficiency
  efficiency.mean                        = VARIABLES$efficiency.mean;
  efficiency.sigma                       = VARIABLES$efficiency.sigma;

  # biomass per individual values
  biomass.per.individual.sources         = VARIABLES$biomass.per.individual.sources;
  number.of.individuals.sources          = VARIABLES$number.of.individuals.sources;

  # end assigning variables from input data
  ########################################




  ########################################
  # write input values out to a file
  output.inputs.filename      = paste(filename.prefix, "_sisus_inputs.txt",    sep="");
    #file.remove(output.inputs.filename);
  write.input(DATA, filename, output.inputs.filename);


  ########################################
  # write execution settings to process log
  write.model.settings(n.mixtures, n.sources, n.isotopes, analysis.name, filename.prefix, mixtures.group.name, SW);

  ########################################
  # set random seed for samples
        p.o = paste("Set random seed =", seed, "\n"); write.progress(p.o, time.start);
  set.seed(seed);

  ########################################
  # prior on p
        # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
  if (SW$prior.include.sw == 1) {
          p.o = paste("Dirichlet prior on p vector", "\n"); write.progress(p.o, time.start);
    alpha.p.dirichlet = prior.on.p(priors.sources, priors.precision, n.sources);
    # if all 1's, then don't resample, otherwise, resample
    resample.dirichlet.sw = 1-as.numeric(identical(TRUE, all.equal(as.vector(alpha.p.dirichlet), rep(1,n.sources))));  # , tol=1e-10
    if (resample.dirichlet.sw == 0) {
            p.o = paste("Dirichlet prior for p is uniform 1s -- don't resample", "\n"); write.progress(p.o, time.start);
    }
    if (resample.dirichlet.sw == 1) {
            p.o = paste("Dirichlet prior for p is not uniform -- resample", "\n"); write.progress(p.o, time.start);
    }
  } else {
    resample.dirichlet.sw = 0;
  } # end if SW$prior.include.sw


  ########################################
  # draw a sample from the joint mvn distribution of isotopes (assume uncorrelated)
        p.o = paste("Draw", n.samples.isotope.mvn, "samples from MVN isotope distribution", "\n"); write.progress(p.o, time.start);
  ISOTOPE.MVN = isotope.mvn.sampling(n.samples.isotope.mvn, isotope.mean, isotope.sigma);
    n.samples.isotope.mvn = ISOTOPE.MVN$n.samples.isotope.mvn;
    isotope.mvn.sample    = ISOTOPE.MVN$isotope.mvn.sample   ;


  ########################################
  # Plot Augmented Convex Hull (http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=61)
  if (SW$plot.convex.hull.include.sw == 1) {
        p.o = paste("Plot Convex Hull", "\n"); write.progress(p.o, time.start);

    if (SW$discrimination.include.sw == 1) { tit.discrimination = "Discrimination Corrected "; } else { tit.discrimination = ""; };
    if (SW$concentration.include.sw  == 1) { tit.concentration  = "Concentration";             } else { tit.concentration  = ""; };
    if (SW$assimeffic.include.sw     == 1) { tit.efficiency     = "Efficiency";                } else { tit.efficiency     = ""; };
    if (SW$concentration.include.sw  == 1 && SW$assimeffic.include.sw == 1) { tit.hyphen = "-";   } else { tit.hyphen = ""; };
    if (SW$concentration.include.sw  == 1 || SW$assimeffic.include.sw == 1) { tit.with = " with "; tit.curves = " Curves"; } else { tit.with = ""; tit.curves = ""; };
    tit.iso = "Isotope Ratios";

    #tit = paste(tit.discrimination, tit.iso, tit.with, tit.concentration, tit.hyphen, tit.efficiency, tit.curves, sep=""); # construct plot title
    tit = tit.iso; # construct plot title

    isotope.sigma.row = matrix(diag(isotope.sigma),nrow=1); # row of standard deviations
    for (i.mixture in seq(1,n.mixtures)) { # one mixture at a time
      ch.n.mixtures = 1; # group convex hull is below
      ch.isotope.mvn.sample.indy.mixture = indy.mixture.isotope.mvn.sample(isotope.mvn.sample, i.mixture, n.mixtures, n.sources, n.isotopes);
      ch.isotope.sigma.indy.mixture      = indy.mixture.isotope.mvn.sample(isotope.sigma.row,  i.mixture, n.mixtures, n.sources, n.isotopes);
      if (n.mixtures == 1) {  # only one mixture
        ch.isotopes.mixtures = isotopes.mixtures;
      } else {                # multiple mixtures
        if (n.isotopes == 1) { # only one isotope
          ch.isotopes.mixtures = isotopes.mixtures[i.mixture];
        } else {                # multiple mixtures
          ch.isotopes.mixtures = isotopes.mixtures[i.mixture,];
        }
      }
      ch.names.mixtures = names.mixtures[i.mixture];
      ch.isotope.mvn.sample = ch.isotope.mvn.sample.indy.mixture;
      ch.isotope.sigma = ch.isotope.sigma.indy.mixture;
      if (n.isotopes > 1) { # skip if only one isotope
          plot.filename = filename.clean(paste("plot_convex_hull_", names.mixtures[i.mixture], sep=""));
        ch.matrix.sw=0;
        s.plot.convex.hull(isotopes.sources, ch.isotopes.mixtures, ch.n.mixtures, n.sources, n.isotopes, names.sources, names.isotopes, ch.names.mixtures, title.mixture, ch.isotope.mvn.sample, concentration.sources, efficiency.sources, ch.isotope.sigma, tit, analysis.name, filename.prefix, plot.filename, plot.format.list, ch.matrix.sw);
      } # if n.isotopes
          plot.filename = filename.clean(paste("plot_convex_hull_matrix_", names.mixtures[i.mixture], sep=""));
        ch.matrix.sw=1;
        s.plot.convex.hull(isotopes.sources, ch.isotopes.mixtures, ch.n.mixtures, n.sources, n.isotopes, names.sources, names.isotopes, ch.names.mixtures, title.mixture, ch.isotope.mvn.sample, concentration.sources, efficiency.sources, ch.isotope.sigma, tit, analysis.name, filename.prefix, plot.filename, plot.format.list, ch.matrix.sw);
    } # end for i.mixture
    if (n.mixtures > 1) { # group convex hull
      ch.n.mixtures = n.mixtures; # individual convex hulls are above
      title.mixture = mixtures.group.name;
      ch.isotopes.mixtures = isotopes.mixtures;
      ch.names.mixtures = names.mixtures;
      ch.isotope.mvn.sample = isotope.mvn.sample;
      ch.isotope.sigma = isotope.sigma.row;
      if (n.isotopes > 1) { # skip if only one isotope
          plot.filename = filename.clean(paste("plot_convex_hull_", title.mixture, sep=""));
        ch.matrix.sw=0;
        s.plot.convex.hull(isotopes.sources, ch.isotopes.mixtures, ch.n.mixtures, n.sources, n.isotopes, names.sources, names.isotopes, ch.names.mixtures, title.mixture, ch.isotope.mvn.sample, concentration.sources, efficiency.sources, ch.isotope.sigma, tit, analysis.name, filename.prefix, plot.filename, plot.format.list, ch.matrix.sw);
      }
          plot.filename = filename.clean(paste("plot_convex_hull_matrix_", title.mixture, sep=""));
        ch.matrix.sw=1;
        s.plot.convex.hull(isotopes.sources, ch.isotopes.mixtures, ch.n.mixtures, n.sources, n.isotopes, names.sources, names.isotopes, ch.names.mixtures, title.mixture, ch.isotope.mvn.sample, concentration.sources, efficiency.sources, ch.isotope.sigma, tit, analysis.name, filename.prefix, plot.filename, plot.format.list, ch.matrix.sw);
    } # end if n.mixtures
  } # end if SW$plot.convex.hull.include.sw

  ################################################################################
  ################################################################################
  # Loop over each mixture.

  for (i.mixture in seq(1,n.mixtures))
  {
 ## i.mixture = 1; # debug
          p.o = paste("\n"); write.out(p.o);
          p.o = paste("==> process Mixture: ", names.mixtures[i.mixture], "\n"); write.out(p.o);

    # use only i.mixture mvn samples
    isotope.mvn.sample.indy.mixture = indy.mixture.isotope.mvn.sample(isotope.mvn.sample, i.mixture, n.mixtures, n.sources, n.isotopes);
    names.mixtures.indy = names.mixtures[i.mixture];


    ########################################
    # Draw from polytope
          p.o = paste("Draw uniform samples from convex solution polytope", "\n"); write.progress(p.o, time.start);
    POLYTOPE.SAMPLE = polytope.multiple.samples(M, n.samples.isotope.mvn, n.isotopes, n.sources, isotope.mvn.sample.indy.mixture, concentration.sources, efficiency.sources, biomass.per.individual.sources, number.of.individuals.sources, lc, skip, burnin, names.sources, SW);
      M.actual           = POLYTOPE.SAMPLE$M.actual          ;
      p.biomass.sam      = POLYTOPE.SAMPLE$sam               ;
      sol.feasible.count = POLYTOPE.SAMPLE$sol.feasible.count;
      n.vertices         = POLYTOPE.SAMPLE$n.vertices        ;
      ##V.sam              = POLYTOPE.SAMPLE$V.sam             ;

      # report when unique solution
      if (M.actual == 1) { p.o = paste("           Unique Solution", "\n"); write.out(p.o);
                  } else { p.o = paste("           Solution Polytope with", n.vertices, "vertices", "\n"); write.out(p.o); }

      # report percentage of feasible solutions
      feasible.percentage = 100*sol.feasible.count/n.samples.isotope.mvn;
          p.o = paste(sol.feasible.count, " of ", n.samples.isotope.mvn, " = ", feasible.percentage, "% of mvn samples lead to feasible solutions", "\n"); write.progress(p.o, time.start);


    ################################################################################
    # continue only if feasible solutions exist
    if (M.actual == 0) {
          p.o = paste("\n"); write.out(p.o);
          p.o = paste("Information:  No feasible solutions.  Check Convex Hull plots.", "\n"); write.out(p.o);
          p.o = paste("\n"); write.out(p.o);
    } else { # solutions exist, continue

      ########################################
      # Resample from dirichlet-weighted samples
      if (resample.dirichlet.sw == 1) {
              p.o = paste("Resample based on Dirichlet prior for p, with parameters alpha=", paste(alpha.p.dirichlet, collapse=" "), "\n"); write.progress(p.o, time.start);
        p.biomass.resam = resample.dirichlet.p(M.actual, p.biomass.sam, alpha.p.dirichlet, n.vertices);
      } # if resample.dirichlet.sw

      ########################################
      # back-solve for isotope contributions
            p.o = paste("Back-solve from biomass to isotopes contributions", "\n"); write.progress(p.o, time.start);
      # p.isotopes has M.actual rows, the columns are groups of sources for each element isotope {iso1(s1, .., sn), .., ison(s1, .., sn)}
        p.isotopes.sam   = model.mass.balance.equation.inverse(M.actual, n.vertices, n.isotopes, n.sources, p.biomass.sam,   concentration.sources, efficiency.sources, biomass.per.individual.sources, number.of.individuals.sources);
      if (resample.dirichlet.sw == 1) {
        p.isotopes.resam = model.mass.balance.equation.inverse(M.actual, n.vertices, n.isotopes, n.sources, p.biomass.resam, concentration.sources, efficiency.sources, biomass.per.individual.sources, number.of.individuals.sources);
      } # if resample.dirichlet.sw

      ########################################
      # display results
            p.o = paste("Plots and Numerical Output", "\n"); write.progress(p.o, time.start);

      # define results stats filename and remove if preexisting
        output.stats.filename = filename.clean(paste(filename.prefix, "_sisus_stats_", names.mixtures.indy, ".txt" , sep=""));
        #file.remove(output.stats.filename);
      # timeseries filename
        output.stats.timeseries.filename = filename.clean(paste(filename.prefix, "_sisus_stats_timeseries_", mixtures.group.name, ".txt" , sep=""));
        ##file.remove(output.stats.timeseries.filename); # if remove file, do outside i.mixture loop
        output.stats.timeseries.filename.csv = filename.clean(paste(filename.prefix, "_sisus_stats_timeseries_", mixtures.group.name, ".csv" , sep=""));
        ##file.remove(output.stats.timeseries.filename.csv); # if remove file, do outside i.mixture loop

      # only output isotope proportions if concentrations or efficiencies were specified
      if ((SW$concentration.include.sw  == 1) || (SW$assimeffic.include.sw == 1)) {
        results.isotopes = seq(1,n.isotopes);   # if concentration or efficiency, then output isotope proportions
      } else {
        results.isotopes = NULL;                # otherwise, just do biomass
      } # end if

      # will do for uniform sampled, then dirichlet resampled if that was requested, in which case resample.dirichlet.sw is 1
      # i.sam=0 is uniform sample, i.sam=1 is dirichlet resample
      for (i.sam in unique(c(0, resample.dirichlet.sw))) {
        if (i.sam == 0) {
          column.names = "sam"; # create column labels for writing out samples
          p.o = paste("   Uniform Sample", "\n"); write.progress(p.o, time.start);
        } # if i.sam
        if (i.sam == 1) {
          column.names = "resam"; # create column labels for writing out samples
          p.o = paste("   Dirichlet Resample", "\n"); write.progress(p.o, time.start);
        } # if i.sam

        for (j.isotope in c(0,results.isotopes)) {  # 0 is p.biomass.resam, 1..k are p.isotopes.resam
          if (j.isotope == 0) {
            p.results.name = "Biomass";
            if (i.sam == 0) { p.results = p.biomass.sam;  } # 0 is p.biomass.sam
            if (i.sam == 1) { p.results = p.biomass.resam;} # 0 is p.biomass.resam
          } else {
            p.results.name = names.isotopes[j.isotope];
            ind.first = 1+(j.isotope-1)*n.sources; ind.last = j.isotope*n.sources;    # defining the columns for the current element isotope
            if (i.sam == 0) { p.results = p.isotopes.sam[,seq(ind.first,ind.last)];}   # 1..k are p.isotopes.sam
            if (i.sam == 1) { p.results = p.isotopes.resam[,seq(ind.first,ind.last)];} # 1..k are p.isotopes.resam
            #if (M.actual==1) {p.results = t(as.matrix(p.results));} # if single solution, still make a matrix
          }

          if (M.actual > 1) {

            #### MCMC was here, moved below 4/23/2008 11:29PM

            ##########################
            # plot marginal histograms
            if (SW$plot.marginal.histogram.include.sw == 1) {
              p.o = paste("     ", p.results.name, "Marginal Histogram", "\n"); write.progress(p.o, time.start);

              if (i.sam == 0) { plot.filename = paste("plot_marginal_hist_unif_sam_",  names.mixtures.indy, "_", p.results.name, sep = ""); plot.filename = filename.clean(plot.filename); }
              if (i.sam == 1) { plot.filename = paste("plot_marginal_hist_dir_resam_", names.mixtures.indy, "_", p.results.name, sep = ""); plot.filename = filename.clean(plot.filename); }
              s.plot.marginal.histogram(p.results, p.results.name, n.sources, n.isotopes, M.actual, hist.bins, names.isotopes, names.mixtures.indy, names.sources, analysis.name, filename.prefix, plot.filename, plot.format.list, SW);
            }; # end if SW$plot.marginal.histogram.include.sw

            ########################################
            # scatterplot matrix
            if (SW$plot.scatterplot.matrix.include.sw == 1) {

              p.o = paste("     ", p.results.name, "Scatterplot Matrix (long time for large samples)", "\n"); write.progress(p.o, time.start);

              if (i.sam == 0) {
                  tit = "Uniform samples from solution polytope";
                  plot.filename = paste("plot_scatterplot_matrix_unif_sam_", names.mixtures.indy, "_", p.results.name, sep = ""); plot.filename = filename.clean(plot.filename);
              }
              if (i.sam == 1) {
                  tit = "Dirichlet resamples from solution polytope";
                plot.filename = paste("plot_scatterplot_matrix_dir_resam_", names.mixtures.indy, "_", p.results.name, sep = ""); plot.filename = filename.clean(plot.filename);
              }
                s.plot.scatterplot.sample(p.results, p.results.name, n.sources, n.isotopes, names.isotopes, names.mixtures.indy, names.sources, analysis.name, filename.prefix, plot.filename, plot.format.list);
               #s.plot.scatterplot.sample(sam,       n.sources, n.isotopes, names.isotopes, names.mixtures, names.sources, tit, analysis.name, filename.prefix, plot.filename, plot.format.list)
               #s.plot.marginal.histogram(p.results, p.results.name, n.sources, n.isotopes, M.actual, hist.bins, names.isotopes, names.mixtures.indy, names.sources, analysis.name, filename.prefix, plot.filename, plot.format.list);
            }; # end if SW$plot.scatterplot.matrix.include.sw
          } # end if M.actual

          #########################
          # numerical summaries
          if (SW$numerical.summaries.include.sw == 1) {
              p.o = paste("     ", p.results.name, "Numerical Summaries", "\n"); write.progress(p.o, time.start);

            if (i.sam == 0) { sample.name = paste("Uniform Sample from Solution Polytope for ",     names.mixtures.indy, sep = "");}
            if (i.sam == 1) { sample.name = paste("Dirichlet Resample from Solution Polytope for ", names.mixtures.indy, sep = "");}
            numerical.summaries(p.results, p.results.name, n.sources, names.sources, names.mixtures.indy, M.actual, n.vertices, analysis.name, sample.name, output.stats.filename, mixtures.group.name, output.stats.timeseries.filename, output.stats.timeseries.filename.csv, i.mixture, SW);
          }; # end if SW$numerical.summaries.include.sw

          # create column labels for writing out samples
          column.names.prefix = p.results.name;
          for (i.sources in seq(1,n.sources)) {
            column.names = c(column.names, paste(column.names.prefix, "_", names.sources[i.sources], sep="", collapse=NULL));
          } # end for i.sources
        } # end for j.isotope

        if (SW$samples.include.sw == 1) {
              p.o = paste("      Write out samples", "\n"); write.progress(p.o, time.start);
          # define samples filenames and remove if preexisting
          output.samples.filename     = filename.clean(paste(filename.prefix, "_sisus_samples_",   names.mixtures.indy,".csv",   sep=""));
            #file.remove(output.samples.filename);
          output.resamples.filename   = filename.clean(paste(filename.prefix, "_sisus_resamples_", names.mixtures.indy, ".csv", sep=""));
            #file.remove(output.resamples.filename);

          ########
          # writing out samples
          if (i.sam == 0) {
            # create table of samples for writing out samples
            #### 12/2/2007 8:40PM bug here in cbind  -FIXED
            # p.biomass.isotopes.sam   = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), p.biomass.sam[(1+n.vertices):(M.actual+n.vertices),]);

            if (M.actual==1) {
              p.biomass.isotopes.sam   = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), t(as.matrix(p.biomass.sam[(1+n.vertices):(M.actual+n.vertices),]))); # if single solution, still make a matrix
              p.isotopes.sam.nvertices = matrix(p.isotopes.sam[(1+n.vertices):(M.actual+n.vertices),],nrow=1); # if single solution, still make a matrix
            } else {
              p.biomass.isotopes.sam   = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), p.biomass.sam[(1+n.vertices):(M.actual+n.vertices),]);
              p.isotopes.sam.nvertices = p.isotopes.sam[(1+n.vertices):(M.actual+n.vertices),];
            }
            if (SW$concentration.include.sw  == 1 || SW$assimeffic.include.sw == 1) { # isotope contributions
              p.biomass.isotopes.sam   = cbind(p.biomass.isotopes.sam, p.isotopes.sam.nvertices);  # p.isotopes.sam is not a matrix -- fixed using *.nvertices 4/22/2008 10:13PM
            }
            colnames(p.biomass.isotopes.sam) = column.names;  # label columns with source names
            write.table(p.biomass.isotopes.sam,   file = output.samples.filename,   append=FALSE, quote=TRUE, sep=",", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE, qmethod="escape");
          };
          # writing out resamples
          if (i.sam == 1) {
            # create table of resamples for writing out resamples
            #### 12/2/2007 8:40PM bug here in cbind  -FIXED
            #p.biomass.isotopes.resam = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), p.biomass.resam[(1+n.vertices):(M.actual+n.vertices),]);
            if (M.actual==1) {
              p.biomass.isotopes.resam = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), t(as.matrix(p.biomass.resam[(1+n.vertices):(M.actual+n.vertices),]))); # if single solution, still make a matrix
              p.isotopes.resam.nvertices = matrix(p.isotopes.resam[(1+n.vertices):(M.actual+n.vertices),],nrow=1); # if single solution, still make a matrix
            } else {
              p.biomass.isotopes.resam = cbind(matrix(data=seq(1,M.actual), nrow = M.actual, ncol = 1), p.biomass.resam[(1+n.vertices):(M.actual+n.vertices),]);
              p.isotopes.resam.nvertices = p.isotopes.resam[(1+n.vertices):(M.actual+n.vertices),];
            }
            if (SW$concentration.include.sw  == 1 || SW$assimeffic.include.sw == 1) { # isotope contributions
              #p.biomass.isotopes.resam = cbind(p.biomass.isotopes.resam, p.isotopes.resam[(1+n.vertices):(M.actual+n.vertices),]);
              p.biomass.isotopes.resam = cbind(p.biomass.isotopes.resam, p.isotopes.resam.nvertices);
            }
            colnames(p.biomass.isotopes.resam) = column.names;  # label columns with source names
            write.table(p.biomass.isotopes.resam, file = output.resamples.filename, append=FALSE, quote=TRUE, sep=",", eol="\n", na="NA", dec=".", row.names=FALSE, col.names=TRUE, qmethod="escape");
          };
        }; # end if SW$samples.include.sw

      }; # end for i.sam

      #########################
      # MCMC Diagnostics
      if (SW$mcmc.diagnostics.sw == 1 && M.actual > 1) {
        i.sam = 0;      # original samples
        j.isotope == 0; # Biomass only
            p.o = paste("      MCMC Diagnostics"); write.progress(p.o, time.start);
        # MCMC diagnostics filename
          output.mcmc.diagnostics.filename = filename.clean(paste(filename.prefix, "_mcmc_diagnostics_", names.mixtures.indy, ".txt" , sep = ""));
        mcmc.diag.status = mcmc.diagnostics(p.biomass.sam, p.results.name, names.mixtures.indy, n.sources, names.sources, n.isotopes, names.isotopes, M.actual, output.mcmc.diagnostics.filename, filename.prefix, analysis.name, plot.format.list);
        if (sum(unlist(mcmc.diag.status)) == 0) { # check whether passed all diagnostics
            p.o = paste(" .. all passed", "\n"); write.out(p.o);
        } else {
            p.o = paste(" .. AT LEAST ONE FAILED, check file:", output.mcmc.diagnostics.filename, "\n"); write.out(p.o);
        }
      } # SW$mcmc.diagnostics.sw

    } # end if, continue only if feasible solutions exist
  } # end for, i.mixture

  #file.remove("Rplots.ps");   # delete useless plot file  (not needed after R 2.7)

        p.o = paste("\n"); write.out(p.o);
        p.o = paste("SISUS Complete", "\n"); write.progress(p.o, time.start);

  file.copy(process.filename, paste(filename.prefix, "_", process.filename, sep=""));  # copy process filename to use prefix

  ### no value returned, many output files produced depending on options specified in workbook

},
  #example
  # # set working directory for many output files with setwd()
  # # see http://statacumen.com/sisus for workbook
  # filename = "http://statacumen.com/sw/sisus/examples/SISUS_v0_09_template.xls";
  # sisus.run(filename)
ex=function(){
  # set working directory for many output files with setwd()
  # see http://statacumen.com/sisus for workbook
  filename = "http://statacumen.com/sw/sisus/examples/SISUS_v0_09_template.xls";
  sisus.run(filename)
}
)
