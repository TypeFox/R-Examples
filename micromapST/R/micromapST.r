#
# Updated Package Version 120828
# Updated Package Version 130426
# Updated Package Version 130506 - V0.94
# Updated Package Version 130510 - V0.95 - fixes.
# Updated Package Version 130511 - V0.96 - attempt to complete - 
# Updated Package Version 130511 - V0.97 (8:00pm) - fixes 
# Updated Package Version 130513 - V0.98 (8:00am) - fixes and testing
# Updated Package Version 130517 - V0.99 - fixes and work with BW.
#                                        - correct ref line color and minor updates.
#                                        - corrected micromapSTDefaults and Arrows errors.
#                                        - label adjustment and fix parameter checking for boxplots
# Updated Package Version 130604 - V1.0.0  - Final Edit and fixes for release.
#                                         - Dynamically defined variables must be globalVariables add.
#                                         - Formal Release of package.
# Updated Package Version 131127 - V1.0.1 - Correct segmented and centered  bars to handle only two data columns
# Updated Package Version 140104 - V1.0.2 - Add diagonal line in scatter plot with equal x and y values.
#                                         - Update NormSeg, Seg, Centered Seg to use variable width bars.
#                                         - Changed method of providing colors and details parameters.
#                                         - Correct median dot in scatter plots
#                                         - Add logic to allow numeric (integer) or column names in col1, col2, col3
#                                         - Correct logic to handle multiple columns in sortVar.
# Updated Package Version 140307 - V1.0.3 - Add Rank Glyph
#                                         - Remove limit on number of time series elements.
#                                         - Plot the median time series data in the panels above and below
#                                           the median row.
#                                         - Adjusted defaults on stacked bar graphs
# Updated Package Version 140712 - V1.0.4 - Correct single and double quote mark usage in examples
# Updated Package Version 150121 - V1.0.5 - Fix panelDesc parameter checking logic. Incorrect error 
#                                           messages generated..
#            
#
#  Update Log and change details by Jim Pearson
#    May 31, 2009 - corrected dates on three column micromap
#        1990-2000 to 2001-5   --> 1996-2000 to 2001-5
#    June 7, 2009 - Added VerStr as a parameter to be able to determine
#        which output files are from this version.
#        - Updated book Micromap-Dot-Arrow-Box plot to use new 
#        data files:
#           WFAgeAdjLungMort2000-4CountyAgeAdj2000.csv
#           WFLungMort19951999AgeAdj2000State.csv
#           WFLungMort20002004AgeAdj2000State.csv
#        and change the titles for the columns in the output to match.
#        - Updated sections to use labels instead of column numbers.
#        - Updated Book micromap to merge two files instead of using
#        one file.  This also changed the column number by +1.
#        Note: future update should look at using column names instead of 
#        numbers.
#        - Updated ARROW chart to plot DOT when difference is zero.
#        - Reduce white space between columns (just a little, cannot be eliminate to
#        maintain readibility.
#    July 22, 2010 - Correct reference value (refVals) code.
#        - add variable for reference value label text (refTexts) per column.
#             panelDesc$refTexts -> vector, one per column.
#        - add variable to color the reference value label test 
#             details$Ref.Text.col
#        - No reference label (legend) is printed if no refTexts for the
#             column is provided.
#    January 30, 2011 - Determine running directory and load
#             panelFunctions.r, panelLayout.Rdata, and micromapST.Rdata 
#             from directory.
#    August 28, 2012 - Cleaned up code and re-packaged it with .onLoad
#        - duplicate variable cleaned up, and unused code removed.
#        - integrated the test/demo code correctly.
#        - made adjustments to handle the micromapST namespace.
#        - changed refVals and refTexts to local variables (lRefVals and lRefTexts) to add clarity.
#        - changed parameter for BoxPlots colMedian to BoxP.Median.col to kill duplication with the colMedian 
#          used on the general graphic
#        - Modified "Details" and "Colors" variable to be unique and
#          re-ordered by subroutine usage.
#    October 5, 2012 - update documentation for review.
#        - deleted second version of panelGroupOutline- in panelFunctions.r
#        - Changed rlStateRefText function to build a legend with a line followed by 
#          the reference text.  Problem was line was on both sides of the label and 
#          in some cases overlaid the text.  This way the line is on the left of the text.
#        - changed default value for reference text from black to mid green to match the line 
#          color.
#    April 26, 2013 - add new panel graphic function - TS and TSConf
#        - added Time Series where each state has a strip within the panel for the line graph.
#        - changed boxPlot argument to panelData to represent more types of auxilary data for the program.
#    May 1-2, 2013  - add new panel graphic functions - ScatDot, StackedBar, and Normalized Bar
#        - add graduated colors to stacked bars and normalized stacked bars.
#        - changed normalized axis labels to percentages.
#        - add Time Series with all plots in one panels (one x-y graph)
#        - change TS confidence band to lighter shade = 10% transparency.
#        - attempted to fix order issues.  On TS series of panels, assume order of the panelData is the 
#          same as the original col1, col2, col3, stateId orders.  When they are re-ordered, Save the 
#          index change to remap back to the old order.  Use this to re-order panelData.
#        - On scatdot and segbar panels, the panelData contains a stateId.  Reordering is 
#          done by using the sorted stateId column in the stateFrame to re-order the panelData frames.
#        - added programing feature to permit adjustments to colsize, left and right margins of a 
#          panel based on the type of panel to be created.  Needed to allow space for the 
#          left axis labels for the time series panels (4).
#    May 4, 2013 - remove prototype strip time series - did not work, code deleted.
#        - Added centered stacked bars.
#        - changed circle size on Scatdot of non-colored dots to 75 smaller.
#        - Changed source of data for "scatdot", "segbar", "normbar", and "ctrbar" from 
#           an extra panelData structure to using columns in the stateFrame call parameters data.frame.
#           Now the col1 and col2 parameters in the panelDesc data.frame indicate which columns or
#           range of columns in the startFrame data.frame to use for the X,Y coordinates or the 
#           set of bar segment values per state.
#    May 6, 2013 - change package name from stateMicromap to micromapST.
#        - updated documentation and added new examples to micromapST.Rd
#    May 8, 2013 - Fixes - change colData to panelData to avoid confusion.
#        - Add parameter value checks to Arrow, Bar, dot, dotSE, dotconf, TS, ScatDot, segbar, normbar, and ctrbar functions.
#        - fix examples 
#    May 9, 2013 - switch the TS Array to be 1=x, 2=y, 3=low-y, 4=high-y.
#    May 10, 2013 - add support for rownames on the time series arrays.
#        - added validation of state ids in boxplots and time series.
#        - added new time series dataset to package.
#        - added panelInBound to generating x and y axis labels.
#    May 11, 2013 - reduced Y axis labels size to get more detail
#        - replaced wflung00cnty data file.
#        - created segbar data file.
#        - fixed problem with saving new time series file - needed names on all dimensions.
#        - fixed problem with at and labels argments on mtext calls.
#        - saved original tests in init/tests directory and replace them 
#          in the micromapST.Rd with the master 6 examples.
#        - cleaned up examples.
#        - added code to try and ensure the min and max values on the y axis 
#          are always printed for the median area (middle).
#        - add code to do Dan's color mixing to get opaque colors in bars.
#    May 17, 2013 - make adjustment for publishing package
#        - adjust grey colors to allow a grey scale color pattern to be used. (based on 
#          ColorBrewer "Greys" for 5 colors.
#        - fixed grey/gray colors issues with dots, etc.  using outline colors.
#        - added circles around dots to make the grey standout more.
#    May 20, 2013 - Added "grays" as an equivalent palette name.
#    May 21, 2013 - Fix ref line color to mid-green, change reftext to black.
#        - check fill color for scat dot, fixed.
#        - changed scat dot median symbol from triangle to dot and filled with blakc.
#        - adjusted box positions on maptail, mapcum, and mapmedian titles.
#        - fixed grays to work with ref lines.
#    May 24, 2013 - finish clean up - fix micromapSTDefaults error during loading.
#        - Final Testing.
#    May 25, 2013 - fixed micromapSTDefaults error on initial load
#        - fixed arror warning by using > .005 as zero.
#        - moved up titles printing to let INTERRUPTED pdf build have titles.
#    May 28, 2013 - fix parameter checking for boxplot list. 
#        - Added names check for box plot,
#        - Added "list" type check for box plot.
#        - Reorganized test to not cause a secondary error.
#        - Added Id.Text.adj parameter to details and rlStateID to adjust text alignment.
#    June 2, 2013 - fix DotSE missing X1 variable - should be x.
#        - Added code to do proper capitalization of state abbreviations and full state names.
#        - Added code to intercept common names for Washington, D. C. and convert to "D.C."
#    June 3, 2013 - Released to CRAN.
#    June 4, 2013 - cran check does not handle automatic variable assignments (around line 3100.)
#          register them with R via globalVariable function to add them to the list for rcmd check.
#          During testing, the variables do not show up as globals and are protected within the 
#          micromapST namespace.  - re-released.
#    Nov. 27, 2013 - Correct the parameter check for segmented and centered bars to permit a 
#          minimum of 2 data columns.
#    Jan 4-9, 2014 - The diagonal line added to the scatter plots must reflect equal x and y values. 
#          Current line is diagonal to the box not the data.
#         - Add option to vary the segment bar width from small to larger from left to right for
#           the NormSeg, SegBar, and Centered SegBar glyphics.
#         - Changed method of setting up details variables within the micromapST namespace.
#           Originally, user had to provide a complete list of all of the details variables.  If
#           one was missing or misspelled, no detection or correction.  New method, starts by 
#           assigning all of the variables from the default values. Then takes the provided details
#           list from the user and merges it into the already declared variables.  If a variable
#           does not exist or is misspelled, it is caught by checking against the default list of names
#           and not processed.  In the future, a similar structure will be used to check the 
#           ranges or types of information to validate the user provided details variable values.
#         - Correct median dot in scatter dot plots to only appear in the 4 and 6 rows (just either side
#           of the median row.
#         - Update logic in sortVar option to correctly handle multiple sort columns.  
#         - Add ability to reference data.frame columns by name in the col1, col2, col3 and sortVar
#           parameters.
#         - Enhanced parameter verification and error checking to help user understand the specific
#           problem and correct it fast.  Don't allow R to abort if possible.
#    March 7, 2014 - Removed limit on the number of points in Time Series
#         - Add code for Rank glyph
#         - The time series line and confidence band are squeezed in the median row space and do not
#           properly show the data.  The median time series data is plotted in the panel above and below
#           to median row to properly present the data using the same aspect ratio as the other data.
#         - Adjusted the defaults for the segbar, ctrbar, and normbar graphics to have no center dot
#           and fixed bar height.
#    July 12, 2014 - Corrected usage of single and double quote marks in examples and test code.
#    December 14, 2014 - The error checking for several panelDesc glyphics and micromapST call 
#           parameters, generating incorrect messages and causing R errors.  
#           The checking logic is rewritten.  The use of is.na and is.null was not appropriate
#           to verify contents and colunms in the panelDesc data.frame.  Changed logic to use names()
#           to check the user provided list names for validity.  In the glyphics, logic was changed 
#           to again use the "names" list to verify the required lists were present for the 
#           specific glyphic.
#           The "colors" parameter is not longer required.  The default color lists are used if
#           an override list is not provided.
#    January 21, 2015 - included citations to JSS article and fixed release NOTES on global variables.
#        .
########

########
#
# Copyrighted 2013, 2014 - by: Dan Carr, GMU and Linda Pickle and Jim Pearson of StatNet Consulting, LLC.
#
########

########
#
#  functions used from RColorBrewer:   brewer.pal
#
#  functions used from graphics:   plot, lines, arrows, polygon, axis, text, mtext, boxplot,
#                                  points, legend, plot.new, plot.default, plot.design, plot.function,
#                                  plot.xy, plot.windows, abline, axTicks, barplot, matplot,
#                                  matpoints, title
#
#  functions used from stats:      qnorm
#
#  functions used from grDevices:  rgb, col2rgb
#
########
#
#  Initial Variables that require setting before running this file:
#
#   current directory <-  location of the three Micromap files
#                   micromapST.r
#                   panelFunctions.r
#                   micromapST.Rbata
#
#   
#
#  The following datasets must be included in the package to provide the boundaries.
#
#   stateNamesFips
#   stateVisBorders
#   stateNationVisBorders
#
#   ONLOAD - micromapST moves these structures into local variables:
#
#    rlStateNamesFips, rlStateVisBorders, rlStateNationVisBorder
#
#   also created is the rlMicromapSTDefaults data.frame with colors and details.
#
#   in the micromapST namespace.
#
######

######
#
#

globalVariables(c("ne","ng","ib","ie",
                "topMar","botMar","botMarLegend","botMardif",
                "leftMarAxis","rowSep","rowSize","groupedRowSize","groupedRowSep",
                
                "Map.width","Id.width","Rank.width",
                
                "sc","pad","padex","padMinus",
                "Title.Line.1.pos","Title.Line.2.pos","Title.Line.3.pos",
                "Title.Line.4.pos","Title.Line.5.pos","lineTiclab",
                "Title.cex",
                
                "Grid.Line.col","Grid.Line.lwd","mgpTop","mgpBottom","padjBottom","mgpLeft",
                
                "Panel.Fill.col","Panel.Outline.col",
                
                "Text.cex",
                
                "Ref.Val.lty","Ref.Val.lwd","Ref.Val.col","Ref.Val.BW.col",
                "Ref.Text.col","Ref.Text.BW.col","Ref.Text.cex",
                
                "Arrow.Head.length","Arrow.lwd","Arrow.cex",
                "Arrow.Shadow.lwd","Arrow.Shadow.col",
             
                "Bar.barht","Bar.Outline.col","Bar.Outline.lwd","Bar.Outline.lty",
             
                "CSNBar.barht",
                "CSNBar.Outline.col","CSNBar.Outline.lwd","CSNBar.Outline.lty",
                "CSNBar.First.barht","CSNBar.Last.barht",
                
                "SNBar.varht","SNBar.two.ended",                
                "SNBar.Middle.Dot","SNBar.MDot.pch","SNBar.MDot.pch.fill","SNBar.MDot.pch.lwd","SNBar.MDot.pch.size",
                "SNBar.MDot.pch.border.col","SNBar.MDot.pch.border.lwd",
                
                "CBar.Zero.Line.col","CBar.Zero.Line.lwd","CBar.Zero.Line.lty",
                "CBar.varht","CBar.two.ended",
          
                "BoxP.thin","BoxP.thick","BoxP.Use.Black",
                "BoxP.Median.Line",
                "BoxP.Median.Dot.col","BoxP.Median.Dot.pch","BoxP.Median.Dot.cex","BoxP.Median.Dot.lwd",
                "BoxP.Median.col",
                "BoxP.Outline.col",
                "BoxP.Outlier.lwd","BoxP.Outlier.cex","BoxP.Outlier.BW.col",
                
                "Dot.pch","Dot.pch.size","Dot.conf","Dot.conf.lwd","Dot.conf.size","Dot.Outline","Dot.Outline.col","Dot.Outline.lwd",
                
                "TS.lwd","TS.Axis.cex","TS.hGrid",
                
                "SCD.Bg.pch","SCD.Bg.pch.lwd","SCD.Bg.pch.size","SCD.Bg.pch.fill",
                "SCD.Fg.pch","SCD.Fg.pch.lwd","SCD.Fg.pch.size",
                "SCD.Median.pch","SCD.Median.pch.lwd","SCD.Median.pch.size","SCD.Median.pch.fill",
                "SCD.Axis.cex",
                "SCD.xsc","SCD.ysc","SCD.hGrid",
                "SCD.DiagLine","SCD.DiagLine.col","SCD.DiagLine.lwd","SCD.DiagLine.lty",
                
                "Id.Dot.pch",
                "Id.Dot.Outline.col","Id.Text.cex","Id.Dot.cex","Id.Text.adj",
                
                "Map.Bg.col","Map.Bg.Line.col","Map.Fg.Line.col","Map.Nation.Line.col",
                "Map.State.Spec.cex","Map.Bg.Line.lwd","Map.Fg.Line.lwd","Map.Nation.Line.lwd",
             
                "MST.Debug",
                
                "stateVisBorders","stateNationVisBorders","stateNamesFips"
                ),
                                     
                "micromapST",add=TRUE)

#
#   Would rather have these variable in the local "micromapST" environment.
#
######

######
#
# Functions 
#
# groupPanelOutline 
#

groupPanelOutline = function (panelGroup, j )
   ## used in micromapST function  - assumes 3 rows in the panels..
{
  for (i in 1:3){
     panelSelect(panelGroup,i,j)  # select a space
     panelScale()               # scale it
     panelOutline()             # outline it.
  }
}   

#
# simpleCap - capitalize each word in a phrase and removes "."s and extra blanks.
#     Not good on vectors - must apply
#

simpleCap <- function (x)
   {
      s <- strsplit(x,"[ ._]")[[1]]
      s1 <- s[s != ""]
      paste(toupper(substring(s1,1,1)),tolower(substring(s1,2)),sep="",collapse=" ")
   }
      
#
# Alternative:
#   gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", name, perl=TRUE)
#

#
#  Subroutine to take the colx vectors, convert numeric to integer, convert character by matching
#   with column names to column numbers.   NA's become "0", Invalid column numbers or names become "0"
#
#.
#  Errors are flaged by the glyphics code.
#
CheckColx <- function(wcol,colname,wnam2,len_wnam) 
   {
     # wcol = col vector from panelDesc
     # colname = character name of col vector for error message.
     # wnam2 = character list of column names and row numbers (in character format)
     # len_wnam = number of original columns.
     #
     # Rules:  "0" means invalid number, out of range number or invalid name.
     #         NAs are converted to "0" values.
     #         Glyphics check for valid values based on need.  
     #
     #
     xwcol <- wcol
     l_wcol <- length(wcol)
     ErrFnd = FALSE
     
     if (is.factor(xwcol))
       { xwcol <- as.character(xwcol) }
       
     if (is.numeric(xwcol))
       {  # have number
       
         rcol <- as.integer(xwcol)  # convert numeric to integer.
         rcol[is.na(rcol)] <- 0     # get rid of NA.  Turn to zeros doesn't et rid of negatives.
      
         if (any(rcol < 0))
           {
             ErrFnd = TRUE
             xmsg <- paste("CCOL-01 The one or more of the column number(s) are negative: ",sep="")
             xmsg <- paste(xmsg,paste(rcol,collapse=", ")," Literal:",wcol)
             warning(xmsg)
           } else {
             if (any(rcol > len_wnam))
               {
                 ErrFnd = TRUE
                 xmsg <- paste("CCOL-02 The one or more of the column number(s) is greater than the columns in the stateFrame data.frame: ",sep="")
                 xmsg <- paste(xmsg,paste(rcol,collapse=", "))
                 warning(xmsg)
               } 
           
           }
         # if ErrFnd = FALSE, the all number in vector are within range.
         # check valid range in glyph  (NA become zeros.) Leave the final check to the glyphics.
       } else {     
         if (is.character(xwcol))
           {  # have character - may be name or number - check each
              # get number for other code, if column name.

              xcol <- match(xwcol,wnam2,nomatch=0)    # match against column names and numbers (as characters)          
              rcol <- ifelse(xcol>len_wnam,xcol-len_wnam,xcol)  # adjust matches to row numbers to real row numbers.       

              # name and character number converted to integer       
              # bad and NA values are "0" and will be caught in the glyphic    
           
           } else {
              # invalid variable type
              ErrFnd = TRUE
              xmsg = paste("CCOL-03 The type of ",colname," panelDesc variable is invalid. ",typeof(xwcol),".  Must be integer or character.",sep="")
              warning(xmsg)
           } 
      }
     if (ErrFnd)
       { 
         return (rep.int(0,l_wcol))
       } else {
         # clean up any NAs in list, set to 0
         rcol[is.na(rcol)]  = 0   # set NA to 0 (invalid)
         return (rcol)
         #print(rcol)
       }
        
   }        


###
#
#  micromapST
#
#  In the "micromapST.Rdata", the micromapST and the 
#  micromapSTSetDefaults functions have been replaced by the following code.
#

micromapST = function(
    stateFrame,
    panelDesc,
    rowNames=c("ab","fips","full")[1],   # default = "ab"
    sortVar=NULL, 
    ascend=TRUE,     
    title=c("",""),
    plotNames=c("ab","full")[2],         # default = "full"
    colors = NULL,
    details = NULL)
{
#
#  Routine:   micromapST
#
#  Created by:  Dr. Dan Carr
#  Updated and Extended by:  Jim Pearson, April 20, 2009
#  Updated and Extended by:  Jim Pearson, August 28, 2012
#  Updated and Extended by:  Jim Pearson, May and June, 2013
#  Updated and Extended by:  Jim Pearson, Nov, 2013
#  Updated and Extended by   Jim Pearson, Jan, 2014
#  Updated and Extended by:  Jim Pearson, March, 2014
#
#  Packaged by: Jim Pearson
#
#  Dependencies:   micromapSTSetDefaults
#  DataSets:
#                  stateNamesFips
#                  stateVisBorders
#                  stateNationVisBorders
#
#           Files: panelFunctions.r
#
#  Call Parameters:
#
#   Defaults List for call simulation
#     stateFrame <- data
#     panelDesc
#     rowNames  <- "ab"
#     sortVar   <- NULL
#     ascend    <- TRUE
#     title     <- c("titles")
#     plotNames <- "full"
#     colors    <- NULL
#     details   <- NULL
#
#
#####
#
# stateFrame  data.frame           # data.frame of state ID and data for micromaps.
#             rownames must be state abbreviations, names, or fips codes
#
#             Used with Dot, DotConf, DotSE, arrows, bars, segbar, ctrbar, and normbar column panels. 
#
#             Not used for boxplots or time series column panels.
#
#             The stateFrame must have the state abbr, state name or fips code as 
#             the rownames of the data.frame.
#     
#             The data.frame must be at least 2 columns for some of the functions
#             in R.  To compensate for possible 1 column data.frames, a column of zero 
#             is appended to the right side of the data.frame to ensure there is always 
#             2 columns.
#
#             An example of the problem:
#               When the structure is ordered xxx[ord,] and then assigned to the working 
#               variable "dat", the dimensions are preserved. 
#               If the data.frame has only one column, the ordering and assigned, 
#               strips the rownames and leaves the dim(dat) = NULL.
#
######
#
# panelDesc   data.frame        # data frame for panel descriptions/definitions               
#             Example
#             panelDesc = data.frame(
#                type=c('mapcum','id','dotconf','dotconf'),                  # manditory column
#                lab1=c('','','White Males','White Females'),                # recommended
#                lab2=c('','','Rate and 95% CI','Rate and 95% CI'),          # optional
#                lab3=c('','','Deaths per 100,000','Deaths per 100,000'),    # optional
#                lab4=c('','','',''),
#                col1=c(NA,NA,2,9),                                          # dependent on "type"
#                col2=c(NA,NA,4,11),                                         # dependent on "type" 
#                col3=c(NA,NA,5,12),                                         # dependent on "type"
#                refVals=c(NA,NA,NA,wflungbUS[,1]),                          # optional
#                refTexts=c(NA,NA,NA,'US Rate'),                             # optional
#                panelData=c('','','','')                                    # required if boxplot or time series used.
#                )
#
#             The first description row describes the first column of panels
#             an so on.  This is a candidate for change since each column
#             describing a column avoids a mental transposition.  
#  
# The type parameter must be present for each panel column.  The other parameters are optionals.
# However, if a parameter is required for any column, it is present for all columns.  
# If not used by a column, the parameter's value for that column should be set to "NA".
#
#  type refers the graphic panel type to be used. The valid types are  
#          "map", "mapcum","maptail","mapmedian",       for maps
#          "id",                                        for state ids
#          "dot", "dotse","dotconf",                    for dot plots
#          "arrow",                                     for arrow plots
#          "bar",                                       for simple bar plots
#          "ts", "tsconf",                              for time series plots
#          "scatdot",                                   for scatter dot plots
#          "normbar","segbar","ctrbar",                 for stacked bar plots
#          "boxplot",                                   for box plot 
#          "rank"                                       for state ranking
#                   
#         For non-highlighted contours:
#             map accumulates states top to bottom
#             maptail accumulates states outside in
#             mapMedian feature above median state above the median and vis versa
#
#         bar  will accept negative values and plot from 0 in that direction.
#
#  col1, col2, col3
#    These values idenfity the column numbers or column names in stateFrame to be 
#       used as data for most of the panel types.  They are used by:
#            "dot", "bar", "dotse", "dotconf", "scatdot", "segbar", "ctrbar", "normbar"
#      ls
#     Panel types using only one column parameter:
#
#       Dot and bar plots require only one column (col1) = value  (height of bar)
#
#     Panel types dotse, arrows, ScatDat, SegBar, CtrBar, NormBar using two column 
#       parameters (col1 and col2):
#  
#       dotse needs:  col1=estimates and col2=standard errors 
#            Plus and minus the SE is draw around the estimates
#
#       arrows needs col1=beginning (older) and col2=ending (newer) values. 
#            The arrow head is on the col2 end of the arrow.
#                  
#       scatdot needs: col1 = x value (horizontal axis), col2 = y value (vertical axis)
#            for each data point (one per state)..
#
#       segbar, ctrbar, and normbar need: col1 is the name or number of the   
#           column in the stateFrame for the first bar segment length values,  
#           col2 is the column name or number of the column in the stateFrame containing 
#           the length of the last bar segment.  The columns between col1 and col2 contain
#           the lengths of the other bar segments in the glyphic.  col1 must preceed col2 in 
#           the stateFrame data.frame.The number of data columns (bar segments) can 
#           range from 2 to 9 columns.
#
#     Panel type dotconf using three column parameters: (col1, col2, col3):
#     
#        dotconf needs: col1=estimate, col2=lower and col3=upper bounds
#
#     Panel following types do not requiring any column parameters:
#
#       boxplots uses the "panelData" vector in panelDesc to provide the name of a saved 
#           boxplot structure.  The boxplot structure is created by saving the 
#           results of aboxplot(...,plot=F) call.
#
#       ts and tsconf use the "panelData" vector in the panelDesc to obtain the name of 
#           a matrix the data for the time series. The name represents a array(51,"x",4).  
#           The first dimension represents the states (51).  The second dimension 
#           represents the number of samples in the time series.  The third dimension 
#           are the "x", "low.y", "y", and "high.y" values for each sample.  
#           For ts glyphics, the "low.y" and "high.y" values are ignored, but required.
#
#  lab1, lab2
#     Two label lines at the top of columns. Use "" for blank, not NA or MULL.
#
#  lab3
#     One label line at the bottom of a each column,
#     typically measurement units
#
#  lab4
#     One label line for used with the Y axis on each panel.  Only used with time series panels.
#
#  refVals           # P-2010/07/23  changed variable from refvals to refVals 
#                    #    to be consistant.
#     name of objects providing a reference values shown
#     as a line down the column 
#
#  refTexts          # JP-2010/07/23 - New 
#     texts to be used as the legend for the reference values.
#     If refTexts for column is NA, then no legend is added.
#
#  colSize           
#     If value > 0 then the calculated column size is overridden by this value.
#     The value is in inches.  If specified, the widths of all columns must be
#     specified. 
#
#  panelData           # (old boxplot column)
#      names a list object with a boxplot data or time series data (x/y or x/yl/ym/yh 
#      data for each state.
#
#      The boxplot list the xxxx$names list must be the abbreviated state id
#      for the entry and the related data in the structure. 
#.
#      Used to link graphic to additional data beyond the 3 data elements 
#      provided in col1, col2, col3 indexes into the stateFrame..
#
#      For boxplot graphics, a list of "boxplot" function values for each state and DC
#        with the names (2 characters) used as the row.names. 
#
#      For time series graphics, the object must be an array(51,"x",4), 
#         where the 1st index is the states (1 to 51), the second index is the number 
#         of time periods ("x") with a minimum of 2 and maximum of 30, and 
#         the third index is the type of variable. The rownames of array must
#         be the associate state id (a 2 character abbreviation).  This 
#         is required so the time series array can be properly associated 
#         with the data in the stateFrame when it's sorted.
#         For time series with no confidence band, column 1 is the x value and 
#             column 2 is the y value.  
#         For time series with a confidence band, column 1 is the x value, 
#             column 2 is the y-low value, column 3 is the y-median value, 
#             and column 4 is the  y-high value.
#                
#      Note:  Some descriptors may be omitted if none of the panel plots need them.
#             often refValues and boxplots can be omitted 
#
#
#####
#
# Individual Parameters:
#
# rowNames: Type of state id used as row.names in stateFrame data.frame. 
#           Acceptable values are: "ab", "full", "fips".
#           The default is "ab" for abbreviation, 
#
# plotNames: State label use in in the plot when an ID column is requested. 
#           The default is the "full" for full name. 
#           Acceptable values are: "ab", "full"
#
# sortVar   The column name or number in the stateFrame to be used as the variable 
#           in sorting.  Can be a vector of column subscripts to break ties.
#           Warning: The sortVar parameter cannot be used to sort a boxplot or 
#           time series, since data is not contained in the stateFrame.
#
# ascend    TRUE default sorts in ascending order.  FALSE indicated descending order.
#
# title     A vector with one or two character strings to use the title.for the page.
#      
#####
#
# List/Control Parameters:  (package default data.frames are used if the colors and 
#      details parameters do not specify an alternate data.frame.  
#      It is strongly recommended to use the default data.frame)
#
# colors   a color palette as a vectors of strings (character-vectors)
#              5 colors for states in a group of 5
#              1 color for the median state
#              1 foreground color for non-highlighted states in the map
#          and 7 matching colors with 20% transparency for time series.
#
#          If a color vector is provided, it's length must = 14.
#
#          If the value of colors is "bw" or "greys", a grey scale is used instead 
#          of the default or user provided colors vector.
#
#      see rlMicromapSTDefaults$colors for more details
#
#
# details   defines the spacing, line widths, colors and many many other details controling the 
#      style and apparence of the generated glyphs.
#
#      see the micromapSTDefaults$details section for more details.
#
#      The function automatically loads the default values into the code when the function 
#      is started.  The user can use the details parameter to override any of the items and values
#      in the micromapST package.  To override a value, create a list as follows:
#
#      details = list(<variable name> = <value>,,,  )
#
#      See the micromapSTSetDefaults function below for a definition of each micromapST 
#      variable and it's default.
#
#####



#______________________Argument Checks______________________
#

micromapSTDefaults <- micromapSTSetDefaults()  # get master list of variables and defaults

#______________stateFrame - data frame______________
#

#  check to see if the stateFrame was provided.

   if (missing(stateFrame) || is.null(stateFrame) || is.na(stateFrame) || !is.data.frame(stateFrame)) 
     { 
       stop("MST-01 First argument (stateFrame) is missing or not a data.frame.")
     }
     
   nr = nrow(stateFrame)
   if (nr!=51)
     {
        stop(paste("MST-02 The first argument (stateFrame) must have 51 rows (states plus DC). It only has",nr,".",sep=""))
     }
     
#
#   JP - Make sure the input data.frame is at least two columns - add one.  A single column data.frame
#        acts differently then a two or more column data.frame under many operations.
#   JP - Dot code (at least) has problems with single column stateFrame structures.
#
#   To protect code and any other areas that may have problems,
#   quick fix is to append "0" column to the right of the provided data.frame.
#   This forces the data.frame to be at least 2 columns.
#

   Ex = rep(0,nr)
   SFrame = cbind(stateFrame,Ex)     # move to SFrame and add Zero column.

# have SFrame and stateFrame put together

#_____________Set up for State names and abbreviation links. (Lazy Data loads)
#

   rlStateNamesFips = stateNamesFips
   rlStateVisBorders = stateVisBorders
   rlStateNationVisBorders = stateNationVisBorders
   
####

   #  Setup for stateId checks
   sortedStateId = sort(stateNamesFips$ab)

   # Get state abbreviation as polygon link
   #  

   fullNames = row.names(rlStateNamesFips)    # List of full state names.

   curnam = row.names(SFrame)                 # Get list of current names in row.names.

   # get proper capitalization of state ab or full names.
   curnam2 = as.vector(sapply(curnam,function(x) simpleCap(x)))
   
   #  Compare against common "DC" names and replace with "D.C."
   DCnames = c("Washington, D. C.",   "Washington D. C.",
               "Washington, D C",   "Washington D C",
               "Washington, Dc",    "Washington Dc",
               "District Columbia", "District Of Columbia",
               "DC","Dc","D C","D. C.")
   
   curnam2[!is.na(match(curnam2,DCnames))] = "D.C."            

#_________ Build the column name list for verification later

   wSFnam      <- names(stateFrame)       # get the column names from data.frame
   len_wSFnam  <- length(wSFnam)          # record the number of "named" rows in list (at front.)
   wSFNameList <- c(wSFnam,seq(from=1, to=len_wSFnam))   # add valid row numbers to the list.

#  wSFNameList now contains a list of the column names and column numbers as character strings.
#  This string will be used to verify any user provided column names or character column numbers.

#
#  headers or US rate rows should not be included in data.format. 
#

#_______________panelDesc structure_______________
#

####### Processed later

#_______________rowNames option___________________
#

   stateId = switch(rowNames,
      # if "ab", use current name
      "ab"=  rlStateNamesFips$ab[match(toupper(curnam), rlStateNamesFips$ab)],
      
      # if "fips", convert to abrv name      
      "fips"= rlStateNamesFips$ab[match(as.integer(curnam), rlStateNamesFips$fips)],
      
      # if "full" state name, convert abrv name
      "full"= rlStateNamesFips$ab[match(curnam2,fullNames)],
      
      #  No match..
      warning("MST-03 Check rownames type, must be 'ab', 'fips', or 'full'.")
   )

   if (any(is.na(stateId)))
     {  # one of the state abrv or full names are not valid
       BadList = paste(curnam[is.na(stateId)],collapse=" ")  # create a list of bad names.
       stop(paste("MST-04 The following row names in the stateFrame data.frame are invalid: ",BadList,sep=""))
     }


#_______________plotNames option__________________
# Get statenames or abbreviations to plot_______________________
   stateNames = switch(plotNames,
          "ab"=stateId,

          "full"= fullNames[match(stateId,rlStateNamesFips$ab)],

          warning("MST-05 Check plotNames type, must be 'ab' or 'full'.")
   )

#_______________title option______________________
#

#  checks missing,, is character, length = 1 or 2.


#_______________sortVar option____________________
#

# sort and store stateFrame, stateid, and stateNames____________
   if (missing(sortVar) || is.null(sortVar) || is.na(sortVar) )
     { 
        # if field omitted (null) sort by state name
        ord = order(stateNames)
     } else  {
        litsortVar = sortVar
        sortVar = CheckColx(litsortVar,'sortVar',wSFNameList,len_wSFnam)
        if (!all(sortVar>0))                     # check to see if all values are good
          {
            warning(paste("MST-06 One of the column names or numbers in the sortVar parameter is out of range or does not exist in the data: ",litsortVar,sep=""))
            ord = order(stateNames)
          } else {  
            if (length(sortVar)>1)
              {
                 # sortVar is a vector of column numbers must do a do.call
                  ord = do.call(order,SFrame[,sortVar])            # if field a numeric and present, sort by specified SFrame column.
              } else {
                  ord = order(SFrame[,sortVar[1]])
              }
         
          }
     }
#_______________ascend option_____________________
#

   if (!(missing(ascend) || is.null(ascend) || is.na(ascend)))
     {
       if (is.logical(ascend))
         {
           if(!ascend)ord = rev(ord)
         } else {
           warning("MST-07 The ascend parameter is not a logical variable.")
         }
     }    


#_______________SORT the data array as requested____________
#
   assign("dat",SFrame[ord,])                       # data fields    "dat" has sorted data frame of the stateFrame
   assign("stateId",stateId[ord])                   # StateID        "stateId" in order of the dat
   assign("stateNames",stateNames[ord])             # StateNames
   assign("datOrder",ord)                           # data order for use with panelData.
   
#

#_______________colors parameter__________________
#

#  See below..

#_______________details overrides_________________
#

#  See below.

#########


#####################
#
# Define panel glyph functions=====================================
#
#    All of these glyph creating functions are internal to the micromapST function.
#

#####
#
# type = 'arrow' =========================================================
#
# rlStateArrow
#
# JP - fixed error when difference is zero.

rlStateArrow = function(j){
  # j = current panel column number
  #  
  #  col1[j] points to the stateFrame column holding the first arrow end point.value
  #  col2[j] points to the startFrame column holding the second arrow end point value
  #
  wnam <- names(dat)
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
    
      xmsg <- paste("ARROW-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (is.na(match('col2',PDUsed))) {
    
      xmsg <- paste("ARROW-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (!ErrFnd) { 
    
      if (is.na(col1[j]) || col1[j] == 0) {  # invalid name or column number
        
           xmsg <- paste("ARROW-01 Specified column name or number in col1 for the first end point is out of range or does not exist: ",litcol1[j]," in stateFrame for column ",j,sep="")
           warning(xmsg)
           ErrFnd = TRUE
      } else {
        if (col1[j] > wdim[2]) {   # invalid name or column number
        
           xmsg <- paste("ARROW-03 Specified column in col1 is too high ",col1[j], " for column ",j,". Literal=",litcol1[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
      if (is.na(col2[j]) || col2[j] == 0) {
        
           xmsg <- paste("ARROW-02 Specified column name or number in col2 for the second end point is out of range or does not exist: ",litcol2[j]," in stateFrame for column ",j,sep="")
           warning(xmsg)
           ErrFnd = TRUE
      } else {
        if (col2[j] > wdim[2]) {
        
           xmsg <- paste("ARROW-04 Specified column in col2 is too high ",col2[j], " for column ",j,". Literal=",litcol2[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
  }  
  if (ErrFnd) return()    # Error warning noted, return from function.
  
  
  x1 = dat[,col1[j]]      # Arrow uses two columns from the state.frame (col1 = arrow start points
  x2 = dat[,col2[j]]      #              col2 = arrow end points.)
  
  refval = lRefVals[j]    # change to lRefVals - JP-2010/07/23   Reference value for column
  reftxt = lRefTexts[j]   # added - JP-2010/07/23                Reference test for column
  
  good1 = !is.na(x1)                   # test to see if both values are present.
  good2 = !is.na(x2)
  good = !is.na(x1+x2)   # used by code to skip bad entries.
  
  if (!all(good1))
    {
       xmsg <- paste("ARROW-05 Missing value in start point data (col1) for column ",j,sep="")
       warning(xmsg)
    }
  if (!all(good2))
    {
       xmsg <- paste("ARROW-06 Missing value in end point data (col2) for column ",j,sep="")
       warning(xmsg)
    }
    
  rx = range(x1,x2,na.rm=T)              # range on of all x1 and x2 values for all states.
  
  rx = sc*diff(rx)*c(-.5,.5)+mean(rx)    # 
                                   #  x-scale extention (sc) = 1.08 *
                                   #  diff of min and max of all * 1/2 + or - to get bracket around mean
                                   #  if range 1 to 25, mean is 13, diff(rx) = 24, --> 0.04 to 25.96 (almost + and - 1)
  ry = c(0,1)                            # Y axis range = 0 to 1.. 

  # ____________labeling and axes_______________

  panelSelect(panels,1,j)               # Select the first panel in the column
  
  panelScale(rx,ry)                     # scale panels for all states, based on above calculations.  rx and ry.
  
  mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)                 # top labels (2)  (above panel # 1)
  mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
 
  atRx <- panelInbounds(rx)
  axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))              # tick labels. 

  panelSelect(panels,ng,j)              # Select the last panel in the column
  
  panelScale(rx,ry)                     # temp set scale to 0 to 1.
  
  # padj in axis needed to make grid line label close
  axis(side=1,mgp=mgpBottom,padj=padjBottom,tck=0,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # bottom pad
  mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)                      # bottom labels.


  #_________________drawing loop__________________
  #  Draw all of the elements - one per state.

  for (i in 1:ng){
     # loop to generate each panel in column
     gsubs = ib[i]:ie[i]       # get range ib to ie (state indexes for this panel) ----  gsubs vector of the indexes for this panel.
     ke = length(gsubs)        # get length  (length = 1 or 5)
     laby = ke:1               # labels 1:1 or 5:1 in most the US state cases.
     
     pen = if(i==6) 6 else 1:ke # if index=6 (?) then pen = 6, else 1:ke (length of line)
     
     panelSelect(panels,i,j)          # select current panel
     panelScale(rx,c(1-pad,ke+pad))   # scale to rx by 1,ke (pad)  (ry = effectively 0.33 to 5.67 (pad = 0.67)
                                      #   Scale = rx by 0.33 to 5.67 with arrows at 1,2,3,4,5...
     panelFill(col=Panel.Fill.col) 
  
     arrLim = max(diff(rx)/par("pin")/1000) * 1.05
  
     axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid lines in panel

     # if a refval is provided then add line.
     if(!is.na(refval))
        {
          lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
        }
      
     panelOutline(col=Panel.Outline.col)     # outline panel 
  
     oldpar = par(lend="butt")        # save old 
  
     for (k in 1:ke){
        # loop through each item in panel (5 or 1)
        m = gsubs[k]     # get index into data array
        if(good[m]){              #  if good values
          # print(paste(k,m,x1[m],x2[m],abs(x1[m]-x2[m])))
          # Getting warning for NON-ZERO length arrows - must be rounding error <> 0.
          #  So, taking liberties to say 0 is .002 and below.  Arrow works in inches??
          #  Alternative is to suppressWarnings...
          if(abs(x1[m]-x2[m])> arrLim){         #  If arrow length is > 1.05/1000 inch do line draw...
             arrows(x1[m],laby[k],x2[m],laby[k],col=colors[pen[k]],
                    length=Arrow.Head.length,lwd=Arrow.lwd)
          } else {
             # length of arrow is zero, so plot a dot..
             points(x1[m],laby[k],pch=20,cex=Dot.pch.size,col=colors[pen[k]])
          }
        }  
     }   
     #  y is from 0 to 6, so the enter line for each arrow is 1,2,3,4,5, etc.
     par(oldpar)
   
   }

  # ____________________________PanelOutline____________________

  groupPanelOutline(panelGroup,j)      # outline full group (column)
  
  # Column done check for reference line.
  
  if(!is.na(refval)) 
             rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/2

}

#####
#
#  type = 'bar' =========================================================
#
#  rlStateBar
#

rlStateBar = function(j){
  # j = current panel column number
  #  
  #  col1[j] points to the stateFrame column holding the bar height from zero.
  #
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
  
      xmsg <- paste("SGLBAR-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
 
  if (!ErrFnd) {
 
      if (is.na(col1[j]) || col1[j] == 0 ) {
         
           xmsg <- paste("SGLBAR-01 Specified column name or number in col1 for the bar height value is out of range, invalid or does not exist: ",litcol1[j]," in stateFrame.",sep="")
           warning(xmsg)
           ErrFnd = TRUE
      } else {
        if (col1[j] > wdim[2] ) {
        
           xmsg <- paste("SGLBAR-02 Column number in col1 is too high: ",col1[j],"  Literal=",litcol1[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
  }
 
  if (ErrFnd) return ()    # error warning found - return.

  py =  Bar.barht*c(-.5,-.5,.5,.5,NA)     #  Bar.barht = 2/3 (0.6667)
  ry = c(0,1)
  refval = lRefVals[j]    # changed to lRefVals - JP-2010/07/23
  reftxt = lRefTexts[j]   # new - JP-2010/07/23

  # ________scale x axis________________________

  x = dat[,col1[j]]             # one column
  
  good = !is.na(x)
  if (!all(good))
    {
      xmsg <- paste("SGLBAR-03 Missing value in bar length data (col1) for column ",j,". Literal=",litcol1[j],sep="")
      warning(xmsg)
    }
  
  rx = range(x,na.rm=T)         # get range of values (min-1, max-2)
  if(rx[2]<=0){                 
      # max < 0..
    rx[2]= 0           # set max to zero
    rx[1] = mean(1,sc)*rx[1]   # adjust min.
  } else if (rx[1] >=0 ) {
       #  min > 0 
       rx[1]= 0      # set min to zero
       rx[2] = rx[2]*(1+sc)/2  # adjust max
    } else {
       # min and max are both > 0 
       rx = sc*diff(rx)*c(-.5,.5)+mean(rx)
    }

  # ____________label axis_______________

  panelSelect(panels,1,j)                         # first panel
  panelScale(rx,ry)                               # scale to match data.
  mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)    # two column top column titles
  mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
  
  atRx = panelInbounds(rx)
  axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # top of column axis labels

  panelSelect(panels,ng,j)                        # last panel
  panelScale(rx,ry)
  # padj in axis needed to make grid line label close
  
  axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))  # both labels
  mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)    # bottom column title

  # _______________drawing loop___________________

  for (i in 1:ng){                       
     gsubs = ib[i]:ie[i]                         # index of elements in panel
     ke = length(gsubs)
     pen = if(i==6)6 else 1:ke                   # Pen indexes.
     laby = ke:1                                 # laby (1 or 1:5)
     
     panelSelect(panels,i,j)                     # select current panel
     panelScale(rx,c(1-pad,ke+pad))              # scale to 1 or 5 entries         
     panelFill(col=Panel.Fill.col)
     
     axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grids
     
     # if a refval is provided then add line.
        if(!is.na(refval))
           {
             lines(rep(refval,2),c(1-padMinus,ke+padMinus), lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
           }          
     panelOutline(col=Panel.Outline.col)                # outline full panel
     
     for (k in 1:ke){
        m = gsubs[k]                             # draw each entry (1 to ke), get index from gsubs
        val = x[m]                               # get value for bar height
        if(good[m]){
           # good value - draw bars are polygons.  (why to polygon)
           polygon(c(0,val,val,0,NA),rep(laby[k],5)+py,col=colors[pen[k]]) 
           polygon(c(0,val,val,0,NA),rep(laby[k],5)+py,
               col=Bar.Outline.col,lwd=Bar.Outline.lwd,density=0)
        }
        lines(c(0,0),c(1-.5*Bar.barht,ke+.5*Bar.barht),col=1) # bar base line  
     }   
  }

  # ____________________________PanelOutline____________________

  groupPanelOutline(panelGroup,j)
 
  # _______Reference Value Legend

  if(!is.na(refval)) 
             rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}

#####
#
#  type = 'boxplot' ======================================================
#
#  rlStateBoxplot
#

rlStateBoxplot = function(j,boxnam){
   
   ErrFnd = FALSE
   
   boxlist = tryCatch(get(boxnam, pos=1),error=function(e) e)
   
   if (inherits(boxlist,"error")) {
        # could not find object named in boxnam.
        warning(paste("BOXP-01 List named:",boxnam," does not exist or is bad.",sep=""))
        ErrFnd = TRUE
     
   } else {
        if (!is.list(boxlist)) {
          
             warning("BOXP-02 Data structure for Boxplots must be a list.")
             ErrFnd = TRUE

        } else {
          
            lnam = names(boxlist)    # names of lists in boxlist
            
            if (is.null(lnam) || is.na(lnam)) {
              
                 warning("BOXP-11 The boxplot data is not valid.")
                 ErrFnd = TRUE

            } else { 
                if (length(lnam) != 6) {
                  
                    # must have at least 1 element and name.
	            warning("BOXP-11 The boxplot data is not a valid structure. Should contain 6 lists.")
	            ErrFnd = TRUE

                } else {
                    nbox = c("stats","n","conf","out","group","names")  # correct list of names for boxplot data.
               
                    if (any(is.na(match(lnam,nbox)))) {
                      
                         # at least one of the list names does not match or is missing.
                         warning("BOXP-04 The boxplot data names do not match the standard boxplot function output list names. Invalid structure.")
                         ErrFnd = TRUE
   
                    } else {
                        nc = dim(boxlist$stat)[2]                # number of rows in boxplot stats data list.
                        if (nc != 51) {
                          
                             warning("BOXP-05 The $stats matrix in the boxplot data must have 51 elements - one for each state and DC.")
                             ErrFnd = TRUE
                        }
   
                        nr = dim(boxlist$stat)[1]                # get number of 
                        if (nr != 5) {
                          
                             warning("BOXP-06 The $stats matrix in the boxplot data does not have 5 values per state/DC.") 
                             ErrFnd = TRUE
                        }
       
                        nn = sort(unique(boxlist$names))          # get list of unique state ids used 
                        if (is.null(nn) || is.na(nn)) {
                          
                             warning("BOXP-03 The list of names is missing from the boxplot data.")
                             ErrFnd = TRUE
                          
                        } else {
                            if (length(nn) != 51) {
                              
                                 warning("BOXP-07 The boxplot list does not contain 51 entries.")
                                 ErrFnd = TRUE
                            }
      
                            tnn = is.na(match(nn,sortedStateId))
                            if (any(tnn)) {   # test to see if any did NOT match
                              
                                 lnn = paste(nn[tnn],collapse=" ")
                                 warning(paste("BOXP-08 The abbreviated state ids found in the boxplot list $names list contain invalid values: ",lnn,sep=""))
                                 ErrFnd = TRUE
                            }
                        }
                    }
                }
            }
        }   
   }
   
   if (ErrFnd) return ()
   
   refval = lRefVals[j]              # get referrence to object, changed 
                                     #    to lRefVals - JP-2010/07/23
   reftxt = lRefTexts[j]             # new - JP-2010/07/23

   #_______________Scaling____________
   
   # y boxplot scaling               # standard - horizontal box - no vertical 
                                     #     (y) dimensions
   py = c(-.5,-.5,.5,.5)
   thiny = BoxP.thin*py
   thicky = BoxP.thick*py 
   medy = BoxP.Median.Line*c(-.5,.5)
  
   ry = c(0,1)                       # used in y scaling for grid lines
  
   #_______________Gather stats and put in State Order______________
  
   # For the moment match on names
   #                     Boxlist = names, stats, out, group, 
   #
   # Boxplot function generates a list value containing:
   #     stats  = matrix - each column is lower, lower hinge, median, upper hinge, upper wicker for plot/group
   #     n      = vector of number of observ in each group
   #     conf   = a matrix which each col contins the low/upper extremes
   #     out    = valies of any data points which lie extremes of whiskers
   #     group  = vector (same length as out) whose elements indicate to which group
   #     names  = vector of names for the groups  (must be 2 char state names)
   #              There must be 51 unique names that match the state abbreviation list.
   #
   
   stats = boxlist$stats       # statistics: 1-low,2-25%,3-median,4-75%,5-high 
                               #   - 5 variables for each state.
   #  indexes to boxplot values.   (pull values into thin and thick)  (set up for "boxes")
   thin = stats[c(1,5,5,1),]   # a column for each state - thin line - outliers (Lower, upper wickers)
                               #   - columns in boxlist (1,5,5,1)
   thick = stats[c(2,4,4,2),]  # a column for each state - thick line - 25% to 75% (lower and upper hinge)
                               #   - columns in boxlist(2,4,4,2)
   med = stats[3,]             # a single value for each state (median)
  
   nam = boxlist$names         # state name.
  
   # conf = boxlist$conf       # matrix of extremes - not used.
   
   outlier = rep(F,length(med))   # build vector of all outliers - set to False
   if(!is.null(boxlist$out))
     {                               # if outliers exist
       out = boxlist$out
       group = boxlist$group
       outlier[unique(group)] = T
       # set to True if we have an outlier to graph.
     }


   #### Need to put in order
   ord = match(stateId,nam)   # ord based on match between boxplot$names and StateIDs.  (Convert XX to index.

   # what about missing values  -  if NA do not plot on that line

   # What about name type inconsistency  
   # I will require use of state name abbreviation
   
   # Fips codes be useful
   #    split() based on first two digits of county fips  
   #    I could stash state fips in stateFrame sorted order

   # For Boxplot median sorting    
   #   Currently the user would need to sort the 
   #   medians in the state frame making sure
   #   the row.names were correct.
   #
   #   JP-no data in col1, col2, or col3 to sort like the other columns... All of the data is in these structures.
   #   
   #   boxlist$stats[3,]   # the median.
   #
   #   at present no re-ordering of the boxplots like the other plots.
   #   JP-if other column is sorted, boxplots will follow that order via the indexes.
   #

   # ___________ scale x axis_______________

   if(is.null(out)) rx = range(stats,na.rm=TRUE) else    # if no outliers - range only on stats
              rx = range(stats,out,na.rm=TRUE)           # if outliers - range on stats and outliers
   
   rx = sc*diff(rx)*c(-.5,.5)+mean(rx)        # min to max range with expansion factors.
   # are these used.
   dx = diff(rx)/200                          # difference / 200 (??? use)
   px= c(-dx,-dx,dx,dx)                       # is this used???

   # ____________titles and labeling axes_______________

   # _____________top of column______
   panelSelect(panels,1,j)   # top panel - add title.
   panelScale(rx,ry)
   mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)             # top column titles
   mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
   
   atRx = panelInbounds(rx)
   axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # top axis labels

   # ____________bottom of column____
   panelSelect(panels,ng,j)  # bottom panel - add sub title and refvals.
   panelScale(rx,ry)
   # padj in axis needed to make grid line label close
   axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # bottom axis labels
   mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)             # bottom column titles

   # _______________drawing loop___________________

   oldpar = par(lend="butt")

   for (i in 1:ng){

      # Cycle through the Row/Groups in the micromap column
      
      gsubs = ib[i]:ie[i]    # get beginning to end row number in group  
      ke = length(gsubs)     # get number of rows in group  
      
      pen = if(i==6) 6 else 1:ke  # if middle group (6), then pen=6, otherwise pen = c(1...x)   
      
      laby = ke:1            # laby = reverse order list for row index.         
      
      panelSelect(panels,i,j)   # select panel for group i in column j)
      panelScale(rx,c(1-pad,ke+pad))   # set scale for panel
      panelFill(col=Panel.Fill.col)           # set fill for panel
      
      axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid lines
  
     # if a refval is provided then add line.
     if(!is.na(refval))
        {
          lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
        }
      
      panelOutline(col=Panel.Outline.col)     # outline panel

      for (k in 1:ke){
         # cycle through row-groups and build each box plot
         
         m = ord[gsubs[k]]   # m is the location of the state in boxlist
         if(is.na(m)) next   #   if no location - skip box plot for state
         kp = pen[k]         # color number
         ht = laby[k]
         
         if(outlier[m]){
            #   plot points for outliers
            vals = out[group==m]
            if (colFull)
              {  # full color do the correct color
                 points(vals,rep(ht,length(vals)),pch=1,
                    col=ifelse(BoxP.Use.Black,"black",colors[kp]),
                    cex=BoxP.Outlier.cex,lwd=BoxP.Outlier.lwd)
              } else {
                 # Greys - do the a grey.
                 points(vals,rep(ht,length(vals)),pch=1,
                    col=BoxP.Outlier.BW.col,
                    cex=BoxP.Outlier.cex,lwd=BoxP.Outlier.lwd)
              }
         }  
 
         # draw thin lower to upper box.
         polygon(thin[,m],rep(ht,4)+ thiny,col=colors[kp],border=NA)
#        polygon(thin[,m],rep(ht,4)+ thiny,col=BoxP.Outline.col,density=0) # don't outline boxes
         # draw middle think box
         polygon(thick[,m],rep(ht,4)+ thicky,col=colors[kp],border=NA)
#        polygon(thick[,m],rep(ht,4)+ thicky,col=BoxP.Outline.col,density=0) # don't outline boxes

#        points(med[m],ht,col=BoxP.Median.Dot.col,pch=BoxP.Median.Dot.pch,cex=BoxP.Median.Dot.cex)  # don't put a dot in.
#        points(med[m],ht,col="black",pch=1,cex=BoxP.Median.Dot.cex)

#        polygon(med[m]+px,ht+BoxP.Median.Line*dy,lwd=1,density=0)   # median line?

#        Lines looked crooked
         segments(med[m],ht+medy[1],med[m],ht+medy[2],         # use segment line.
               col=BoxP.Median.col,lwd=BoxP.Median.Dot.lwd)
#        lines(rep(med[m],2),ht+medy,col=BoxP.Median.col,lwd=BoxP.Median.Dot.lwd)
      }   
   }
   par(oldpar)
   # ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
         rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23


}


#####
#
# type = 'dot'   =====================================================
#
# rlStateDot
#

rlStateDot = function(j){

  # Single Dot, no extra line or interval
  #
  # j = current panel column number
  #  
  #  col1[j] points to the stateFrame column holding the first arrow end point.value
  #
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
    
      xmsg <- paste("SGLDOT-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (!ErrFnd) {  
    
      if (is.na(col1[j]) || col1[j] == 0) {
        
           xmsg<-paste("SGLDOT-01 Specified column name or number in col1 for the dot value is out of range or does not exist: ",litcol1[j]," in stateFrame in column ",j,sep="")
           warning(xmsg)
           ErrFnd = TRUE
      } else {
        if (col1[j] > wdim[2]) {
        
           xmsg <- paste("SGLDOT-02 Specified column name or number in col1 is too high: ",col1[j]," for column ",j,". Literal=",litcol1[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
  } 
  if (ErrFnd)  return ()    # error warning found - return

  #  JB - add "as.double(as.vector(" to handle variation in how objects are converted.
  
  xdat = as.double(as.vector(dat[,col1[j]]))   # one value - the dot.
 
  good = !is.na(xdat)
  if (!all(good))
    {
       xmsg<-paste("SGLDOT-03 Missing value in dot data (col1) for column ",j,". Literal=",litcol1[j],sep="")
       warning(xmsg)
    }
  refval = lRefVals[j]    # get reference value for this column, changed 
                         #   to lRefVals - JP-2010/07/23
  reftxt = lRefTexts[j]   # new - JP-2010/07/23

  ry = c(0,1)

  #____________scale x axis______________________
  rx = range(xdat,na.rm=TRUE)
  rx = sc*diff(rx)*c(-.5,.5)+mean(rx)

  # ____________labeling axis_______________
  panelSelect(panels,1,j)
  panelScale(rx,ry)
  mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)   # top column titles
  mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)

  atRx = panelInbounds(rx)
  axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))  # top axis labels
 
  panelSelect(panels,ng,j)
  panelScale(rx,ry)
  # padj in axis needed to make grid line label close
  axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))  # bottom axis labels
  mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)   # bottom column titles

  # _______________drawing loop___________________
  for (i in 1:ng){
     gsubs = ib[i]:ie[i]
     ke = length(gsubs)
     pen = if(i==6) 6 else 1:ke
     laby = ke:1 
     panelSelect(panels,i,j)
     panelScale(rx,c(1-pad,ke+pad))
     panelFill(col=Panel.Fill.col)
     axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid
     
     # if a refval is provided then add line.
     if(!is.na(refval))
        {
          lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
        }
     
     panelOutline(Panel.Outline.col) 
     for (k in 1:ke){
        # step through values for this panel
        m=gsubs[k]
        if(good[m]){    # if good - plot dot.
           if (doDotOutline) 
             {
               points(xdat[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,lwd=Dot.Outline.lwd, col=Dot.Outline.col,bg=colors[pen[k]])         
             } else {
               points(xdat[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,col=NA, bg=colors[pen[k]])
             }
        }
        
     }
  }

  # ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
             rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}

#####
#
#  type = 'dotconf' ====================================================
#
#  flStateDotConf
#

rlStateDotConf = function(j){
  #
  #  j is the current panel column index
  #
  #   col1 indicates the column number for the dot value in the stateFrame.
  #   col2 indicates the column number for the lower confidence value in the stateFrame.
  #   col3 indicates the column number for the upper confidence value in the stateFrame.
  
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
      xmsg <- paste("DOTCONF-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (is.na(match('col2',PDUsed))) {
      xmsg <- paste("DOTCONF-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (is.na(match('col3',PDUsed))) {
      xmsg <- paste("DOTCONF-12 'col3' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (!ErrFnd) {  
    
      if (is.na(col1[j]) || col1[j] == 0 ) {
         
           warning(paste("DOTCONF-01 Specified column name or number in col1 for dot values is out of range or does not exist: ",litcol1[j]," in stateFrame for column ",j,sep=""))
           ErrFnd = TRUE
      } else {
        if (col1[j] > wdim[2] ) {
        
           xmsg<-paste("DOTCONF-04 Specified column number is too high for col1: ",col1[j],". Literal=",litcol1[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
      if (is.na(col2[j]) || col2[j] == 0 ) {
         
           warning(paste("DOTCONF-02 Specified column name or number in col2 for lower confidence values is out of range or does not exist: ",litcol2[j]," in stateFrame for column ",j,sep=""))
           ErrFnd = TRUE
      } else {
        if (col2[j] > wdim[2] ) {
        
           xmsg<-paste("DOTCONF-05 Specified column number is too high for col2: ",col2[j],". Literal=",litcol2[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
      if (is.na(col3[j]) || col3[j] == 0 ) {
      
           warning(paste("DOTCONF-03 Specified column name or number in col3 for upper confidence values is out of range or does not exist: ",litcol3[j]," in stateFrame for column ",j,sep=""))
           ErrFnd = TRUE
      } else {
        if (col3[j] > wdim[2] ) {
        
           xmsg<-paste("DOTCONF-06 Specified column number is too high for col3: ",col3[j],". Literal=",litcol3[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
 
  }
  
  if (ErrFnd) return ()        # error warning found - return
 
  # get data and verify it
  
  x = dat[,col1[j]]              # Col 1 = DOT - median/mean
  
  lower = dat[,col2[j]]          # Col 2 = lower
  upper = dat[,col3[j]]          # Col 3 = upper
 
  good1  = !is.na(x)              # Good Col1 values (dot value)
  good2L = !is.na(lower)         # Good col2 values (lower)
  good2U = !is.na(upper)         # Good col3 values (upper)
  good2  = !is.na(upper+lower)
  
  if (!all(good1))
    {
      xmsg <- paste("DOTCONF-07 Missing value in median for column ",j,sep="")
      warning(xmsg)
    }
  if (!all(good2L))
    {
      xmsg <- paste("DOTCONF-08 Missing value in lower confidence interval for column ",j,sep="")
      warning(xmsg)
    }
  if (!all(good2U))
    {
      xmsg <- paste("DOTCONF-09 Missing value in upper confidence interval for column ",j,sep="")
      warning(xmsg)
    }
 
  refval = lRefVals[j]           # changed to lRefVals, JP-2010/07/23
  reftxt = lRefTexts[j]          # new - JP-2010/07/23

  ry = c(0,1)

  #_____________scale x axis________________
  rx = range(upper,lower,x,na.rm=TRUE)
                        #  x may not be needed???
  rx = sc*diff(rx)*c(-.5,.5)+mean(rx)

  # ____________labeling axes_______________
  
  # panel 1 has line 1 and line 2, top Axis + later image.

  panelSelect(panels,1,j)     # labels (line 1, line 2 and top axis)
  panelScale(rx,ry)
  mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)    # top column titles
  mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
  
  atRx = panelInbounds(rx)
  axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))  # top axis labels
 
  # panel ng has bottom axis and line 3 + later image.
  panelSelect(panels,ng,j)    # labels (bottom -> axis and line 3)
  panelScale(rx,ry)
  # padj in axis needed to make grid line label close
  axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # bottom axis labels
  mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)  # bottom column title
 
  # line 4 is added if refval is present

  #_____________drawing loop___________________
  
  for (i in 1:ng){
     gsubs = ib[i]:ie[i]
     ke = length(gsubs)
     pen = if(i==6) 6 else 1:ke
     laby = ke:1
     panelSelect(panels,i,j)   
     panelScale(rx,c(1-pad,ke+pad))   # Adjusted scale for interior
     panelFill(col=Panel.Fill.col)
     axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # vertical grid lines
     
     # if a refval is provided then add line.
     if(!is.na(refval))
        {
          lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
        }
     
     panelOutline(col=Panel.Outline.col)     # outline scaled image.
     for (k in 1:ke){ 
        m = gsubs[k]
        
        #
        #  Not checking good1 - change code and test.
        #

        if(good2[m])  # if valid upper value.
          {
             lines(c(lower[m],upper[m]),rep(laby[k],2),
                   col=colors[pen[k]],lwd=Dot.conf.lwd)
          }
        if (doDotOutline) 
          {
             points(x[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,lwd=Dot.Outline.lwd,col=Dot.Outline.col,bg=colors[pen[k]])         
          } else {
             points(x[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,col=NA, bg=colors[pen[k]])
          }
     }   

     #   segments(lower[gsubs],laby,upper[gsubs],laby,col=color[pen],lwd=Dot.conf.lwd)
     #   points(x[gsubs],laby,pch=pch,cex=dotCex,col=colors[pen])
     #   points(x[gsubs],laby,pch=1,cex=Dot.pch.size,col=Dot.Outline.col,
     #         lwd=Dot.Outline.lwd)
   }

  # ____________________________PanelOutline____________________

  groupPanelOutline(panelGroup,j)
  
  #  Put legend at end of column for reference value.
  if(!is.na(refval)) 
             rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}

#####
#
# type = 'dotse' =======================================================
#
# rlStateDotSe
#

rlStateDotSe = function(j){
  #   j = current panel column
  #
  #   col1 indicates the column number for the dot value in the stamicroteFrame.
  #   col2 indicates the column number for the SE value in the stateFrame.
  
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
    
      xmsg <- paste("DOTSE-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (is.na(match('col2',PDUsed))) {
  
      xmsg <- paste("DOTSE-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
    }
  if (!ErrFnd) {
  
      if (is.na(col1[j]) || col1[j] == 0 ) {
      
           warning(paste("DOTSE-01 Specified column name or number in col1 for dot values is out of range, invalid, or does not exist: ",litcol1[j]," in stateFrame.",sep=""))
           ErrFnd = TRUE
      } else {
        if (col1[j] > wdim[2] ) {
        
           xmsg<-paste("DOTSE-03 Column number for col1 data is too high: ",col1[j]," for column ",j,". Literal=",litcol1[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
      if (is.na(col2[j]) || col2[j] == 0)  {
      
           warning(paste("DOTSE-02 Specified column name or number in col2 for SE values is out of range, invalid, or does not exist: ",litcol2[j]," in stateFrame.",sep=""))
           ErrFnd = TRUE
      } else {
         if (col2[j] > wdim[2] ) {
        
           xmsg<-paste("DOTSE-04 Column number for col1 data is too high: ",col2[j]," for column ",j,". Literal=",litcol2[j],sep="")
           warning(xmsg)
           ErrFnd = TRUE
        }
      }
  }
  
  if (ErrFnd) return ()   # error warning found - return
  
  x = dat[,col1[j]]
  zval = qnorm(.5+Dot.conf/200)
  inc = zval*dat[,col2[j]]
  upper = x+inc
  lower = x-inc
 
  good1 = !is.na(x)
  if (!all(good1))
    {
      xmsg <- paste("DOTSE-05 Missing data in the dot median data (col1) for column ",j,sep="")
      warning(xmsg)
    }
  good2 = !is.na(upper)
  if (!all(good2))
     { 
       xmsg <- paste("DOTSE-06 Missing data in the SE data (col2) for column ",j,sep="")
       warning(xmsg)
     }
  
  if(any(is.na(inc)))
    {
       warning(paste("DOTSE-07 Missing Value in Standard Errors for column ",j,sep=""))
    }
  refval = lRefVals[j]          # changed to lRefVals, JP-2010/07/23
  reftxt = lRefTexts[j]         # new - JP-2010/07/23

  ry=c(0,1)

#_______________scale x axis__________________
  rx = range(upper,lower,x,na.rm=TRUE)  # use upper, lower and x to find "range" of x
          # x may not be needed at all. But best to leave.
  rx = sc*diff(rx)*c(-.5,.5)+mean(rx)

# ____________labeling axes_______________

  panelSelect(panels,1,j)
  panelScale(rx,ry)
  mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)   # top column titles
  mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
  
  atRx = panelInbounds(rx)
  axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))  # top axis labels


  panelSelect(panels,ng,j)
  panelScale(rx,ry)

  # padj in axis needed to make grid line label close
  axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # bottom axis labels
  mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)   # bottom column titles

#__________________drawing loop________________

  for (i in 1:ng){
     gsubs = ib[i]:ie[i]
     ke = length(gsubs)

     pen = if(i==6)6 else 1:ke

     laby = ke:1 

     panelSelect(panels,i,j)
     panelScale(rx,c(1-pad,ke+pad))

     panelFill(col=Panel.Fill.col)

     axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid
     
     # if a refval is provided then add line.
     if(!is.na(refval))
        {
          lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
        }
     
     panelOutline(Panel.Outline.col)
     for (k in 1:ke){
        m = gsubs[k]
        if(good2[m])  # if upper value good.
          {
             lines(c(lower[m],upper[m]),rep(laby[k],2),
                   col=colors[pen[k]],lwd=Dot.conf.lwd)
          }
        if(good1[m])
          {
            if (doDotOutline) 
             {
               points(x[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,lwd=Dot.Outline.lwd,col=Dot.Outline.col,bg=colors[pen[k]])         
             } else {
               points(x[m],laby[k],pch=Dot.pch,cex=Dot.pch.size,col=NA, bg=colors[pen[k]])
             }
          }  
     }   
   }

# ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
             rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

 }

#####
#
# type = 'id' =======================================================
#
# rlStateId
#

rlStateId = function(j){
  #  j = panel column number
  
  py = Bar.barht*c(-.5,-.5,.5,.5,NA)
  px = c(.04,.095,.095,.04,NA)
  idstart = .137

#_____________ Scaling ______________________ 
 
  rx = c(0,diff(panels$coltabs[j+1,])) # width in inches
  ry = c(0,1)

#______________________panel labels_____________

  panelSelect(panels,1,j)      # start at I = 1, but j= is the current column.
  panelScale(rx,ry)
  
  mtext('U.S.',side=3,line=Title.Line.1.pos,cex=Text.cex)
  mtext('States',side=3,line=Title.Line.2.pos,cex=Text.cex)


# Cycle thought the GROUPS (ng)
  for (i in 1:ng){
     gsubs = ib[i]:ie[i]           # first element of group to last element of group.
     ke = length(gsubs)            # number of elements.
     laby = ke:1
     pen = if(i==6)6 else 1:ke
     panelSelect(panels,i,j)
     npad = ifelse(i==6,.57,pad)
     panelScale(rx,c(1-npad,ke+npad))
     gnams = stateNames[gsubs]
     polygon(rep(px,ke),rep(laby,rep(5,ke)) + py,col=colors[pen])
     polygon(rep(px,ke),rep(laby,rep(5,ke)) + py,col=Id.Dot.Outline.col,density=0)
     text(rep(idstart,ke), laby-Id.Text.adj, gnams, adj=c(0,0), cex=Text.cex,xpd=T)
  }

  # No reference values for this type of column
}


#####
#
#  State Rank Number ================================================================
#
#  rlStateRank   # based ID dot.
#      display the sorted rank.
#

rlStateRank = function(j){
  #  j = panel column number

  #________________ Scaling _______________

  rx = c(0,1)
  ry = c(0,1)
  rankstart = .137
 
  #______________________panel labels_____________

  panelSelect(panels,1,j)
  panelScale(rx,ry)
  mtext('Rank',side=3,line=Title.Line.1.pos,cex=Text.cex)
  # mtext('States',side=3,line=Title.Line.2.pos,cex=Text.cex)
 
  for (i in 1:ng){
     gsubs = ib[i]:ie[i]
     ke = length(gsubs)
     laby = ke:1
     pen = if(i==6)6 else 1:ke
     panelSelect(panels,i,j)
     panelScale(rx,c(1-pad,ke+pad))
     #gnams = stateNames[gsubs]
     #points(dotstart,laby,pch=Id.Dot.pch,col=colors[pen],cex=Dot.pch.size)
     #points(dotstart,laby,pch=1,col=Dot.Outline.col,cex=Dot.pch.size)
     Fgsubs <- formatC(gsubs,format="f",width=3,digits=0)
     text(rep(rankstart,ke),laby+.1,Fgsubs,adj=0,cex=Text.cex)
  }

  #  No reference values for this type of column.
}



#####
#
# type = 'map'  =========================================================
#
# rlStateMap
#

rlStateMap = function(j){

  # Works using state abbreviations
  # bnd.ord gives abbreviations in the
  #           the boundary are stored.
  # stateId give the abbreviations in the order plotted 

  bnd.ord = rlStateVisBorders$st[is.na(rlStateVisBorders$x)] # State abbrev
  rxpoly = range(rlStateVisBorders$x,na.rm=TRUE)
  rypoly = range(rlStateVisBorders$y,na.rm=TRUE)

  # ____________labeling and axes_______________
  
  panelSelect(panels,1,j)
  panelScale()
  par(xpd=T)
  
  mtext("Highlighted",side=3,line=Title.Line.1.pos,cex=Text.cex)
  mtext("States",side=3,line=Title.Line.2.pos,cex=Text.cex)

  # Drawing Loop

  for (i in 1:ng){

    if(i==6){                   # line break in maps.   Group 6 - middle group of 11.
      panelSelect(panels,6,j)
      panelScale()
      panelFill (col=Panel.Fill.col)
      
      panelOutline()
      text (.5,.55,'Median For Sorted Panel',cex= Text.cex*0.8)
      next  # skip to next FOR item
    }
    
    panelSelect(panels,i,j)     # Do map in - Panels by group...
    panelScale(rxpoly,rypoly)
    gsubs = ib[i]:ie[i]

    #  Add median state coloring to the row above and below the median line.
    if(i==5) gsubs = c(gsubs,26)  # slot 5 - add 26 to this group
    if(i==7) gsubs = c(gsubs,26)  # slot 7 - add 26 to this group
    #    26 (median) must always be the 6th name on the list to get painted black.
  
    gnams = stateId[gsubs]    # index to state id's (translation)
  
    # now find the state regions to plot
    
    back = is.na(match(rlStateVisBorders$st,gnams))   # list of states not involved.
    
    if(any(back)){
      polygon(rlStateVisBorders$x[back], rlStateVisBorders$y[back],
          density=-1, col=Map.Bg.col,  border=FALSE)         # fill in states
      polygon(rlStateVisBorders$x[back], rlStateVisBorders$y[back],
          density=0,  col=Map.Bg.Line.col, lwd=Map.Bg.Line.lwd)       # outline states
    }

    fore = !back    # reverse to list of states involved in this row. (by definition this is 5 states.
    
    pen = match(bnd.ord,gnams,nomatch=0)      
    pen = pen[pen>0]
    
    polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
          density=-1, col=colors[pen], border=FALSE)        # fill in states (1 to 5, 6)
    
    polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
          density=0,  col=Map.Fg.Line.col, lwd=Map.Fg.Line.lwd)     # outline states
  
    polygon(rlStateNationVisBorders$x, rlStateNationVisBorders$y,
          density=0, col=Map.Nation.Line.col, lwd=Map.Nation.Line.lwd)      # outside US boundary
  
    # might be made a function
    if (i==1)
      {
        text(135,31,'DC',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(22, 17,'AK',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(47, 8, 'HI',cex=Map.State.Spec.cex,adj=.5, col=1)
      }
  
  }  # i loop
  
  # no reference values for this type of column. If present - ignor.
}

#####
#
# type = 'mapcum'   ========================================================
#
# rlStateMapCum
#

rlStateMapCum = function(j){

  # Works using state abbreviations
  # bnd.ord gives abbreviations in the order the boundary are stored.
  # stateId give the abbreviations in the order plotted 

  bnd.ord = rlStateVisBorders$st[is.na(rlStateVisBorders$x)]   # State abbrev for states with boundaries
  rxpoly = range(rlStateVisBorders$x,na.rm=TRUE)
  rypoly = range(rlStateVisBorders$y,na.rm=TRUE)

  # ____________labeling and axes_______________

  panelSelect(panels,1,j)
  panelScale()
  
  
  box.x = rep(c(.14,.14,.208,.208,NA),2)-.04     
  par(xpd=T)
  y.ht = c(.05,.172)
  y.sep = .19*legfactor + 0.05            #  .185
  box.y = 1.025*legfactor +c(y.ht,rev(y.ht),NA) + 0.07 - 0.06  ## down 0.06

  polygon(box.x,c(box.y,box.y+y.sep),col=c(Map.Bg.col,colors[7]))
  polygon(box.x,c(box.y,box.y+y.sep),col=1,density=0)
  
  mtext("Cumulative Maps",side=3,line=Title.Line.1.pos,cex=Text.cex)
  mtext('States Featured Above',side=3,line=Title.Line.2.pos,at=.20,cex=Text.cex,adj=0)
  mtext('States Featured Below',side=3,line=lineTiclab,at=.20,cex=Text.cex,adj=0)

  # Drawing Loop

  for (i in 1:ng){

     if(i==6){
        panelSelect(panels,6,j)
        panelScale()
        panelFill (col=Panel.Fill.col)
        panelOutline()
        text (.5,.55,'Median For Sorted Panel',cex=Text.cex)
        next
     }
     panelSelect(panels,i,j)
     panelScale(rxpoly,rypoly)
     gsubs = ib[i]:ie[i]
     
     #  cont is list of states up to and including this row.
     if(i < 5)  cont = stateId[1:(5*i)] else cont = stateId[1:(5*i-4)]
                 # i = 1,   cont = stateId[1:5]
                 # i = 2,   cont = stateId[1:10]
                 # i = 7,   cont = stateId[1:31]  (5*7=35-4=31) (-4 compensates for the median row.
                 # i = 11,  cont = stateId[1:51]
                 
     if(i == 5) {gsubs = c(gsubs,26); cont = stateId[1:26]}
     if(i == 7) gsubs = c(gsubs,26) 

     gnams = stateId[gsubs]    # translate from sequence number to sorted order of states (abbrev)
                       # list of states in this row.

     # now find the state regions to plot
 
     back = is.na(match(rlStateVisBorders$st,cont))   # list of states not reference or used yet...
              # back = TRUE is no match between boundary file and cont - needing background color
 
     if (any(back))
       {
         polygon(rlStateVisBorders$x[back],rlStateVisBorders$y[back],
               density=-1, col=Map.Bg.col, border=F)        # fill in states with general fill.
         polygon(rlStateVisBorders$x[back], rlStateVisBorders$y[back],
               density=0,  col=Map.Bg.Line.col,lwd=Map.Bg.Line.lwd)   # outline states
       }

     fore = !back                                     # fore is a list of active states from other rows and this one..
     
     pen = match(bnd.ord,gnams,nomatch=0)             # match current states in row to boundary list, if found = TRUE, otherwize = 0 (no boundary).
                                                      # list is 51 long.
                                                      
     pen = ifelse(pen==0 & match(bnd.ord,cont,nomatch=0)>0, 7, pen)   # colors 1-6 or 7.
           # if pen=0 (not in this row) and in cont list (used), then use color 7.
     
     pen = pen[pen>0]                                 # colors - Categories (5), Black (1), used (1) = 7 total
     
     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=-1, col=colors[pen],   border=F)           # fill in states
     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=0,  col=Map.Fg.Line.col, lwd=Map.Fg.Line.lwd)      # outline states

     
     polygon(rlStateNationVisBorders$x,rlStateNationVisBorders$y,
        col=Map.Nation.Line.col,density=0,lwd=Map.Nation.Line.lwd)      # US outside boundary

     # might be made a function
     if(i==1){
        text(135,31,'DC',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(22,17,'AK',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(47,8,'HI',cex=Map.State.Spec.cex,adj=.5, col=1)
     }

   }  # i loop

  # no reference values for this type of column. If present - ignor.

}

#####
#
# type = 'mapmedian'  =================================================
#
# rlStateMapMedian
#

rlStateMapMedian = function(j){

  # Works using state abbreviations
  # bnd.ord gives abbreviations in the
  #           the boundary are stored.
  # stateId give the abbreviations in the order plotted
  # This MapMedian cream colors all states above and below the median state. 

  bnd.ord = rlStateVisBorders$st[is.na(rlStateVisBorders$x)] # State abbrev
  rxpoly = range(rlStateVisBorders$x,na.rm=TRUE)
  rypoly = range(rlStateVisBorders$y,na.rm=TRUE)

  # ____________labeling and axes_______________

  panelSelect(panels,1,j)
  panelScale()
  box.x = rep(c(.14,.14,.208,.208,NA),2)+.02   
  par(xpd=T)
  y.ht = c(.05,.172)
  y.sep = .19*legfactor + 0.05    # .185
  box.y = 1.025*legfactor +c(y.ht,rev(y.ht),NA) + 0.07 - 0.03
  
  polygon(box.x,c(box.y,box.y+y.sep),col=c(Map.Bg.col,colors[7]))
  polygon(box.x,c(box.y,box.y+y.sep),col=1,density=0)

  mtext("Median Based Contours",side=3,line=Title.Line.1.pos,cex=Text.cex)
  mtext('States In This Half',side=3,line=Title.Line.2.pos,at=.26,cex=Text.cex,adj=0)
  mtext('States In Other Half',side=3,line=lineTiclab,at=.26,cex=Text.cex,adj=0)

  # Drawing Loop

  for (i in 1:ng){

     if (i==6)
       {
         panelSelect(panels,6,j)
         panelScale()
         panelFill (col=Panel.Fill.col)
         panelOutline()
         text (.5,.58,'Median For Sorted Panel',cex=Text.cex)
         next   # exit for loop  
       }
     panelSelect(panels,i,j)
     panelScale(rxpoly,rypoly)
     gsubs = ib[i]:ie[i]

     if(i <= 5) cont = stateId[1:26] else cont = stateId[26:51]

     if(i == 5) gsubs = c(gsubs,26)
     if(i == 7) gsubs = c(gsubs,26) 

     #  gsubs = current state list
     #  cont  = state list to be colored cream.

     gnams = stateId[gsubs]

     # now find the state regions to plot
     back = is.na(match(rlStateVisBorders$st,cont))   # every state not active.
   
     if(any(back)){
          polygon(rlStateVisBorders$x[back],rlStateVisBorders$y[back],
             density=-1,col=Map.Bg.col,border=F)          # fill in states
          polygon(rlStateVisBorders$x[back], rlStateVisBorders$y[back],
             density=0, col=Map.Bg.Line.col,lwd=Map.Bg.Line.lwd) # outline states

     }

     fore = !back     # 
     
     pen = match(bnd.ord,gnams,nomatch=0)
     
     pen = ifelse(pen==0 & match(bnd.ord,cont,nomatch=0)>0, 7, pen)
     
     pen = pen[pen>0]

     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=-1,col=colors[pen], border=F)       # fill in states
     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=0, col=Map.Fg.Line.col, lwd=Map.Fg.Line.lwd) # outline states

     polygon(rlStateNationVisBorders$x,rlStateNationVisBorders$y,
        density=0, col=Map.Nation.Line.col,lwd=Map.Nation.Line.lwd)  # outside US boundary

     if(i==1){
        text(135,31,'DC',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(22,17,'AK',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(47,8,'HI',cex=Map.State.Spec.cex,adj=.5, col=1)
     }

   }   # i loop

  # no reference values for this type of column. If present - ignor.

}

#####
#
# type = 'maptail' ====================================================
#
# rlStateMapTail
#

rlStateMapTail = function(j){

  # Works using state abbreviations
  # bnd.ord gives abbreviations in the
  #           the boundary are stored.
  # stateId give the abbreviations in the order plotted
  # MapTail shows current states in a group as colored and
  # a tail of states (in cream color) from the outside inward.  
  # 

  bnd.ord = rlStateVisBorders$st[is.na(rlStateVisBorders$x)] # State abbrev
  rxpoly = range(rlStateVisBorders$x,na.rm=TRUE)
  rypoly = range(rlStateVisBorders$y,na.rm=TRUE)

  # ____________labeling and axes_______________

  # column header titles and "box"
  panelSelect(panels,1,j)    #  Line 1 and Line 2 - panel 1
  panelScale()
  # JP - as column labels move around or are repositioned,
  #    the associated BOX below does not follow it.
  box.x = c(.14,.14,.208,.208,NA)-.00 
  par(xpd=T)
  y.ht = c(.05,.172)
  y.sep = .185*legfactor
  
  box.y = 1.025*legfactor +c(y.ht,rev(y.ht),NA) + .07 - 0.01  # down 0.01
  polygon(box.x,box.y+y.sep,col=colors[7])
  polygon(box.x,box.y+y.sep,col=1,density=0)

  mtext("Two Ended Cumulative Maps",side=3,line=Title.Line.1.pos,cex=Text.cex)
  mtext('States Highlighted',side=3,line=Title.Line.2.pos,at=.25,cex=Text.cex,adj=0)
  
  #  JP - removed - temp
  #  mtext('Further From Median',side=3,line=lineTiclab,at=.15,cex=Text.cex,adj=0)

  # Drawing Loop

  for (i in 1:ng){

     if(i==6){
        panelSelect(panels,6,j)
        panelScale()
        panelFill (col=Panel.Fill.col)
        panelOutline()
        text (.5,.58,'Median For Sorted Panel',cex=Text.cex)
        next
     }
     panelSelect(panels,i,j)  
     panelScale(rxpoly,rypoly)
     # get list of states in this group.
     gsubs = ib[i]:ie[i]
     
     if(i < 5) cont = stateId[1:(5*i)]
     
     if(i==5)  {gsubs = c(gsubs,26); cont = stateId[1:26]}
     if(i==7)  {gsubs = c(gsubs,26); cont = stateId[26:51]} 
     
     if(i > 7) cont = stateId[(5*i-8):51]
     
     # get list of group state names 
     gnams = stateId[gsubs]

     # now find the state regions to plot
     #   plot states with cream filling (reported states)
     
     back = is.na(match(rlStateVisBorders$st,cont))
     
     if(any(back)){
         # paint fill
         polygon(rlStateVisBorders$x[back],rlStateVisBorders$y[back],
              density=-1,col=Map.Bg.col,border=F)          # fill in states
         # paint lines
         polygon(rlStateVisBorders$x[back], rlStateVisBorders$y[back],
              density=0, col=Map.Bg.Line.col,lwd=Map.Bg.Line.lwd) # outline states
     }

     fore = !back
     #  current 5 states with colors.
     
     pen = match(bnd.ord,gnams,nomatch=0)
     pen = ifelse(pen==0 & match(bnd.ord,cont,nomatch=0)>0,7, pen)
     pen = pen[pen>0]

     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=-1,col=colors[pen], border=F)        # fill in states
     polygon(rlStateVisBorders$x[fore], rlStateVisBorders$y[fore],
        density=0, col=Map.Fg.Line.col, lwd=Map.Fg.Line.lwd) # outline states

     #  The US border.
     polygon(rlStateNationVisBorders$x,rlStateNationVisBorders$y,
        density=0, col=Map.Nation.Line.col, lwd=Map.Nation.Line.lwd) # outside US boundary

     if(i==1){
        text(135,31,'DC',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(22,17,'AK',cex=Map.State.Spec.cex,adj=.5, col=1)
        text(47,8,'HI',cex=Map.State.Spec.cex,adj=.5, col=1)
     }

   }   #  i loop

  # no reference values for this type of column. If present - ignor.

}  

#
#

###################################################
#
#  For TS, and TSConf I could not find a way to  use to have stateIDs as the names of 
#  each state matrix, in list or data.frame.   So, the out at this time is
#  to assume the original panelData array is in the order of the original stateFrame data.frame.
#  When stateFrame is re-ordered, I have captured the re-ordering. Using the "order" index
#  the raw panelData is used via the order index to associate the line on the micromap to the data.
#   
#  Boxplot uses $names to look up to find out the record and link the Boxplot list to the 
#  stateFrame data.
#
#
#####


#####
#
# type = TS and TSConf   =====================================================
#
# rlStateTSConf  (Time Series with and without confidence interval in panel groups)
#
#     Plot all data for panel's states as one graph in panel.
#


rlStateTSConf = function(j,dataNam,conf=TRUE){
   #
   #  j = panel column number
   #
   #  dataNam = Name of large data array containing the x, y (or y low, med and high) values 
   #     for each time period and state.  Data element is three dimensions (State, sample, value)
   #     The state index is limited to 1:51.  The value index is limited ot 1:4.  
   #     The sample index is not limited, but a practical limit is around 200-250 samples.
   #
   #  conf = logical.  
   #    If TRUE, do the confidence band using y-low, y-med, and y-high values (columns 2, 3, 4)
   #    If FALSE, only plot the Y value (column 2)
   #
   
   ErrFnd = FALSE
   MsgLabel <- "TS"
   if (conf) MsgLabel <- "TSCONF"
   
   
   # Check data
   
   DataList = tryCatch(get(dataNam,pos=1),error=function(e) e)      # get name of array data object list.
   
   if (inherits(DataList,"error"))
     {
     
        warning(paste(MsgLabel,"-01 data.frame ",dataNam," does not exist or is not valid.",sep=""))
        ErrFnd = TRUE
        
     } else {
     
        # data.frame exists - can do other checks
        workDArr <- DataList
        wDArrNames <- rownames(workDArr)  # get rownames
  
        if (!is.array(workDArr))
          {
            warning(paste(MsgLabel,"-02 The data structured passed in the panelData field is not an array. Structure name = ",dataNam,sep=""))
            ErrFnd = TRUE
          }
   
        dimDArr <- dim(workDArr)
   
        if (dimDArr[1] != 51)
          {
            warning(paste(MsgLabel,"-03 The time serial array\'s 1st dimension is not 51 (states and DC). It is ",dimDArr[1],".",sep=""))
            ErrFnd = TRUE
          }
        #if (dimDArr[2] < 2 || dimDArr[2] > 31)
        if (dimDArr[2] < 2 )
          {
            #warning(paste(MsgLabel,"-04 The time serial array\'s 2nd dimension is not between 2 and 30 (number of time periods).  It is ",dimDArr[2],".",sep=""))
            warning(paste(MsgLabel,"-04 The time serial array\'s 2nd dimension must have at least 2 points / time periods.  It is ",dimDArr[2],".",sep=""))
            ErrFnd = TRUE
          }
  
        if (conf)   # TSCONF option.  
          {
            # Time Series with Confidence Bands
            if (dimDArr[3] !=4)
              {
                warning(paste(MsgLabel,"-05 The time series array\'s 3rd dimension is not 4.  It is ",dimDArr[3],",",sep=""))
                ErrFnd = TRUE
              }
       
          } else {
            # Time Series without Confidence Bands
            
            if (dimDArr[3] < 2)
              {
                warning(paste(MsgLabel,"-10 The time series array\'s 3nd dimension must be at least 2.  It is ",dimDArr[3],".",sep=""))
                ErrFnd = TRUE
              }
            
            if (dimDArr[3] != 2 && dimDArr[3] != 4)
              {  # accept confidence data - don't stop run.
                warning(paste(MsgLabel,"-06 The time series array\'s 3rd dimension must be 2 or 4. It is ",dimDArr[3],".",sep=""))
                ErrFnd = TRUE
              }
          }
  
        if (is.null(wDArrNames))
          {  # names are not present
               warning(paste(MsgLabel,"-11 The time series array does not have rownames assigned to the 1st dimension. Data cannot be paired up with stated.",sep=""))   
               ErrFnd=TRUE 
          } else {
            tnn <- is.na(match(wDArrNames,stateId))
            if (any(tnn))    # non-match found.
              {
                 lnn <- paste(wDArrNames[tnn],collapse=" ")
                 warning(paste(MsgLabel,"-07 Rownames on array do not match state ID list. The bad state IDs are:",lnn,sep=""))
                 ErrFnd = TRUE
              }
         } 
     }

   if (ErrFnd) return ()
  
   refval = lRefVals[j]              # get referrence to object, changed 
                                     #    to lRefVals - JP-2010/07/23
   reftxt = lRefTexts[j]             # new - JP-2010/07/23

   # structure of dataArr
   #     dataList is a 3 dim array :
   #          a * b * c, where: 
   #          a is the state index number (1 to 51) (area)
   #          b is the time period index (2 to "n" range) (Limited only by R and memory)
   #          c is the type of value (1=x, 2=low, 3=mid, 4=high) or (1=x, 2=y)
   #
      
   
   #_______________Scaling of TS Axis____________
   
   # x scaling
   rx <- range(workDArr[,,1],na.rm=TRUE)           # x range from all values in vector
   rx <- sc*diff(rx)*c(-.5,.5)+mean(rx)        # min to max range with expansion factors.
   
   # y scaling                  
   if (conf)
     {
        # range of line, high and low.
        ry <- range(workDArr[,,c(-1)],na.rm=TRUE)       # range of all Y values
     } else {
        # range of line.
        ry <- range(workDArr[,,2],na.rm=TRUE)           # range for the one Y value
     }
   ry = sc*diff(ry)*c(-.5,.5)+mean(ry)        # min to max range with expansion factors.
  
   
   #_______________Find range/min/max of median row line/high/low.____________
   
  
   #_______________Gather stats and put in State Order______________
  
   #
   #   JP-no data in col1, col2, or col3 to sort like the other columns... 
   #     All of the data is in these structures.
   #   
   #   at present no re-ordering of the time series like the other plots.
   #   JP-if other column is sorted, time series will follow that order via the indexes.
   #

   ####
   
   # ____________column titles and axis_______________

   # ___________  x Axis Labels ____________
   
   atRx = panelInbounds(rx)    # get labels for X axis
   
   # _____________top of column panel______ 
   
   panelSelect(panels,1,j)   # top panel - add title.
   panelScale(rx,ry)
   
   # _____________top titles________________
   
   mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)   # top column titles
   mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
  
   # ___________ top of column X axis ________________
   
   axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # top axis labels

   ####
   
   # ____________bottom of column panel____
   
   panelSelect(panels,ng,j)  # bottom panel - add sub title and refvals.
   panelScale(rx,ry)
   
   # ____________bottom of column X axis________________
   
   # padj in axis needed to make grid line label close
   
   axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx)) # bottom axis labels
   
   # ____________bottom title_________________
   
   mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex)  # bottom column titles

   ####
   
   # _______________drawing loop (panels 1->11)___________________

   oldpar = par(lend="butt")

   for (i in 1:ng)      #   1,2,3,4,5,    6,     7,8,9,10,11    ng=11
    {

      # Cycle through the Row/Groups in the micromap column
      
      gsubs = ib[i]:ie[i]               # get beginning to end index row number in group  
    
      if (i == 5) 
        {   # panel before the median row
           gsubs <- c(gsubs,ib[i+1]:ie[i+1]) # extend one more to get median row
        }
      if (i == 7)
        {   # panel after the median row
           gsubs <- c(gsubs,ib[i-1]:ie[i-1]) # extend to include at end of the list
        }
    
      ke = length(gsubs)                # get number of rows in group  (5 or 1)  

      gnams = stateId[gsubs]            # get list of state ids for data group of data.

      # adjust if middle group      
      pen = if(i==6) 6 else 1:ke        # if middle group (6), then pen=6 (Black), otherwise pen = c(1...5) or c(1...6)   
      
      #  do panel - 
      panelSelect(panels,i,j)           # select panel for group i in column j)
      panelScale(rx,ry)                 # set scale for panel  (should this be ry * 5 or 1?)
                                        # scale x and y to the shape of the panel (6 - median is squeezed.)
      panelFill(col=Panel.Fill.col)     # set fill for panel
      
      # draw grid lines in panel - vertical (x axis)
      axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid lines (x axis)
      
      if (i==6)  # median panel
        {
        
          # median panel (# 6)
          atRy = c(saveAtRy[1],saveAtRy[length(saveAtRy)])    # median panel range (Get copy of first and last number)  
        
        } else {
        
          # all other panels 
          atRy = panelInbounds(ry)   # get labels for y-axis
     
        }
      if (TS.hGrid)   # horizontal grids on Y axis
        {
           axis(side=2,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd, at=atRy) # Grid lines
        }
     
      # Y axis values and labels
      axis(side=2,tick=F,mgp=mgpLeft,cex.axis= TS.Axis.cex*.75 , at=atRy, labels=as.character(atRy)) # Y axis labels
      mtext(lab4[j],side=2,line=Title.Line.5.pos,cex=TS.Axis.cex)  # Y axis title
      
      #if(!is.na(refval))
      #   # add vertical line for reference.
      #   lines(rep(refval,2),c(1-padMinus,ke+padMinus), 
      #         lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
      
      panelOutline(col=Panel.Outline.col)     # outline panel
   
      #####
      # Issue with median row - line drawing.  The y axis is squeezed
      # to about 1/5 of the scale used in the other rows.  This distorts
      # the line graph and any confidence band.
      #####
      
      for (k in 1:ke)                  # Process each slot of panel - step 1 to 5 or 1 to 1

       {
         # cycle through row-groups and build each time series
         
         #m = datOrder[gsubs[k]]

         #m = ord[gsubs[k]]   # m is the location of the state in DataList$DArr matrix

         #if(is.na(m)) next   #     if no location - skip box plot for state

         kp = pen[k]         # color number
  
         wDArr <- workDArr[gnams[k],,]
         
         wX   <- wDArr[,1]    # get X values for line and polygon plots
	 wLine = wDArr[,2]    #  Get Y values for mid line 
         
         if (conf)
           {
             #  build polygon of confidence band to fill (y-low to y-high) and draw first.
            
             wArr <- wDArr[,c(3,4)]  # adjust data for plotting in the right slot.
             wPoly <- cbind(c(wX[1],wX,rev(wX)),c(wArr[1,1],wArr[,2],rev(wArr[,1])))
             # draw the polygon for the band with filling.
             #    filling is always current color + 7 -> transparent colors.
             polygon(wPoly,col=colors[kp+7],border=NA)
        
           }
    
         #  Plot mid Line
         lines(wX,wLine,col=colors[kp],lwd=TS.lwd)
       }
       
      saveAtRy = atRy
    }
   par(oldpar)
   # ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
         rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}

####
#
# type = 'segbar' and 'normbar'  ====================================
#
#  rlStateSegBar   (Segmented Bar chart)
#
#  Segmented bars is actually a stacked bar chart. Each segment is the length of one value.
#  The total length is the sum of the lengths of all segments.
#  The x scale of the column panels will be set to the "max" length of any bar.
#
#  In the normalized mode, the total for the segments is divided into value of each 
#  segment to get a percentage (0 to 100%).  The segments are then plotted as stacked
#  bars using the percentage.  The complete bar will be drawn from the left to right edge of 
#  the panel.
#
#  The data structure can have between 2 to 9 values per state.
#  Each state must have the same number of values. This limitation may be removed in the future.
#
#  Feature added to make each segment a different thickness. 1/4/2014
#
#  panelData => data.frame where each row is a state with the stateId as the row.name.
#     The columns are the bar segment values.

#

rlStateSegBar = function(j, SBnorm=FALSE) {
   #  j = the panel column number
   #  SBnorm  (FALSE = stacked,  TRUE = normalized)

   #   col1 indicates the starting or first column in the stateFrame data for bar segment values.
   #   col2 indicates the ending or last column in the stateFrame data.
   #
   #   The bar segment values are in the stateFrame for each state in columns "col1" to "col2".
   #
   wdim <- dim(dat)
   ErrFnd = FALSE
   
   if (is.na(match('col1',PDUsed))) {
    
      xmsg <- paste("SEGBAR-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
   }
   if (is.na(match('col2',PDUsed))) {
    
      xmsg <- paste("SEGBAR-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
   }
   if (!ErrFnd) {
  
     if (is.na(col1[j]) || col1[j] == 0 ) {
         
         warning(paste("SEGBAR-01 Specified column name or number in col1 for the first segment bar values is out of range, invalid, or does not exist: ",litcol1[j]," in stateFrame.",sep=""))
         ErrFnd = TRUE
     }
     if (is.na(col2[j]) || col2[j] == 0 ) {
         
         warning(paste("SEGBAR-02 Specified column name or number in col2 for the last segment bar values is out of range, invalid, or does not exist: ",litcol2[j]," in stateFrame.",sep=""))
         ErrFnd = TRUE
     }
   
   } 
   if (!ErrFnd) {
     
     if (col1[j] >= col2[j]) {
          
         warning(paste("SEGBAR-05 The first column number of bar values cannot be => the last column number.",litcol1[j],"::",litcol2[j],sep=""))
         ErrFnd = TRUE
     } else {
   
         wD =  ( col2[j] - col1[j] + 1 )   # corrected to calculate the number of data columns
         if (wD < 2 || wD > 9) {
              
             warning(paste("SEGBAR-06 The number of bar values for segmented/normalized bar plots must be between 2 and 9.",wD))
             ErrFnd = TRUE
         }
     }
   }

   if (ErrFnd) return ()    # error warning found - return
 
   workSB = dat[,c(col1[j]:col2[j])]  # bar segment values
   
   refval = lRefVals[j]              # get referrence to object, changed 
                                     #    to lRefVals - JP-2010/07/23
   reftxt = lRefTexts[j]             # new - JP-2010/07/23

   good = !is.na(rowSums(workSB))    # good values.
   
   #
   #
   # Colors - added transparency from x in steps of number of Segments up to 100%
   #   so 2 step = 50, 100
   #      3 step = 33.3, 66.6, 100
   #      4 step = 25, 50, 75, 100
   #      5 step = 20, 40, 60, 80, 100
   #      6 step = 16.6, 33.3, 50, 66,6, 83.3, 100
   #    etc.
   #    1/(NumSegs)*step = transparency
   #
   #   Dan's addition ==> 
   #    as the colors are generated from the base color
   #
   #    pInc = 1 / NumSegs
   #
   #    cSteps = cumsum(rep(pInc,NumSegs))^1.35
   #
   #    thickness = constant  vs.  very based on 2 to 9th segment
   #
   
   #_______________Gather stats and put in State Order______________
  
   #  Sorting has already been done - by stateId or value.
   #  The stateId list has therefore been re-ordered accordingly.  
   #  Reorder the DataList to match.  The assumption was that the input data order for the panelData 
   #  matched the order of the original data in the stateFrame.
   #
   
   workMatSB <- as.matrix(workSB)
   
   SBLen    = apply(workMatSB,1,length)  # get length of each row.
   SBLRange = range(SBLen,na.rm=TRUE)

   NumSegs  = SBLRange[2]                # number of segments (Max Length)
 
   SBBarPt <- cbind(rep(0,51),workMatSB)
   SBBarPt <- t(apply(SBBarPt,1,cumsum))
 
   #_______________Scaling____________
   
   # x scaling
 
   rMax  <- max(SBBarPt)
   if (SBnorm)
     {
       rx <- c(0,100)
     } else {
       rx <- c(0,rMax*1.02)
     }

   ry = c(0,1)
   
   pyPat = c(-0.5,-0.5,0.5,0.5,NA)
   py =  CSNBar.barht * pyPat     #  SNBar.barht = 2/3 (0.6667) (fixed)
       # py <- c( -1/3, -1/3, +1/3, +1/3, NA)
   
   # variable bar height calculations
   wYPdelta <- (CSNBar.Last.barht - CSNBar.First.barht)/(NumSegs-1)  # increment
   wYP1 <- CSNBar.First.barht - wYPdelta
      
   # _____________ Color Patterns _______________
   
   #  Build color patterns for all bar charts
   baseColors  <- t(col2rgb(colors[1:7]))
   bgColors    <- t(col2rgb("white"))
   
   #  New Way with lighter colors - but opaque 
   x1 = cumsum(rep(1/NumSegs,NumSegs))
   x2 = x1 ^ 1.9
   pInc = (x2 * .6) + .4
 
   # baseColors -- base 255...
   baseCol2 <- baseColors/255
  
   # baseCol2[Colors,RGB]
   baseCol3 <- sapply(pInc,function(x) baseCol2 * x)  # colors(1-7),  segment(1-5) for (Rgb=RED)
                                                      # colors(8-14), segment(1-5) for (Rgb=GREEN)
                                                      # colors(15-21),segment(1-5) for (Rbg=BLUE)

   # baseCol3[(Colors-Red,Colors-Grn,Colors-Blu),Segments]
   
   baseColMod <- array(baseCol3,c(7,3,NumSegs))
                         #   [x,,]   x = color (1-7)
                         #   [,,y]   y = segment (1-5)
                         #   [,z,]   z = RGB 1=RED, 2=GREEN, 3=BLUE
                         #
                         #   [1,2,3]   1 fills first, 2 fills next, 3 fills last.
                         #  
   
   pIncM <- 1-pInc
   bgCol2 <- bgColors/255
   bgCol3 <- sapply(pIncM,function(x) bgCol2 * x)   # [rgb,segment]
   bgColMod <- t(bgCol3)                            # [segment, rgb]
   #  bgColMod[Segments,RGB]   (Segment =5 ==> 0)p
   #   NumSegs, RGB value
   
   baseColRgb <- matrix(rep(0,7*NumSegs),nrow=7,ncol=NumSegs)
   #  baseColRgb[Colors, Segment]
      
   for (isg in 1:NumSegs)  # [,,isg]   Level
      {
        for (icl in 1:7)   # colors   [icl,,]
          {
            wC <- baseColMod[icl,,isg] + bgColMod[isg,]          
            baseColRgb[icl,isg] <- rgb(wC[1],wC[2],wC[3])
          }
      }
   #
   #  Resulting colors are in baseColRgb[color,segment]
   #
   #  Now I have a matrix of colors - [x,y] where
   #   x is the color base - 1 to 7 (we use 1 to 6).
   #   y is the level based on the number of segments = 1 : NumSegs
   #
   #   rows - color ID
   #   columns - segment 1:x
   
   # ____________titles and labeling axes_______________
   
   # _____________top of column______
   
   panelSelect(panels,1,j)   # top panel - add title.
   panelScale(rx,ry)

   mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)
   mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
   
   atRx = panelInbounds(rx)
   axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))
   
   # ____________bottom of column____
   
   panelSelect(panels,ng,j)  # bottom panel - add sub title and refvals.
   panelScale(rx,ry)
   
   # padj in axis needed to make grid line label close
   axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))
  
   mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex) 

   # ___________________drawing loop_____________________

   oldpar = par(lend="butt")
   
   #  build each panel for each stacked bar set.
   
     for (i in 1:ng)
      {
        gsubs = ib[i]:ie[i]               # get beginning to end index row number in this group  
        ke = length(gsubs)                # get number of rows in group  (5 or 1)  
        # adjust if median group      
        pen = if(i==6) 6 else 1:ke        # if median group (6)(black), then pen=6, otherwise pen = c(1...x)   
        laby = ke:1 
   
        panelSelect(panels,i,j)
        panelScale(rx,c(1-pad,ke+pad)) #   1 to 5 are the y values for each bar.
        panelFill(col=Panel.Fill.col)
 
        axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid
        
        # if a refval is provided then add line, before segments
        if(!is.na(refval))
           {
             lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
           }
        
        #
        #  Not checking "good" values provided.
        #
        
        #
        #  Process each state's line. 
        #
        for (k in 1:ke)
         {
           # cycle through row-groups and assign colors to associated states dots.
           m = gsubs[k]
           wX = SBBarPt[m,]            # Get Row of data.
         
           if (SBnorm) 
             {
                wX = wX / wX[NumSegs+1] * 100   # last segment value is in NumSegs + 1 to get last column (end point)
             }
   
           #wYP <- rep(laby[k],5)+py   # height of segment (laby[k] => center line of segbar)
           wYP <- rep(laby[k],5)       # height of segment (laby[k] => center line of segbar)
      
           # calculate box for each segment
           wYPht = wYP1
      
           for (ik in 1:NumSegs)
            {
              if (SNBar.varht)
                {
                  # variable height bar segments
                  
                  wYPht <- wYPht + wYPdelta
                  wYP2 <- wYP + (pyPat * wYPht)
                  #print(paste("Seg:",ik,"  wYP2:",wYP2,sep=""))
                  
                } else {
                  # fixed height bar segments
                  wYP2 <- wYP + py
                }
              val0 = wX[ik]     # start
              val1 = wX[ik+1]   # end position
              wXP <- c(val0,val1,val1,val0,NA)
              # good value - draw bars are polygons.  (why to polygon)
     
              polygon(wXP,wYP2,col=baseColRgb[pen[k],ik],lwd=CSNBar.Outline.lwd,border=CSNBar.Outline.col,lty=CSNBar.Outline.lty) 
              #polygon(wXP,wYP2,col=CSNBar.Outline.col,density=0)
            } # end of ik loop (plotting Segments)
            #
            if (SNBar.Middle.Dot)    # do we graph a middle dot on the row?
              {
                mY <- laby[k]      # get Y position
                # put dot on boundary if even number of segments or in middle of middle segment if odd.
                if ((NumSegs %% 2)==1)
                  {
                     # put dot in middle of middle segment.                 
                    mSeg <- NumSegs %/% 2 + 1
                    mX <- (wX[mSeg] + wX[mSeg+1])/2   # middle of segment
                  } else {
                    # put dot on border between two middle segments.                 
                    mSeg <- NumSegs %/% 2
                    mX <- wX[mSeg+1]
                  }
                if (SNBar.MDot.pch >= 21 && SNBar.MDot.pch <= 25)
                  #  treat filled and non-filled symbols differently - get close to same results.
                  #  with filled, fill is bg, col and lwd deal with border
                  #  with non-filled, fill is col, lwd deals with border using col.
                  {  # filled symbol
                     points(mX,mY,pch=SNBar.MDot.pch, cex=SNBar.MDot.pch.size, 
                          bg=SNBar.MDot.pch.fill,      # fill color  
                          col = SNBar.MDot.pch.border.col,    # border color 
                          lwd = SNBar.MDot.pch.border.lwd)
                  } else {
                     # non filled symbol
                     points(mX,mY,pch=SNBar.MDot.pch, cex=SNBar.MDot.pch.size, 
                          col = SNBar.MDot.pch.fill,   # fill and border color 
                          lwd = SNBar.MDot.pch.border.lwd)
                  }
              }  # end of Middle Dot drawing.
           
           
         }  # end of k loop   
        # finish up panel
        
        panelOutline(Panel.Outline.col)
 
      } # end of i loop
  
  par(oldpar)
  # ____________________________PanelOutline____________________

  groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
         rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}


####
#
# type = 'ctrbar'  ====================================
#
#  rlStateCtrBar   (Centered Bar chart)
#
#  The centered bars is a stacked bar chart with the middle segment centered on the "0" value
#  of the chart and extending 1/2 it's value in both directions (plus and minus).
#  The other bar segments are plotted to it's left and right as appropriate.
#
#
#  The data structure can have between 2 to 9 data values per state.
#  Each state must have the same number of values. This limitation may be removed in the future.
#
#  panelData => data.frame where each row is a state with the stateId as the row.name.
#     The columns are the bar segment values.

#

rlStateCtrBar = function(j) {
   #  j = the panel column number
   #  
   #   col1 and col2 indentify the starting column and ending column number in the stateFrame
   #   that contains the bar values for each state.
   #
   wdim <- dim(dat)
   ErrFnd = FALSE
   
   if (is.na(match('col1',PDUsed))) {
     
       xmsg <- paste("SEGBAR-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
       warning(xmsg)
       ErrFnd = TRUE
   }
   if (is.na(match('col2',PDUsed))) {
     
       xmsg <- paste("SEGBAR-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
       warning(xmsg)
       ErrFnd = TRUE
   }
   if (!ErrFnd) {
     
       if (is.na(col1[j]) || col1[j] == 0 ) {
         
            warning(paste("CTRBAR-01 Specified column name or number in col1 for the first segment bar values is out of range, invalid, or does not exist: ",litcol1[j]," in stateFrame.",sep=""))
            ErrFnd = TRUE
       }
       if (is.na(col2[j]) || col2[j] == 0 ) {
        
            warning(paste("CTRBAR-02 Specified column name or number in col2 for the last segment bar values is out of range, invalid, or does not exist: ",litcol2[j]," in stateFrame.",sep=""))
            ErrFnd = TRUE
       }
   } 
   if (!ErrFnd) {
    
       if (col1[j] >= col2[j]) {
          
          warning(paste("CTRBAR-05 The first column name or number of bar values must proceed the last column name or number.",litcol1[j],"::",litcol2[j],sep=""))
          ErrFnd = TRUE
       } else {
   
          wD = ( col2[j] - col1[j] + 1 )  # corrected to properly calculate the number of data columns.
          if (wD < 2 || wD > 9) {
               
              warning(paste("CTRBAR-06 The number of bar values for centered bar plots must be between 2 and 9.",wD))
              ErrFnd = TRUE
          }
       }
   }

   if (ErrFnd) return ()

   workCB = dat[,c(col1[j]:col2[j])]  # get bar segment data from the stateFrame.
   
   refval = lRefVals[j]              # get referrence to object, changed 
                                     #    to lRefVals - JP-2010/07/23
   reftxt = lRefTexts[j]             # new - JP-2010/07/23
   
   good = !is.na(rowSums(workCB))    # good values

   #
   # Colors - series of lighter colors of the base colors for each bar.
   #   Use an adjusted list of percentages based on the Number of Segments.
   #      2 step = 50, 100
   #      3 step = 33.3, 66.6, 100
   #      4 step = 25, 50, 75, 100
   #      5 step = 20, 40, 60, 80, 100
   #      6 step = 16.6, 33.3, 50, 66,6, 83.3, 100
   #    etc.
   #    1/(NumSegs)*step = transparency or lightness level  (100% full)
   
   #   Dan's addition ==> 
   #    as the colors are generated from the base color
   #
   #    pInc = 1 / NumSegs
   #
   #    cSteps = cumsum(rep(pInc,NumSegs))^1.35
   #
   #    thickness = constant  vs.  very based on 2 to 9th segment
   #

   #_______________Gather stats and put in State Order______________
  
   #  Sorting has already been done - by stateId or value.
   #  The stateId list has therefore been re-ordered accordingly.  
   #  Reorder the DataList to match.  The assumption was that the input data order for the panelData 
   #  matched the order of the original data in the stateFrame.
   #
   
   workMatCB = as.matrix(workCB)

   CBLen    = apply(workMatCB,1,length)  # get length of each row.
   CBLRange = range(CBLen,na.rm=TRUE)

   NumSegs  = CBLRange[2]                # number of segments

   CBBarPt  = cbind(rep(0,51),workMatCB)
   CBBarPt  = t(apply(CBBarPt,1,cumsum))
   
   #_____________Centering_____________

   CtrSeg   = as.integer(NumSegs/2) + 1  # center segment

   if ((NumSegs %% 2) != 0)
     {   # old number of segments
        CtrPt = workMatCB[,CtrSeg]/2 + CBBarPt[,CtrSeg]
     } else {
         # even number of segments
        CtrPt = CBBarPt[,CtrSeg]
     }

   CBPlotPts = CBBarPt - CtrPt

   #_______________Scaling____________
   
   # x scaling
   
   rx = range(CBPlotPts,na.rm=TRUE)
   rx = sc*diff(rx)*c(-.5,.5)+mean(rx)
 
   ry = c(0,1)
 
   pyPat = c(-.5,-.5,.5,.5,NA) 
   py =  CSNBar.barht * pyPat            #  CBar.barht = 2/3 (0.6667) (fixed)
       # py <- c( -1/3, -1/3, +1/3, +1/3, NA)
   
   # variable bar height calculations
   wYPdelta <- (CSNBar.Last.barht - CSNBar.First.barht)/(NumSegs-1)  # increment
   wYP1 <- CSNBar.First.barht - wYPdelta
      
   # _____________ Color Patterns _______________
   
   #  Build color patterns for all bar charts
   baseColors <- t(col2rgb(colors[1:7]))
   bgColors   <- t(col2rgb("white"))
   
   #  New Way with lighter colors - but opaque 

   x1 = cumsum(rep(1/NumSegs,NumSegs))
   x2 = x1 ^ 1.9
   pInc = (x2 * .6) + .4

   # baseColors -- base 255...
   baseCol2 <- baseColors/255

   # baseCol2[Colors,RGB]
   baseCol3 <- sapply(pInc,function(x) baseCol2 * x)  # colors(1-7),  segment(1-5) for (Rgb=RED)
                                                      # colors(8-14), segment(1-5) for (Rgb=GREEN)
                                                      # colors(15-21),segment(1-5) for (Rbg=BLUE)

   # baseCol3[(Colors-Red,Colors-Grn,Colors-Blu),Segments]
   
   baseColMod <- array(baseCol3,c(7,3,NumSegs))
                         #   [x,,]   x = color (1-7)
                         #   [,,y]   y = segment (1-5)
                         #   [,z,]   z = RGB 1=RED, 2=GREEN, 3=BLUE
                         #
                         #   [1,2,3]   1 fills first, 2 fills next, 3 fills last.
                         #  
   
   pIncM <- 1-pInc
   bgCol2 <- bgColors/255
   bgCol3 <- sapply(pIncM,function(x) bgCol2 * x)   # [rgb,segment]
   bgColMod <- t(bgCol3)                            # [segment, rgb]
   #  bgColMod[Segments,RGB]   (Segment =5 ==> 0)p
   #   NumSegs, RGB value
   
   baseColRgb <- matrix(rep(0,7*NumSegs),nrow=7,ncol=NumSegs)
   #  baseColRgb[Colors, Segment]
      
   for (isg in 1:NumSegs)  # [,,isg]   Level
      {
        for (icl in 1:7)   # colors   [icl,,]
          {
            wC <- baseColMod[icl,,isg] + bgColMod[isg,]          
            baseColRgb[icl,isg] <- rgb(wC[1],wC[2],wC[3])
          }
      }
   #
   #  Resulting colors are in baseColRgb[color,segment]
   #
   #  Now I have a matrix of colors - [x,y] where
   #   x is the color base - 1 to 7 (we use 1 to 6).
   #   y is the level based on the number of segments = 1 : NumSegs
   #
   #   rows - color ID
   #   columns - segment 1:x

   # ____________titles and labeling axes_______________
   
   # _____________top of column______
   
   panelSelect(panels,1,j)   # top panel - add title.
   panelScale(rx,ry)
   
   mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)
   mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
 
   atRx = panelInbounds(rx)
   atRxAbs = abs(atRx)
   #
   # Special Note - for the centered segmented bars, the values to the left and right are offsets
   #   and listed as positive -  NOT negative to left and positive to right.
   #
   axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRxAbs))
   
   # ____________bottom of column____
   
   panelSelect(panels,ng,j)  # bottom panel - add sub title and refvals.
   panelScale(rx,ry)
   
   # padj in axis needed to make grid line label close
   axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRxAbs))
  
   mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex) 

   # ___________________drawing loop_____________________

   oldpar = par(lend="butt")
 
   #  build each panel for each stacked bar set.
   
     for (i in 1:ng)
      {
        gsubs = ib[i]:ie[i]               # get beginning to end index row number in this group  
        ke = length(gsubs)                # get number of rows in group  (5 or 1)  
        # adjust if median group      
        pen = if(i==6) 6 else 1:ke        # if median group (6)(black), then pen=6, otherwise pen = c(1...x)   
        laby = ke:1 
   
        panelSelect(panels,i,j)
        panelScale(rx,c(1-pad,ke+pad)) #   1 to 5 are the y values for each bar.
        panelFill(col=Panel.Fill.col)
 
        axis(side=1,at=atRx, tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid
        
        # if a refval is provided then add line.
        if(!is.na(refval))
          {
            lines(rep(refval,2),c(1-padMinus,ke+padMinus),lty=Ref.Val.lty,lwd=Ref.Val.lwd,col=iRef.Val.col)
          }
        
        #
        #  Not checking "good" values provided.
        #
        
        #
        #  Process each state's line. 
        #
        for (k in 1:ke)
         {
           # cycle through row-groups and assign colors to associated states dots.
           m = gsubs[k]
           wX = CBPlotPts[m,]            # Get Row of data.
      
           #wYP <- rep(laby[k],5)+py
           wYP <- rep(laby[k],5)
      
           # calculate box for each segment
           wYPht = wYP1
      
           for (ik in 1:NumSegs)
            {
              if (CBar.varht)
                {
                  # variable height bar segments
                  
                  wYPht <- wYPht + wYPdelta
                  wYP2 <- wYP + (pyPat * wYPht)
                } else {
                  # fixed height bar segments
                  wYP2 <- wYP + py
                }
              val0 = wX[ik]     # start
              val1 = wX[ik+1]   # end position
              wXP <- c(val0,val1,val1,val0,NA)
              # good value - draw bars are polygons.  (why to polygon)
              
              polygon(wXP,wYP2,col=baseColRgb[pen[k],ik],lwd=CSNBar.Outline.lwd,border=CSNBar.Outline.col,lty=CSNBar.Outline.lty) 
              #polygon(wXP,wYP2,col=CSNBar.Outline.col,density=0)
            }
         }   
       # finish up panel
       
       # draw vertical line at zero.
       lines(rep(0,2),c(1-padMinus,ke+padMinus),
               lty=CBar.Zero.Line.lty,lwd=CBar.Zero.Line.lwd,col=CBar.Zero.Line.col)
       
       panelOutline(Panel.Outline.col)
 
      }
  
   par(oldpar)
   # ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
         rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

   
}

#####
#
# type = 'ScatDot'   =====================================================
#
# rlStateScatDot  (Scattered Plot Dots)
#

rlStateScatDot = function(j){
  #
  #  j = panel column number
  #
  #  col1 and col2 point to the X and Y data values in the stateFrame data.frame (known here as "dat").
  # 
  #
  wdim <- dim(dat)
  ErrFnd = FALSE
  
  if (is.na(match('col1',PDUsed))) {
     
      xmsg <- paste("SEGBAR-10 'col1' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (is.na(match('col2',PDUsed))) {
    
      xmsg <- paste("SEGBAR-11 'col2' vector is missing from the panelDesc data.frame.",sep="")
      warning(xmsg)
      ErrFnd = TRUE
  }
  if (!ErrFnd) {
    
      if (is.na(col1[j]) || col1[j] == 0 ) {
        
         warning(paste("SCATDOT-01 Specified column name or number in col1 for X values is out of range, invalid, or does not exist: ",litcol1[j]," in the data frame.",sep=""))
         ErrFnd = TRUE
      }
      if (is.na(col2[j]) || col2[j] == 0 ) {
         
         warning(paste("SCATDOT-02 Specified column name or number in col2 for Y values is out of range, invalid, or does not exist: ",litcol2[j]," in the data frame.",sep=""))
         ErrFnd = TRUE
      }
   }

   if (ErrFnd) return ()
  
   workSCD = dat[,c(col1[j],col2[j])]     # get x and y data from the stateFrame.
   colnames(workSCD) <- c("x","y")
   
   refval = lRefVals[j]              # get referrence to object, changed 
                                     #    to lRefVals - JP-2010/07/23
   reftxt = lRefTexts[j]             # new - JP-2010/07/23

   #_______________Gather stats and put in State Order______________
  
   #  Sorting has already been done of the stateFrame (dat) by stateId or value 
   #     in the function startup.
   
   #_______________Scaling____________
   
   # x scaling
   rx = range(workSCD$x,na.rm=TRUE)       # range of X values
   rx = SCD.xsc*diff(rx)*c(-.5,.5)+mean(rx)     # min to max range with expansion factors.
   
   # y scaling                  
   ry = range(workSCD$y,na.rm=TRUE)       # range of Y values
   ry = SCD.ysc*diff(ry)*c(-.5,.5)+mean(ry)
   
   # diagonal end points
   dx <- max(rx[1],ry[1])
   
   diagr <- c(max(rx[1],ry[1]), min(rx[2],ry[2]))
  
   # ____________titles and labeling axes_______________

 
   # _____________top of column______
   
   panelSelect(panels,1,j)   # top panel - add title - above.
   panelScale(rx,ry)         # scaled for data.
   
   mtext(lab1[j],side=3,line=Title.Line.1.pos,cex=Text.cex)
   mtext(lab2[j],side=3,line=Title.Line.2.pos,cex=Text.cex)
 
   atRx = panelInbounds(rx)
   axis(side=3,mgp=mgpTop,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))

   # ____________bottom of column____
   
   panelSelect(panels,ng,j)  # bottom panel - add sub title and refvals.
   panelScale(rx,ry)
   
   # padj in axis needed to make grid line label close
   axis(side=1,mgp=mgpBottom,padj=padjBottom,tick=F,cex.axis=Text.cex,at=atRx,labels=as.character(atRx))
   
   mtext(side=1,lab3[j],line=Title.Line.3.pos,cex=Text.cex) 

   # ___________________drawing loop_____________________

   # in the ordered list, the median should be 26 of 51 items.
     
   oldpar = par(lend="butt")

   #  build each panel for scatter plot dots
   
   for (i in 1:ng)   # groups from 1 to 5, 6, 7 to 11   ##  6 is the median group.
    {
      # Cycle through the Row/Groups in the micromap column

      # Set defaults for dots for this panel
     
      workSCD$bg  <- SCD.Bg.pch.fill    # default color   - was SCD.Bg.pch.fill
      workSCD$cex <- SCD.Bg.pch.size    # default size, except median
      workSCD$lwd <- SCD.Bg.pch.lwd     # default line weight
      workSCD$pch <- SCD.Bg.pch         # default pch code.
      
      if (i >= 5 && i <= 7) {    # force median dot to be highlighted in 5, 6 and 7 rows. 
        
          workSCD$bg[26]  <- SCD.Median.pch.fill
          workSCD$cex[26] <- SCD.Median.pch.size
          workSCD$lwd[26] <- SCD.Median.pch.lwd
          workSCD$pch[26] <- SCD.Median.pch
        }  
      gsubs = ib[i]:ie[i]               # get beginning to end index row number in this group  
      ke = length(gsubs)                # get number of rows in group  (5 or 1)  

      # adjust if median group      
      pen = if(i==6) 6 else 1:ke        # if median group (6)(black), then pen=6, otherwise pen = c(1...x)   
      
      panelSelect(panels,i,j)           # select panel for group i in column j)
      panelScale(rx,ry)                 # set scale for panel  (should this be ry * 5 or 1?)
      panelFill(col=Panel.Fill.col)            # set fill for panel
      
      # vertical grid lines.
      axis(side=1,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd) # grid lines
      
      # y axis labels
      if (i==6)  # median panel 
        {
          atRy = c(saveAtRy[1],saveAtRy[length(saveAtRy)])   # for margin panel, print the lowest and highest.   
        } else {
          atRy = panelInbounds(ry)                           # prettyprint a range.
        }
      #  horizontal grid.
      if (SCD.hGrid)
        {
           axis(side=2,tck=1,labels=F,col=Grid.Line.col,lwd=Grid.Line.lwd, at=atRy) # Grid lines
        }
      
      #  
      axis(side=2,tick=F,cex.axis=TS.Axis.cex*.75,mgp=mgpLeft,at=atRy,labels=as.character(atRy))
      mtext(lab4[j],side=2,line=Title.Line.5.pos,cex=TS.Axis.cex)
      
      panelOutline(col=Panel.Outline.col)     # outline panel

      # dv <- c(gsubs[1:ke],26)
      #
      # draw diagonal line of symetry from c(min (x, y),min(x,y)) to 
      #     c(max(x,y), max(x,y)), all point have x=y.
      #
      if ((diagr[1] < diagr[2]) && SCD.DiagLine)   
        {  # draw symetric line if within box range.
          dx = c(diagr[1],diagr[2])
          dy = c(diagr[1],diagr[2])
          lines(dx,dy, col=SCD.DiagLine.col, lwd=SCD.DiagLine.lwd, lty=SCD.DiagLine.lty)  # place a diagonal line on plot.
          #  print out the statistics for the line
          if (MST.Debug == 1)
            {
               print("line:")
               print(c(dx,dy))
               print("usr:")
               print(par("usr"))
               print("pin:")
               print(par("pin"))
               MST.Debug = 0 # turn off.
            }
        }      
    
      if (i == 6)
        {
          wS <- workSCD[gsubs[1],]      # get one entry - the median.
        
        } else {
          for (k in 1:ke)                  # Process each slot of panel - step 1 to 5 or 1 to 1
            {
               # cycle through row-groups and assign colors to associated state's dots.
               m = gsubs[k]
               workSCD$bg[m]   <- colors[pen[k]]       # set approvate color to circle fill.
               workSCD$cex[m]  <- SCD.Fg.pch.size
               workSCD$lwd[m]  <- SCD.Fg.pch.lwd
               workSCD$pch[m]  <- SCD.Fg.pch
            }
          wS <- workSCD[order(workSCD$cex,decreasing=FALSE),]
          # plot all points by size, others first, colored and median last.   
   
        }
       points(wS$x,wS$y,pch=wS$pch, col="black",bg=wS$bg, cex=wS$cex, lwd=wS$lwd)  # removed 
       saveAtRy <- atRy  # save for possible use on median panel.
    }
   par(oldpar)
   # ____________________________PanelOutline____________________

   groupPanelOutline(panelGroup,j)

  if(!is.na(refval)) 
         rlStateRefText(j,reftxt)  # added reftxt field - JP-2010/07/23

}

############################################


####
#
# ============================================================
# rlStateRefText  function inserts text lines at bottom of column as
#   legend.   Reference line has already been placed in each panel.
#
#   Changed-10/05/12- legend line overprints the texts.
# 

rlStateRefText= function(j,reftext){
   par(xpd=T)
   panelSelect(panels,ng,j)  # select one beyond the last panel (legend)
   panelScale()
   
   if (!is.null(reftext))
       if (!is.na(reftext)) {
   
          xt = (strwidth(reftext,cex=Text.cex) + .2)/2  # length of line and text
          xl = par("pin")[2]/2  # x axis
          xp = xl - xt
          # graphic space is 0 to 1 from left to right.
          # reference line is from 0.24->.37 and 0.63->.76 (x) at Title.Line.4.pos = y; 
          #   text inserted from .37->,63 (???)
        
          # way to find graphic length of string --> sw <- strwidth(reftext,cex=Text.cex)

          # print reference line for legend.
          #lines(c(.24,.37,NA,.63,.76),rep(Title.Line.4.pos,5),
          #      lty=Ref.Val.lty,col=iRef.Val.col,lwd=Ref.Val.lwd)
          #  10/10/12-changed above to print section of line on left of reference legend label.
          
          #text(.50,Title.Line.4.pos+.01,reftext,cex=Ref.Text.cex,col=iRef.Text.col,adj=.5)
          #  10/10/12-changed above to print reference label to right of dash line. Not centered.   
      
          # add text definition for legend.   (5/21/13 - added color to line)
    
          lines(c(xp, xp+.15),rep(Title.Line.4.pos,2), lty=Ref.Val.lty, lwd=Ref.Val.lwd, col=iRef.Val.col)
          text(xp + .2,Title.Line.4.pos+.01,reftext,cex=Text.cex,col=iRef.Text.col,adj=0)
       }
   # later add code to center line and text dependent on length of text (pixels)
}

#### end of micromap functions


#####  end of functions  #####


##############################################


#_________________ Check panel description content and formats _____________
#

if(!is.data.frame(panelDesc))
    stop("MST-10 Panel descriptor argument (2nd argument) must be a data.frame")

#______________Check for panelDesc$type validity______________

valid = c("map","mapcum","maptail","mapmedian",
          "rank","id","arrow","bar",
          "dot","dotse","dotconf",
          "ts","tsconf",
          "scatdot",
          "segbar","normbar","ctrbar",
          "boxplot")              # idDot and rank are not currently implemented
          
#____________________ List of expected and valid parameters in the panelDesc

PDParms <- c('type',
             'lab1','lab2','lab3','lab4',
             'col1','col2','col3',
             'refVals','refTexts',
             'panelData'
            )
      
PDUsed <- names(panelDesc)

PDPmatch <- match(PDUsed,PDParms)    # is if all entries in panelDesc are valid

if (any(is.na(PDPmatch))) {
   #  one of more panelDesc parameters are bad
   PDErrorList <- paste0(PDUsed[is.na(PDPmatch)],collapse=" ")
   xmsg <- paste0("PD-01 The following parameters in panelDesc are not valid:",PDErrorList,".  Execution stopped.")
   stop(xmsg)
}

#___________________the panelDesc parameters (column names) are good ___________

#________________type parameter
#

if (is.na(match('type',PDUsed))) {
   # Error 'type' parameter is not present
   xmsg <- 'PD-02 The required "type" parameter is not in the panelDesc data.frame.'
   stop(xmsg)
}

# get type vector

type = as.character(panelDesc$type)

# test contents of type vector for validity
subs = match(type,valid)

if (any(is.na(subs))) {
    PDErrorList <- paste0(type[is.na(subs)],collapse=" ")
    stop(paste("PD-03 The panelDesc data.frame has an incorrect panel type: ",PDErrorList))

}

#  Set up number of glyphics columns

numCol = nrow(panelDesc)

#
#_________________panelDesc$labx____________________
#

blank = rep('',numCol)  # empty vector for labels

# a NULL column cannot exist in a data.frame.  If the name is present, it exist!

if (is.na(match('lab1',PDUsed))) lab1 = blank else
              lab1 = as.character(panelDesc$lab1)

if (is.na(match('lab2',PDUsed))) lab2 = blank else
              lab2 = as.character(panelDesc$lab2)

if (is.na(match('lab3',PDUsed))) lab3 = blank else
              lab3 = as.character(panelDesc$lab3)

if (is.na(match('lab4',PDUsed))) lab4 = blank else
              lab4 = as.character(panelDesc$lab4)


#_________Save panelDesc Parameters in to namespace____________
#
   assign('lab1',lab1)
   assign('lab2',lab2)
   assign('lab3',lab3)
   assign('lab4',lab4)


# more panelDesc checks and setups after the function definitions.

#_______________________panelDesc$colx_____________________
#

NAList <- rep(NA,numCol)

# number of columns based on the presence of Descriptions for Column

  if (!is.na(match('col1',PDUsed))) {
   
     litcol1 <- panelDesc$col1
     col1 <- CheckColx(litcol1,"col1",wSFNameList,len_wSFnam)
   
  } else {
   
     litcol1 <- NAList
     col1 <- NAList
  
  }
  
  if (!is.na(match('col2',PDUsed))) {
    
    litcol2 <- panelDesc$col2
    col2 <- CheckColx(litcol2,"col2",wSFNameList,len_wSFnam)
    
  } else {
 
    litcol2 <- NAList
    col2 <- NAList
    
  }
  if(!is.na(match('col3',PDUsed))) {
   
    litcol3 <- panelDesc$col3
    col3 <- CheckColx(litcol3,"col3",wSFNameList,len_wSFnam)
 
  } else {
 
    litcol3 <- NAList
    col3 <- NAList

  }

#_____________panelDesc$refxxx________________
#

   if (!is.na(match('refVals',PDUsed)))
      assign('lRefVals',panelDesc$refVals)  else
      assign('lRefVals',NAList)

   if (!is.na(match('refTexts',PDUsed)))
      assign('lRefTexts',panelDesc$refTexts) else
      assign('lRefTexts',NAList)
      
      
#_____________panelDesc$panelData_______________
#

#  if present is the typeof correct ?

   if (!is.na(match('panelData',PDUsed))) { 
  
      wPanelData = panelDesc#panelData             # save pointer to panelD
  
   } else {
  
      wPanelData = NAList
  
   }
  
   assign('panelData',wPanelData)
   rm(wPanelData)
       
   #  


# ________________Detail defaults_______________________________
#

#  Process defaults into the local variables as before.

   wDetails <- micromapSTDefaults$details

   # dynamic assignment of defaults to individual variables in "micromapST"
   #  namespace.
   defNam = names(wDetails)
   for (i in 1:length(wDetails))
     {
       assign(defNam[i],wDetails[[i]])    # assign default values.
     }
   
   # All details names must be in the globalVariable call to be visible to CRAN checks.
   #
   # The defaults have been moved to the individual variables.
   # Keep the list of names around to be to verify user supplied names.
   #

   
# Now overlay with any values provided by user.

   #
   # dynamic assignment of detail data.frame to individual variables in the 
   #  "micromapST' namespace..
   #
   # Should I add code to verify names provided?
   #
   if (!(missing(details) || is.null(details) || is.na(details)))
     {
       if (typeof(details) == "list")
         {
           nam = names(details)                 # parse the details list into variable that can be
           nam_match = match(nam,defNam)
       
           for (i in 1:length(details))         #  referenced using the list's name.
             {
               if(is.na(nam_match[i]))
                 {
                    # invalid variable name in details
                    xmsg <- paste("MST-12 Invalid details variable name: ",nam[i], " in the details list. IGNORED.",sep="")
                    warning(xmsg)
                 } else {
                    # valid name
                    assign(nam[i],details[[i]])
                 }
             }
         } else {
             stop("MST-13 The details parameter is not a list.")
         }
     }
   
   # 
   # This is the code the rcmd check could not detect the scope of the detail$ variables.
   #


#_________________colors _______________________________________
#

#  Must do after completing the details list processing
#
#  Verify "colors=" argument
#
#  Second purpose is to set the graphics colors not in the "colors" vector to grays or colors.
#  
#

colFull = TRUE                  # control logical = TRUE doing Color, FALSE doing Greys
NoErrs = TRUE
doDotOutline = Dot.Outline

if (missing(colors) || is.null(colors) || is.na(colors))
  {
     colors = micromapSTDefaults$colors
  } 

if (typeof(colors) == "character")
  {
    if (length(colors) != 14)
      {
        if (length(colors) == 7)    # check for the basic colors.
          {
            # we have the basic 7 colors. Expand to the list of 14.
            colorlab <- names(colors)
            TransColors <- adjustcolor(colors,0.2)
            colors <- c(colors, TransColors)
            if (!is.null(colorlab))
              { names(colors) <- c(colorlab,paste("l_",colorlab,sep="")) }
        
          } else {
      
           if (length(colors) == 1)
             {
               if (colors == "bw" || colors == "greys" || colors == "grays")
                 {
                   #  set up the colors for printing in BW or Gray tones
                   
                   #  Get the main greys for the 5 colors (the middle 3-7 grays in the RColorBrewer scale.
                   #    and add the black for the median and a grey for the state highlight color.
                   xbw <- brewer.pal(name="Greys",8)
                   greyColors <- c(xbw[c(3:7)],"#000000","#E8E8E8")
                   
                   #  Build the transparent colors for the segmented bar charts.
                   TransColors <- adjustcolor(greyColors,0.2)
                  
                   #  Set up the grey color vector as requested.
                   colors <- c(greyColors,TransColors)
                   
                   #  Set up running parameters.
                   colFull = FALSE
                   doDotOutline = TRUE  # outline dots in dot glyphs.
                  
                 } else {
                   
                   warning("MST-15 A single value is provided for the colors. It''s not 'BW', 'greys', or 'grays' and no other valid is supported.  The default colors vector will be used.")
                   colors = micromapSTDefaults$colors
                
                }
             } else {
               warning("MST-16 The colors vector has the incorrect number of elements.  It must have 1 or 14 entries.")
             }
          }
      } else {
        # have 14 values in vector
       
      }
  } else {
    warning("MST-17 The colors vector does not contain character elements, Vector type is invalid.  The default colors are used.")
    colors = micromapSTDefaults$colors
  }

   assign("colors",colors)

#____ end of color check and adjustments.___



#___panelDesc_________User specificed column width processing and checking

###  Add check of column type to table of miniumal or statics column widths.
###  Must have details lists processed to do this.

plotWidth = par("pin")[1]   #  width in inches.

IdW = Id.width[1]

if(is.null(panelDesc$colSize)) {

     #  no colSize provided by User - create the default version.
     
     colSize = rep(0,length(type))             # set vector to zeros. Length equal the number of columns requested.

     # check for "map..." type panel columns
     loc = substring(type,1,3)=='map'
     # was "map" the start of the type name for column?
     if(any(loc))  colSize[loc] = Map.width  # set size for map for columns doing maps.

     # check for "id" type panel columns
     loc = type=='id'   # is column type = "id"
     if(any(loc))
       {
         sub = ifelse(plotNames=="full",1,2)  # yes, set size for ID (ab or full)
         colSize[loc] = Id.width[sub]
         IdW = Id.width[sub]
       }
       
     # check for "rank" type of panel columns
     loc = type == 'rank'  # is column type = 'rank'
     if(any(loc)) colSize[loc] = Rank.width
     
     # Get plot width and calculate size of each remaining column.
     #   Assume equal width for each non-id or non-map column.
     #
     equalWidth= (plotWidth-sum(colSize))/sum(colSize==0)
     colSize = ifelse(colSize==0,equalWidth,colSize)
  } else {
     colSize = panelDesc$colSize
  }

#
#if (sum(colSize) >= plotWidth)
#  {
#     warning("The sum of colSize vector provided in the panelDesc argument is greater then the plotting area width.")
#  }
# did not work - try again later..
#

#  
#  Make adjustments for color or grays
#

if (colFull)
  {  
    # set color values to work variables
    iRef.Val.col = Ref.Val.col
    iRef.Text.col = Ref.Text.col
    
  } else {
  
    # set gray values to work variables
    iRef.Val.col = Ref.Val.BW.col
    iRef.Text.col = Ref.Text.BW.col
  }


  
# __________________________layout

  cparm <- data.frame(cSize=numeric(0),lSep=numeric(0),rSep=numeric(0))
  
  ncol = length(type)
  
  for (j in 1:ncol)
    {
      # Test type of column to be built and call build routine.
    cparm2 =  switch(type[j],
                  #  colSize, col left sep, col right sep
         "map"=      c(Map.width,0,0),
         "mapcum"=   c(Map.width,0,0),
         "maptail"=  c(Map.width,0,0),
         "mapmedian"=c(Map.width,0,0),
         "id"=       c(IdW,0,0),
         "dot"=      c(0,0,0),
         "dotse"=    c(0,0,0),
         "dotconf"=  c(0,0,0),
         "arrow"=    c(0,0,0),
         "bar"=      c(0,0,0),
         "segbar" =  c(0,0,0),
         "boxplot"=  c(0,0,0),
         "ts" =      c(0,.2,0),
         "tsconf" =  c(0,.2,0),
         "scatdot" = c(0,.2,0),
         "segbar"  = c(0,0,0),
         "normbar" = c(0,0,0),
         "ctrbar"  = c(0,0,0),
         "rank"    = c(Rank.width,0,0),
         "nomatch" = c(0,0,0)
      )
     cparm <- rbind(cparm,cparm2)
    }
   colnames(cparm) <- c("cSize","lSep","rSep")

   #topMar = details[["topMar"]]
   #botMar = details[["botMar"]]
   legfactor=1

   # add space if reference values provided.
   # JP-2010/07/23 0 change to refVals to be consistent.
   if(!is.null(panelDesc$refVals)){
      # if not null field.
      if(any(!is.na(panelDesc$refVals))){
         # value provided - not NA
         #botMar = details[["botMarLegend"]]
         # revisit calculation below to be more precise
         #legfactor= 9/(9-details[['botMardif']])
      
         botMar = botMarLegend
         # revisit calculation below to be more precise
         legfactor= 9/(9-botMardif)
      }      
   }
   
   
   
   assign('legfactor',legfactor,sys.frame(which = -1))  # set legfactor in environment -1 (caller's space.)

   ncol   = length(type)
   
   
   
   ladj    = c(0,cparm[2:ncol,2],0)          #  0, x2, x2, x2, x2, ... , x2, 0   for ncol + 2
   radj    = c(0,cparm[1:ncol-1,3],0)        #  0, y
   
   colSep = c(0,rep(.1,ncol-1),0)  + ladj + radj  # column separators = 0.1 
                                        # in 4 columns minimum -> c(0,.1,.1,.1,0) Side don't get sep.
                                        
   
   # build panels from panelLayout and pieces of rlStateDefaults$Details
       # nrow = 11 -> 5,5,5,5,5,1,5,5,5,5,5 states = 11 groups
   # individual panels (rows(11) and columns)
   assign("panels",panelLayout(nrow=11,ncol=ncol,
                        topMargin=topMar,                    # 0.95
                        leftMargin=0,                      
                        bottomMargin=botMar,                 # 0.5
                        #rowSep=details[["rowSep"]],         # c(0,0,0,0,0,.1,.1,0,0,0,0,0)
                        #rowSize=details[["rowSize"]],       # c(7,7,7,7,7,1.65,7,7,7,7,7)
                        rowSep=rowSep,                       # c(0,0,0,0,0,.1,.1,0,0,0,0,0)
                        rowSize=rowSize,                     # c(7,7,7,7,7,1.65,7,7,7,7,7)
                        colSize=colSize,                     # calculated colsizes (???)
                        colSep=colSep))                      # c(.1,.1,.1) for 3

   #grounpedRowSize = details[["groupedRowSize"]]            # c(35,1.65,35)
   #groupedRowSep   = details[["groupedRowSep"]]             # c(0,.1,.1,0)

   # Major panel group  title-top, panels, title-bottom  by columns (overlays panels)
   # section of panels (top(25), median(1), bottom(25) and "N" columns wide.
   assign("panelGroup",panelLayout(nrow=3,ncol=ncol,
                        topMargin=topMar,
                        leftMargin=0,
                        bottomMargin=botMar,
                        rowSize=groupedRowSize,
                        rowSep=groupedRowSep,
                        colSize=colSize,
                        colSep=colSep))

   # Page panel - top, middle, bottom - 1 across.                   (overlays panels)
   #   One column wide for each group (top, median, bottom)
   assign("panelOne",panelLayout(nrow=3,ncol=1,
                        topMargin=topMar,
                        leftMargin=0,
                        bottomMargin=botMar,
                        rowSize=groupedRowSize,
                        rowSep=groupedRowSep))

#####
# ____________________Main loop______________________________
#####

#  Build images of each column

   for (j in 1:ncol)
    {
      # Test type of column to be built and call build routine.
      switch(type[j],
         "map"=      rlStateMap(j),
         "mapcum"=   rlStateMapCum(j),
         "maptail"=  rlStateMapTail(j),
         "mapmedian"=rlStateMapMedian(j),
         "id"=       rlStateId(j),
         "dot"=      rlStateDot(j),
         "dotse"=    rlStateDotSe(j),
         "dotconf"=  rlStateDotConf(j),
         "arrow"=    rlStateArrow(j),
         "bar"=      rlStateBar(j),
         "boxplot"=  rlStateBoxplot(j,  as.character(panelDesc$panelData[j]) ),
         "ts" =      rlStateTSConf(j,   as.character(panelDesc$panelData[j]),conf=FALSE),
         "tsconf" =  rlStateTSConf(j,   as.character(panelDesc$panelData[j]),conf=TRUE),
         "scatdot" = rlStateScatDot(j),
         "segbar"  = rlStateSegBar(j),
         "normbar" = rlStateSegBar(j,SBnorm=TRUE),
         "ctrbar"  = rlStateCtrBar(j),
         "rank"    = rlStateRank(j),
         "nomatch"
      )
    }

 
      # All columns are built and sitting in the panel.
      panelSelect(panelOne,margin="top")    # full page top label area.
      panelScale()
 
      if(length(title)==1){
         text(.5,.77,title,cex=Title.cex)
      } else {
         text(.5, .9,title[1],cex=Title.cex)
         text(.5,.65,title[2],cex=Title.cex)
   }
   
} # end of micromapST Function

###  End of micromapST


###############################################################

###
#
#  micromapSTSetDefaults function
#
#  Must be run once to generate the default lists.
#  If you customize - then make a copy and change the copy.
#
#  Call by .onload at package load.  Reference is exported to globlal space for user's access.
#
###

micromapSTSetDefaults = function()
   {
 
#
#  build micromapSTDefaults data.frame so it can be exported.
#

# Candidate colors________________________________________
colorsRefRgb = matrix(c(
 1.00,1.00,1.00,  # white
  .92, .92, .92,  # lighter gray            # changed from .90
  .78, .78, .78,  # light gray              # changed from .80
  .50, .50, .50,  # middle gray
  .30, .30, .30,  # dark gray
  .00, .00, .00,  # black
 
  .93,1.00, .93,  # light green
  .00, .50, .00,  # mid green
 1.00,1.00, .84,  # light yellow foreground  
  .90, .80,1.00,  # bright yellow foreground 
  .80, .90,1.00,  # light green blue
  .60, .70, .85), # mid green blue
  ncol=3,byrow=TRUE)

colorsRef = grDevices::rgb(colorsRefRgb[,1],colorsRefRgb[,2],colorsRefRgb[,3])
names(colorsRef) = c("white","lighter gray","light gray",
                     "mid gray","dark gray", "black",

                     "light green","mid green",
                     "light yellow","bright yellow",
                     "light green blue","mid green blue")           

# Region colors________________________________________________

colorsRgb = matrix(c(                       # the basic 7 colors.
 1.00, .15, .15,  #region 1: red
  .90, .55, .00,  #region 2: orange
  .00, .65, .00,  #region 3: green
  .20, .50,1.00,  #region 4: greenish blue
  .50, .20, .70,  #region 5: lavendar 
  .00, .00, .00,  #region 6: black for median
 1.00,1.00, .80), #non-highlighted foreground
  ncol=3,byrow=TRUE)

colors = c( grDevices::rgb(colorsRgb[,1],colorsRgb[,2],colorsRgb[,3]),            # solid colors
            grDevices::rgb(colorsRgb[,1],colorsRgb[,2],colorsRgb[,3],.2))         # translucent colors.

names(colors) =c("red","orange","green","greenish blue", "purple","black","light yellow",
                 "l_red","l_orange","l_green","l_greenish blue", "l_purple","l_black","l_light yellow")       


# Details variable list _________________________________________

## JP added temp variables so function would read in in R 2.7
#      cannot use values within the details list since it's not really built yet.

tempne      = 5                           # number of states per panel
tempGrid.Line.col = colorsRef["white"]          # grid line color

tempcolFill = colorsRef["lighter gray"]   # panel and default fill color
tempText.cex = .7


details = list(

# panel layout grouping 
    ne = tempne,                   # number of item per group 1
    ng = ceiling(51/tempne),       # number of groups of states 2 
    ib =  c(1, 6,11,16,21,26,27,32,37,42,47), #group lower index 3
    ie =  c(5,10,15,20,25,26,31,36,41,46,51), #group upper index 4

# panel layout margin allocation
    # JP - changed median row size to 1.5.
    topMar       = 0.95,               # margin panel height (inches) 5
    botMar       = 0.5,                # no legend bottom margin 6
    botMarLegend = 0.5,                #                         7
    botMardif    = 0.2,                # maybe not needed 8
    leftMarAxis  = 0.2,                # left margin adjustment when Y axis is printed 9.
    #                1 2 3 4 5   6   7 8 9 10 11
    rowSep       = c(0,0,0,0,0,0.1,0.1,0,0,0,0,0),     # 10 (inches)   
    #                1 2 3 4 5 6 7 8 9 10 11           # row seperators (inches)
    rowSize      = c(7, 7, 7, 7, 7, 1.65, 7, 7, 7, 7, 7),  # JP change 1.5 to 1.65 on Median strip 11.
    #                 1-5 6 7-11                       # working units
    groupedRowSize = c(35, 1.65, 35),                  # JP changed 1.5 to 1.65 to give median a little more room  12
    #               1-5 5-6 6-7 7-11                   # working units.
    groupedRowSep  = c(0,0.1,0.1,0),                   # 13 - rowGroup separators (inches)

# panel column width allocation
             ## JP changed map width to 1.4
    Map.width    = 1.4,                # map width should be set portionally to the height of the panel (working units)
    Id.width     = c(0.9,0.30),        # full and ab 15 column width (working units)
    Rank.width   = 0.25,               # rank width of column   (working units)
    YAxis.width  = 0.2,                # width for Y axis labels when used. (working units? or inches?)
    
# panel scaling
    sc       = 1.08,                   # x and y axis scale expansion factor               16
    pad      = 0.67,                   # y axis padding for integer plotting locates       17
                                       # ry = c(1-pad,ke+pad),ke = no. items in panel
    padex    = 0.34,                   # total panel padding                               18
                                       # (.67-.5)=.17 padding at top and bottom of panel
    padMinus = 0.63,                   # .67 - .04 # keep reference line off panel edge    19

# mtext line placement (Titles)
    ##  JP adjusted placement of lines (titles)
    Title.Line.1.pos  = 1.75,          # top panel 1st line placement       20
    Title.Line.2.pos  = 1.05,          # top panel 2nd line placement       21
    Title.Line.3.pos  = 0.65,          # bottom panel line placement        22
    Title.Line.4.pos  = -0.7,          # reference line (below panel)       23
    Title.Line.5.pos  = 0.40,          # Y axis titles for ScatDot and TS.  24
    Title.cex         = 1.0,           #                                    25
    lineTiclab        = 0.2,           # lowest line for map legend text    26
  
# grid line parameters
    Grid.Line.col     = tempGrid.Line.col,  # grid line color               27
    Grid.Line.lwd     = 1,             # weight of grid line                28
    mgpTop            = c(2,0.1,0),    # gridline (tick) placement          29
    mgpBottom         = c(2,0,0),      # gridline (tick) placement          30
    padjBottom        = -0.7,          # gridline (tick  placement          31
    mgpLeft           = c(0.75,0.1,0), # left axis labels                   32.

# panels
    Panel.Fill.col    = tempcolFill,   # panel fill color                   33
    Panel.Outline.col = colorsRef["black"],     # panel outline color       34

# Title and Text - cex for character size
    Text.cex          = tempText.cex,  ## JP decreased text size.  Used almost everywhere. 35

# refVals parameters

    # see padMinus above for other parameters 
    Ref.Val.lty       = "dashed",                     # line type for Ref Value (dashed)  36
    Ref.Val.lwd       = 1.5,                          # line weight for Ref Value         37
    Ref.Val.col       = colorsRef["mid green"],       # line color                        38
    Ref.Val.BW.col    = colorsRef["black"],           # line color when "grays"           39

# refText parameters
    Ref.Text.cex      = tempText.cex,                 # Ref Text Size                     40
    Ref.Text.col      = colorsRef["black"],           # JP 10/10/12-changed from black to mid green.  # 41
                                          #  5/21/13 - changed back to black.
    Ref.Text.BW.col   = colorsRef["black"],           # Ref Text color when "grays"       42

#__________________________________________________________ 
# working parameters for each panel graphing subfunction within micromapST

# arrow plot parameters
    Arrow.lwd         = 2.5,                          ## JP decrease arrow width.  44
    Arrow.cex         = .08,                          # Not Used   45
    Arrow.Head.length = .08,                          #  Length of arrow head in inches.  43
    Arrow.Shadow.col  = colorsRef["black"],           # Not Used.  46
    Arrow.Shadow.lwd  = 4.0,                          # Arrows shadow when border needed. (Not used) 47
 
# bar plot parameters
    Bar.barht                  = 2/3,                 # fraction of line spacing  48
    Bar.Outline.col            = colorsRef["black"],  #  49
    Bar.Outline.lwd            = .5,                  #  50
    Bar.Outline.lty            = "solid",             #  51

# segmented bar parameters - segbar and normbar only
    SNBar.varht                = FALSE,               #  59   (default fixed height)
    SNBar.two.ended            = FALSE,               #  62   (not implemented)
    SNBar.Middle.Dot           = FALSE,               #  draw dot in middle point of segmented bars (default - no mid-poing dot)
    SNBar.MDot.pch             = 21,                  #  middle point symbol  64
    SNBar.MDot.pch.fill        = colorsRef["white"],  # middle point symbol fill/color    65
    SNBar.MDot.pch.size        = 0.6,                 # middle point symbol size          66
    SNBar.MDot.pch.border.lwd  = 0.6,                 # middle point symbol border lwd    67
    SNBar.MDot.pch.border.col  = 'black',             # middle point symbol.border.col with using filled symbols 68

# segmented bar parameters - centered bar only
    CBar.Zero.Line.col         = colorsRef["white"],  #  52
    CBar.Zero.Line.lwd         = 1,                   #  53
    CBar.Zero.Line.lty         = "dotted",            #  54
    CBar.varht                 = FALSE,               #  69   (default = fixed height)
    CBar.two.ended             = FALSE,               #  70   (not implemented)
    #  CtrSeg uses the rest of the segbar parameters.

# segmented bar parameters for all (segbar, normbar and ctrbar)
    CSNBar.barht               =  2/3,                #  bar heights (percentage of row) 55
    CSNBar.Outline.col         = colorsRef["black"],  #  bar outline border color 56
    CSNBar.Outline.lwd         = .75,                 #  bar outline border width 57
    CSNBar.Outline.lty         = "solid",             #  bar outline border type  58
                                                      # parameters when variable height is requested.
    CSNBar.First.barht         = 0.3333,              #  60
    CSNBar.Last.barht          = 0.80,                #  61

    
# box plot parameters
    BoxP.thin                  =.2,                   # was .29     ## JP decreased line width   71
    BoxP.thick                 =.60,                  # was .58                                  72
   
    BoxP.Use.Black             = FALSE,               # FALSE = Use the Color for outliners;  TRUE = use black  73
  
    BoxP.Median.Line           = 0.80,                #                                  74
  
    BoxP.Median.Dot.col        = colorsRef["white"],  ## JP changed from colDotMedian for clarity  75
    BoxP.Median.Dot.pch        = 19,                  #  76
    BoxP.Median.Dot.cex        = 0.95,                #  77
    BoxP.Median.Dot.lwd        = 2,                   #  78
    BoxP.Median.col            = colorsRef["black"],  ## JP changed to BoxP.Median.col from colMedian - was duplicate - set to black  79.
  
    BoxP.Outline.col           = colorsRef["dark gray"],  # color for outline when using colors 80
  
    BoxP.Outlier.lwd           = 0.4,                 ## JP decreased dot border line width   * 81
    BoxP.Outlier.cex           = 0.7,                 # see Dot.pch.size  ## JP decreased dot size  (was .6)   * 82
    BoxP.Outlier.BW.col        = colorsRef["dark gray"], # color for outline when using BW mode 83
  
# dot plot parameters (dot, dotconf, dotse)
    Dot.pch                    = 21,                  # plotting character  (1 circle, 16 dot, 21 filled circle)  # 84
    Dot.pch.size               = 0.9,                 # dot size            ## JP adjusted dot size.              # 85
                                                      # use default lwd, lty, and border color (.col = black)     
    Dot.Outline                = FALSE,               ## JP added option to control Dot outline.     # 89
    Dot.Outline.col            = colorsRef["black"],                                                     # 90
    Dot.Outline.lwd            = 0.5,                                                                    # 91
 
# dot conf parameters
    Dot.conf                   = 95,                  # % confidence interval                        # 86
    Dot.conf.lwd               = 2,                                                                  # 87
    Dot.conf.size              = 0.55,                # Not Used                                     # 88
                                                      # use default lty, and border color (.col = black)   
 
# ts and tsconf parameters
    TS.lwd                 = 1.1,                     # TS Line weight                                       # 92
    TS.Axis.cex            = tempText.cex * .7,                                                      # 93
    TS.hGrid               = FALSE,                                                                  # 94
    
# scatdot parameters 
    SCD.Bg.pch             =  21,                     # Background symbol pch                            95
    SCD.Bg.pch.lwd         =  0.6,                    # Background symbol border line weight             96
    SCD.Bg.pch.size        =  0.75,                   # Background symbol size                           97
    SCD.Bg.pch.fill        =  'transparent',          # Background symbol fill (bg) color (19:25)               98
    SCD.Fg.pch             =  21,                     # Foreground symbol pch                            99
    SCD.Fg.pch.lwd         =  0.6,                    # Foreground symbol border line weight            100
    SCD.Fg.pch.size        =  1,                      # Foreground symbol size                          101
    SCD.Median.pch         =  21,                     # median symbol PCH value (21 = filled circle)    102
    SCD.Median.pch.lwd     =  0.6,                    # median symbol border line weight                103
    SCD.Median.pch.size    =  1,                      # median symbol border size (cex)                 104
    SCD.Median.pch.fill    = colorsRef["black"],      # median median symbol fill color.                105
    SCD.Axis.cex           = tempText.cex * .7,       # not used                                        106
    SCD.xsc                = 1.08,                    # fudge for margins to try and not clip circles.(not used) 107
    SCD.ysc                = 1.12,                    # fudge for margins to try and not clip circles.(not used) 108
    SCD.hGrid              = FALSE,                   # draw horizontal grid.                            109
 
    SCD.DiagLine           = TRUE,                    #    110
    SCD.DiagLine.col       = tempGrid.Line.col,       #    111
    SCD.DiagLine.lwd       = 1.25,                    #    112
    SCD.DiagLine.lty       = "solid",                 #    113

# id State Dot parameters (link - State Lab and Dot)
    Id.Text.cex            = .9,                      ## JP decreased ID text size.     114
    Id.Text.adj            = 1/3,                     # subtract to lower text baseline for state names.  115
    Id.Dot.pch             = 21,                      # rlStateIdDot              116
    Id.Dot.cex             = .6,                      # Not Used.                      117
    Id.Dot.Outline.col     = colorsRef["dark gray"],                  # 118

# map parameters  
    Map.Bg.col             = tempcolFill,             # map/state background fill color    119
    Map.Bg.Line.col        = tempGrid.Line.col,       #   120
    Map.Bg.Line.lwd        = 1,                       #   121
    Map.Fg.Line.col        = colorsRef["black"],      #   122
    Map.Fg.Line.lwd        = 1,                       #   123
    Map.Nation.Line.col    = colorsRef["black"],      #   124
    Map.Nation.Line.lwd    = 1,                       #   125
    Map.State.Spec.cex     = .32,                     # label size for AK, HI, DC in top map.   #  126 

# debug parameter
    MST.Debug              = 0                        # debug switch - for use by developers only. (default=0)

  )

#
# When something is added or deleted from this structure, change the 
# globalVariables call at the start of this module.
#

#
#  Set up variable in the micromapST namespace - used by micromapST and micromapSTSetDefaults functions.
#

micromapSTDefaults = list(colors=colors,details=details)


     return(micromapSTDefaults)
   }

####
#
#   .onLoad function - executed when the package is loaded initially.
#      builds a non-changeable micromapSTDefault data.frame for use
#      as the default when colors and/or details are not specified.
#
#    Added by JP - Oct, 2012 - Setup permanent micromapSTDefault data.frame for 
#          use as the default values on the call.
#
####

.onLoad = function (libraryName, pkgName)

   { 
     #packageStartupMessage(".onLoad")
     #packageStartupMessage(libraryName)
     #packageStartupMessage(pkgName)
     # generate default data.frame for micromapST.
     #rlMicromapSTDefaults <- micromapSTSetDefaults()
     #micromapSTDefaults <<- rlMicromapSTDefaults
  
    }

#  
####  
#
# End of load and variable initialization
#
####

