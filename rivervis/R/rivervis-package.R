

#' Ballinderry River Basin Dataset
#' 
#' A ballinderry river basin dataset for demonstrating purposes.
#' 
#' This contains the following 8 datasets and 2 character vectors.  
#' \describe{
#' 
#' \item{B.elevation}{ This dataset can be used to plot elevation
#' profile of Ballinderry Rivers. It is a data frame with 90 observations on
#' the following 3 variables.  
#' \describe{ 
#' \item{River}{Rivers on which the elevation sites are located}
#' \item{Distance}{The along-the-river distance between the elevation sites and the mouth of the river.} 
#' \item{Elevation}{A numeric vector of elevation values} } } 
#' 
#' \item{B.reach}{ Selected Ballinderry River reaches. It is a data
#' frame with 1 observation on the following 6 variables.  
#' \describe{
#' \item{Reach}{Reach names.} 
#' \item{River}{Rivers on which the monitoring sites are located.} 
#' \item{From}{A numeric vector of starting points of reaches.} 
#' \item{To}{A numeric vector of ending points of reaches.} 
#' \item{Group}{A vector of reach group names. This
#' indicates to which group the reaches belong.} 
#' \item{Style}{A vector
#' of reach styles and the location of reach lines.} } }
#' 
#' \item{B.river}{ Main rivers and tributaries in Ballinderry Basin for
#' \code{RiverMap} and \code{RiverLayout}. It is a data frame with 8
#' observations on the following 5 variables.  
#' \describe{
#' \item{River}{River names.} 
#' \item{Length}{Length of rivers.}
#' \item{Position}{Relative relations between rivers and their parent rivers.} 
#' \item{Parent}{Parent rivers.}
#' \item{Distance}{Distance between the mouths of each river and the
#' mouths of each river's parent.} } } 
#' 
#' \item{B.siteaspt}{
#' ASPT scores measured in the Ballinderry River Basin in Spring and Autumn, 2009. It is a
#' data frame with 15 observations on the following 5 variables.  
#' \describe{
#' \item{Site}{Site ID} 
#' \item{River}{Rivers on which the sites are located} 
#' \item{Distance"}{The along-the-river distance between the 
#' sites and the mouth of the river.}
#' \item{ASPT_Spring}{ASPT measure in Spring.} 
#' \item{ASPT_Autumn}{ASPT measure in Autumn.} 
#' } }
#' 
#' \item{B.sitehm}{ Selected hydromorphological
#' results from RHAT. The hydromorphological variables are ordinary factors
#' which have five grades: High, Good, Moderate, Poor and Bad. It is a data
#' frame with 17 observations on the following 9 variables.  
#' \describe{
#' \item{Site}{Monitorings sites.} 
#' \item{River}{Rivers on which the monitoring sites are located.} 
#' \item{Distance}{The along-the-river distance between the sites and the mouth of the river.}
#' \item{ChanVeg}{Channel vegetation condition.}
#' \item{ChanFlow}{Channel flow condition.}
#' \item{BankVegLeft}{Left bank vegetation condition.}
#' \item{BankVegRight}{Right bank vegetation condition.}
#' \item{RipLULeft}{Left riparian land-use condition.} 
#' \item{RipLURight}{Right riparian land-use condition.} } }
#' 
#' \item{B.sitenh4n}{
#' NH4-N values measured in the Ballinderry River Basin in Spring and Autumn, 2009. It is a
#' data frame with 17 observations on the following 5 variables.  
#' \describe{
#' \item{Site}{Site ID} 
#' \item{River}{Rivers on which the sites are located} 
#' \item{Distance"}{The along-the-river distance between the 
#' sites and the mouth of the river.}
#' \item{NH4N_Spring}{NH4-N measure in Spring.} 
#' \item{NH4N_Autumn}{NH4-N measure in Autumn.} 
#' } }
#' 
#' \item{B.soi}{ This dataset provides sites of interest in the Ballinderry River Basin. 
#' The sites have three types: towns, conjunctions of left and right tributaries. 
#' It is a data frame with 5 observations on the following 4 variables.  
#' \describe{
#' \item{SOI}{Sites of interest.} 
#' \item{River}{Rivers on which the sites are located.} 
#' \item{Distance}{The along-the-river distance
#' between the sites and the mouth of the river.} 
#' \item{Group}{Groups of the sites.} } } 
#' 
#' \item{B.town}{ This dataset provides 2 main towns in the Ballinderry River Basin. 
#' It has the following 4 variables.  
#' \describe{
#' \item{Town}{Town names.} 
#' \item{River}{Rivers on which the sites are located.} 
#' \item{Distance}{The along-the-river distance between the sites and the mouth of the river.} 
#' \item{Group}{Groups of the sites.} } }
#' 
#' \item{fivegrades}{This vector contains five grades, which are High, Good, Moderate,
#' Poor and Bad}
#' 
#' \item{fivecolours}{This vector contains five colours representing the five grades.
#' The five colours are blue(#5381FFFF), green(#7BE859FF), yellow(#FFC944FF), 
#' orange(#E87539FF) and red(#FF3931FF).} }
#' 
#' @name Ballinderry
#' @aliases B.elevation B.reach B.river B.siteaspt B.sitehm B.sitenh4n B.soi B.town fivecolours fivegrades
#' @docType data
#' @source North Ireland Environment Agency
#' @keywords datasets
NULL





#' rivervis Package: River Visualisation Tool
#' 
#' The \bold{rivervis} package is designed to visualise river ecosystem data.
#' 
#' In general, the \bold{rivervis} package draws two types of diagrams - river
#' charts and river block charts. River charts can present points, lines, bars
#' and blocks in relation to the topological structure of the river network.
#' River block charts show qualitative data without the river network
#' structure. It is recommended to run the examples below and in each function
#' manual. The \bold{rivervis} package contains 15 functions in total.
#' \describe{ \item{RiverLayout}{ This calculates best fit plotting
#' coordinates for rivers to be shown on river charts. The output is a list,
#' which is used when plotting the river chart and the information on the river
#' chart. It provides an opportunity to change the coordinates and other
#' plotting parameters before actually plotting.  } \item{RiverDraw}{
#' This plots the river charts according to the output list of
#' \code{RiverLayout}.  } \item{RiverMap}{ This can be understood as a
#' combination of \code{RiverLayout} and \code{RiverDraw}. It not only
#' calculates best fit plotting coordinates for rivers to be shown on river
#' charts, but also plots the river charts according to the calculated
#' coordinates. This implies that the coordinates cannot be changed before
#' river chart plotting.  } \item{RiverFrame}{ This plots river frames,
#' lead lines and archor points.  } \item{RiverPoint}{ This plots
#' points or broken lines on the river chart.  } \item{RiverBar}{ This
#' plots bars for quantitative data on the river chart.  }
#' \item{RiverBlock}{ This plots blocks for qualitative data on the
#' river chart.  } \item{RiverSite}{ This plots sites of interest on
#' the river chart.  } \item{RiverLabel}{ This adds the name labels to
#' the plotted rivers.  } \item{RiverTM}{ This adds tick marks to the
#' river chart.  } \item{RiverAxisLabel}{ This adds left or right axis
#' labels to the river chart.  } \item{RiverReach}{ This highlights
#' river reaches on the river chart.  } \item{RiverDirection}{ This
#' adds a flow direction arrow on the river chart.  }
#' \item{RiverScale}{ This adds a plotting scale on the river chart.  }
#' \item{RiverBlockChart}{ This function plots a river block chart for
#' qualitative data without the topological structure of the river network. The
#' function does not require the output list from \code{RiverLayout} or
#' \code{RiverMap}.  } }
#' 
#' @name rivervis-package
#' @aliases rivervis
#' @docType package
#' @author Feng Mao, Yichuan Shi, and Keith Richards
#' @keywords hplot
#' @examples
#' 
#' 
#' 
#' data(Ballinderry)
#' 
#' riverlayout <- RiverLayout(B.river$River, B.river$Length, B.river$Parent, 
#'                            B.river$Position, B.river$Distance, direction = -1)
#' 
#' # Example Figure 1
#' 
#' RiverDraw(riverlayout)
#' RiverLabel(riverlayout, offset = -1, corner = "lt", srt = 0, adj = c(0, -0.7))
#' 
#' RiverBar(B.siteaspt$Site, B.siteaspt$River, B.siteaspt$Distance, 
#'          B.siteaspt[4:5], riverlayout, range = c(0,8), 
#'          bar.col = c("#5381FFFF", "#FF3931FF"),lbl.adj = c(0.5,1.3))
#' 
#' RiverPoint(B.sitenh4n$Site, B.sitenh4n$River, B.sitenh4n$Distance, 
#'            B.sitenh4n$NH4N_Spring, riverlayout, type = "o", 
#'            pt.col = "#5381FFFF", pt.pch = 21, pt.bg = "lightblue")
#' RiverPoint(B.sitenh4n$Site, B.sitenh4n$River, B.sitenh4n$Distance, 
#'            B.sitenh4n$NH4N_Autumn, riverlayout, type = "o", 
#'            pt.col = "#FF3931FF", pt.pch = 21, pt.bg = "pink")
#' 
#' RiverSite(B.town$Town, B.town$River, B.town$Distance, B.town$Group, 
#'           riverlayout, pt.pch = 22, lbl.shw = FALSE, 
#'           pt.bg = "orange", pt.col = "black")
#' 
#' RiverSite(B.soi$SOI, B.soi$River, B.soi$Distance, B.soi$Group, riverlayout, 
#'           pt.pch = c(25, 24, NA), lbl.shw = FALSE, pt.bg = NA, pt.col = "black")
#' 
#' RiverTM(c(0,2,4,6,8,10), B.siteaspt[4:5], riverlayout, pos=-1, side = "L", 
#'         range = c(0,8), label = c(0,2,4,6,8))
#' 
#' RiverTM(c(0,0.04,0.08,0.12), B.sitenh4n[4:5], riverlayout, pos=-1, side = "R", 
#'         range = c(0,0.15), label = c(0,0.04,0.08,0.12))
#' 
#' RiverAxisLabel("ASPT score", riverlayout, adj = c(0.5, -3))
#' 
#' RiverAxisLabel(expression(paste("N ",H[4],"-N (mg/L)")), 
#'                riverlayout, side = "R", 
#'                srt = 270, adj = c(0.5, -3))
#' 
#' legend(0.8, 0.43, inset=0.05, title = "Legend", 
#'        c("ASPT Spring", "ASPT Autumn", 
#'          expression(paste(NH[4],"-N Spring")), 
#'          expression(paste(NH[4],"-N Autumn")),
#'          "Town", "Unshown left tribs",
#'          "Unshown right tribs"), 
#'        lty = c(-1,-1,1,1,-1,-1,-1), 
#'        pch = c(22,22, 21,21, 22, 25, 24),
#'        col= c("black", "black", "#5381FFFF", "#FF3931FF",
#'               "black", "black", "black"),
#'        pt.bg = c("#5381FFFF", "#FF3931FF", "lightblue", 
#'                  "pink", "orange", NA, NA),
#'        pt.cex = c(2, 2, 1, 1, 1,1,1),
#'        cex = 0.8)
#' 
#' RiverScale(2, "2 km", riverlayout, loc = c(0.6, 0.10),lbl.cex = 0.8)
#' 
#' RiverDirection(riverlayout, arw.length = 0.03, 
#'                loc = c(0.6, 0.05), lbl.cex = 0.8)
#' 
#' # Example Figure 2
#' 
#' RiverDraw(riverlayout)
#' RiverLabel(riverlayout, offset = -1, corner = "lt", 
#'            srt = 0, adj = c(0, -0.7))
#' 
#' RiverReach(B.reach$Reach, B.reach$River, B.reach$From, B.reach$To, 
#'            B.reach$Group, B.reach$Style, riverlayout, rea.lwd = 4, 
#'            rea.lty = 3,rea.col = "#51B0A8FF")
#' 
#' RiverPoint(NA,B.elevation$River, B.elevation$Distance, 
#'            B.elevation$Elevation, riverlayout)
#' 
#' RiverTM(c(0, 100, 200, 300, 400, 500), B.elevation[3], riverlayout, 
#'         pos=-1, side = "R", range = c(0,500), 
#'         label = c(0, 100, 200, 300, 400, 500))
#' 
#' RiverAxisLabel("Elevation (m)", riverlayout, side = "R", 
#'                srt = 270, adj = c(0.5, -4))
#' 
#' RiverBlock(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'            B.sitehm[4:9], riverlayout, c(1,1,2,2), 
#'            block.col = fivecolours, lbl.adj = c(0.5,1.3), 
#'            par.txt = c("ChanVeg", "ChanFlow", "BankVegLeft", 
#'                        "Right", "RipLULeft", "Right"))
#' 
#' legend(0.8, 0.43, inset=0.05,title = "Legend",
#'        c("High", "Good", "Moderate", "Poor", "Bad",
#'          "Elevation", "Upper Ballinderry SAC"),
#'        pch = c(22,22,22,22,22,NA,22),
#'        pt.bg = c(fivecolours, "black","#51B0A8FF"),
#'        pt.cex = c(2,2,2,2,2,NA,3),
#'        lty = c(NA,NA,NA,NA,NA,1,NA),
#'        cex = 0.8)
#' 
#' RiverScale(2, "2 km", riverlayout, loc = c(0.6, 0.10),lbl.cex = 0.8)
#' 
#' RiverDirection(riverlayout, arw.length = 0.03, 
#'                loc = c(0.6, 0.05), lbl.cex = 0.8)
#' 
#' # Figure 3
#' 
#' RiverBlockChart(B.sitehm$Site, B.sitehm$River, B.sitehm$Distance, 
#'                 B.sitehm[4:9],  c(1,1,2,2), mar = 0.15, 
#'                 site.ofs = 1, site.cex = 0.7, 
#'                 site.order = "R", block.col = fivecolours)
#' 
#' legend(0.88, 0.6, inset=0.05, title = "Legend", 
#'        c("High", "Good", "Moderate", "Poor", "Bad"), 
#'        border = rep("black", 5), 
#'        fill = fivecolours,
#'        cex = 0.8)
#'        
#'        
NULL



