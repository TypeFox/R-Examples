# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"data.checkout"<-
    function(obj=NULL,
             datafile = ".ask.",
             hlin = -99,
             dotcol = "black",
             dotpch = 16,
             dotcex = 1,
             idlab="ID",
             csv=NULL,
             main="Default",
             ...) {

        if (!is.null(obj) &&
            missing(dotcol) &&
            missing(dotpch) &&
            missing(dotcex) ){
            dotcol = obj@Prefs@Graph.prefs$bwdotcol
            dotpch = obj@Prefs@Graph.prefs$bwdotpch
            dotcex = obj@Prefs@Graph.prefs$bwdotcex
        }

        if (datafile==".ask.") {
            cat("Please enter the name of the data file you wish to check:\n")
            datafile <- readline()
        }

        if(is.null(csv)){
            cat("Is the data in CSV format? (y/n):\n")
            csv <- readline()

            if ((csv == "y") ||
                (csv == "Y") ||
                (csv == "yes") ||
                (csv == "YES") ||
                (csv == "Yes")) {
                csv <- TRUE
            } else {
                csv <- FALSE
            }
        }

        if(!is.readable.file(datafile)) {
            cat("There is no file with that name in the current directory!\n")
            return()
        }

        if (hlin == -99) {
            cat("Please enter the line number in the data file\n")
            cat("containing the column headers. Note that if you have\n")
            cat("commented them out, with a '#', for example, you may get\n")
            cat("unpredictable results.) (1): ")
            hlin <- readline()
        }
        if(hlin == "")
            hlin <- 1
        hlin <- as.numeric(hlin)
        if(hlin < 0) {
            cat("The line number has to be > 0!\n"
                )
            invisible()
            return()
        }

        hlin <- hlin - 1

        if (csv) {
            data <- read.csv(datafile, skip = hlin, header = T)
        } else {
            data <- read.table(datafile, skip = hlin, header = T)
        }

        delcol <- grep("XXXX*", names(data))

        if(length(delcol) != 0) {
            data <- data[,  - delcol]
        }



        tit <- paste("Data set checkout of", datafile)

        ## create list for plots
        number.of.plots <- dim(data)[2] - 1
        plotList <- vector("list",number.of.plots)
        plot.num <- 0 # initialize plot number

        idcol <- match(idlab,names(data))
        if(is.na(idcol)) idcol <- 1
        #ids <- unique(data[, idcol])

        for(j in 1:dim(data)[2]) {

            if (j==idcol) next
            xlb <- names(data)[j]
            ylb <- names(data)[idcol]

            xplot <- dotplot(data[, idcol] ~ data[, j],
                             main = NULL,
                             scales = list(x = list(cex = 0.7),
                             y = list(cex = 0.45)),
                             xlab = xlb,
                             ylab = ylb,
                             col = dotcol,
                             cex = dotcex,
                             pch = dotpch
                             )

            plot.num <- plot.num+1
            plotList[[plot.num]] <- xplot
        }

        default.plot.title <- tit
        plotTitle <- xpose.multiple.plot.title(object=obj,
                                               plot.text = default.plot.title,
                                               main=main,
                                               ...)

        obj <- xpose.multiple.plot(plotList,plotTitle,...)
        return(obj)

    }
