#############################################################################################
## File: PKoutput.R
## Author: Xiaoyong Sun
## Date: 10/12/2009
## Goal: PKoutput
## Notes:
##      -
#############################################################################################


## user should not use this function
## match lattice figure against ggplot figure to output both figure parallel
#  input: - lattice.note: figure note
#         - ggplot.note: figure note
#  output: vector index for ggplot.note that matches lattice.note
match.index <- function(lattice.note, ggplot.note)
{
    tmp.lattice.list <- strsplit(lattice.note, "lattice")
    tmp.ggplot.list <- strsplit(ggplot.note, "ggplot2")

    tmp.lattice <- unlist(tmp.lattice.list)
    rm.index <- which(tmp.lattice==")")
    tmp.lattice <- tmp.lattice[-rm.index]

    tmp.ggplot <- unlist(tmp.ggplot.list)
    rm.index <- which(tmp.ggplot==")")
    tmp.ggplot <- tmp.ggplot[-rm.index]

    ggplot.index <- match(tmp.lattice, tmp.ggplot)

    return(ggplot.index)
}

# Create fig file with html directory from input and fig format from input
#   - input: 1. the html directory to keep figures
#            2. the save format for figures, in order to find the matching figs in the folder
#   - output: a html file with figures
writeHTML <- function(htmldir, save.format)
{
        # web browser compatibility
        html.start <- "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
                              \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
        html.start <- paste(html.start, "<html><header><link rel=\"stylesheet\"
                              type=\"text/css\" href=\"",
                              #system.file("doc", package="PKreport"),
                              "http://129.186.62.12/PKreport.css",
                              "\"/>", sep="")

        html.start <- paste(html.start, "<title>PKoutput", htmldir, "Figures</title></header><body>",
                              sep=" ")
        html.start <- paste(html.start, "<hr><h3><p class=\"subtitle\">", htmldir,
                        "Figures</p></h3>", sep=" ")
        
        html.content <- NULL

        fig.files <- dir(htmldir, pattern=save.format)
        fig.files <- paste(htmldir, fig.files, sep="/")

        txt <- paste("<img src=", fig.files, sep="")
        txt <- paste(txt, ">", sep="")

        # lattice package
        if (.pkplot$getConfig("package")==0)
        {
            lattice.txt <- txt[grep("lattice", txt)]
            lattice.no <- sapply(1:length(lattice.txt), function(i) strsplit(lattice.txt[i], "_")[[1]][2])
            lattice.note <- .pkplot$getAllPKCodeNote()[c(as.numeric(lattice.no))]
            lattice.note <- unlist(lattice.note)

            # add content
            # html.content: Notes for figures
            # html.body: figures
            # html code: 1. html.content; 2. html.body
            html.tmp <- paste("<a href=\"#", lattice.no, "\">", lattice.note, "</a>", sep="") 
            html.tmp <- paste("<td>", html.tmp, "</td>", sep="")

            ## every 4 contents in one row
            # add <tr>
            tr.index <- which((((1:length(lattice.note))-1)%%4) == 0)
            html.tmp[tr.index] <- paste("<tr>", html.tmp[tr.index], sep="")
            # add </tr>
            tr.index <- which(((1:length(lattice.note))%%4) == 0)
            html.tmp[tr.index] <- paste(html.tmp[tr.index], "</tr>", sep="")
             
            html.content <- paste(html.content, html.tmp, collapse="")
            html.content <- paste("<a name=\"top\"><h3><table cellpadding=5 align=center>", html.content, "</table></h3></a>", sep="") 

            title.body <- "<p><center><h2>Lattice package</h2></center></p>"
            fig.tmp <- paste("<a name=\"", lattice.no, "\">", lattice.no, "</a>", sep="")
            fig.body <- paste("<p>Figure ID: ", fig.tmp, "</p><a href=\"#top\">Top</a>", sep="")
            fig.body <- paste(fig.body, "<p> Note: ", lattice.note, "</p>", sep="")
            fig.body <- paste(fig.body, lattice.txt, collapse="")
            html.body <- paste(title.body, fig.body, sep="")
        }
     
        ## ggplot package
        if (.pkplot$getConfig("package")==1)
        {
            ggplot.txt <- txt[grep("ggplot", txt)]
            if (length(ggplot.txt) > 0)
            {
                ggplot.no <- sapply(1:length(ggplot.txt), function(i) strsplit(ggplot.txt[i], "_")[[1]][2])
                ggplot.note <- .pkplot$getAllPKCodeNote()[c(as.numeric(ggplot.no))]
                ggplot.note <- unlist(ggplot.note)               

                # add content
                html.tmp <- paste("<a href=\"#", ggplot.no, "\">", ggplot.note, "</a>", sep="") 
                html.tmp <- paste("<td>", html.tmp, "</td>", sep="")

                ## every 4 contents in one row
                # add <tr>
                tr.index <- which((((1:length(ggplot.note))-1)%%4) == 0)
                html.tmp[tr.index] <- paste("<tr>", html.tmp[tr.index], sep="")
                # add </tr>
                tr.index <- which(((1:length(ggplot.note))%%4) == 0)
                html.tmp[tr.index] <- paste(html.tmp[tr.index], "</tr>", sep="")
                 
                html.content <- paste(html.content, html.tmp, collapse="")
                html.content <- paste("<a name=\"top\"><h3><table cellpadding=5 align=center>", html.content, "</table></h3></a>", sep="") 

              
                title.body <- "<p><center><h2>ggplot package</h2></center></p>"
                fig.tmp <- paste("<a name=\"", ggplot.no, "\">", ggplot.no, "</a>", sep="")
                fig.body <- paste("<p>Figure ID: ", fig.tmp, "</p><a href=\"#top\">Top</a>", sep="")
                fig.body <- paste(fig.body, "<p> Note: ", ggplot.note, "</p>", sep="")
                fig.body <- paste(fig.body, ggplot.txt, collapse="")
                html.body <- paste(title.body, fig.body, sep="")
               
            }
        }

        if (.pkplot$getConfig("package")==2)
        {
            lattice.txt <- txt[grep("lattice", txt)]
            ggplot.txt <- txt[grep("ggplot", txt)]
            
            if (length(lattice.txt)==0) stop("There is no lattice figure!")
            if (length(ggplot.txt)==0) stop("There is no ggplot figure!")

            lattice.no <- sapply(1:length(lattice.txt), function(i) strsplit(lattice.txt[i], "_")[[1]][2])
            lattice.note <- .pkplot$getAllPKCodeNote()[c(as.numeric(lattice.no))]
            lattice.note <- unlist(lattice.note)
            
            ggplot.no <- sapply(1:length(ggplot.txt), function(i) strsplit(ggplot.txt[i], "_")[[1]][2])
            ggplot.note <- .pkplot$getAllPKCodeNote()[c(as.numeric(ggplot.no))]
            ggplot.note <- unlist(ggplot.note)
            
            ggplot.index <- match.index(lattice.note, ggplot.note)
            ggplot.txt <- ggplot.txt[ggplot.index]
            ggplot.note <- ggplot.note[ggplot.index]
            ggplot.no <- ggplot.no[ggplot.index]

            # add content
            # html.content: Notes table for figures
            # html.body: figures
            # html code: 1. html.content; 2. html.body
            lattice.note.1 <- strsplit(lattice.note, " \\(")
            lattice.note.2 <- unlist(lattice.note.1)
            tmp.index <- which(lattice.note.2 == "lattice)")
            lattice.note.3 <- lattice.note.2[-tmp.index]
            
            html.tmp <- paste("<a href=\"#", lattice.no, "\">", lattice.note.3, "</a>", sep="")
            html.tmp <- paste("<td>", html.tmp, "</td>", sep="")

            ## every 4 contents in one row
            # add <tr>
            trs.index <- which((((1:length(lattice.note))-1)%%4) == 0)
            html.tmp[trs.index] <- paste("<tr>", html.tmp[trs.index], sep="")
            # add </tr>
            tre.index <- which(((1:length(lattice.note))%%4) == 0)
            html.tmp[tre.index] <- paste(html.tmp[tre.index], "</tr>", sep="")

            html.content <- paste(html.content, html.tmp, collapse="")
            # control last component
            if (length(trs.index) != length(tre.index)) html.content <- paste(html.content, "</tr>", sep="")
            html.content <- paste("<a name=\"top\"><h3><table class=\"figcenter\">", html.content, "</table></h3></a>", sep="")

            # setup 'figure id'
            fig.tmp.lattice <- paste("<a name=\"", lattice.no, "\">", lattice.no, "</a>", sep="")
            fig.id.lattice <- paste("<td>Figure ID: ", fig.tmp.lattice, "<a href=\"#top\"> Top </a></td>", sep="")
            fig.id.ggplot <- paste("<td>Figure ID: ", ggplot.no, "</td>", sep="")
            fig.id <- paste("<tr>", fig.id.lattice, fig.id.ggplot, "</tr>", sep="")

            # setup 'figure.note'
            fig.note.lattice <- paste("<td> Note: ", lattice.note, "</td>", sep="")
            fig.note.ggplot <- paste("<td> Note: ", ggplot.note, "</td>", sep="")
            fig.note <- paste("<tr>", fig.note.lattice, fig.note.ggplot, "</tr>", sep="")

            # setup 'figure'
            fig.fig.lattice <- paste("<td>", lattice.txt, "</td>", sep="")
            fig.fig.ggplot <- paste("<td>", ggplot.txt, "</td>", sep="")
            fig.fig <- paste("<tr>", fig.fig.lattice, fig.fig.ggplot, "</tr>", sep="")
            
            fig.body <- paste(fig.id, fig.note, fig.fig, sep="")
            fig.body <- paste(fig.body, collapse="")
            fig.body <- paste("<table class=\"center\">", fig.body, "</table>", sep="")

            title.body <- "<p><center><h2>lattice and ggplot package</h2></center></p>"
            html.body <- paste(title.body, fig.body, sep="")
        }

        html.end <- paste("<p></p><p></p><hr><p align=\"right\" class=\"bottomFont\">PKreport", packageDescription("PKreport")$Version,
                          date(), "</p></body></html>", sep=" ")                          
        filename <- paste("PK", htmldir, ".html", sep="")
        PK.html <- file(filename, "w")  # open an output file connection
        html.all <- paste(html.start, html.content, html.body, html.end, sep="\n")
        writeLines(html.all, con = PK.html, sep = "\n")
        close(PK.html)
}

# Create file with file name from input and content from input
#   - input: 1. the html filename you want to creat
#            2. the title of html file
#            2. the content as a table
#   - output: a html file

writeTable <- function(filename, mytitle, mydata)
{
    
        # web browser compatibility
    html.start <- "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
                              \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
    html.start <- paste(html.start, "<html><header><link rel=\"stylesheet\"
                              type=\"text/css\" href=\"",
                              #system.file("doc", package="PKreport"),
                              "http://129.186.62.12/PKreport.css",
                              "\"/>", sep="")
  
    html.end <- paste("<p></p><p></p><hr><p align=\"right\" class=\"bottomFont\">PKreport", packageDescription("PKreport")$Version,
                          date(), "</p></body></html>", sep=" ")
                            
    # generate data file

    html.start<- paste(html.start, "<title>", mytitle, 
                          "</title></header><body>", sep=" ")
    html.start <- paste(html.start, "<hr><h3><p class=\"subtitle\">", mytitle,
                        "Table</p></h3>", sep=" ")
    ## Add description
    html.content <- paste("<h4>Color Description (scaled in column direction)</h4>", sep="")
    html.content <- paste(html.content, 
        "<table><tr><td BGCOLOR=	#33FFFF>25%</td><td BGCOLOR=	#FFFF99>50%</td><td BGCOLOR=	#FFCC00>75%</td><td BGCOLOR=	#FF6600>100%</td></tr></table>", sep="")                           
    html.content <- paste(html.content, "<h4>Data</h4><table>", sep="")

    html.colname <- paste("<td BGCOLOR= #66FF66>", colnames(mydata), "</td>", sep="", collapse="")
    html.content <- paste(html.content, "<tr>", html.colname, "</tr>", sep="")
        
    html.data <- NULL
    
    for (i in 1:ncol(mydata))
    {
        tmp <- mydata[,i]

        if (length(unique(tmp))>=5 && is.numeric(tmp))
        {
            order.index <- order(as.numeric(tmp))
            
            # TODO: have more color choice?
            quantile.index <- round(quantile(1:length(order.index), c(0.25,0.5,0.75)))
            
            # q25, 50, 75, 100
            q25.index <- c(1:quantile.index[1])
            if (quantile.index[1]!=quantile.index[2]) q50.index <- c((quantile.index[1]+1):quantile.index[2])
            if (quantile.index[2]!=quantile.index[3]) q75.index <- c((quantile.index[2]+1):quantile.index[3])
            if (quantile.index[3] < length(order.index)) q100.index <- c((quantile.index[3]+1):(length(order.index))) 
            
            tmp[order.index[q25.index]] <- paste("<td BGCOLOR=	#33FFFF>", tmp[order.index[q25.index]], "</td>", sep="")
            tmp[order.index[q50.index]] <- paste("<td BGCOLOR=	#FFFF99>", tmp[order.index[q50.index]], "</td>", sep="")
            tmp[order.index[q75.index]] <- paste("<td BGCOLOR=	#FFCC00>", tmp[order.index[q75.index]], "</td>", sep="")
            tmp[order.index[q100.index]] <- paste("<td BGCOLOR=	#FF6600>", tmp[order.index[q100.index]], "</td>", sep="")                        

        }
        else 
        {
            if (length(unique(tmp)) >=3  && is.numeric(tmp))
            {
              #min
              min.index <- which(tmp==min(tmp))
              max.index <- which(tmp==max(tmp))
              
              tmp[min.index] <- paste("<td BGCOLOR=	#33FFFF>", tmp[min.index], "</td>", sep="")
              tmp[max.index] <- paste("<td BGCOLOR=	#FF6600>", tmp[max.index], "</td>", sep="")
              tmp[-c(min.index, max.index)] <- paste("<td>", tmp[-c(min.index, max.index)], "</td>", sep="")
            }
            else
            {
              tmp <- paste("<td>", tmp, "</td>", sep="")
            }
        }
        
        mydata[,i] <- tmp
                
    } 
       
    for (i in 1:nrow(mydata))
    {
         html.data <- paste(mydata[i,], sep="", collapse="")
         html.content <- paste(html.content, "<tr>", html.data, "</tr>", sep="")         
    }

        
    html.content <- paste(html.content, "</table>", sep="") 
        

    PK.html <- file(filename, "w")  # open an output file connection
    html.all <- paste(html.start, html.content, html.end, sep="\n")
    writeLines(html.all, con = PK.html, sep = "\n")
    close(PK.html)         
                     
}

# Create file with file name from input and content from input
#   - input: 1. the html filename you want to creat
#            2. the content as a string
#   - output: a html file
writeLst <- function(filename, lstString)
{
        # web browser compatibility
    html.start <- "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
                              \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
    html.start <- paste(html.start, "<html><header><link rel=\"stylesheet\"
                              type=\"text/css\" href=\"",
                              #system.file("doc", package="PKreport"),
                              "http://129.186.62.12/PKreport.css",
                              "\"/>", sep="")
  
    html.end <- paste("<p></p><p></p><hr><p align=\"right\" class=\"bottomFont\">PKreport", packageDescription("PKreport")$Version,
                          date(), "</p></body></html>", sep=" ")
                            
    # generate data file
    #browser()
    html.start<- paste(html.start, "<title>lst file",
                          "</title></header><body>", sep=" ")
    html.start <- paste(html.start, "<hr><h3><p class=\"subtitle\">lst file</p></h3>", sep=" ")
    html.content <- paste(lstString, "<br>", sep="")
    html.content <- paste(html.content, collapse="")
    html.content <- paste("<p><pre>", html.content, "</pre></p>", sep="")
    
    PK.html <- file(filename, "w")  # open an output file connection
    html.all <- paste(html.start, html.content, html.end, sep="\n")
    writeLines(html.all, con = PK.html, sep = "\n")
    close(PK.html) 
        
}

writeLst.tab <- function(filename, lstString)
{
        # web browser compatibility
    html.start <- "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
                              \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
    html.start <- paste(html.start, "<html><header><link rel=\"stylesheet\"
                              type=\"text/css\" href=\"",
                              #system.file("doc", package="PKreport"),
                              "http://129.186.62.12/PKreport.css",
                              "\"/>", sep="")
  
    html.end <- paste("<p></p><p></p><hr><p align=\"right\" class=\"bottomFont\">PKreport", packageDescription("PKreport")$Version,
                          date(), "</p></body></html>", sep=" ")
                            
    # generate data file
    html.start<- paste(html.start, "<title>lst file",
                          "</title></header><body>", sep=" ") 
    
          
    html.content <- "<table>"
    for (i in 1:length(lstString))
    {
        tmp <- strsplit(lstString[i], " ")[[1]]
        tmp.html <- paste(tmp, collapse="</td><td>")
        html.content <- paste(html.content, "<tr><td>", tmp.html, "</td></tr>", sep="")
    }
    html.content <- paste(html.content, "</table>", sep="")
    
    
    PK.html <- file(filename, "w")  # open an output file connection
    html.all <- paste(html.start, html.content, html.end, sep="\n")
    writeLines(html.all, con = PK.html, sep = "\n")
    close(PK.html) 
        
}

PKoutput <- function(nonmemObj=NULL, table.Rowv=FALSE, table.Colv=FALSE)
{
    general.list <- .pkplot$getGlobalConfig()

        # web browser compatibility
    html.start <- "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\"
                              \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">"
    html.start <- paste(html.start, "<html><header><link rel=\"stylesheet\"
                              type=\"text/css\" href=\"",
                              #system.file("doc", "PKreport.css", package="PKreport"),
                              "http://129.186.62.12/PKreport.css",
                              #"PKreport.css",
                              "\"/>", sep="")
    html.index <- paste(html.start, "<title>PKreport</title></header><body>", sep="")

    ## add symbol image
    html.index <- paste(html.index, "<h1><p class=\"title\">PKreport</p></h1>", sep="")

    ## general output
    if ((!is.null(nonmemObj)) && class(nonmemObj) == "nonmem")
    {
        # start table
        html.index <- paste(html.index,
            "<h3><p class=\"subtitle\">NONMEM results</p></h3><table class=\"center\" border=0>",
            sep="")
        # lst file
        if (length(nonmemObj@file.lst) != 0)
        {
            writeLst("lst.html", nonmemObj@file.lst)
            html.index <- paste(html.index, "<tr><td><a href=lst.html>lst file</a></td></tr>", sep="")
        }

        # tab file, match that in "nonmem" class    
        if (nrow(nonmemObj@tabdata) != 0)
        {
            order.data <- nonmemObj@tabdata
            if (table.Colv || table.Rowv) 
            {
                tmp.data <- order.data[,sapply(order.data, is.numeric)]
                tmp.heat <- heatmap(as.matrix(tmp.data))
                if (table.Rowv) order.data <- order.data[rev(tmp.heat$rowInd),]
                if (table.Colv) order.data <- order.data[, tmp.heat$colInd]
            }
            
            writeTable("tab.html", "tab file", order.data)
            html.index <- paste(html.index, "<tr><td><a href=tab.html>Tab file</a></td></tr>", sep="")
        }
            
        if (nrow(nonmemObj@file.cov$data) != 0)
        {
            order.data <- nonmemObj@file.cov$data
            if (table.Colv || table.Rowv) 
            {
                tmp.data <- order.data[,sapply(order.data, is.numeric)]
                tmp.heat <- heatmap(as.matrix(tmp.data))
                if (table.Rowv) order.data <- order.data[rev(tmp.heat$rowInd),]
                if (table.Colv) order.data <- order.data[, tmp.heat$colInd]
            }        
            writeTable("cov.html", "cov file",order.data)
            html.index <- paste(html.index, "<tr><td><a href=cov.html>Cov file</a></td></tr>", sep="")
        }
        
        if (nrow(nonmemObj@file.cor$data) != 0)
        {
            order.data <- nonmemObj@file.cor$data
            if (table.Colv || table.Rowv) 
            {
                tmp.data <- order.data[,sapply(order.data, is.numeric)]
                tmp.heat <- heatmap(as.matrix(tmp.data))
                if (table.Rowv) order.data <- order.data[rev(tmp.heat$rowInd),]
                if (table.Colv) order.data <- order.data[, tmp.heat$colInd]
            }    
            writeTable("cor.html", "cor file",order.data)
            html.index <- paste(html.index, "<tr><td><a href=cor.html>Cor file</a></td></tr>", sep="")
        }

        if (nrow(nonmemObj@file.coi$data) != 0)
        {
            order.data <- nonmemObj@file.coi$data
            if (table.Colv || table.Rowv) 
            {
                tmp.data <- order.data[,sapply(order.data, is.numeric)]
                tmp.heat <- heatmap(as.matrix(tmp.data))
                if (table.Rowv) order.data <- order.data[rev(tmp.heat$rowInd),]
                if (table.Colv) order.data <- order.data[, tmp.heat$colInd]
            }            
            writeTable("coi.html", "coi file",order.data)
            html.index <- paste(html.index, "<tr><td><a href=coi.html>Coi file</a></td></tr>", sep="")
        }
          
        if (nrow(nonmemObj@file.phi$data) != 0)
        {
            order.data <- nonmemObj@file.phi$data
            if (table.Colv || table.Rowv) 
            {
                tmp.data <- order.data[,sapply(order.data, is.numeric)]
                tmp.heat <- heatmap(as.matrix(tmp.data))
                if (table.Rowv) order.data <- order.data[rev(tmp.heat$rowInd),]
                if (table.Colv) order.data <- order.data[, tmp.heat$colInd]
            }         
            writeTable("phi.html", "phi file",order.data)
            html.index <- paste(html.index, "<tr><td><a href=phi.html>phi file</a></td></tr>", sep="")
        }                       
    }
    
    #########################
    ## add separator
    #########################
    html.index <- paste(html.index, "</table><hr><h3><p class=\"subtitle\">Diagnostics</p></h3><table class=\"center\" border=0>", sep="")
    
    ## univariate figures
    if (file.exists(general.list$univar.dir))
    {
        writeHTML(general.list$univar.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKunivar.html>Univariate Figures</a></td></tr>", sep="")
        
    }

    ## bivariate figures
    if (file.exists(general.list$bivar.dir))
    {
        writeHTML(general.list$bivar.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKbivar.html>Bivariate Figures</a></td></tr>", sep="")
    }

    ## bivariate figures
    if (file.exists(general.list$ind.dir))
    {
        writeHTML(general.list$ind.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKind.html>Individual Figures</a></td></tr>", sep="")
    }

    ## gof figures
    if (file.exists(general.list$gof.dir))
    {
        writeHTML(general.list$gof.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKgof.html>Figures for Goodness of Fit</a></td></tr>", sep="")
    }

    if (file.exists(general.list$struct.dir))
    {
        writeHTML(general.list$struct.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKstruct.html>Figures for Structural Model Dagnostics</a></td></tr>", sep="")
    }
    if (file.exists(general.list$resid.dir))
    {
        writeHTML(general.list$resid.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKresid.html>Figures for Residual Model Dagnostics</a></td></tr>", sep="")
    }

    if (file.exists(general.list$para.dir))
    {
        writeHTML(general.list$para.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKpara.html>Figures for Parameters Dagnostics</a></td></tr>", sep="")
    }

    if (file.exists(general.list$cov.dir))
    {
        writeHTML(general.list$cov.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKcov.html>Figures for Covariate Model Dagnostics</a></td></tr>", sep="")
    }
    
    if (file.exists(general.list$eta.dir))
    {
        writeHTML(general.list$eta.dir,general.list$save.format)
        html.index <- paste(html.index, "<tr><td><a href=PKeta.html>Figures for Random Effects Dagnostics</a></td></tr>", sep="")
    }
    
    if (file.exists("PKcode.txt"))
    {
        html.index <- paste(html.index, "<tr><td><a href=PKcode.txt>R Code for Figures</a></td></tr>", sep="")
    }

    # end table
    html.index <- paste(html.index, "</table>", sep="")
    
    html.end <- paste("<p></p><p></p><hr><p align=\"right\" class=\"bottomFont\">PKreport", packageDescription("PKreport")$Version,
                          date(), "</p></body></html>", sep=" ")
    html.index <- paste(html.index, html.end, sep="")
    PK.html <- file("PKindex.html", "w")  # open an output file connection
    writeLines(html.index, con = PK.html, sep = "\n")
    close(PK.html)
}



