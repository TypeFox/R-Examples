## genlatex.R --- 
## Author          : Claus Dethlefsen
## Created On      : Tue May 07 10:10:39 2002
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Nov 03 13:34:12 2005
## Update Count    : 48
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

genlatex <-  function(nwl,
                      outdir="pic/",
                      prefix="scoretable",
                      picdir="",
                      picpre="pic",
                      ncol=5,
                      nrow=7,
                      width=12/ncol,
                      vadjust=-1.8) {
    ## Create latex-table of pictex figures with references to the
    ## generated pictex-files.
    ##
    ## nwl: networkfamily
    ## outdir: where the file is stored
    ## prefix: the filename minus extension (which is .tex)
    ## picdir: where to find the picfiles (the path is inserted in the
    ##                                     latex files)
    ## picpre: the filenames of the picfiles are 'picpre'xx.pic, where
    ##                                  xx is the index of the network
    ## ncol: the number of columns in the table
    ## nrow: the number of rows in the table
    ## width: the width of each cell in the table
    ## vadjust: Vertical adjustment
    
    ## uses: fmt+findexponent defined locally
    ## and network-attributes: score, relscore

    findexponent <- function(x) {
        ## find exponent:
        n <- 0
        y <- x
        while (floor(y)==0) {
            n <- n+1
            y <- y*10
        }
        n
    }
    fmt <- function(x,digits=2) {
        ## format a number to a LaTeX string in scientific notation
        ## Used by: genlatex
        
        if (x==1) return("\\footnotesize{$1$}")
        if (x==0) return("\\footnotesize{$0$}")
        
        n <- findexponent(x)
        y <- x*10^n
        yy<- signif(y,digits)
        
        y <- as.character(signif(y,digits))
        
        if ( (yy*10)%%10==0) 
            y <- paste(y,".0",sep="")
        
        fod  <- paste("\\footnotesize{$",y)
        expo <- ifelse(n==0, "",paste("\\cdot 10^{-",n,"}"))
        paste(fod,expo,"$}")
    }
    
    
    dir.create(outdir)
    
    ff <- file(paste(outdir,prefix,".tex",sep=""),"w") ## output filename
    ## filename of picfile i
    pf <- function(i)  paste(picdir,picpre,i,".tex",sep="")
    
    ## how to include one picfile as a minipage with score and relscore
    putfig <- function(i) paste("\\vspace{",vadjust/2,"cm}",
                                "\\begin{minipage}[t]{",width,"cm}\n",
                                "\\input{",pf(i),"}\n",
                                "\\vspace{",vadjust,"cm}",
                                fmt(nwl[[i]]$score),"\\\\\n",
                                fmt(nwl[[i]]$relscore),"\n",
                                "\\end{minipage}\n",sep="")
    
    finished <- FALSE
    
    cat("%% generated automatically",date(),"- Don't edit by hand\n",file=ff)
    cat("%% A master file:\n",file=ff)
    cat("%% \\documentclass{article}\n",file=ff)
    cat("%% \\usepackage{array,pictex}\n",file=ff)
    cat("%% \\begin{document}\n",file=ff)
    cat("%% \\input{scoretable}\n",file=ff)
    cat("%% \\end{document}\n",file=ff)
    fig <- 1
    
    while (!finished) {
        ## header
        if (fig %% (ncol*nrow) == 1 || fig == 1) {
            cat("\\begin{tabular}{",file=ff)
            for (i in 1:ncol)
                cat("|m{",width,"cm}",sep="",file=ff)
            cat("|}\\hline\n",file=ff)
        }
        
        ## figs
        for (i in 1:ncol) {
            if (fig==length(nwl)) {
                cat(putfig(fig),"\\\\ \n\\hline",file=ff)
                finished <- TRUE; break }
            if (i %% ncol == 0) {
                cat(putfig(fig),"\\\\",file=ff)
                if (fig %% (ncol*nrow) != 0) cat("[-9mm]",sep="",file=ff)
                cat("\n\\hline",file=ff)
            }
            else 
                cat(putfig(fig),"&\n",sep="",file=ff)
            fig <- fig + 1
        }
        
        ## footer
        if (fig %% (ncol*nrow) == 1 || finished) 
            cat("\\end{tabular}\\clearpage\n",file=ff)
        
    }
    close(ff)
    invisible()
}

genpicfile <- function(nwl,outdir="pic/",prefix="pic",w=1.6,h=1.6,bigscale=3) {
    ## Create latex-table of pictex figures with references to the
    ## generated pictex-files.
    ##
    ## nwl: networkfamily
    ## outdir: where the files are stored
    ## prefix: the filename prefix of all files
    ## w: width of pictex object
    ## h: height of pictex object
    ## bigscale: scaling of the best network, which is output in 'nice.tex'
    
    ## uses: plot.network
    

    cat("\nGenerating pic-files...")

    dir.create(outdir)
    
    ## the best
    pictex(
           paste(outdir,prefix,"nice.tex",sep=""),
           width=w*bigscale,height=h*bigscale
           )
    plot(nwl[[1]])
    dev.off()
    
    ## the rest
    for (i in 1:length(nwl)) {
        name <- paste(outdir,prefix,i,".tex",sep="")
        pictex(name,width=w,height=h)
        plot(nwl[[i]],cexscale=3,arrowlength=0.05,notext=TRUE)
        dev.off()
    }
    cat("complete\n")
}
