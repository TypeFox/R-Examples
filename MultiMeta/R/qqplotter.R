#'Quantile-quantile plot for observed vs. expected p-values
#'
#'\code{qqplotter} returns a pdf file containing the QQ-plot for a given set of p-values
#'
#'The function is useful for comparing expected and observed distributions of statistics arising from genome-wide testing. 
#'Two different types of input are allowed: 
#'\itemize{ 
#'\item an object such as vector or data-frame column, containing p-values; 
#'\item a file name and the column number for retrieving p-values from file. Note that separator is set to white space,
#'can be changed accordingly if needed. The special character "*" is allowed for selecting more than on file (e.g. in case there is
#'one file for each chromosome) but the p-values column position must be the same for all files.
#'}
#'
#'@param p_values A vector or a column in a data-set containing p-values to plot. This is an alternative to \code{file}.
#'@param res.file File name of the \code{file} containing p-values to retrieve. This alternative does not require loading the results and is an alternative to
#'the \code{p_values} parameter.
#'@param p.num.col Number indicating which column of the \code{file} contains p-values. Specify only if you are using the \code{file} option.
#'@param sep File separator to be used when reading \code{file}. Default is white space. 
#'@param title Title to appear on the plot.
#'@param out.file File name for the output plot.
#'@return The output will be a pdf file containing the QQplot. 
#'
#'@export



qqplotter <-
  function(p_values=NULL,res.file=NULL,p.num.col=NULL,sep=" ",out.file="QQplot.pdf",title=""){
    #checking the input
    if(length(res.file)!=0 & length(p.num.col)==0) stop("Please indicate which column contains p-values by changing the argument p.num.col")
    if(length(p_values)!=0 & length(res.file)!=0) stop("Only one parameter among p_values and file is required. See help(qqplotter) for details.")
    if(length(res.file)!=0 & length(p.num.col)!=0){ #read only the p-value column from file
      res=read.table(pipe(paste("cut -d \"",sep,"\" -f ",p.num.col," ",res.file,sep="")), header=T, stringsAsFactors=F)
      p_values=res[,1]	
      p_values=suppressWarnings(as.numeric(p_values))
    }
    p_values <- p_values[which(!is.na(p_values))]
    if (max(p_values) <= 1) {
      chisq <- qchisq(p_values, 1, lower.tail = FALSE)
    }
    chisq <- sort(chisq)
    ppoi <- ppoints(p_values)
    ppoi <- sort(qchisq(1 - ppoi, 1))
    s=median(chisq,na.rm=T)/qchisq(0.5,1)
    rm(chisq)
    ppoi=pchisq(ppoi,1,lower.tail=F)
    p_values=p_values[order(p_values)]
    ppoi=ppoi[order(ppoi)]
    #lighten the data
    if(length(p_values>500000)){
      good=which(p_values<1e-2)
      regood=sample(which(p_values>=1e-2),5000)
      good=c(good,regood)
      p_values=p_values[good]
      ppoi=ppoi[good]
    }
    
    pdf(file=out.file)
    plot(-log10(ppoi), -log10(p_values), xlab = "Expected", ylab = "Observed",pch=19,main=paste(title,"Lambda",round(s,6)))
    abline(a = 0, b = 1)
    garbage<-dev.off() 
  }
