require(lattice) # for use of parallelplot ;

#' Perform correlation and linear regression for a data.frame automatically
#'
#' @param mycor Object to mycor
mycor=function(x,...,digits) UseMethod("mycor")

#' @describeIn mycor for class data.frame
#' @param x A data.frame.
#' @param ... further arguments to be passed to \code{\link{cor.test}}.
#' @param digits integer indicating the number of decimal places (round) or
#'     significant digits (signif) to be used.
#' @return mycor returns as object of class "mycor"
#'
#'     The function summary is used to print a summary of the result. The function
#'     plot is used to plot the results using \code{\link{pairs}} and \code{\link[lattice]{parallelplot}}.
#'
#'     An object of class "mycor:" is a list containing at least following components:
#'     \describe{
#'        \item{df}{a data.frame}
#'        \item{select}{logical vectors returns if columns of df is.numeric}
#'        \item{out}{a list of class "htest" from \code{\link{cor.test}}
#'           between the last paired samples in a data.frame.}
#'        \item{r}{a matrix consist of r values from \code{\link{cor.test}}
#'           between all pairs of numeric data from a data.frame}
#'        \item{p}{a matrix consist of p values from \code{\link{cor.test}}
#'           between all pairs of numeric data from a data.frame}
#'        \item{slope}{a matrix consist of slope values from \code{\link{lm}}
#'           between all pairs of numeric data from a data.frame}
#'        \item{intercept}{a matrix consist of intercept values from \code{\link{lm}}
#'           between all pairs of numeric data from a data.frame}
#'     }
#'  @examples
#'  out=mycor(iris)
#'  plot(out)
#'  plot(out, groups=Species)
#'  plot(out,type=2,groups=species)
#'  plot(out,type=4,groups=species)
#'  out1=mycor(~mpg+disp+wt+hp,data=mtcars,alternative="greater",methods="kendall",
#'             conf.level=0.95)
#'  plot(out1,type=3)
#'  plot(out1,type=4,groups=cyl)
mycor.default=function(x,...,digits=3){
    # select numeric data ony
    select<-(lapply(x,function(x) is.numeric(x))==TRUE)
    num_data=x[select]
    y<-names(num_data)
    ncol=length(num_data)
    # initialize data with matrix filled with zero
    r.value<-matrix(0,ncol,ncol)
    colnames(r.value)<-rownames(r.value)<-y
    p.value<-slope<-intercept<-r.value

    for(i in 1:length(y)){
        for(j in 1:length(y)) {

            out=mylm(num_data[[j]],num_data[[i]],...,digits=digits)
            r.value[i,j]=out$result[1]
            p.value[i,j]=out$result[2]
            slope[j,i]=out$result[3]
            intercept[j,i]=out$result[4]
        }
    }
    result<-list(df=x,select=select,out=out$out,r=r.value,p=p.value,slope=slope,intercept=intercept)
    class(result)<-c("mycor")
    result
}

#' @describeIn mycor for class "formula"
#' @param formula a formula of the form ~ u + v, where each of u and v are
#'       numeric variables giving the data values for one sample. The samples
#'       must be of the same length.
#' @param data A data.frame
mycor.formula=function(formula,data,...,digits=3){
    f=formula
    myt=terms(f,data=data)
    x=labels(myt)
    result=mycor(data[x],...,digits=digits)
    result$df=data
    result
}


#' Correlation and Fitting linear model function for function "mycor"
#'
#' @param y numeric vectors of data values
#' @param x numeric vectors of data values
#' @param ... further arguments to be passed to or from methods.
#' @param digits integer indicating the number of decimal places (round) or
#'     significant digits (signif) to be used.
#' @return mylm returns a list of following components
#'
#'     \describe{
#'        \item{out}{a list of class "htest" from \code{\link{cor.test}}
#'           between the last paired samples in a data.frame.}
#'        \item{result}{a numeric vector of length 4, consist of r and p values
#'              from \code{\link{cor.test}},slope and intercept values from
#'              \code{\link{lm}} between numeric vector y and x}
#'     }
mylm=function(y,x,...,digits=3){
    # performing cor.test
    out1=cor.test(y,x,...)
    my.r.value= round(out1$estimate,digits)
    my.p.value= round(out1$p.value,digits)
    # performing lm to get slope and intercept
    out=lm(y~x)
    result=c(my.r.value,my.p.value,
             round(out$coef[2],max(2,digits-1)),
             round(out$coef[1],max(2,digits-1)))

    # Return list consist of output of cor.test
    # as weel as r, p, slope, intercept
    list(out=out1,result=result)
}

#' Summarizing function for class "mycor"
#'
#' @param object an object of class "mycor", a result of a call to \code{\link{mycor}}.
#' @param ... further arguments to be passed to or from methods.
#' @examples
#' out=mycor(iris)
#' summary(out)
summary.mycor=function(object,...){
    cat("\n")
    cat("$ r value by",object$out$method,"\n\n")
    print(object$r)
    cat("\n$ p value (",object$out$alternative,")\n\n")
    print(object$p)
    cat("\n$ slope \n\n")
    print(object$slope)
    cat("\n$ intercept \n\n")
    print(object$intercept)
}

#' Print function for class "mycor"
#'
#' @param x an object of class "mycor", a result of a call to \code{\link{mycor}}.
#' @param ... further arguments to be passed to or from methods.
#' @examples
#' out=mycor(iris)
#' print(out)
print.mycor=function(x,...) {
    cat("\n")
    cat("$ r value by",x$out$method,"\n\n")
    print(x$r)
    cat("\n$ p value (",x$out$alternative,")\n\n")
    print(x$p)
}

#' Plot for an object of class "mycor"
#' @param x an object of class "mycor"
#' @param ... further arguments to be passed to \code{\link[graphics]{pairs}} or
#'     \code{\link[lattice]{parallelplot}}(in case of "type" argument is 4).
#' @param groups a variable to be evaluated in a data.frame x$df, expected to
#'      act as a grouping variable within each panel, typically used to
#'      distinguish different groups by varying graphical parameters like color and line type.
#' @param type specify the type of plot
#'    \describe{
#'       \item{1}{makes plot with \code{\link[graphics]{pairs}}}
#'       \item{2}{makes plot with \code{\link[graphics]{pairs}} using
#'             \code{\link{panel.hist}} as a diagonal panel}
#'       \item{3}{makes plot with \code{\link[graphics]{pairs}} using
#'             \code{\link{panel.cor}} as a upper panel}
#'       \item{4}{makes plot with \code{\link[lattice]{parallelplot}} using
#'             \code{\link{panel.cor}} as a upper panel}
#'  }
#'  @examples
#'  out=mycor(iris)
#'  plot(out)
#'  plot(out, groups=Species)
#'  plot(out,type=2,groups=species)
#'  out1=mycor(mtcars[1:5],alternative="greater",methods="kendall",
#'             conf.level=0.95)
#'  plot(out1,type=3)
#'  plot(out1,type=4,groups=cyl)
plot.mycor=function(x,...,groups=-1,type=1) {
    # select subset of dataframe
    df=x$df[names(x$select[x$select==TRUE])]
    # in case of type 4, use parallelplot

    name=deparse(substitute(groups))
    # in case no groups specified
    if(name=="-1") {
        if(type==1) pairs(df,...)
        else if(type==2) pairs(df,...,panel=panel.smooth,
                               cex=1, pch=21,bg="light blue",
                               diag.panel=panel.hist,cex.labels=1.5,font.labels=1.5)
        else if(type==3) pairs(df,...,row1attop=FALSE, gap=2,
                               lower.panel=panel.smooth,upper.panel=panel.cor)
        else if(type==4) lattice::parallelplot(df,...)
    }
    #  in case groups specified
    else{
        result=which(grepl(name,colnames(x$df),ignore.case=TRUE))
        if(length(result)<1) {
            cat("no matched column :",name,"\n")
            return()
        }
        else if(length(result)>1) result=result[1]
        focus=x$df[[result]]
        if(type==4){
            mydf=data.frame(df,x$df[result])
        }
        if(length(levels(focus))<=3) mybg=c("red","green3","blue")[factor(focus)]
        else mybg=factor(focus)
        if(type==1) pairs(df,...,pch=21,bg=mybg)

        else if(type==2) pairs(df,...,panel=panel.smooth,
                               cex=1, pch=21,bg=mybg,
                               diag.panel=panel.hist,cex.labels=1.5,font.labels=1.5)
        else if(type==3) pairs(df,...,row1attop=FALSE, gap=2,
                               lower.panel=panel.smooth,upper.panel=panel.cor)
        else if(type==4) lattice::parallelplot(~mydf,...,data=x$df,groups=x$df[[result]],
                          auto.key=TRUE)

    }
}

#' Make plot with histogram for plot of class "mycor"
#' @param x a numeric vector
#' @param ... further arguments to be passed to or from methods.
panel.hist<-function(x,...){
    usr<-par("usr");on.exit(par(usr))
    par(usr=c(usr[1:2],0,1.5))
    h<-hist(x,plot=FALSE)
    breaks<-h$breaks;nB<-length(breaks)
    y<-h$counts;y<-y/max(y)
    rect(breaks[-nB],0,breaks[-1],y,col="cyan",...)
}

#' Make correlation plot for plot of class "mycor"
#' @param x a numeric vector
#' @param y a numeric vector
#' @param digits integer indicating the number of decimal places (round) or
#'     significant digits (signif) to be used.
#' @param prefix a character vector
#' @param cex.cor a numeric variable
panel.cor<-function(x,y,digits=2,prefix="",cex.cor)
{
    usr<-par("usr");on.exit(par(usr))
    par(usr=c(0,1,0,1))
    r<-abs(cor(x,y))
    txt<-format(c(r,0.123456789),digits=digits)[1]
    txt<-paste(prefix,txt,sep="")
    if(missing(cex.cor)) cex<-0.8/strwidth(txt)
    text(0.5,0.5,txt,cex=cex*r)
}
