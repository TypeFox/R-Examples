roundup=function (value, digits) {
    format(round(value, digits), nsmall=digits, scientific=FALSE) 
}
formatInt=function (x, digits, fill="0", ...) {
    formatC(x, format="d", flag=fill, width=digits) 
}
formatDouble=roundup

# don't have to transpose x
mywrite=function(x, ...){
    if (is.list(x)) x=fill.jagged.array(x)
    if (is.null(ncol(x))) i=length(x)
    else i=ncol(x)
    write (t(x), ncolumns=i, ...)
}

# print a matrix/table or a list of them to a latex file as xtable
# note file.name can not have space in it
# e.g. mytex(matrix(0,2,2));
# e.g. mytex(matrix(0,2,2), digits=4);
# e.g. mytex(list(matrix(0,2,2), c(1,1))); 
# default arguments: file.name="temp"; digits=NULL; display=NULL; align="r"; append=FALSE; preamble=""; keep.row.names=TRUE
mytex.begin=function(file.name,preamble=""){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%+%".tex"
    cat ("\\documentclass{article}\n", file=file.name, append=FALSE)
    cat (preamble, file=file.name, append=TRUE)
    cat("\n\\usepackage{geometry}\n", file=file.name, append=TRUE)    
    cat("\n\\begin{document}\n", file=file.name, append=TRUE)    
}
mytex.end=function(file.name){
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    if (!endsWith(file.name,".tex")) file.name=file.name%+%".tex"
    cat ("\n\\end{document}", file=file.name, append=TRUE)
}

# adapted from print.xtable.R
sanitize.text <- function(str) {
    result <- str
    result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
    result <- gsub("$", "\\$", result, fixed = TRUE)
    result <- gsub(">", "$>$", result, fixed = TRUE)
    result <- gsub("<", "$<$", result, fixed = TRUE)
    result <- gsub("|", "$|$", result, fixed = TRUE)
    result <- gsub("{", "\\{", result, fixed = TRUE)
    result <- gsub("}", "\\}", result, fixed = TRUE)
    result <- gsub("%", "\\%", result, fixed = TRUE)
    result <- gsub("&", "\\&", result, fixed = TRUE)
    result <- gsub("_", "\\_", result, fixed = TRUE)
    result <- gsub("#", "\\#", result, fixed = TRUE)
    result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
    result <- gsub("~", "$\\sim$", result, fixed = TRUE) # this is changed by Y.F. 
    result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                   result, fixed = TRUE)
    return(result)
}
sanitize.numbers <- function(x) {
    result <- x
#    if ( math.style.negative ) {
        ## Jake Bowers <jwbowers@illinois.edu> in e-mail
        ## from 2008-08-20 suggested disabling this feature to avoid
        ## problems with LaTeX's dcolumn package.
        ## by Florian Wickelmaier <florian.wickelmaier@uni-tuebingen.de>
        ## in e-mail from 2008-10-03 requested the ability to use the
        ## old behavior.
        for(i in 1:length(x)) {
            result[i] <- gsub("-", "$-$", result[i], fixed = TRUE)
#        }
    }
    return(result)
}

        
mytex=function(dat=NULL, file.name="temp", 
    digits=NULL, display=NULL, align="r", 
    include.rownames=TRUE, include.dup.rownames=FALSE, include.colnames=TRUE,
    col.headers=NULL,
    comment=FALSE, floating=FALSE, 
    lines=TRUE, hline.after=NULL, 
    add.to.row=NULL, 
    sanitize.text.function = NULL, #function(x) x,
    append=FALSE, preamble="", stand.alone=TRUE,
...) {
        
#    if(exists("tablePath") && file.exists(tablePath)) {
#        file.name=tablePath%+%"/"%+%file.name
#    } else {
#        file.name=file.name
#    }    
    
    if (endsWith(file.name,".tex")) file.name=substr(file.name, 1, nchar(file.name)-4)
    if (stand.alone) {
        # create two files, one stand alone and one not, to facilitate debugging latex code
        file.name.2=file.name%+%".tex"
        tmp=strsplit(file.name, split="/")[[1]]
        if (length(tmp)==1) path="./" else path=concatList(tmp[-length(tmp)], "/")
        foldername=path%+%"/stdaln"; 
        if(!file.exists(foldername)) dir.create(foldername) 
        file.name=foldername%+%"/"%+%tmp[length(tmp)]%+%"_stdaln"
    }
    file.name=file.name%+%".tex"
    
    if (include.dup.rownames) include.rownames=F
    
    if(is.data.frame(dat)) dat=list(dat)
    if (!is.list(dat)) dat=list(dat)
    
    if (!append) { #start a new file
        #document tag, preamble etc
        mytex.begin(file.name, preamble)
        if (stand.alone) {
            #empty file
            cat ("", file=file.name.2, append=FALSE)
        }
    } 
    
    if (length(dat)>0) {
        names(dat)=gsub("_"," ",names(dat))
        for (i in 1:length(dat)) {
            dat1 = dat[[i]]        
            .ncol=ncol(dat1)
            if (is.null(.ncol)) {
                if (is.null(nrow(dat1))) .ncol=1
                else .ncol=nrow(dat1)
            }
            
            top.1=concatList("& \\multicolumn{1}{c}{"%+%sanitize.text(colnames(dat1))%+%"} ") %+% "\\\\ \n"%+% # center aligned column titles
                "\\hline\n" # insert at the beginning of table, "\n" is added so that there is no need to keep it in col.title
            if(!include.rownames) top.1=substr(top.1, 2,10000)
            top=if(!include.colnames)  "" else top.1
                
            if (include.colnames & is.null(hline.after)) hline.after=c(nrow(dat1)) # cannot use default due to add.to.row    
            include.colnames=FALSE
            if (!is.null(col.headers)) top=col.headers%+%top else top="\\hline  "%+%top
            
            if (is.null(add.to.row)) {
                add.to.row=list(list(0), top)
            } else {
                add.to.row=list(c(list(0), add.to.row[[1]]), c(top, add.to.row[[1]]))
            }
        
            if (!is.matrix(dat1) & is.character(dat1)) {
                cat (dat1%+%"\n\n\n", file=file.name, append=TRUE)
                if(stand.alone) cat (dat1%+%"\n\n\n", file=file.name.2, append=TRUE)
            } else {        
                if (is.vector(dat1)) dat1=as.matrix(dat1)
                
                if (stand.alone & length(dat)>1) cat (names(dat)[i]%+%"\n\n", file=file.name, append=TRUE)
                if (!is.null(dat1)) {
                    if (!is.null(attr(dat1,"caption"))) caption=attr(dat1,"caption") else caption=NULL
                    
                    if (include.dup.rownames & !is.null(rownames(dat1))) {
                        tmp=suppressWarnings(data.frame(rownames(dat1),data.frame(dat1))) # warnings about duplicate row names
                        if (!is.null(colnames(dat1))) colnames(tmp)[-1]=colnames(dat1)
                        dat1=tmp
                        .ncol=.ncol+1
                    }
                    if (is.null(hline.after)) {
                        if (lines) hline.after=c(-1,0,nrow(dat1)) else hline.after=c(nrow(dat1))
                    }
                    if (length(align)==.ncol+1) {
                        # no need to do anything
                    } else {
                        if (length(align)==1) align=rep(align,.ncol+1)
                        if (include.dup.rownames) align[2]="l" else align[1]="l" # col names
                    }
                    print(..., xtable::xtable(dat1, 
                            digits=(if(is.null(digits)) rep(3, .ncol+1) else digits), # cannot use ifelse here!!!
                            display=(if(is.null(display)) rep("f", .ncol+1) else display), # or here
                            align=align, caption=caption, ...), 
                        hline.after=hline.after, type = "latex", file = file.name, append = TRUE, floating = floating, 
                        include.rownames=include.rownames, include.colnames=include.colnames, comment=comment, 
                        add.to.row=add.to.row, sanitize.text.functio =sanitize.text.function )
                    if (stand.alone) print(..., xtable::xtable(dat1, 
                            digits=(if(is.null(digits)) rep(3, .ncol+1) else digits), # cannot use ifelse here!!!
                            display=(if(is.null(display)) rep("f", .ncol+1) else display), # or here
                            align=align, caption=caption, ...), 
                        hline.after=hline.after, type = "latex", file = file.name.2, append = TRUE, floating = floating, 
                        include.rownames=include.rownames, include.colnames=include.colnames, comment=comment, 
                        add.to.row=add.to.row, sanitize.text.functio =sanitize.text.function )
                }
                cat ("\n", file=file.name, append=TRUE)
                if(stand.alone) cat ("\n", file=file.name.2, append=TRUE)
            }
        
        }
    }
    
    if(!append) {
        mytex.end(file.name)
    }
    cat ("Saving data to "%+%getwd()%+%"/"%+%file.name%+%"\n")
}
#x=matrix(0,2,2)
#attr(x,"caption")="cap"
#mytex(x, floating=TRUE)


# write a table that contains mean and sd to temp.tex in the current working directory, getwd()
# models can be a list of models, or a single model
make.latex.coef.table = function (models, model.names=NULL, row.major=FALSE, round.digits=NULL) {
# e.g.: models=list(gam1, gam2); round.digits= c(3,3,3,3,3); model.names=c("gam1", "gam2");  row.major=TRUE   
    if (! ("list" %in% class (models) ) ) {models=list(models)}
    
    numParams = nrow (getFixedEf(models[[1]]))
    numModels = length (models)
    
    if (is.null (model.names)) {model.names=rep("",numModels)}
    if (is.null(round.digits)) round.digits=rep(3,numParams)    
    
    coef.table = mysapply (1:numModels, function (i.model) {
        temp = getFixedEf(models[[i.model]]) [,1:2,drop=FALSE]
        for (i.param in 1:numParams) {
            temp[i.param,] = round (temp[i.param,], round.digits[i.param])
        }
        temp2 = paste (format(temp[,1]), "(", format(temp[,2]), ")")
        names (temp2) = dimnames(temp)[[1]]
        temp2
    })
    dimnames (coef.table)[[1]] = model.names
    
    if (row.major) mytex ( coef.table, align="r" ) 
    else mytex (t(coef.table), align="r") 
}


# default row.names to FALSE
# file name needs no file extension
mywrite.csv = function(x, file="tmp", row.names=FALSE, digits=NULL, ...) {  
    if (!is.null(digits)) {
        if(length(digits)==1) {
            x=round(x,digits)
        } else {
            for (i in 1:ncol(x)) {
                x[,i]=round(x[,i], digits[i])
            }                
        }
    }
    print("Writing csv file to "%+%getwd()%+%"/"%+%file%+%".csv", quote=FALSE)
    write.csv(x, file=file%+%".csv", row.names=row.names, ...)
}


myprint <- function(object, ...) UseMethod("myprint") 

# this function is placed at the bottom of the file because it contains "\""), which makes all the following line be miss-interpreted as being in quotes
myprint.default = function (..., newline=TRUE, digits=3) {   
    digits.save=getOption("digits")
    options(digits=digits)
    object <- as.list(substitute(list(...)))[-1]
    x=list(...)
    for (i in 1:length(x)) {
        if (is(x[[i]],"formula")) {cat(as.character(x[[i]]), "; "); next}
        tmpname <- deparse(object[[i]])[1]
        #str(tmpname)
        #str(gsub("\\\\","\\",gsub("\"", "", tmpname)))
        #str(x[[i]])
        #if (gsub("\\\\","\\",gsub("\"", "", tmpname))!=x[[i]]) {
        if (contain(tmpname, "\"") | contain(tmpname, "\\")) {
            for (a in x[[i]]) cat(a)
        } else {
            cat (tmpname %+% " = ")
            for (a in x[[i]]) cat(a,"") # by putting "" there, a space is introduced b/c cat prints a sep
            if (i!=length(x)) cat ("; ")
        }
    }
    if (newline)  cat("\n")
    options(digits=digits.save)
}
#a="hello"; b="you"; myprint (a, b); myprint ("test"); myprint.default ("\t")
