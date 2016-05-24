#' Place two or more ztables or figures side by side in Latex or HTML format
#'
#' Place two or more ztables or figures side by side in Latex or HTML format.
#' Requires Latex "boxedminipage" package in preamble.
#' The ztable for this purpose should be made by function ztable with tabular="TRUE".
#' @param width a numeric vector specifies the width to which the tables or
#'          figures should be scaled
#' @param listTables a list consists of object of "ztable" or valid figure name
#' @param type Type of table to produce. Possible values for type are "latex" or
#'        "html". Default value is "latex".
#' @examples
#' require(ztable)
#' z=ztable(head(mtcars[1:3]),tabular=TRUE)
#' parallelTables(c(0.4,0.3),list(z,z))
#' parallelTables(c(0.5,0.5),list(z,z))
#' parallelTables(c(0.5,0.5),list(z,z,type="html"))
#' z1=ztable(head(iris[1:3]),turn=TRUE,angle=10,zebra=1)
#' z2=ztable(head(iris[1:3]),turn=TRUE,angle=-10,zebra=2)
#' parallelTables(c(0.5,0.5),list(z1,z2))
parallelTables=function(width,listTables,type="latex"){
    a=length(width)
    if(class(listTables)!="list"){
        cat("\nThe 2nd parameter listTables sholud be a list\n")
        return(invisible())
    }
    b=length(listTables)
    if(a!=b) {
        cat("\nLengths of width and tables are different\n")
        cat(paste("length of width=",a,",length of tables=",b,sep=""))
        return(invisible())
    }
    if(type=="html") parallelTablesHTML(width,listTables)
    else parallelTablesLatex(width,listTables)
}

#' Place two or more ztables or figures side by side in Latex format
#'
#' Place two or more ztables or figures side by side in HTML format.
#' The ztable for this purpose should be made by function ztable with tabular="TRUE".
#' @param width a numeric vector specifies the width to which the tables or
#'          figures should be scaled
#' @param listTables a list consists of object of "ztable" or valid figure name
parallelTablesLatex=function(width,listTables){
    a=length(width)
    cat("\\begin{table}[!htb]\n")
    for(i in 1:a){
        cat(paste("\\begin{minipage}{",width[i],"\\linewidth}\n\\centering\n",
                  sep=""))
        if(class(listTables[[i]])=="ztable") print(listTables[[i]],type="latex")
        else if(class(listTables[[i]])=="character") {
            cat(paste("\\includegraphics[width=1\\linewidth]{",
                      listTables[[i]],"}\n",sep=""))
        }
        cat("\\end{minipage}\n")
    }
    cat("\\end{table}")
}

#' Place two or more ztables or figures side by side in HTML format
#'
#' Place two or more ztables or figures side by side in HTML format.
#' The ztable for this purpose should be made by function ztable with tabular="TRUE".
#' @param width a numeric vector specifies the width to which the tables or
#'          figures should be scaled
#' @param listTables a list consists of object of "ztable" or valid figure name
parallelTablesHTML=function(width,listTables){
    a=length(width)
    cat("<table width=\"100%\" cellspacing=\"5px\" cellpadding=\"5px\" border=\"0\">
        \n<colgroup>\n")
    for(i in 1:a) cat(paste("<col width=",
                            ifelse(width[i]<=1,width[i]*100,width[i]),
                            "%>\n",sep=""))
    cat("</colgroup>\n<tr>")
    for(i in 1 :a){
        cat("<td>")
        if(class(listTables[[i]])=="ztable") print(listTables[[i]],type="html")
        else if(class(listTables[[i]])=="character")
            cat(paste("<img src=\"",listTables[[i]],"\">",sep=""))
        cat("</td>\n")
    }
    cat("</tr>\n</table>\n")
}
