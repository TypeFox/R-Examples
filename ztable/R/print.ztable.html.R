#' Find rgb value from color name
#'
#'@param name a valid color name
#'@return rgb value
name2rgb=function(name){
    number=grep(paste("^",name,sep=""),ztable::zcolors$name)
    if(length(number)<1) result="white"
    else{
        rgb=ztable::zcolors[number[1],2]
        result=paste("#",rgb,sep="")
    }

    result
}

#' Delete first components of align
#'
#' @param align A character for define the align of column in Latex format
align2nd=function(align){
    if(substr(align,1,1)=="|") {
        result=substr(align,2,nchar(align))
        result=align2nd(result)
    } else result=substr(align,2,nchar(align))
    result
}

#' Count the number of align
#'
#' @param align A character for define the align of column in Latex format
alignCount=function(align){
    result=unlist(strsplit(align,"|",fixed=TRUE))
    temp=c()
    for(i in 1:length(result)) temp=paste(temp,result[i],sep="")
    nchar(temp)
}


#' Check the validity of align
#'
#' @param align A character for define the align of column in Latex format
#' @param ncount An integer equals of ncol function
#' @param addrow An integer
alignCheck=function(align,ncount,addrow){
    count=alignCount(align)
    #cat("align=",align,"count=",count,"\n")
    while(count != (ncount+addrow)){
        if(count< (ncount+addrow)) align=paste(align,"c",sep="")
        else if(count > (ncount+addrow)) align=align2nd(align)
        count=alignCount(align)
        #cat("align=",align,"count=",count,"\n")
    }
    result=align
    result
}


#' Convert the align in Latex format to html format
#'
#' @param align A character of align in Latex format
align2html=function(align){
    result=c()
    for(i in 1:nchar(align)){
        temp=substr(align,i,i)
        if(temp=="|") next
        temp=ifelse(temp=="l","left",ifelse(temp=="r","right","center"))
        result=c(result,temp)
    }
    result
}

#' Extract column position information only(without vertical line specifier)
#'
#' @param align A character string indicating align for latex table
extractAlign=function(align){
    result=c()
    for(i in 1:nchar(align)){
        temp=substr(align,i,i)
        if(temp=="|") next
        result=c(result,temp)
    }
    temp=result[1]
    if(length(result)>1)
        for(i in 2:length(result)) {
           temp=paste(temp,result[i],sep="")
        }
    temp
}

#' Add or delete vertical lines in a ztable
#'
#' @param z An object of ztable
#' @param type An integer or one of c("none","all")
#' @param add An integer vector indicating columns where the width of vertical lines added
#' @param del An integer vector indicating columns where the width of vertical lines subtracted
vlines=function(z,type=NULL,add=NULL,del=NULL){

    if(is.null(type) & is.null(add) & is.null(del)) {
        cat("\nvlines : add or delete vertical lines to a ztable\n
Usage: type must be one of these or NULL: 0-1 or \"none\",\"all\"\n
       add and del: An integer vector indicating position to add or delete vertical line(s)\n")

        return(z)
    }
    align=extractAlign(z$align)
    vlines=align2lines(z$align)
    colcount=colGroupCount(z)
    addrow=ifelse(z$include.rownames,1,0)
    #align=alignCheck(align,ncol(z$x),addrow)
    count=nchar(align)

    if(!is.null(type)) {
        vltype=NULL
        if(!is.numeric(type)) {
            if(toupper(type) == "NONE") vltype=0
            else if(toupper(type) == "ALL") vltype=1
            else return(z)
        }
        if((type>=0) & (type<=1)) vltype=type
        if(vltype==0) vlines=rep(0,count+1)
        else vlines=rep(1,count+1) #vltype=1

    }
    if(!is.null(add)){
        if(is.numeric(add)){
            for(i in 1:length(add)) {
                if(add[i]<1 | add[i]>(count+1)) next
                vlines[add[i]]=vlines[add[i]]+1
            }
        }
    }
    if(!is.null(del)){
        if(is.numeric(del)){
            for(i in 1:length(del)){
                if(del[i]<1 | del[i]>(count+1)) next
                if(vlines[del[i]]>0) vlines[del[i]]=vlines[del[i]]-1
            }
        }
    }
    newalign=vline2align(align,vlines)
    z$align=newalign
    z
}


#' Add or delete horizontal lines in a ztable
#'
#' @param z An object of ztable
#' @param type An integer or one of c("none","all")
#' @param add An integer vector indicating rows where the horzontal lines added
#' @param del An integer vector indicating rows where the horizontal lines deleted
hlines=function(z,type=NULL,add=NULL,del=NULL){

    if(is.null(type) & is.null(add) & is.null(del)) {
        cat("\nhlines : add or delete horizontal lines to a ztable\n
            Usage: type must be one of these or NULL: 0-1 or \"none\",\"all\"\n
            add and del: An integer vector indicating position to add or delete horizontal line(s)\n")

        return(z)
    }

    count=nrow(z$x)
    if(!is.null(z$hline.after)) result=z$hline.after
    else result=c(-1,0,count)

    if(!is.null(type)) {
        if(!is.numeric(type)) {
            if(toupper(type) == "NONE") hltype=0
            else if(toupper(type) == "ALL") hltype=1
            else return(z)
        }
        if((type>=0) & (type<=1)) hltype=type
        if(hltype==0) result=c(-1,0,count)
        else result=c(-1,0,1:count)

    }
    if(!is.null(add)){
        if(is.numeric(add)){
            for(i in 1:length(add)) {
                result=c(result,add)
            }
        }
    }
    if(!is.null(del)){
        if(is.numeric(del)){
            result1=c()
            for(i in 1:length(result)){
                if(!(result[i] %in% del)) result1=c(result1,result[i])
            }
            result=result1
        }
    }
    z$hline.after=result
    z
}

#' Make a latex "align" from a string and vertical line specifier
#'
#' @param align A character string indicating align of latex table
#' @param vlines An integer vector indicating vertical line position
vline2align=function(align,vlines){
    newalign=c()
    for(i in 1:nchar(align)) {
        if(vlines[i]>0) for(j in 1:vlines[i]) newalign=c(newalign,"|")
        newalign=c(newalign,substr(align,i,i))
    }
    last=vlines[length(vlines)]
    if(last>0) for(j in 1:last) newalign=c(newalign,"|")
    temp=newalign[1]
    if(length(newalign)>1)
        for(i in 2:length(newalign)) {
            temp=paste(temp,newalign[i],sep="")
        }
    temp
}

#' count the vertical column lines from align of Latex format
#'
#' @param align A string of align Latex format
#' @return a numeric vector consists of vertical lines of each column
align2lines=function(align){
    result=c()
    length=nchar(align)
    count=0
    number=alignCount(align)
    for(i in 1:length){
        temp=substr(align,1,1)
        if(temp=="|") {
            count=count+1
            if(i==length) result=c(result,count)
        }
        else{
            result=c(result,count)
            count=0
        }
        align=substr(align,2,nchar(align))
    }
    if(length(result)==number) result=c(result,0)
    result
}

#' Make a charater string indicating the alignment of components of table.
#'
#' @param z An object of ztable
getNewAlign=function(z){
    #cat("z$align=",z$align,"\n")
    if(is.null(z$cgroup)) return(z$align)
    lines=align2lines(z$align)
    exAlign=extractAlign(z$align)
    ncount=ncol(z$x)
    addrow=ifelse(z$include.rownames,1,0)
    colCount=colGroupCount(z)
    result=c()
    start=2-addrow
    # Add column group align "c" if lines
    for(i in 1:length(colCount)){
        #cat("start=",start,"stop=",colCount[i]+addrow,",")
        result=paste(result,substr(exAlign,start=start,stop=(colCount[i]+1)),sep="")
        #cat("i=",i,",start=",start,"stop=",(colCount[i]+1),",result=",result)
        start=colCount[i]+2
        #cat(",line[start]=",start,"\n")
        if(lines[start]==0) result=paste(result,"c",sep="")
        #cat("result=",result,"\n")
    }
    result
    if(colCount[length(colCount)]<ncount)
        result=paste(result,substr(exAlign,start=start,stop=nchar(z$align)),sep="")
    result
    newlines=c()
    for(i in 1:length(lines)){
        if(i==1) newlines=lines[1]
        else newlines=c(newlines,lines[i])
        if((i-1) %in% colCount[-length(colCount)])
            if(lines[i+1]==0) newlines=c(newlines,0)
    }
    temp=c()
    for(i in 1:length(newlines)){
        if(newlines[i]>0) for(j in 1:newlines[i]) temp=paste(temp,"|",sep="")
        if(i>nchar(result)) break
        temp=paste(temp,substr(result,start=i,stop=i),sep="")
    }
    #temp=paste(temp,"c",sep="")
    temp
}


#' print html style
myhtmlStyle=function(){
    cat("<head>")
    cat("<style>
        table {
              font-family: serif;
              text-align: right;}
        th {
              padding: 1px 1px 5px 5px;
	        }
        td {
             padding: 1px 1px 5px 5px; }
      </style>")
    cat("</head>")
}

#' Print HTML head if ztable object a has a colgroup
#'
#' @param z An object of ztable
printHTMLHead=function(z){
    if(is.null(z$cgroup)) return
    if(is.null(z$n.cgroup)) return
    #colCount=colGroupCount(z)
    ncount=ncol(z$x)
    addrow=ifelse(z$include.rownames,1,0)
    cGroupSpan=cGroupSpan(z)
    totalCol=totalCol(z)

    vlines=align2lines(z$align)

    for(i in 1:nrow(z$cgroup)){
        cat("<tr>\n")
        if(z$include.rownames) {
            cat("<td style=\"")
            if(i==1) cat("border-top: 2px solid gray; border-bottom: hidden;")
            cat(paste(" border-left: ",vlines[1],"px solid black;",sep=""))
            if(z$cgroupcolor[i,1]!="white")
                cat(paste("background-color: ",name2rgb(z$cgroupcolor[i,1]),sep=""))
            cat("\"> </td>\n")
        }
        colSum=1
        for(j in 1:ncol(z$cgroup)) {
            if(is.na(z$cgroup[i,j])) {
                cat("<td colspan=\"",cGroupSpan[i,j],"\" align=\"center\" ")
                cat("style=\"")
                if(i==1) cat("border-top: 2px solid gray;")
                cat("border-bottom: hidden;")
                cat(paste(" border-left: ",vlines[colSum+1],"px solid black;",sep=""))
                colSum=colSum+cGroupSpan[i,j]
                #if(colSum==ncol(z$x)+1)
                cat(paste("border-right:",vlines[colSum+1],"px solid black;",sep=""))
                if(z$cgroupcolor[i,j+1]!="white")
                    cat(paste("background-color: ",name2rgb(z$cgroupcolor[i,j+1]),";",sep=""))
                cat(paste("\"></td>\n",sep=""))
            } else {
                cat("<td colspan=\"",cGroupSpan[i,j],"\" align=\"center\" ")
                if(z$colnames.bold) cat("style=\"font-weight: bold;")
                else cat("style=\"font-weight: normal;")
                if(i==1) cat("border-top: 2px solid gray;")
                if(z$cgroup[i,j]!="") cat(" border-bottom: 1px solid gray;")
                else cat(" border-bottom: hidden;")
                cat(paste(" border-left: ",vlines[colSum+1],"px solid black;",sep=""))
                colSum=colSum+cGroupSpan[i,j]
                if(colSum==ncol(z$x)+1)
                cat(paste("border-right:",vlines[colSum+1],"px solid black;",sep=""))
                if(z$cgroupcolor[i,j+1]!="white")
                    cat(paste("background-color: ",name2rgb(z$cgroupcolor[i,j+1]),";",sep=""))
                cat(paste("\">",z$cgroup[i,j],"</td>\n",sep=""))
            }
            #if((j < ncol(z$cgroup)) & ((colSum+j-1)<totalCol)) {
            if(j < ncol(z$cgroup)) {
                result=colSum+1
                if(result<=length(vlines)) {
                    if(vlines[result]==0){
                        cat("<td style=\"")
                        if(i==1) cat("border-top: 2px solid gray;")
                        cat("border-bottom: hidden\">&nbsp;</td>\n")
                    }
                }
            }
        }
        cat("</tr>\n")
    }
}


#' Print an object of class "ztable" to html table
#'
#' @param z An object of class "ztable"
#' @param xdata A formatted data.frame
ztable2html=function(z,xdata){
    ncount=ncol(z$x)
    addrow=ifelse(z$include.rownames,1,0)
     # caption position
    if(z$caption.position=="r") cposition="right"
    else if(z$caption.position=="l") cposition="left"
    else cposition="center"
    fontsize=ifelse(z$size>=5,11+(z$size-5)*2,10-(4-z$size))
    headingsize=fontsize-2

    rgroupcount=0
    printrgroup=1
    if(!is.null(z$n.rgroup)){
        if(length(z$n.rgroup)>1) {
            for(i in 2:length(z$n.rgroup)) {
                printrgroup=c(printrgroup,printrgroup[length(printrgroup)]+z$n.rgroup[i-1])
            }
        }
        rgroupcount=1
    }

    NewAlign=getNewAlign(z)
    totalCol=totalCol(z)
    colCount=colGroupCount(z)

    rgroupcount=0
    printrgroup=1
    if(!is.null(z$n.rgroup)){
        if(length(z$n.rgroup)>1) {
            for(i in 2:length(z$n.rgroup)) {
                printrgroup=c(printrgroup,printrgroup[length(printrgroup)]+z$n.rgroup[i-1])
            }
        }
        rgroupcount=1
    }

    # table position
    if(z$position=="flushleft") tposition="left"
    else if(z$position=="flushright") tposition="right"
    else tposition="center"
    #cat("<table class='gmisc_table'")
    myhtmlStyle()
    cat("<table ")
    cat(paste("align=\"",tposition,"\" style=\"border-collapse: collapse; caption-side:",
              z$caption.placement,"; font-size:",as.integer(fontsize),"pt;\">",sep=""))
    cat(paste("<caption style=\"text-align:",cposition,";",sep=""))
    if(z$caption.bold) cat("font-weight: bold")
    cat(paste("\">",z$caption,"</caption>",sep=""))
    if((z$show.heading==TRUE) & (!is.null(attr(z$x,"heading")))) {
        head=attr(z$x,"heading")
        for(i in 1:length(head)) {
            if(nchar(head[i])<1) next
            cat(paste("<tr>\n<td style=\"border-top: hidden; font-size: ",
                      as.integer(headingsize),"pt; padding: 0px 0px;\" colspan=\"",ncount+addrow,
                      "\"  align=\"left\" >",head[i],sep=""))
            cat("</td>\n</tr>\n")

        }
    }
    vlines=align2lines(z$align)
    printtop=1
    if(!is.null(z$cgroup)) {
        printHTMLHead(z)
        printtop=0
    }
    if(z$include.colnames) {
        cat("<tr>\n")
        subcolnames=ifelse(is.null(z$subcolnames),0,1)
        if(z$include.rownames) {
            result=1
            if(!is.null(isspanCol(z,1,1)))
                cat(paste("<th colspan=\"",isspanCol(z,1,1),"\"",sep=""))
            else if(!is.null(isspanRow(z,1,1))){
                result=isspanRow(z,1,1)
                if(result>0) cat(paste("<th rowspan=\"",result,"\"",sep=""))
            } else cat("<th ")
            cat(paste("style=\"border-left: ",vlines[1],
                                  "px solid black;",
                                  "background-color: ",name2rgb(z$cellcolor[1,1]),";",sep=""))
            if(printtop) cat("border-top: 2px solid gray;")
            if(subcolnames==0) cat("border-bottom: 1px solid gray;")
            else cat("border-bottom: hidden;")
            cat(paste("\">&nbsp;</th>\n",sep=""))
        }
        colpos=align2html(z$align)
        for(i in 1:ncol(z$x)) {
            result=1
            if(!is.null(isspanCol(z,1,(i+1)))){
                result=isspanCol(z,1,(i+1))
                if(result>0) cat(paste("<th colspan=\"",result,"\"",sep=""))
                else if(result==0) next
            } else if(!is.null(isspanRow(z,1,(i+1)))){
                result=isspanRow(z,1,(i+1))
                if(result>0) cat(paste("<th rowspan=\"",isspanRow(z,1,(i+1)),"\"",sep=""))
                else cat("<th")
            } else cat("<th ")
            if(result!=0){
                 cat("<th ")
                 drawbottom=0
                 if((subcolnames==1)) {
                     if(is.na(z$subcolnames[i])){
                         cat("rowspan=\"2\" ")
                         drawbottom=1
                     }
                 }
                 cat(paste("align=\"center\" ",sep=""))
                 if(z$colnames.bold) cat("style=\"font-weight: bold;")
                 else cat("style=\"font-weight: normal;")
                 cat(paste("border-left: ",vlines[i+1],"px solid black;",sep=""))
                 if((i==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                     cat(paste("border-right:",vlines[i+2],"px solid black;",sep=""))
                 if((subcolnames==0) | (subcolnames+drawbottom==2))
                     cat("border-bottom: 1px solid gray;")
                 else cat("border-bottom: hidden;")
                 if(printtop) cat("border-top: 2px solid gray;")
                 if(z$cellcolor[1,i+1]!="white")
                     cat(paste("background-color: ",name2rgb(z$cellcolor[1,i+1]),";",sep=""))
                                  cat(paste("\">",colnames(z$x)[i],"</th>\n",sep=""))
                 if(i %in% colCount[-length(colCount)]) {
                     if(vlines[i+2]==0){
                        if(subcolnames==0) cat("<th style=\"border-bottom: 1px solid gray;")
                        else cat("<th style=\"border-bottom: hidden;")
                        if(printtop) cat("border-top: 2px solid gray; ")
                        if((z$cellcolor[1,i+1]!="white") & (z$cellcolor[1,i+1]==z$cellcolor[1,i+2]))
                            cat("background-color: ",name2rgb(z$cellcolor[1,i+1]),";")
                        cat("\">&nbsp;</th>\n")
                     }
                 }
            }
        }
        cat("</tr>\n")
        printtop=0
        if(subcolnames){
            cat("<tr>\n")
            if(addrow) {
                cat(paste("<th style=\"border-left: ",vlines[1],
                          "px solid black;","border-bottom: 1px solid gray;",
                          "background-color: ",name2rgb(z$cellcolor[1,1]),";",sep=""))
                cat(paste("\">&nbsp;</th>\n",sep=""))
            }
            for(i in 1:length(z$subcolnames)){
                if(is.na(z$subcolnames[i])) {
                    if(vlines[i+2]==0){
                        if(i!=length(z$subcolnames)){
                            cat("<th style=\"border-bottom: 1px solid gray;")
                            #if(printtop) cat("border-top: 2px solid gray;")
                            if((z$cellcolor[1,i+1]!="white") & (z$cellcolor[1,i+1]==z$cellcolor[1,i+2]))
                                cat("background-color: ",name2rgb(z$cellcolor[1,i+1]),";")
                            cat("\">&nbsp;</th>\n")
                        }
                    }
                    next
                }
                cat("<th align=\"center\" ")
                if(z$colnames.bold) cat("style=\"font-weight: bold;")
                else cat("style=\"font-weight: normal;")
                cat(paste("border-left: ",vlines[i+1],"px solid black;",sep=""))
                if((i==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                    cat(paste("border-right:",vlines[i+2],"px solid black;",sep=""))
                cat("border-bottom: 1px solid gray;")
                if(z$cellcolor[1,i+1]!="white")
                    cat(paste("background-color: ",name2rgb(z$cellcolor[1,i+1]),";",sep=""))
                cat(paste("\">",z$subcolnames[i],"</th>\n",sep=""))
                if(i %in% colCount[-length(colCount)]) {
                    if(vlines[i+2]==0){
                        cat("<th style=\"border-bottom: 1px solid gray;")
                        #if(printtop) cat("border-top: 2px solid gray;")
                        if((z$cellcolor[1,i+1]!="white") & (z$cellcolor[1,i+1]==z$cellcolor[1,i+2]))
                            cat("background-color: ",name2rgb(z$cellcolor[1,i+1]),";")
                        cat("\">&nbsp;</th>\n")
                    }
                }
            }
            cat("</tr>\n")
        }
    }
    colpos=align2html(z$align)
    rgroupprinted=0
    for(i in 1:nrow(z$x)){
        if(rgroupcount>0) {

            if(i %in% printrgroup) {
                rgroupprinted=1
                if(is.null(z$cspan.rgroup)){
                    temp=paste("<tr>\n<td colspan=\"",totalCol,
                               "\"  align=\"left\""," style=\"font-weight: bold;",sep="")
                    if(z$colcolor[1]!="white")
                        temp=paste(temp,"background-color:",name2rgb(z$colcolor[1]),";",sep="")
                    temp=paste(temp," border-left: ",vlines[1],"px solid black; ",sep="")
                    temp=paste(temp,"border-right:",vlines[ncol(z$x)+2],"px solid black;",sep="")
                    temp=paste(temp,"border-bottom: 1px solid black;",sep="")
                    temp=paste(temp,"border-top: 1px solid black;",sep="")
                    temp=paste(temp,"\">",z$rgroup[rgroupcount],"</td>\n",sep="")
                }
                else {
                    if(z$cspan.rgroup==1) {
                        temp=paste("<tr>\n<td align=\"left\""," style=\"font-weight: bold;",sep="")
                        if(z$colcolor[1]!="white")
                            temp=paste(temp,"background-color:",name2rgb(z$colcolor[1]),";",sep="")
                        temp=paste(temp," border-left: ",vlines[1],"px solid black; ",sep="")
                        #temp=paste(temp,"border-bottom: 1px solid black;",sep="")
                        if(i!=1) temp=paste(temp,"border-top: hidden; ",sep="")
                        if(!is.null(z$hline.after)){
                            if((i-1) %in% z$hline.after)
                                temp=paste(temp,"border-top: 1px solid black;")
                        }
                        temp=paste(temp,"\">",z$rgroup[rgroupcount],"</td>\n",sep="")
                        for(j in 1:(ncount+addrow-1)){
                            temp1=paste("<td style=\"border-left: ",
                                        vlines[j+1],"px solid black; ",sep="")
                            if(!is.null(z$hline.after)){
                                if((i-1) %in% z$hline.after)
                                    temp1=paste(temp1,"border-top: 1px solid black;")
                            }
                            else if(i!=1) temp1=paste(temp1,"border-top: hidden; ",sep="")
                            if((j==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                                temp1=paste(temp1,"border-right:",vlines[j+2],"px solid black;",sep="")
                            if(!is.null(z$colcolor)) {
                                if(z$colcolor[j+1]!="white")
                                    temp1=paste(temp1,"background-color:",
                                            name2rgb(z$colcolor[j+1])," ",sep="")
                            }
                            temp1=paste(temp1,"\"></td>\n",sep="")
                            if(is.null(isspanRow(z,i+1,j+1))) temp=paste(temp,temp1,sep="")
                            else if(isspanRow(z,i+1,j+1)>0) temp=paste(temp,temp1,sep="")

                            if(!is.null(colCount)){
                                if(j %in% colCount[-length(colCount)]) {
                                    if(vlines[j+2]==0){
                                        #if((z$cellcolor[i+1,j+1]!="white")&(z$cellcolor[i+1,j+1]==z$cellcolor[i+1,j+2]))
                                        #    temp=paste(temp,"<td style=\"background-color: ",
                                        #           name2rgb(z$cellcolor[i+1,j+1]),"\"></td>\n",
                                        #           sep="")
                                        #else temp=paste(temp,"<td></td>\n",sep="")

                                        temp=paste(temp,"<td",sep="")
                                        if(i!=1) temp=paste(temp,"style=\"border-top: hidden;\"")
                                        temp=paste(temp,"></td>\n",sep="")

                                    }
                                }
                            }
                        }
                    } else {
                        if(z$cspan.rgroup<1 | z$cspan.rgroup>(ncount+addrow))
                            z$cspan.rgroup=ncount+addrow

                        temp=paste("<tr>\n<td colspan=\"",z$cspan.rgroup,
                                   "\"  align=\"left\""," style=\"font-weight: bold;",sep="")
                        if(z$colcolor[1]!="white")
                            temp=paste(temp,"background-color:",name2rgb(z$colcolor[1]),";",sep="")
                        temp=paste(temp," border-left: ",vlines[1],"px solid black; ",sep="")
                        temp=paste(temp,"border-bottom: 1px solid black;",sep="")
                        temp=paste(temp,"border-top: 1px solid black;",sep="")
                        if(!is.null(z$hline.after)){
                            if((i-1) %in% z$hline.after)
                                temp=paste(temp,"border-top: 1px solid black;")
                        }
                        temp=paste(temp,"\">",z$rgroup[rgroupcount],"</td>\n",sep="")

                        if(z$cspan.rgroup<(ncount+addrow)) {
                            for(j in (z$cspan.rgroup):(ncount+addrow-1)) {
                                temp1=paste("<td style=\"border-left: ",
                                            vlines[j+1],"px solid black; ",sep="")
                                if((j==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                                    temp1=paste(temp1,"border-right:",vlines[j+2],"px solid black;",sep="")
                                #temp1=paste(temp1,"border-bottom: 1px solid black;",sep="")
                                #temp1=paste(temp1,"border-top: 1px solid black;",sep="")
                                if(!is.null(z$hline.after)){
                                    if((i-1) %in% z$hline.after)
                                        temp1=paste(temp1,"border-top: 1px solid black;")
                                }
                                else if(i!=1) temp1=paste(temp1,"border-top: hidden; ",sep="")
                                if(!is.null(z$colcolor)) {
                                    if(z$colcolor[j+1]!="white")
                                        temp1=paste(temp1,"background-color:",
                                                name2rgb(z$colcolor[j+1])," ",sep="")
                                }
                                temp1=paste(temp1,"\"></td>\n",sep="")
                                if(is.null(isspanRow(z,i+1,j+1))) temp=paste(temp,temp1,sep="")
                                else if(isspanRow(z,i+1,j+1)>0) temp=paste(temp,temp1,sep="")

                                if(!is.null(colCount)){
                                    if(j %in% colCount[-length(colCount)]) {
                                        if(vlines[j+2]==0) {
                                            #if((z$cellcolor[i+1,j+1]!="white")&(z$cellcolor[i+1,j+1]==z$cellcolor[i+1,j+2]))
                                            #    temp=paste(temp,"<td style=\"background-color: ",
                                            #           name2rgb(z$cellcolor[i+1,j+1]),"\"></td>\n",
                                            #           sep="")
                                            #else temp=paste(temp,"<td></td>\n",sep="")
                                            if(i!=1) temp=paste(temp,"<td style=\"border-top: hidden;\"",sep="")
                                            else temp=paste(temp,"<td",sep="")
                                            temp=paste(temp,"></td>\n")
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                cat(temp,"</tr>\n")
                rgroupcount=rgroupcount+1
            }
        }
        bcolor="white"
        #if(i %in% z$prefix.rows)
        #    if(is.numeric(z$zebra)) bcolor=z$zebra.color[i]
        #        cat("<tr style=\"background-color:",name2rgb(bcolor),"\">")
        cat("<tr>\n")
        if(z$include.rownames) {
            result=1
            if(!is.null(isspanCol(z,(i+1),1)))
                cat(paste("<td colspan=\"",isspanCol(z,i+1,1),"\"",sep=""))
            else if(!is.null(isspanRow(z,(i+1),1))){
                result=isspanRow(z,(i+1),1)
                if(result>0) cat(paste("<td rowspan=\"",result,"\"",sep=""))

            } else cat("<td ")
            if(result!=0){
                cat(paste(" style=\"border-left: ",vlines[1],"px solid black; ",sep=""))
                if(i==1 & printtop) cat("border-top: 2px solid gray;")
                else if(i!=1 | rgroupprinted) cat("border-top: hidden;")
                if(!is.null(z$hline.after)){
                    if((i-1) %in% z$hline.after)
                        if(!(i %in% printrgroup)) cat("border-top: 1px solid black;")
                }
                if(z$cellcolor[i+1,1]!="white")
                    cat(paste("background-color: ",name2rgb(z$cellcolor[i+1,1]),"; ",sep=""))
                cat(paste("\">",rownames(z$x)[i],"</td>\n",sep=""))
            }

        }
        for(j in 1:ncount) {
            if(is.null(isspanCol(z,(i+1),(j+1)))){
                if(is.null(isspanRow(z,(i+1),(j+1)))){
                    result=-1
                    cat("<td ")
                } else {
                    result=isspanRow(z,(i+1),(j+1))
                    if(result > 0) {
                        cat("<td rowspan=\"",result,"\" ")
                    }
                }
                if((result==-1)|(result>1)){
                    cat(paste("align=\"",colpos[j+1],"\" style=\"border-left: ",
                              vlines[j+1],"px solid black;",sep=""))
                    if((j==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                        cat(paste("border-right:",vlines[j+2],"px solid black;",sep=""))
                    if(i==1 & printtop) cat("border-top: 2px solid gray;")
                    else if(i!=1 | rgroupprinted) cat("border-top: hidden;")
                    if(!is.null(z$hline.after)){
                        if((i-1) %in% z$hline.after)
                            if(!(i %in% printrgroup)) cat("border-top: 1px solid black;")
                    }
                    if(z$cellcolor[i+1,j+1]!="white")
                        cat(paste("background-color: ",name2rgb(z$cellcolor[i+1,j+1]),";",sep=""))
                    cat("\">")
                    cat(paste(xdata[i,j],"</td>\n",sep=""))
                }
                if(j %in% colCount[-length(colCount)]) {
                    if(vlines[j+2]==0) {
                        backcolor=NULL
                        if(!is.null(z$rowcolor)){
                            if(z$rowcolor[i+1]!="white") backcolor=z$rowcolor[i+1]
                        }
                        if(is.null(backcolor)){
                            if((z$cellcolor[i+1,j+1]!="white")&(z$cellcolor[i+1,j+1]==z$cellcolor[i+1,j+2]))
                                backcolor=z$cellcolor[i+1,j+1]
                        }
                        cat("<td style=\"")
                        if(i==1 & printtop) cat("border-top: 2px solid gray;")
                        else if(i!=1 | rgroupprinted) cat("border-top: hidden;")

                        if(!is.null(backcolor)) cat(" background-color: ",name2rgb(backcolor),";")
                        cat("\"></td>\n")

                    }
                }
            } else {
                result=isspanCol(z,(i+1),(j+1))
                if(result>0) {
                    width=spanColWidth(z,(i+1),(j+1))
                    cat(paste("<td colspan=\"",result,"\" align=\"",colpos[j+1],"\" style=\"border-left: ",
                              vlines[j+1],"px solid black;",sep=""))
                    #if((j==ncol(z$x)) & (length(vlines)>ncol(z$x)+1))
                    cat(paste("border-right:",vlines[j+width+1],"px solid black;",sep=""))
                    if(i==1 & printtop) cat("border-top: 2px solid gray;")
                    else if(i!=1 | rgroupprinted) cat("border-top: hidden;")
                    if(!is.null(z$hline.after)){
                        if((i-1) %in% z$hline.after)
                            if(!(i %in% printrgroup)) cat("border-top: 1px solid black;")
                    }
                    if(z$cellcolor[i+1,j+1]!="white")
                        cat(paste("background-color: ",name2rgb(z$cellcolor[i+1,j+1]),";",sep=""))
                    cat("\">")
                    cat(paste(xdata[i,j],"</td>\n",sep=""))
                    if(isGroupCol(j,result,colCount)) {
                        if(vlines[j+width+1]==0) {

                            cat("<td style=\"")
                            if(i==1 & printtop) cat("border-top: 2px solid gray;")
                            else if(i!=1 | rgroupprinted) cat("border-top: hidden;")

                            if(!is.null(backcolor)) cat(" background-color: ",name2rgb(backcolor),";")
                            cat("\"></td>\n")
                        }
                    }
                }
            }

        }
        cat("</tr>\n")
    }
    if((z$show.footer!=TRUE) | (is.null(attr(z$x,"footer")))) footer=""
    else footer=attr(z$x,"footer")
    cat("<tr>\n")
    cat(paste("<td colspan=\"",totalCol,
              "\" align=\"left\" style=\"font-size:",as.integer(headingsize),
              "pt ;border-top: 1px solid black; border-bottom: hidden;\">",footer,"</td>\n",sep=""))
    cat("</tr>\n")
    cat("</table>\n")
}

#' Print an object of ztable via rstudio::viewer
#'
#' @param z An object of ztable
ztable2viewer=function(z){
    temp.f=tempfile(fileext=".html")
    sink(temp.f)
    cat(paste("<html>",
              "<head>",
              "<meta http-equiv=\"Content-type\" content=\"text/html;charset=UTF-8\">",
              "</head>",
              "<body>",
              "<div style=\"margin: 0 auto; display: table; margin-top: 1em;\">",
              sep="\n"))
    print(z,type="html")
    cat(paste("</div>","</body>","</html>",sep="\n"))
    sink()

    viewer <- getOption("viewer")
    if (!is.null(viewer) &&
            is.function(viewer)){
        # (code to write some content to the file)
        viewer(temp.f)
    }else{
        utils::browseURL(temp.f)
    }
}
