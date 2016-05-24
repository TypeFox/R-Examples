#' Exporting a R object to an object of class "ztable"
#'
#' Exporting a R object to an object of class "ztable"
#' @param x An R object, mainly data.frame
#'@param digits Numeric vector of length equal to one (in which case it will be
#'       replicated as necessary) or to the number of columns of the resulting table
#' @param ... arguments to be passed to \code{\link{ztable_sub}}
ztable=function(x,digits=NULL,...)  UseMethod("ztable")

#'@describeIn ztable
#'
ztable.default=function(x,digits=NULL,...){
  cat(paste("\n Sorry, Currently function ztable() cannot handle",
      " the object of class ",class(x),"!\n",sep=""))
  invisible()
}

#'@describeIn ztable
#'
ztable.data.frame=function(x,digits=NULL,...){
    z=ztable_sub(x,digits=digits,...)
    class(z) <-c("ztable")
    z
}

#' Exporting "data.frame" to an object of class "ztable"
#'
#' Exporting "data.frame" to an object of class "ztable"
#'@param x A data.frame
#'@param size An integer from 1 to 10 indicating font size= c("tiny","scriptsize",
#'       "footnotesize","small","normalsize","large","Large","LARGE","huge","Huge")
#'       respectively. Defaulting is 5(= "normalsize").
#'@param color A character indicating color of ztable
#'@param tablewidth A numeric value indicating desired table width as a ratio to linewidth.
#'       This value is only useful when caption is longer than table length.
#'       Default value is 0.3.
#'@param type character indicating formats of ztable, either "html" or "latex".
#'       Default value is "latex"
#'@param include.rownames A logical value whether or not include rownames in the table
#'       Default value is TRUE.
#'@param placement The table will have placement given by placement where placement
#'       must be NULL or contain only elements of {"h","t","b","p","!","H"}.
#'       Default value is "!hbtp".
#'@param position The table will be have placed at the center of the paper
#'        if position is "center" or "c", and at the left side of the paper
#'        if it equals "left" or "l", and at the right side of the paper
#'        if it equals "right" or "r". The position is translated to specificed
#'        latex environments such as "flushright" or "flushleft" or "center"
#'        (provided as a character vector) will enclose the tabular environment.
#'        Default value is "center".
#'@param show.heading A logical value whether or not include headings in the table.
#'        Default value is TRUE.
#'@param show.footer A logical value whether or not include headings in the table.
#'        Default value is TRUE.
#'@param caption A character
#'@param caption.placement The caption will be have placed at the top of the table
#'        if caption.placement is "top" and at the bottom of the table
#'        if it equals "bottom". Default value is "top".
#'@param caption.position The caption will be have placed at the center of the table
#'        if caption.position is "center" or "c", and at the left side of the table
#'        if it equals "left" or "l", and at the right side of the table
#'        if it equals "right" or "r". Default value is "center".
#'@param caption.bold whether or not use bold font for caption
#'@param align Character vector : nchar equal to the number of columns of the
#'       resulting table indicating the alignment of the corresponding columns.
#'@param digits Numeric vector of length equal to one (in which case it will be
#'       replicated as necessary) or to the number of columns of the resulting table
#'@param display Character vector of length equal to the number of columns of the
#'       resulting table indicating the format for the corresponding columns.
#'       Since the row names are printed in the first column, the length of display
#'       is one greater than ncol(x) if x is a data.frame. These values are passed
#'       to the formatC function. Use "d" (for integers), "f", "e", "E", "g", "G",
#'       "fg" (for reals), or "s" (for strings). "f" gives numbers in the usual
#'       xxx.xxx format; "e" and "E" give n.ddde+nn or n.dddE+nn (scientific format);
#'       "g" and "G" put x[i] into scientific format only if it saves space to do so.
#'       "fg" uses fixed format as "f", but digits as number of significant digits.
#'       Note that this can lead to quite long result strings. Default value is NULL.
#'       the class of x.
#'@param sidewaystable Logical value whether or not set the tabular environment=
#'       "sidewaystable". Requires Latex "rotating" package in preamble.
#'       Default value is FALSE.
#'@param longtable Logical value whether or not set the tabular environment=
#'       "longtable". Requires Latex "longtable" package in preamble.
#'       Default value is FALSE.
#'@param wraptable Logical value whether or not set the tabular environment=
#'       "wraptable". Requires Latex "wrapfig" package in preamble.
#'       Default value is FALSE.
#'@param rotate Logical value whether or not set the tabular environment=
#'       "rotate". No special arrangement is made to find space for the resut.
#'       Requires Latex "rotating" package in preamble.
#'       If TRUE, requires the rotate angle(counterclockwise).
#'       Default value is FALSE.
#'@param turn Logical value whether or not set the tabular environment=
#'       "turn". In this environment, Latex leaves space for the rotated table.
#'       Requires Latex "rotating" package in preamble.
#'       If TRUE, requires the rotate angle.
#'       Default value is FALSE.
#'@param angle An integer indicate the angle to rotate(degree); range -180 to 180.
#'       Default value is 0.
#'@param wraptablewidth A integer indicate wraptable width in centimeter. Default=12.
#'@param tabular Logical value whether or not set the tabular environment.
#'       If TRUE, no tabular environment is set. Default value is FALSE.
#'@param label Character vector of length 1 containing the LaTeX label or HTML anchor.
#'       Set to NULL to suppress the label. Default value is NULL.
#'@param hline.after A vector of numbers between -1 and "nrow(x)", inclusive,
#'       indicating the rows after which a horizontal line should appear.
#'       If NULL is used no lines are produced. Default value is c(-1,0,nrow(x))
#'       which means draw a line before and after the columns names and at the
#'       end of the table. Repeated values are allowed.
#'@param booktabs Logical value. If TRUE, the toprule, midrule and bottomrule tags
#'       from the LaTex "booktabs" package are used rather than hline for the
#'       horizontal line tags. Requires Latex "booktabs" package in preamble.
#'       Default value is TRUE.
#'@param prefix.rows A numeric vector contains the position of rows on which
#'       extra Latex commands should be added as a prefix.
#'@param commands A character vector of the length 1 or same length of the nrow of
#'       data.frame which contains the command that should be added as a prefix at
#'       the specified rows. Default value is NULL, i.e. do not add commands.
#'@param top.command A character vector of the length 1 which contains the command
#'       that should be added as a prefix at the colnames row.
#'@param zebra Null or an integer of 0 or 1 or 2 or 3. The arguments zebra and zebra.color are
#'       used to make a Zebra striping table(table with alternating background colors)
#'       easly. A value of 1 sets background color of all odd rows/columns with specified with
#'       zebra.color. A value of 2 sets all even rows/columns. A value of 0 sets
#'       background colors of all rows/columns with colors specified with zebra.color.
#'       When zebra is 1 or 2, the parameters of prefix.rows and commands ignored.
#'       When zebra=3, the background colors can be defined by addRowColor, addColColor
#'       and addCellColor functions.
#'       Default is NULL.
#'@param zebra.color A color name or a numeric value indicating pre-defined color.
#'       When parameter zebra is 0 or 1 or 2 and zebra.color is NULL, then zerba.color
#'       is set to "platinum". Numeric values between 1 to 13 is converted to
#'       predefined color names. The predefined color names are c("peach","peach-orange",
#'       "peachpuff","peach-yellow","pear","pearl","peridot","periwinkle","pastelred",
#'       "pastelgray"). Default is NULL.
#'@param zebra.colnames whether or not use background colors in column names row,
#'       Default value is FALSE
#'@param zebra.rownames whether or not use background colors in row names column,
#'       Default value is TRUE
#'@param zebra.type An integer of 0 or 1 or 2 or 3. A value of 1 sets background colors by row.
#'       A value of 2 sets background colors by column. A value of 0 sets background colors of all cells.
#'       A value of 3 sets background colors of cells specified with zebra.list.
#'       Default value is 1.
#'@param zebra.list A list consists of y,x,color. zebra.list is used only when zebra.type=3.
#'       zebra.list sets the cells specified with rows of vector "y" and columns of vector "x" with "color".
#'       The y and x are integer vactor indicating rows and columns. NA value of y or x indicating all columns or rows.
#'       The color is an character vector consists of names of color.
#'@param colnames.bold whether or not use bold font for column names, Default value is FALSE
#'@param include.colnames Logical. If TRUE the column names is printed. Default value is TRUE.
#'@param cgroup A character vector or matrix indicating names of column group. Default value is NULL
#'@param n.cgroup A integer vector or matrix indicating the numbers of columns included in each cgroup
#'       Dafault value is NULL
#'@param rgroup A character vector indicating names of row group. Default value is NULL
#'@param n.rgroup A integer vector indicating the numbers of rows included in each rgroup
#'       Dafault value is NULL
#'@param cspan.rgroup The number of columns that an rgroup should span. It spans by default all
#'       columns but you may want to limit this if you have column colors that you want to retain.
#'@examples
#' require(ztable)
#' x=head(iris)
#' ztable(x)
#' ztable(x,size=3,caption="Table 1. ztable Test")
#' ztable(x,size=7,caption="Table 1. ztable Test",caption.position="l")
#' ztable(x,size=7,caption="Table 1. ztable Test",caption.placement="bottom",
#'       caption.position="l")
#' fit=lm(mpg~.,data=mtcars)
#' ztable(fit)
#' data(USArrests)
#' pr1 <- prcomp(USArrests)
#' ztable(pr1)
#' ztable(summary(pr1))
#' require(survival)
#' data(colon)
#' attach(colon)
#' out <- glm(status ~ rx+obstruct+adhere+nodes+extent, data=colon, family=binomial)
#' ztable(out)
#' colon$TS = Surv(time,status==1)
#' out1=coxph(TS~rx+obstruct+adhere+differ+extent+surg+node4,data=colon)
#' ztable(out1)
#' ztable(head(mtcars),zebra=1)
#' ztable(head(mtcars),zebra=1,zebra.type=2)
ztable_sub=function(x,
                    size=5, # normal size, range 1-10
                    color=getOption("ztable.color","black"),
                    tablewidth=0.3,
                    type=getOption("ztable.type","latex"),
                    include.rownames=getOption("ztable.include.rownames",TRUE),
                    placement="!hbtp",position="c",
                    show.heading=getOption("ztable.show.heading",TRUE),
                    show.footer=getOption("ztable.show.footer",TRUE),
                    caption=NULL,
                    caption.placement=getOption("ztable.caption.placement","top"),
                    caption.position=getOption("ztable.caption.position","c"),
                    caption.bold=getOption("ztable.caption.bold",FALSE),
                    align=NULL,digits=NULL,display=NULL,
                    sidewaystable=FALSE,
                    longtable=FALSE,
                    rotate=FALSE,
                    turn=FALSE,
                    angle=0,
                    wraptable=FALSE,wraptablewidth=12,
                    tabular=FALSE,
                    label=NULL,hline.after=NULL,
                    booktabs=getOption("ztable.booktabs",TRUE),
                    prefix.rows=NULL,commands=NULL,top.command=NULL,
                    zebra=getOption("ztable.zebra",NULL),
                    zebra.color=getOption("ztable.zebra.color",NULL),
                    zebra.type=getOption("ztable.zebra.type",1),
                    zebra.colnames=getOption("ztable.zebra.colnames",FALSE),
                    zebra.rownames=getOption("ztable.zebra.rownames",TRUE),
                    zebra.list=NULL,
                    colnames.bold=getOption("ztable.colnames.bold",FALSE),
                    include.colnames=getOption("ztable.include.colnames",TRUE),
                    cgroup=NULL,n.cgroup=NULL,
                    rgroup=NULL,n.rgroup=NULL,cspan.rgroup=NULL){

    ncount=ncol(x)
    nrow=nrow(x)
    cn=colnames(x)
    if(identical(caption.placement,"bottom") | identical(caption.placement,"b"))
        caption.placement="bottom"
    else caption.placement="top"
    if(identical(caption.position,"left")|identical(caption.position,"l"))
        caption.position="l"
    else if(identical(caption.position,"right")|identical(caption.position,"r"))
        caption.position="r"
    else caption.position="c"

    if(identical(position,"left")|identical(position,"l"))
        position="flushleft"
    else if(identical(position,"right")|identical(position,"r"))
        position="flushright"
    else position="center"

    addrow=ifelse(include.rownames,1,0)
    logicals <- sapply(x, is.logical)

    x[logicals] <- lapply(x[logicals], as.character)

    characters <- sapply(x, is.character)
    factors <- sapply(x, is.factor)
    ints <- sapply(x, is.integer)

    if(is.null(align)){
        y <- c("r", c("r","l")[(characters | factors) + 1])
        if(include.rownames) for(i in 1:length(y)) align=paste(align,y[i],sep="")
        else for(i in 2:length(y)) align=paste(align,y[i],sep="")
    }
    if(is.null(digits)) digits=c(0,rep(2,ncol(x)))
    if(length(digits)==1) digits=rep(digits,ncount+1)
    if (is.null(display)) {
        display <- rep("f", ncol(x))
        display[ints] <- "d"
        display[characters | factors] <- "s"
        display <- c("s", display)
    }
    if(!is.null(zebra)) {
        if(zebra==0){
            prefix.rows=1:nrow(x)
            if(is.null(zebra.color)) zebra.color=1:10
        } else if(zebra==1) {
            prefix.rows=seq(2,nrow(x),by=2)
            if(is.null(zebra.color)) zebra.color=2 #peach-orange
        } else {
            zebra=2
            prefix.rows=seq(1,nrow(x),by=2)
            if(is.null(zebra.color)) zebra.color=2 #peach-orange
        }
        mycolor=c("peach","peach-orange","peachpuff","peach-yellow","pear",
                  "pearl","peridot","periwinkle","pastelred","pastelgray")
        if(zebra!=0) {
            zebra.color[1]=validColor(zebra.color[1],mycolor)
            zebra.color=rep(zebra.color[1],nrow)
        }
        else {     # zebra==0; all rows
            result=c()
            for(i in 1:length(zebra.color)){
                result=c(result,validColor(zebra.color[i],mycolor))
            }
            zebra.color=result
            if(length(zebra.color)<nrow)
                zebra.color=rep(zebra.color,1+(nrow/length(zebra.color)))
        }
    }
    cellcolor=make.cell.color(x,zebra,zebra.color,zebra.type,zebra.list,
                              zebra.colnames,zebra.rownames)
    if(!is.null(prefix.rows) & (length(commands)==1))
        commands=rep(commands,nrow)
    if((0 %in% prefix.rows) & is.null(top.command) &(length(commands)>0))
        top.command=commands[1]
    if(!is.numeric(size)) size=5
    else if(size<0 | size>10) size=5

    result=list(x=x,
                cellcolor=cellcolor,
                size=size,
                color=color,
                tablewidth=tablewidth,
                type=type,
                include.rownames=include.rownames,
                placement=placement,
                position=position,
                show.heading=show.heading,
                show.footer=show.footer,
                caption=caption,
                caption.placement=caption.placement,
                caption.position=caption.position,
                caption.bold=caption.bold,
                align=align,
                digits=digits,
                display=display,
                sidewaystable=sidewaystable,
                longtable=longtable,
                wraptable=wraptable,
                wraptablewidth=wraptablewidth,
                tabular=tabular,
                rotate=rotate,
                turn=turn,
                angle=angle,
                label=label,
                hline.after=hline.after,
                booktabs=booktabs,
                prefix.rows=prefix.rows,
                commands=commands,
                top.command=top.command,
                zebra=zebra,
                zebra.color=zebra.color,
                zebra.type=zebra.type,
                zebra.list=zebra.list,
                zebra.colnames=zebra.colnames,
                zebra.rownames=zebra.rownames,
                include.colnames=include.colnames,
                colnames.bold=colnames.bold,
                cgroup=cgroup,
                n.cgroup=n.cgroup,
                rgroup=rgroup,
                n.rgroup=n.rgroup,
                cspan.rgroup=cspan.rgroup
    )
    class(result) <-c("ztable")
    result
}

#' Make a data.frame named "cellcolor" from ztable call
#'
#'@param x a data.frame
#'@param zebra Null or an integer of 0 or 1 or 2. The arguments zebra and zebra.color are
#'       used to make a Zebra striping table(table with alternating background colors)
#'       easly. A value of 1 sets background color of all odd rows/columns with specified with
#'       zebra.color. A value of 2 sets all even rows/columns. A value of 0 sets
#'       background colors of all rows/columns with colors specified with zebra.color.
#'       When zebra is 1 or 2, the parameters of prefix.rows and commands ignored.
#'       Default is NULL.
#'@param zebra.color A color name or a numeric value indicating pre-defined color.
#'       When parameter zebra is 0 or 1 or 2 and zebra.color is NULL, then zerba.color
#'       is set to "platinum". Numeric values between 1 to 13 is converted to
#'       predefined color names. The predefined color names are c("peach","peach-orange",
#'       "peachpuff","peach-yellow","pear","pearl","peridot","periwinkle","pastelred",
#'       "pastelgray"). Default is NULL.
#'@param zebra.type An integer of 0 or 1 or 2 or 3. A value of 1 sets background colors by row.
#'       A value of 2 sets background colors by column. A value of 0 sets background colors of all cells.
#'       A value of 3 sets background colors of cells specified with zebra.list.
#'       Default value is 1.
#'@param zebra.list A list consists of y,x,color. zebra.list is used only when zebra.type=3.
#'       zebra.list sets the cells specified with rows of vector "y" and columns of vector "x" with "color".
#'       The y and x are integer vactor indicating rows and columns. NA value of y or x indicating all columns or rows.
#'       The color is an character vector consists of names of color.
#'@param zebra.colnames whether or not use background colors in column names row,
#'       Default value is FALSE
#'@param zebra.rownames whether or not use background colors in row names column,
#'       Default value is TRUE
make.cell.color=function(x,zebra,zebra.color,zebra.type,zebra.list,
                         zebra.colnames,zebra.rownames){

    temp=rep("white",nrow(x)+1)
    cellcolor=c()
    for(i in 1:(ncol(x)+1)) cellcolor=cbind(cellcolor,temp)
    colnames(cellcolor)=c(" ",colnames(x))
    rownames(cellcolor)=c(" ",rownames(x))
    if(!is.null(zebra)){

        if(is.null(zebra.color)) {
            if(zebra==0) color=1:10
            else color=2
        }
        else color=zebra.color

        mycolor=c("peach","peach-orange","peachpuff","peach-yellow","pear",
                  "pearl","peridot","periwinkle","pastelred","pastelgray")
        result=c()
        for(i in 1:length(color)){
                result=c(result,validColor(color[i],mycolor))
        }
        color=result
        if(zebra==0){
            if(zebra.type==1){
               cellcolor=apply(cellcolor,2,repColor,color=color)
            }
            else if(zebra.type==2){
                cellcolor=apply(cellcolor,1,repColor,color=color)
                cellcolor=t(cellcolor)
            }
            else if(zebra.type==0){
                temp=rep(color,1+length(cellcolor)/length(color))
                for(i in 1:nrow(cellcolor))
                    for(j in 1:ncol(cellcolor))
                        cellcolor[i,j]=temp[(i-1)*ncol(cellcolor)+j]
            }
        }
        else if(zebra==1){
            if(zebra.type==1){
                select=seq(2,nrow(cellcolor),by=2)
                for(i in select)
                    for(j in 1:ncol(cellcolor)) cellcolor[i,j]=color[1]
            }
            else if(zebra.type==2){
                select=seq(1,ncol(cellcolor),by=2)
                for(i in 1:nrow(cellcolor))
                    for(j in select)
                    cellcolor[i,j]=color[1]
            }
        }
        else if(zebra==2){
            if(zebra.type==1){
                select=seq(1,nrow(cellcolor),by=2)
                for(i in select)
                    for(j in 1:ncol(cellcolor)) cellcolor[i,j]=color[1]
            }
            else if(zebra.type==2){
                select=seq(2,ncol(cellcolor),by=2)
                for(i in 1:nrow(cellcolor))
                    for(j in select)
                        cellcolor[i,j]=color[1]
            }
        }
    }
    if(zebra.colnames==FALSE) {
        cellcolor[1,1:ncol(cellcolor)]="white"
    }
    if(zebra.rownames==FALSE) {
        cellcolor[1:nrow(cellcolor),1]="white"
    }
    if(zebra.type==3){
        if(!is.null(zebra.list)){
           ylength=length(zebra.list$y)
           if(length(zebra.list$color)<ylength){
                zebra.list$color=rep(zebra.list$color,1+ylength/length(zebra.list$color))
            }
            while(length(zebra.list$x)<ylength){
                zebra.list$x=c(zebra.list$x,NA)
            }
            for(i in 1:ylength){
                if(is.na(zebra.list$y[i])) {
                    if(is.na(zebra.list$x[i])) next
                    for(j in 1:nrow(cellcolor)) cellcolor[j,zebra.list$x[i]]=zebra.list$color[i]
                }
                else if(is.na(zebra.list$x[i])){
                    for(j in 1:ncol(cellcolor)) cellcolor[zebra.list$y[i],j]=zebra.list$color[i]
                }
                else cellcolor[zebra.list$y[i],zebra.list$x[i]]=zebra.list$color[i]
            }
        }

    }
    cellcolor
}

#' Make vector x from vector color
#'
#' Internal function of make.cell.color
#' @param x A destination vector
#' @param color A charactor vector consists of color names
repColor=function(x,color){
    #cat("color=",color,"\n")
    temp=rep(color,1+length(x)/length(color))
    for(i in 1:length(x)) x[i]=temp[i]
    x
}

#'Update ztable before print
#'
#'Update options of ztable before print
#'@param z An object of class "ztable"
#'@param size An integer from 1 to 10 indicating font size= c("tiny","scriptsize",
#'       "footnotesize","small","normalsize","large","Large","LARGE","huge","Huge")
#'       respectively.
#'@param color A character indicating color of ztable
#'@param tablewidth A numeric indicating desired table width as a ratio to linewidth.
#'       Default value is 0.3.
#'@param type character indicating formats of ztable, either "html" or "latex".
#'@param include.rownames A logical value whether or not include rownames in the table
#'@param placement The table will have placement given by placement where placement
#'       must be NULL or contain only elements of {"h","t","b","p","!","H"}.
#'@param position The table will be have placed at the center of the paper
#'        if position is "center" or "c", and at the left side of the paper
#'        if it equals "left" or "l", and at the right side of the paper
#'        if it equals "right" or "r". The position is translated to specificed
#'        latex environments such as "flushright" or "flushleft" or "center"
#'        (provided as a character vector) will enclose the tabular environment.
#'@param show.heading A logical value whether or not include headings in the table.
#'@param show.footer A logical value whether or not include headings in the table.
#'@param caption A character
#'@param caption.placement The caption will be have placed at the top of the table
#'        if caption.placement is "top" and at the bottom of the table
#'        if it equals "bottom".
#'@param caption.position The caption will be have placed at the center of the table
#'        if caption.position is "center" or "c", and at the left side of the table
#'        if it equals "left" or "l", and at the right side of the table
#'        if it equals "right" or "r".
#'@param caption.bold whether or not use bold font for caption
#'@param align Character vector : nchar equal to the number of columns of the
#'       resulting table indicating the alignment of the corresponding columns.
#'@param digits Numeric vector of length equal to one (in which case it will be
#'       replicated as necessary) or to the number of columns of the resulting table
#'@param display Character vector of length equal to the number of columns of the
#'       resulting table indicating the format for the corresponding columns.
#'       Since the row names are printed in the first column, the length of display
#'       is one greater than ncol(x) if x is a data.frame. These values are passed
#'       to the formatC function. Use "d" (for integers), "f", "e", "E", "g", "G",
#'       "fg" (for reals), or "s" (for strings). "f" gives numbers in the usual
#'       xxx.xxx format; "e" and "E" give n.ddde+nn or n.dddE+nn (scientific format);
#'       "g" and "G" put x[i] into scientific format only if it saves space to do so.
#'       "fg" uses fixed format as "f", but digits as number of significant digits.
#'       Note that this can lead to quite long result strings.
#'@param sidewaystable Logical value whether or not set the tabular environment=
#'       "sidewaystable". Requires Latex "rotating" package in preamble.
#'@param longtable Logical value whether or not set the tabular environment=
#'       "longtable". Requires Latex "longtable" package in preamble.
#'@param wraptable Logical value whether or not set the tabular environment=
#'       "wraptable". Requires Latex "wrapfig" package in preamble.
#'@param rotate Logical value whether or not set the tabular environment=
#'       "rotate". No special arrangement is made to find space for the resut.
#'       Requires Latex "rotating" package in preamble.
#'       If TRUE, requires the rotate angle(counterclockwise).
#'@param turn Logical value whether or not set the tabular environment=
#'       "turn". In this environment, Latex leaves space for the rotated table.
#'       Requires Latex "rotating" package in preamble.
#'       If TRUE, requires the rotate angle.
#'@param angle An integer indicate the angle to rotate(degree); range -180 to 180.
#'@param wraptablewidth A integer indicate wraptable width in centimeter.
#'@param tabular Logical value whether or not set the tabular environment.
#'       If TRUE, no tabular environment is set.
#'@param label Character vector of length 1 containing the LaTeX label or HTML anchor.
#'       Set to NULL to suppress the label.
#'@param hline.after A vector of numbers between -1 and "nrow(x)", inclusive,
#'       indicating the rows after which a horizontal line should appear.
#'       If NULL is used no lines are produced. Default value is c(-1,0,nrow(x))
#'       which means draw a line before and after the columns names and at the
#'       end of the table. Repeated values are allowed.
#'@param booktabs Logical value. If TRUE, the toprule, midrule and bottomrule tags
#'       from the LaTex "booktabs" package are used rather than hline for the
#'       horizontal line tags. Requires Latex "booktabs" package in preamble.
#'@param prefix.rows A numeric vector contains the position of rows on which
#'       extra Latex commands should be added as a prefix.
#'@param commands A character vector of the length 1 or same length of the nrow of
#'       data.frame which contains the command that should be added as a prefix at
#'       the specified rows.
#'@param top.command A character vector of the length 1 which contains the command
#'       that should be added as a prefix at the colnames row.
#'@param zebra Null or a integer of 1 or 2. The arguments zebra and zebra.color are
#'       used to make a Zebra striping table(table with alternating background colors)
#'       easly. A value of 1 sets background color of all odd rows with specified with
#'       zebra.color. A value of 2 sets all even rows. when zebra is 1 or 2,
#'       the parameters of prefix.rows and commands ignored.
#'@param zebra.color A color name or a numeric value indicating pre-defined color.
#'       When parameter zebra is 0 or 1 or 2 and zebra.color is NULL, then zerba.color
#'       is set to "platinum". Numeric values between 1 to 13 is converted to
#'       predefined color names. The predefined color names are c("peach","peach-orange",
#'       "peachpuff","peach-yellow","pear","pearl","peridot","periwinkle","pastelred",
#'       "pastelgray").
#'@param zebra.type An integer of 0 or 1 or 2 or 3. A value of 1 sets background colors by row.
#'       A value of 2 sets background colors by column. A value of 0 sets background colors of all cells.
#'       A value of 3 sets background colors of cells specified with zebra.list.
#'       Default value is 1.
#'@param zebra.list A list consists of y,x,color. zebra.list is used only when zebra.type=3.
#'       zebra.list sets the cells specified with cells[y,x] with "color". The y and x are
#'       integer indicating rows and columns. NA value of y or x indicating all columns or rows.
#'@param zebra.colnames whether or not use background colors in column names row,
#'       Default value is FALSE
#'@param zebra.rownames whether or not use background colors in row names column,
#'       Default value is TRUE
#'@param colnames.bold whether or not use bold font for column names.
#'@param include.colnames Logical. If TRUE the column names is printed.
#'@param cgroup A character vector or matrix indicating names of column group. Default value is NULL
#'@param n.cgroup A integer vector or matrix indicating the numbers of columns included in each cgroup
#'       Dafault value is NULL
#'@param rgroup A character vector indicating names of row group. Default value is NULL
#'@param n.rgroup A integer vector indicating the numbers of rows included in each rgroup
#'       Dafault value is NULL
#'@param cspan.rgroup The number of columns that an rgroup should span. It spans by default all
#'       columns but you may want to limit this if you have column colors that you want to retain.
update_ztable=function(z,
                       size=NULL, # normal size, range 1-10
                       color=NULL,
                       tablewidth=NULL,
                       type=NULL,
                       include.rownames=NULL,
                       placement=NULL,position=NULL,
                       show.heading=NULL,
                       show.footer=NULL,
                       caption=NULL,
                       caption.placement=NULL,
                       caption.position=NULL,
                       caption.bold=NULL,
                       align=NULL,digits=NULL,display=NULL,
                       sidewaystable=NULL,
                       longtable=NULL,
                       rotate=NULL,
                       turn=NULL,
                       angle=NULL,
                       wraptable=NULL,wraptablewidth=NULL,
                       tabular=NULL,
                       label=NULL,hline.after=NULL,
                       booktabs=NULL,
                       prefix.rows=NULL,commands=NULL,top.command=NULL,
                       zebra=NULL,
                       zebra.color=NULL,
                       zebra.type=NULL,
                       zebra.list=NULL,
                       zebra.colnames=NULL,
                       zebra.rownames=NULL,
                       colnames.bold=NULL,
                       include.colnames=NULL,
                       cgroup=NULL,
                       n.cgroup=NULL,
                       rgroup=NULL,
                       n.rgroup=NULL,
                       cspan.rgroup=NULL){

     if(!is.null(size)) z$size=size
     if(!is.null(color)) z$color=color
     if(!is.null(tablewidth)) z$tablewidth=tablewidth
     if(!is.null(type)) z$type=type
     if(!is.null(include.rownames)) z$include.rownames=include.rownames
     if(!is.null(placement)) z$placement=placement
     if(!is.null(position)) z$position=position
     if(!is.null(show.heading)) z$show.heading=show.heading
     if(!is.null(show.footer)) z$show.footer=show.footer
     if(!is.null(caption)) z$caption=caption
     if(!is.null(caption.placement)) z$caption.placement=caption.placement
     if(!is.null(caption.position)) z$caption.position=caption.position
     if(!is.null(caption.bold)) z$caption.bold=caption.bold
     if(!is.null(align)) z$align=align
     if(!is.null(digits)) {
         if(is.null(digits)) digits=c(0,rep(2,ncol(z$x)))
         if(length(digits)==1) digits=rep(digits,ncol(z$x)+1)
         z$digits=digits
     }
     if(!is.null(display)) z$display=display
     if(!is.null(sidewaystable)) z$sidewaystable=sidewaystable
     if(!is.null(longtable)) z$longtable=longtable
     if(!is.null(rotate)) z$rotate=rotate
     if(!is.null(turn)) z$turn=turn
     if(!is.null(angle)) z$angle=angle
     if(!is.null(wraptable)) z$wraptable=wraptable
     if(!is.null(wraptablewidth)) z$wraptablewidth=wraptablewidth
     if(!is.null(tabular)) z$tabular=tabular
     if(!is.null(label)) z$label=label
     if(!is.null(hline.after)) z$hline.after=hline.after
     if(!is.null(booktabs)) z$booktabs=booktabs
     if(!is.null(prefix.rows)) z$prefix.rows=prefix.rows
     if(!is.null(commands)) z$commands=commands
     if(!is.null(top.command)) z$top.command=top.command
     if(!is.null(zebra)) z$zebra=zebra
     if(!is.null(zebra.color)) z$zebra.color=zebra.color
     if(!is.null(zebra.type)) z$zebra.type=zebra.type
     if(!is.null(zebra.list)) z$zebra.list=zebra.list
     if(!is.null(zebra.colnames)) z$zebra.colnames=zebra.colnames
     if(!is.null(zebra.rownames)) z$zebra.rownames=zebra.rownames
     if(!is.null(colnames.bold)) z$colnames.bold=colnames.bold
     if(!is.null(include.colnames)) z$include.colnames=include.colnames
     if(!is.null(cgroup)) z$cgroup=cgroup
     if(!is.null(n.cgroup)) z$n.cgroup=n.cgroup
     if(!is.null(rgroup)) z$rgroup=rgroup
     if(!is.null(n.rgroup)) z$n.rgroup=n.rgroup
     if(!is.null(cspan.rgroup)) z$cspan.rgroup=cspan.rgroup
     if(!is.null(z$zebra)) { if(z$zebra!=3) z$cellcolor=make.cell.color(x=z$x,zebra=z$zebra,zebra.color=z$zebra.color,
                                 zebra.type=z$zebra.type,
                                 zebra.list=z$zebra.list,
                                 zebra.colnames=z$zebra.colnames,
                                 zebra.rownames=z$zebra.rownames)
     }
     z
}

#' Print an object of class "ztable"
#'
#' @param x An object of class "ztable"
#' @param ... further argument passed to other function
print.ztable=function(x,...){
    z=update_ztable(z=x,...)
    print_ztable(z)
}


#' Print an object of class "ztable"
#'
#' @param z An object of class "ztable"
print_ztable=function(z){
    xdata=data2table(z)
    if(z$type=="latex") ztable2latex(z,xdata)
    else if(z$type=="viewer") ztable2viewer(z)
    else ztable2html(z,xdata)
}

#'Subfunction used in ztable2latex
#'
#' @param string a character vector
tr=function(string) {
    string=gsub("%","\\%",string,fixed=TRUE)
    string=gsub(" -","\\hspace{0.5cm}",string,fixed=TRUE)
    string
}


#'Subfunction used in ztable2html
#'
#' @param string a character vector
tr2=function(string) {
    string=gsub(" -","&nbsp;&nbsp;&nbsp;",string,fixed=TRUE)
    string=gsub("     ","",string,fixed=TRUE)
    string
}


#' Convert data to formatted data for table
#'
#' @param z An object of class "ztable"
tableLength=function(z){
    xdata=data2table(z)
    a=apply(xdata,2,function(x) max(nchar(x)))
    if(z$include.colnames){
        b=nchar(colnames(xdata))
        l=c()
        for(i in 1:ncol(xdata)){
            l=c(l,max(a[i],b[i]))
        }
        length=sum(l)
    }
    else length=sum(a)
    result=length+ncol(xdata)-1
    if(z$include.rownames) result=result+max(nchar(rownames(xdata)))+1
    result
}

#' Convert data to formatted data for table
#'
#' @param z An object of class "ztable"
data2table=function(z){
    data<-z$x
    ncount=ncol(data)
    nrow=nrow(data)

    select=sapply(data,is.factor)
    data[select]=lapply(data[select],as.character)
    #data

    for(i in 1:nrow){
        for(j in 1:ncount) {
            if(z$display[j+1]=="s"){
                temp=data[i,j]
                if(z$type=="latex") temp<-tr(temp)
                if(z$type=="html") temp<-tr2(temp)
            }
            else{
                if(is.na(z$x[i,j])) {
                    temp<-""
                } else{
                    temp<-formatC(z$x[i,j],digits=z$digits[j+1],
                                 format=z$display[j+1])
                }

            }
            data[i,j]<-temp
        }
    }
    data
}

#' Convert long caption to minipage
#'
#' @param z An object of ztable
#' @param caption A character vector to convert
caption2minipage=function(z,caption){
    tlength=tableLength(z)
    if(nchar(caption)>tlength){
        tablewidth=max(z$tablewidth,tlength/85)
        mycaption=paste("\\begin{minipage}[c]{",tablewidth,"\\linewidth}",
                    caption,"\\end{minipage}",sep="")
    }
    else mycaption=caption
    mycaption
}

#' Print an object of class "ztable" to Latex table
#'
#' @param z An object of class "ztable"
#' @param xdata A formatted data.frame
ztable2latex=function(z,xdata){
    ncount=ncol(z$x)
    nrow=nrow(z$x)
    cn=colnames(z$x)
    addrow=ifelse(z$include.rownames,1,0)

    NewAlign=getNewAlign(z)
    #NewAlign=z$align
    totalCol=totalCol(z)
    colCount=colGroupCount(z)

    vlines=align2lines(z$align)

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
    Fontsize=c("tiny","scriptsize","footnotesize","small","normalsize",
           "large","Large","LARGE","huge","Huge")


    if(z$tabular) sort="tabular"
    else if(z$sidewaystable) sort="sidewaystable"
    else if(z$wraptable) sort="wraptable"
    else if(z$rotate) sort="rotate"
    else if(z$turn) sort="turn"
    else sort="table"
    headingsize=ifelse(z$size>3,z$size-2,1)
    define_colors(z$cellcolor)
    if(!is.null(z$cgroupcolor)) define_colors(z$cgroupcolor)
    align=alignCheck(z$align,ncount,addrow)
    if(z$longtable){
        cat(paste("\\color{",z$color,"}\n",sep=""))
        cat(paste("\\begin{",Fontsize[z$size],"}\n",sep=""))
        cat(paste("\\begin{longtable}{",NewAlign,"}\n",sep=""))

    } else {
        if(z$wraptable) {
            if(z$position=="flushright") wrapposition<-"r"
            else wrapposition<-"l"
            cat(paste("\\begin{wraptable}{",wrapposition,"}[10pt]{",
                      z$wraptablewidth,"cm}\n",sep=""))

        } else if((sort=="rotate") | (sort=="turn")){
            cat(paste("\\begin{",sort,"}{",z$angle,"}\n",sep=""))
        } else if(sort!="tabular"){      # sidewaystable or table
            cat(paste("\\begin{",sort,"}[",z$placement,"]\n",sep=""))
            cat(paste("\\begin{",z$position,"}\n",sep=""))
        }
        cat(paste("\\begin{",Fontsize[z$size],"}\n",sep=""))
        cat(paste("\\color{",z$color,"}\n",sep=""))
        cat(paste("\\begin{tabular}{",NewAlign,"}\n",sep=""))
    }
    if(!is.null(z$caption) & z$caption.placement=="top"){
        mycaption=caption2minipage(z,z$caption)
        cat(paste("\\multicolumn{",totalCol,"}{",
                  z$caption.position,"}{",sep=""))
        if(z$caption.bold) cat(paste("\\textbf{",mycaption,"}",sep=""))
        else cat(mycaption)
        cat("}\\\\ \n")
    }
    if((z$show.heading==TRUE) & (!is.null(attr(z$x,"heading")))) {
        head=attr(z$x,"heading")
        for(i in 1:length(head)) {
            h1=gsub("~","$\\sim$",head[i],fixed=TRUE)
            if(nchar(head[i])<1) next
            cat(paste("\\multicolumn{",totalCol,"}{l}{\\",Fontsize[headingsize],
                      "{",h1,"}}\\\\ \n",sep=""))
        }
    }
    if(is.null(z$hline.after)) cat(ifelse(z$booktabs,"\\toprule[1.2pt]\n","\\hline\n"))
    else if(-1 %in% z$hline.after) cat(ifelse(z$booktabs,"\\toprule[1.2pt]\n","\\hline\n"))
    if(!is.null(z$cgroup)) printLatexHead(z)
    subcolnames=ifelse(is.null(z$subcolnames),0,1)

        if(subcolnames) {
            if(is.na(z$subcolnames[1])) firstcn=paste("\\multirow{2}{*}{}",sep="")
            else firstcn=cn[1]
        }
        else firstcn=cn[1]
        if(z$colnames.bold) firstcn=paste("\\textbf{",firstcn,"}",sep="")

        if(z$cellcolor[1,2]!="white") firstcn=paste("\\cellcolor{",z$cellcolor[1,2],"}",firstcn,sep="")
        if(z$include.rownames) {

            result=1
            if(!is.null(isspanCol(z,1,1)))
                first=paste("\\multicolumn{",isspanCol(z,1,1),"}{c}{}",sep="")
            else if(!is.null(isspanRow(z,1,1))){
                 result=isspanRow(z,1,1)
                 if(result>0) first=paste("\\multirow{",result,"}{*}{}",sep="")
            } else first=""
            if(z$cellcolor[1,1]!="white")
                first=paste("\\cellcolor{",z$cellcolor[1,1],first,"}",sep="")
            firstrow=paste(first,"&",firstcn,sep="")

        }
        else firstrow=firstcn

        if(ncount>1) {
            for(i in 2:ncount) {
                firstrow=paste(firstrow,"&",sep="")
                if((i==2)&(!is.null(colCount))){
                    if(1 %in% colCount[-length(colCount)]) {
                        if(vlines[1+2]==0) firstrow=paste(firstrow,"&",sep="")
                    }
                }
                if(z$cellcolor[1,i+1]!="white")
                    firstrow=paste(firstrow,"\\cellcolor{",z$cellcolor[1,i+1],"}",sep="")
                if(z$colnames.bold) boldcn=paste("\\textbf{",cn[i],"}",sep="")
                else boldcn=cn[i]
                result=1
                if(!is.null(isspanCol(z,1,(i+1)))){
                    result=isspanCol(z,1,(i+1))
                    if(result>0) boldcn=paste("\\multicolumn{",result,"}{c}{",boldcn,"}",sep="")
                    else if(result==0) next
                } else if(!is.null(isspanRow(z,1,(i+1)))){
                    boldcn=paste("\\multirow{",isspanRow(z,1,(i+1)),"}{*}{",boldcn,"}",sep="")
                }
                if((subcolnames==1)) {
                    if(is.na(z$subcolnames[i])){
                       # boldcn=paste("\\multirow{2}{*}{",boldcn,"}",sep="")
                        boldcn=""
                    }
                }
                firstrow=paste(firstrow,boldcn,sep="")
                if(!is.null(colCount)){
                    if(i %in% colCount[-length(colCount)]) {
                        if(vlines[i+2]==0) {
                            #if(z$cellcolor[1,i+1]!="white")
                            #    firstrow=paste(firstrow,"&\\cellcolor{",z$cellcolor[1,i+1],"}",sep="")
                            #else firstrow=paste(firstrow,"&",sep="")
                            firstrow=paste(firstrow,"&",sep="")
                        }
                    }
                }
            }
        }

    if((0 %in% z$prefix.rows) & !is.null(z$top.command)) cat(z$top.command)
    if(z$include.colnames) {
        cat(paste(firstrow,"\\\\ \n",sep=""))
        if(subcolnames){
            if(z$include.rownames) {
                if(z$cellcolor[1,1]!="white")
                    cat(paste("\\cellcolor{",z$cellcolor[1,1],"} &",sep=""))
                else cat("&")

            }
            for(i in 1:length(z$subcolnames)){
                if(is.na(z$subcolnames[i])) {
                    temp=paste("\\multirow{-2}{*}{",colnames(z$x)[i],"}",sep="")
                    if(!is.null(z$colcolor)){
                        if(z$cellcolor[1,i+1]!="white")
                            temp=paste("\\cellcolor{",z$cellcolor[1,i+1],"}",temp,sep="")
                    }
                    cat(temp)
                    if(i!=length(z$subcolnames)) cat("&")
                    if(i %in% colCount[-length(colCount)]) {
                        if(vlines[i+2]==0){
                            if((z$cellcolor[1,i+1]!="white") & (z$cellcolor[1,i+1]==z$cellcolor[1,i+2]))
                                cat(paste("\\cellcolor{",z$cellcolor[1,i+1],"}&",sep=""))
                            else cat("&")
                        }
                    }
                    next
                }
                if(z$colnames.bold) boldcn=paste("\\textbf{",z$subcolnames[i],"}",sep="")
                else boldcn=z$subcolnames[i]

                if(z$cellcolor[1,i+1]!="white")
                    cat(paste("\\cellcolor{",z$cellcolor[1,i+1],"}",boldcn,"&",sep=""))
                else cat(paste(boldcn,"&",sep=""))
                if(i %in% colCount[-length(colCount)]) {
                    if(vlines[i+2]==0){
                    if((z$cellcolor[1,i+1]!="white") & (z$cellcolor[1,i+1]==z$cellcolor[1,i+2]))
                        cat(paste("\\cellcolor{",z$cellcolor[1,i+1],"}&",sep=""))
                    else cat("&")
                    }
                }
            }
            cat("\\\\ \n")
        }
        if(is.null(z$hline.after)) cat(ifelse(z$booktabs,"\\midrule\n","\\hline\n"))
        else if(0 %in% z$hline.after) cat(ifelse(z$booktabs,"\\midrule\n","\\hline\n"))
    }

    for(i in 1:nrow){
        printcline=0
        if(rgroupcount>0) {
            if(i %in% printrgroup) {
                for(k in 1:length(printrgroup)){
                    if(i == printrgroup[k]){
                        if(is.na(z$rgroup[k])) break
                        if(z$rgroup[k]=="") break
                        printRowGroup(z,i)
                        break
                    }
                }

            }
        }
        if(i %in% z$prefix.rows) {
            #if(is.numeric(z$zebra))
            #   cat(paste("\\rowcolor{",z$zebra.color[i],"}",sep=""))
            if(!is.null(z$commands[i])) cat(z$commands[i])
        }
        tempo=NULL
        if(z$include.rownames) {
            if(z$cellcolor[i+1,1]=="white") tempo=rownames(z$x)[i]
            else tempo=paste("\\cellcolor{",z$cellcolor[i+1,1],"}",
                            rownames(z$x)[i],sep="")

            if(!is.null(isspanCol(z,(i+1),1)))
                tempo=paste("\\multicolumn{",isspanCol(z,i+1,1),"}{c}{",tempo,"}",sep="")
            else if(!is.null(isspanRow(z,(i+1),1))){
                result=isspanRow(z,(i+1),1)
                if(result<0) tempo=paste("\\multirow{",result,"}{*}{",tempo,"}",sep="")
            }
            cat(tempo)
        }

        for(j in 1:ncount) {
            skip=0
            if(z$cellcolor[i+1,j+1]=="white") temp1=xdata[i,j]
            else temp1=paste("\\cellcolor{",z$cellcolor[i+1,j+1],"}",
                             xdata[i,j],sep="")

            if(is.null(isspanCol(z,(i+1),(j+1)))){
                if(is.null(isspanRow(z,(i+1),(j+1)))){
                    result=1

                } else {
                    result=isspanRow(z,(i+1),(j+1))
                    if(result < 0) {
                        k=getspanRowData(z,i+1,j+1)
                        if(z$cellcolor[i+1,j+1]=="white") temp2=xdata[k+1,j]
                        else temp2=paste("\\cellcolor{",z$cellcolor[i+1,j+1],"}",
                                         xdata[k-1,j],sep="")
                        temp1=paste("\\multirow{",result,"}{*}{",temp2,"}",sep="")
                    }
                    else {
                        skip=1
                        result=0 #
                        if(z$cellcolor[i+1,j+1]=="white") skipcolor=""
                        else skipcolor=paste("\\cellcolor{",z$cellcolor[i+1,j+1],"}",sep="")
                    }

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
                        if(is.null(backcolor)) temp1=paste(temp1,"&",sep="")
                        else temp1=paste(temp1,"&\\cellcolor{",backcolor,"}",sep="")
                        #temp1=paste(temp1,"&",sep="")
                    }
                }

            } else {
                result=isspanCol(z,(i+1),(j+1))
                if(result>0) {
                    width=spanColWidth(z,(i+1),(j+1))
                    mcalign="c"
                    mclinecount=vlines[j+width+1]
                    if(mclinecount > 0) {
                        for(k in 1:mclinecount)
                            mcalign=paste(mcalign,"|",sep="")
                    }
                    temp1=paste("\\multicolumn{",result,"}{",mcalign,"}{",temp1,"}",sep="")
                    if(isGroupCol(j,result,colCount))
                        if(vlines[j+width+1]==0)
                            #if((j+result)<ncol(z$x))
                                temp1=paste(temp1,"&",sep="")
                    #if((j+result-1) %in% colCount[-length(colCount)])
                    #    if(vlines[j+result+1]==0) temp1=paste(temp1,"&",sep="")
                }
                else next
            }
            #browser()
            if(is.null(tempo)) {
                cat(temp1)
                tempo=temp1
            }
            else {
                if(result!=0) cat(paste("&",temp1,sep=""))
                else if(skip) cat(paste("&",skipcolor,sep=""))
            }

            if(!is.null(colCount)){
                count=j
                if(!is.null(isspanCol(z,i+1,j+1))){
                    result=isspanCol(z,(i+1),(j+1))
                    if(result>0) count=count+result
                }
                #if(count %in% colCount[-length(colCount)]) {
                    #if(vlines[count+2]==0) cat("&")
                    #if(z$cellcolor[i+1,j+1]=="white") cat("&")
                    #else cat(paste("&\\cellcolor{",z$cellcolor[i+1,j+1],"}",sep=""))
                #}
            }
        }
        cat(paste("\\\\ \n",sep=""))
        if(i %in% z$hline.after)
            cat(ifelse(z$booktabs,ifelse(i==nrow,"\\bottomrule[1.2pt]\n","\\midrule"),"\\hline\n"))
    }
    if(is.null(z$hline.after)) cat(ifelse(z$booktabs,"\\bottomrule[1.2pt]\n","\\hline\n"))

    footer=attr(z$x,"footer")
    if(!is.null(footer) & (z$show.footer)){
        myfooter=caption2minipage(z,footer)
        myfooter=gsub("~","$\\sim$",myfooter,fixed=TRUE)
        cat(paste("\\multicolumn{",totalCol,"}{l}{\\",Fontsize[headingsize],
                  "{",myfooter,"}}\\\\ \n",sep=""))
    }
    if(!is.null(z$caption) & z$caption.placement=="bottom"){
        mycaption=caption2minipage(z,z$caption)
        if(z$caption.bold) cat(paste("\\multicolumn{",totalCol,"}{",
                                     z$caption.position,"}{\\textbf{",mycaption,"}}\\\\ \n",sep=""))
        else cat(paste("\\multicolumn{",totalCol,"}{",
                       z$caption.position,"}{",mycaption,"}\\\\ \n",sep=""))
    }

    if(z$longtable) {
        if(!is.null(z$label)) cat(paste("\\label{",z$label,"}\n",sep=""))
        cat("\\end{longtable}\n")
        cat(paste("\\end{",Fontsize[z$size],"}\n",sep=""))
    } else {
        cat("\\end{tabular}\n")
        cat(paste("\\end{",Fontsize[z$size],"}\n",sep=""))
        if(!is.null(z$label)) cat(paste("\\label{",z$label,"}\n",sep=""))
        if(sort!="tabular") {
            if((sort=="table") | (sort=="sidewaystable"))
               cat(paste("\\end{",z$position,"}\n",sep=""))
            cat(paste("\\end{",sort,"}\n",sep=""))
        }
    }
    cat("\\color{black}\n")
}


#' Print Row Groups in a latex table
#'
#' @param z An object of class ztable
#' @param i An integer indicating row
printRowGroup=function(z,i){

    ncount=ncol(z$x)
    nrow=nrow(z$x)
    cn=colnames(z$x)
    addrow=ifelse(z$include.rownames,1,0)

    NewAlign=getNewAlign(z)
    totalCol=totalCol(z)
    colCount=colGroupCount(z)

    vlines=align2lines(z$align)

    printrgroup=1
    if(!is.null(z$n.rgroup)){
        if(length(z$n.rgroup)>1) {
            for(j in 2:length(z$n.rgroup)) {
                printrgroup=c(printrgroup,printrgroup[length(printrgroup)]+z$n.rgroup[j-1])
            }
        }
    }
    printrgroup
    printcline=0
    rgroupcount=0

    for(k in 1:length(printrgroup)){
        if(i == printrgroup[k]){
            rgroupcount=k
            break
        }
    }

if(i %in% printrgroup) {
    if(is.null(z$cspan.rgroup)){
        if(i>1) cat(paste("\\cline{1-",totalCol,"}\n",sep=""))
        vlines=align2lines(NewAlign)
        #mcalign=substr(extractAlign(NewAlign),start=1,stop=1)
        mcalign="l"
        if(vlines[1]>0)
            for(k in 1:vlines[1]) mcalign=paste("|",mcalign,sep="")
        if(vlines[totalCol+1]>0)
            for(k in 1:vlines[totalCol+1]) mcalign=paste(mcalign,"|",sep="")
        temp=paste("\\multicolumn{",totalCol,"}{",mcalign,"}{",sep="")
        if(z$colcolor[1]!="white")
            temp=paste(temp,"\\cellcolor{",z$colcolor[1],"}",sep="")
        temp=paste(temp,"\\textbf{",z$rgroup[rgroupcount],"}}",sep="")
        printcline=totalCol
    }
    else {
        if(z$cspan.rgroup==1) {
            if(z$colcolor[1]!="white")
                temp=paste("\\cellcolor{",z$colcolor[1],"}",sep="")
            else temp=""
            temp=paste(temp,"\\textbf{",z$rgroup[rgroupcount],"}",sep="")
            for(j in 1:(ncount+addrow-1)){
                temp1=""
                if(z$colcolor[j+1]!="white")
                    temp1=paste("\\cellcolor{",z$colcolor[j+1],"}",sep="")
                else {
                    if(!is.null(isspanRow(z,i+1,j+1))){
                        #cat("i=",i,",j=",j,"isspanRow(z,i,j+1)=",isspanRow(z,i+1,j+1),"\n")
                        if(isspanRow(z,i+1,j+1)<=0) {
                            #for(k in 1:nrow(z$spanRow)) {
                            #    if(z$spanRow[k,1]!=j+1) next
                            #    if(z$spanRow[k,2]>=i+1) next
                            #    if(z$spanRow[k,3]==i+1) break
                            #}
                            temp1=paste("\\cellcolor{",z$cellcolor[i+1,j+1],"}",sep="")
                        }
                    }
                    else temp1=""
                }
                temp=paste(temp,temp1,sep="&")
                if(!is.null(colCount)){
                    if(j %in% colCount[-length(colCount)]) {
                        if(vlines[j+2]==0) {
                            #if(z$colcolor[j+1]!="white")
                            #    temp=paste(temp,"&\\cellcolor{",z$colcolor[j+1],"}",sep="")
                            #else temp=paste(temp,"&",sep="")
                            temp=paste(temp,"&",sep="")
                        }
                    }
                }
            }
        } else {
            if(z$cspan.rgroup<1 | z$cspan.rgroup>(ncount+addrow))
                z$cspan.rgroup=ncount+addrow
            printcline=z$cspan.rgroup
            nvlines=align2lines(NewAlign)
            #mcalign=substr(extractAlign(NewAlign),start=1,stop=1)
            mcalign="l"
            if(nvlines[1]>0)
                for(k in 1:vlines[1]) mcalign=paste("|",mcalign,sep="")
            if(nvlines[printcline+1]>0)
                for(k in 1:nvlines[printcline+1]) mcalign=paste(mcalign,"|",sep="")
            temp=paste("\\multicolumn{",z$cspan.rgroup,"}{",mcalign,"}{\\textbf{",
                       "\\cellcolor{",z$colcolor[1],"}",
                       z$rgroup[rgroupcount],"}}",sep="")
            #temp=paste("\\cellcolor{",z$colcolor[1],"}",temp,sep="")
            if(z$cspan.rgroup<(ncount+addrow)) {
                for(j in (z$cspan.rgroup):(ncount+addrow-1)) {
                    if(z$colcolor[j+1]!="white")
                        temp=paste(temp,"&\\cellcolor{",z$colcolor[j+1],"}",sep="")
                    else {
                        if(!is.null(isspanRow(z,i+1,j+1))){
                            if(isspanRow(z,i+1,j+1)<=0) {
                                #for(k in 1:nrow(z$spanRow)) {
                                #    if(z$spanRow[k,1]!=j+1) next
                                #    if(z$spanRow[k,2]>=i+1) next
                                #    if(z$spanRow[k,3]==i+1) break
                                #}
                                temp=paste(temp,"&\\cellcolor{",z$cellcolor[i+1,j+1],"}",sep="")
                            }
                            else temp=paste(temp,"&",sep="")
                        }
                        else temp=paste(temp,"&",sep="")
                    }
                    if(!is.null(colCount)){
                        if(j %in% colCount[-length(colCount)]) {
                            if(vlines[j+2]==0){
                                #if(z$colcolor[z$cspan.rgroup+j]!="white")
                                #    temp=paste(temp,"&\\cellcolor{",z$colcolor[z$cspan.rgroup+j],"}",sep="")
                                #else temp=paste(temp,"&",sep="")
                                temp=paste(temp,"&",sep="")

                            }
                        }
                    }
                }
            }
        }
    }

    cat(paste(temp,"\\\\ \n",sep=""))
    if(printcline>0) cat(paste("\\cline{1-",printcline,"}\n",sep=""))
    rgroupcount=rgroupcount+1
}
}


#' Find valid color name
#'
#' @param a An integer or a character
#' @param mycolor predefined color names
#' @return a valid Latex color name
validColor=function(a,mycolor){
    if(is.numeric(a)) {
        if(a>0 && a <11)
            a=mycolor[a]
        else a="peach"
    } else {
        a=validColor2(a)
    }
    a
}

#' Find valid color name
#'
#' @param a An integer or a character
#' @return a valid Latex color name
validColor2=function(a){
    if(!is.character(a)) a="peach"
    else{
        result=grep(paste("^",a,sep=""),ztable::zcolors$name,ignore.case=TRUE)
        if(length(result)>0) a=ztable::zcolors$name[result[1]]
        else a="peach"
    }
    a
}

#' Define colors
#'
#' Define colors of mycolors
#' @param mycolors chracters vectors of color names
define_colors=function(mycolors) {
    if(is.null(mycolors)) return
    uniquecolors=unique(as.vector(unique(mycolors)))
    for(i in 1:length(uniquecolors)) {
        if(uniquecolors[i]=="white") next
        number=grep(paste("^",uniquecolors[i],sep=""),ztable::zcolors$name)
        if(length(number)<1) next
        else{
            definition=ztable::zcolors[number[1],3]
            cat(definition,"\n")
        }
    }
}

