wdTable <-
    function (data,
              caption = "",
              caption.pos="below",
              bookmark = NULL,
              pointsize = 9,
              padding = 5,
              autoformat = 1,
              row.names=TRUE,
              align = if(row.names) c("l", rep("r", ncol(data))) else c(rep("r",ncol(data))),
              hlines=NULL,
              wdapp = .R2wd)
{
    if (autoformat<0) stop("inadmissible autoformat")
    if (!is.null(hlines) && length(hlines)>nrow(data)+1) stop("length of hlines must be equal to the number of rows in the table + 1")
    wdsel <- try(wdapp[["Selection"]])
    if (is.null(wdsel)||class(wdsel)=="try-error") stop("Word not connected. Run wdGet() first")
    wdopt<-wdapp[["Options"]]
    wddoc <- wdapp[["ActiveDocument"]]
    wdsel$TypeParagraph()
    wdInsertBookmark("R2wdEndmark")
    bookmarkcounter <- wddoc[["Bookmarks"]][["Count"]]
    wdsel$MoveUp()
    nr <- nrow(data)
    nc <- ncol(data)
    if (row.names){
    out <- matrix("", nrow = nr + 1, ncol = nc + 1)
    out[1 + (1:nr), 1 + (1:nc)] <- as.matrix(data)
    out[1, 1 + (1:nc)] <- colnames(data)
    out[1 + (1:nr), 1] <- row.names(data)
     } else{
    out <- matrix("", nrow = nr + 1, ncol = nc )
    out[1 + (1:nr),(1:nc)] <- as.matrix(data)
    out[1, (1:nc)] <- colnames(data)
     }
    tt<-paste(apply(out,1,paste,collapse="\t"),collapse="\n")
    wdsel[['Text']]<-tt
    tab<-wdsel[['Range']]$ConvertToTable(1,nr+1,nc+ifelse(row.names,1,0))
tryout<-try({
    tab$AutoFormat(autoformat)
    if (as.numeric(.R2wd[['Version']])>10)
    {
        tabrows<-tab[["Rows"]]
        try(tabrows[["Height"]] <- pointsize + padding,silent=TRUE)
        try(tabrows[["HeightRule"]] <- 2,silent=TRUE)
        tabcells<-tab[["Range"]][["Cells"]]
        try(tabcells[["VerticalAlignment"]] <- 1,silent=TRUE)
    }
    tab$AutoFitBehavior(1)
    tab[["Range"]]$Select()
    wdsel[["Font"]][["Size"]] <- pointsize
    ## now deal with alignment and vertical lines
    ## check if the table starts with a vertical line
    if (align[1]=="|") {
            tab[["Columns"]]$Item(1)$Select()
            tmp<-wdsel[["Borders"]]$Item(-2)
            tmp[["LineStyle"]]<-wdopt[["DefaultBorderLineStyle"]]
            tmp[["LineWidth"]]<-wdopt[["DefaultBorderLineWidth"]]
            tmp[["Color"]]<-wdopt[["DefaultBorderColor"]]
            align<-align[-1]
        }
    ## now treat the rest
    ii<-0
    for (i in 1:length(align)) {
        if (align[i]=="|"){
            tab[["Columns"]]$Item(ii)$Select()
            tmp<-wdsel[["Borders"]]$Item(-4)
            tmp[["LineStyle"]]<-wdopt[["DefaultBorderLineStyle"]]
            tmp[["LineWidth"]]<-wdopt[["DefaultBorderLineWidth"]]
            tmp[["Color"]]<-wdopt[["DefaultBorderColor"]]
        } else {
            ii<-ii+1
            tab[["Columns"]]$Item(ii)$Select()
            wdselpar<-wdsel[["ParagraphFormat"]]
            wdselpar[["Alignment"]] <- c(l = 0,c = 1, r = 2)[align[i]]
        }
    }
    ## now add row lines if asked for
    if (!is.null(hlines)){
        for (i in 1:length(hlines)){
            if(hlines[i]!="n"){
                tab[["Rows"]]$Item(i)$Select()
                ## add bottom lines
                if(hlines[i]%in%c("b","bt")){
                    tmp<-wdsel[["Borders"]]$Item(-3)
                    tmp[["LineStyle"]]<-wdopt[["DefaultBorderLineStyle"]]
                    tmp[["LineWidth"]]<-wdopt[["DefaultBorderLineWidth"]]
                    tmp[["Color"]]<-wdopt[["DefaultBorderColor"]]
                }
                ## add top lines
                if(hlines[i]%in%c("t","bt")){
                    tmp<-wdsel[["Borders"]]$Item(-1)
                    tmp[["LineStyle"]]<-wdopt[["DefaultBorderLineStyle"]]
                    tmp[["LineWidth"]]<-wdopt[["DefaultBorderLineWidth"]]
                    tmp[["Color"]]<-wdopt[["DefaultBorderColor"]]
                }
            }
        }
    }
    tab[["Range"]]$Select()
    caption <- paste(" ", caption, sep = "")
    wdsel$InsertCaption("Table", caption, "",
                        ifelse(caption.pos=="above",0,1), 0)

    if (is.null(bookmark))
        bookmark <- paste("Table", bookmarkcounter + 1, sep = "")
    tab[["Range"]]$Select()
    wdInsertBookmark(bookmark)
    wddoc[["Bookmarks"]]$Item(bookmarkcounter)$Select()
    })
    if (class(tryout)=="try-error") {
        warning("Error in table construction, removing")
        tab$Delete()
    }
    wdGoToBookmark("R2wdEndmark")
    return()
}
