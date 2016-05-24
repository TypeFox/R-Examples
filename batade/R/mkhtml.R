mkhtml <-
function(filename, data, foot=TRUE, charset="CP932", lang="JP"){
	require(hwriter)
	p <- openPage(filename, charset=charset, lang=lang)
  if(is.data.frame(data))data <- as.matrix(data)
  hwrite(hmakeTag("style","
  .LL{
  font-size: 300%;
  color: #FFFFFF;
  background-color: #009933;
  padding: 4px;
  }

  .L{
  font-size:200%;
  border: solid;
  border-width: 0px 0px 2px 0px;
  border-color: #009933;  
  padding: 4px;
  }

  .M{
  font-size: 100%;
  border: solid;
  border-width: 0px 0px 0px 5px;
  border-color: #009933;
  padding: 4px;
  }

  .S{
  padding: 4px;  
  }

  .foot{
  text-align: right;
  font-size: 100%;
  border: solid;
  border-width: 2px 0px 0px 0px;
  border-color: #009933;
  padding: 4px;  
  }

  "), p)
  hwrite(data[1,1], p, br=TRUE, contenteditable="TRUE", class="LL", div=TRUE)
  targ <- ifelse(foot, nrow(data)-1, nrow(data))
  for(i in 2:targ){
    ROW <- data[i,]
    content <- unlist(strsplit(split="\\|", ROW[1]))
     if(length(content)>1){
         if(isTRUE(as.logical(grep(".*\\.png$|.*\\.jpg$|.*\\.jpeg$|.*\\.gif$|.*\\.tiff$", content, perl=TRUE))) && nchar(content) < 15){
         type <- "DI"
         }else if(isTRUE(as.logical(grep(".*\\.htm$|.*\\.html$", content, perl=TRUE))) && nchar(content) < 15){
         type <- "DH"
         }       
     }else{
       if(isTRUE(as.logical(grep(".*\\.png$|.*\\.jpg$|.*\\.jpeg$|.*\\.gif$|.*\\.tiff$", content, perl=TRUE))) && nchar(ROW[1]) < 15){
         type <- "I"
         }else if(isTRUE(as.logical(grep(".*\\.htm$|.*\\.html$", content, perl=TRUE)))){
         type <- "H"
         }else{
         type <- "T"
         }
     }
     switch(type,
       "T" = hwrite(content, p, br=TRUE, contenteditable="TRUE", class=ROW[2], div=TRUE),
       "I" = hwriteImage(content, p, br=TRUE, div=TRUE, center=TRUE, border=0),
       "DI" = hwriteImage(content, p, br=TRUE, div=TRUE, center=TRUE, border=0),
       "H" = hwrite(paste('<iframe src=', content, ' frameborder="0" width="1200" height="600" scrolling="no"></iframe>', sep=""), p, center=TRUE),
       "DH" = hwrite(c(paste('<iframe src=', content, ' frameborder="0" width="600" height="400" scrolling="no"></iframe>', sep=""),paste('<iframe src=', content, ' frameborder="0" width="600" height="400" scrolling="no"></iframe>', sep="")), p, center=TRUE),
       stop(message = "rayout error")
       )
    }
  if(foot){
    hwrite(data[nrow(data),1], p, br=TRUE, contenteditable="TRUE", class="foot", div=TRUE)
    closePage(p, splash=FALSE)
    }else{
    closePage(p, splash=FALSE)
    }
  }

