`aln2html` <-
function(aln, file = "alignment.html",
                     Entropy = 0.5,
                     append = TRUE,
                     caption.css = "color: gray; font-size: 9pt",
                     caption = "Produced by <a href=http://thegrantlab.org/bio3d/>Bio3D</a>",
                     fontsize = "11pt",
                     bgcolor = TRUE,
                     colorscheme="clustal") {

  if(is.list(aln)) {
    x=aln$ali
    id=aln$id
  } else {
    x=aln
    id=dimnames(x)[[1]]
  }

  if(is.null(id))
    stop("No $id list component or rownames for the alignment object")
  
  back <- ""; bold <- "" #bold <- "font-weight: bold"
  if (bgcolor) { back <- "background-"; bold <- "" }


  head <- paste("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n<html>\n  <head>\n
    <STYLE type=\"text/css\">
      SPAN.A { ", back, "color:  #80b3e6; ", bold, " } /* Small and hydrophobic
*/
      SPAN.V { ", back, "color:  #80b3e6; ", bold, " }
      SPAN.L { ", back, "color:  #80b3e6; ", bold, " }
      SPAN.I { ", back, "color:  #80b3e6; ", bold, " }
      SPAN.M { ", back, "color:  #80b3e6; ", bold, " }
      SPAN.F { ", back, "color:  #80b3e6; ", bold, " }
      SPAN.W { ", back, "color:  #80b3e6; ", bold, " }

      SPAN.S { ", back, "color:  #1acc1a; ", bold, " } /* Hydroxyl and amine amino acids */
      SPAN.T { ", back, "color:  #1acc1a; ", bold, " }
      SPAN.N { ", back, "color:  #1acc1a; ", bold, " }
      SPAN.Q { ", back, "color:  #1acc1a; ", bold, " }

      SPAN.D { ", back, "color:  #cc4dcc; ", bold, " } /* Charged */
      SPAN.E { ", back, "color:  #cc4dcc; ", bold, " }

      SPAN.R { ", back, "color:  #e6331a; ", bold, " } /* Charged */
      SPAN.K { ", back, "color:  #e6331a; ", bold, " }

      SPAN.H { ", back, "color:  #1ab3b3; ", bold, " } /* Histidine and tyrosine */
      SPAN.Y { ", back, "color:  #1ab3b3; ", bold, " }

      SPAN.G { ", back, "color:  #e6994d; ", bold, " }

      SPAN.P { ", back, "color:  #cccc00; ", bold, " }

      SPAN.C { ", back, "color:  #80b3e6; ", bold, " }

      SPAN.boxres { border-width: 1; border: solid; text-align: center}
      SPAN.bot  { ", back, "color:  #ccffff; ", bold, " }
      SPAN.low  { ", back, "color:  #66ccff; ", bold, " }
      SPAN.mid  { ", back, "color:  #6699ff; ", bold, " }
      SPAN.high { ", back, "color:  #6666ff; ", bold, " }
      #helix { background-color: #ffcccc; }
      #strand { background-color: green }

    </STYLE>\n  </head>\n", collapse="",sep="" )


  body <- paste("  <body>
  <div style=\"font-size: ",fontsize, "; font-family: mono, courier; font-weight: normal; margin: 10px; white-space:nowrap\">\n")


##  id <- paste("\n<br /><a>",dimnames(x)[[1]],"&nbsp;&nbsp;&nbsp;&nbsp;</a>",sep="")

  #- id justification
#  id <- dimnames(x)[[1]]
  len <- nchar(id, type="chars")
  pad <- NULL; if(!all(max(len)==len)) {
    for(i in 1:length(id)) {
      pad <- c(pad, paste( rep("&nbsp;", max(len)-len[i]),collapse="" ))
    }
  }
  id <- paste("\n<br /><a>",id,pad,"&nbsp;&nbsp;&nbsp;</a>",sep="")


  if(colorscheme=="clustal") {
    # Clustal coloring
    al <- matrix( paste("<SPAN class=\"",x,"\">",x,"</SPAN>",sep=""), nrow=nrow(x))
  } else {
    #- Color by entropy score
    he <- entropy(x)
    score <- he$H.10.norm; score[ which(he$freq[c("-"),]>0.6) ] = 0

    rn <- cbind( (score > 0.4), (score > 0.575), (score > 0.75), (score > 0.9) )
    rn[ rn[,2], 1] = FALSE; rn[ rn[,3], 2] = FALSE; rn[ rn[,4], 3] = FALSE

    al=x
    b <- matrix( paste("<SPAN class=\"bot\">",x,"</SPAN>",sep=""), nrow=nrow(x))
    l <- matrix( paste("<SPAN class=\"low\">",x,"</SPAN>",sep=""), nrow=nrow(x))
    m <- matrix( paste("<SPAN class=\"mid\">",x,"</SPAN>",sep=""), nrow=nrow(x))
    h <- matrix( paste("<SPAN class=\"high\">",x,"</SPAN>",sep=""), nrow=nrow(x))

    al[ , which(rn[,1])] = b[ , which(rn[,1])]
    al[ , which(rn[,2])] = l[ , which(rn[,2])]
    al[ , which(rn[,3])] = m[ , which(rn[,3])]
    al[ , which(rn[,4])] = h[ , which(rn[,4])]
  }

  #- Dont color unconserved positions
  if(Entropy > 0) {
    if(colorscheme=="clustal")
      he <- entropy(x)
    execlude <- unique( c(which(he$H.10.norm < Entropy),
                          which(he$freq[c("-"),]>0.6)) )
    al[,execlude] = x[,execlude]
  }
  #- Dont color gaps
  ind<-which(x=="-",arr.ind=TRUE); al[ind]=x[ind]


  cat(head, body, file=file, append=append)
  for(i in 1:length(id))
    cat( id[i], al[i,], sep="", file=file, append=TRUE)
  cat(paste("\n  </div>\n",
            "<div style=\"", caption.css,"\">", caption ,"</div>\n",
            "    </body>\n</html>\n"), file=file, append=TRUE)


}

