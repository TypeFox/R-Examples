butdoc<-function(tag, doc, NEW=FALSE)
  {

    ###  button documentation
###############  given a vector of tags and vector docs (explanations)
    ####  click and display select documentations
###    tag = c(tag)
###      doc   = c(doc)


    pcols = pastel.colors(length(tag), seed=1)

    if(missing(NEW)) NEW=FALSE
    
    if(NEW) dev.new()


    upar = par(  no.readonly = TRUE)
    nf <- layout(matrix(c(1,2), 2,1, byrow=TRUE), widths = c(1,1), heights = c(3,1), respect=FALSE)


    P3 = NULL

    while(TRUE)
      {
        par(mai=c(0,0,.4,0)) 
        pkachoo =  chooser( tag, ncol=5, nsel=NA, newdev=FALSE, STAY=TRUE,  cols=pcols)

        pkachoo =  unique(pkachoo)

        mPP = match(pkachoo, tag)

        PP =    paste(sep=" = " , tag[match( pkachoo , tag)  ] , doc[match( pkachoo , tag)  ] )

        acols = pcols[mPP[!is.na(mPP)]    ]
        
        PP = c(PP, "MORE")
        acols = c(acols, "white")

        
        
        par(mai=c(0,0,0,0))
        
        P2 = chooser( PP, ncol=1, nsel=1, newdev=FALSE, STAY=TRUE, main="",  cols=acols )

        if(length(P2)<1) break
        if(is.na(P2) ) break
        if(any(P2!="MORE")) break

        P3 = P2
      }

    par( upar ) 


    invisible(P3)

    
  }

