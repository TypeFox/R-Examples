UNCkeytron<-function(g, fout, LAB="KEY")
  {
    if(missing(LAB)) { LAB="EXAM KEY" }
    ##   source("/Users/lees/SCANTRON/keytron.R")
    
    ##  keytron(0, "sc.eps")
######  scantron answer sheet key
    N = length(g)
    
    cat(file=fout, "%!PS-Adobe-3.0 EPSF-3.0", sep="\n")

    cat(file=fout, "%whole page use: boundingbox: 0 0 612 792", sep="\n",append=TRUE)
    cat(file=fout, "%%Creator: Geotouch(Jonathan Lees)", sep="\n",append=TRUE)
    cat(file=fout, "%%Title: ", sep="\n",append=TRUE)
    cat(file=fout, "%%Creationdate: Mon Jul 22 11:48:40 EDT 1996", sep="\n",append=TRUE)
    cat(file=fout, "%%Pages: 1", sep="\n",append=TRUE)
    cat(file=fout, "%%DocumentFonts: (atend)", sep="\n",append=TRUE)
    cat(file=fout, "%%EndComments", sep="\n",append=TRUE)
    cat(file=fout, "%    POSTSCRIPT CREATED BY   J. M. Lees ", sep="\n",append=TRUE)
    cat(file=fout, "%             Thu Dec 10 21:38:01 EST 2009", sep="\n",append=TRUE)
    
    cat(file=fout, "%  ", sep="\n",append=TRUE)
    cat(file=fout, "%   the next 2 lines create a dictionary which needs to be ", sep="\n",append=TRUE)
    cat(file=fout, "%         finished off with the end command before showpage", sep="\n",append=TRUE)
    cat(file=fout, "%         these are needed primarily for older, postscript I printers", sep="\n",append=TRUE)


    cat(file=fout, "/cm { 28.35 mul } def", sep="\n",append=TRUE)
###############
    cat(file=fout, "/Times-Roman findfont 12 scalefont setfont 0 setgray", sep="\n",append=TRUE)

    cat(file=fout, "/HEL12 { /Helvetica findfont 12 scalefont setfont } def", sep="\n",append=TRUE)
    cat(file=fout, "/HEL10 { /Helvetica findfont 10 scalefont setfont } def", sep="\n",append=TRUE)
    cat(file=fout, "/HEL8 { /Helvetica findfont 8 scalefont setfont } def", sep="\n",append=TRUE)
###############
    cat(file=fout, "/censhow {", sep="\n",append=TRUE)
    
    cat(file=fout, " dup stringwidth", sep="\n",append=TRUE)
    cat(file=fout, "-.5 mul exch -.5 mul exch", sep="\n",append=TRUE)
    cat(file=fout, " rmoveto show } def", sep="\n",append=TRUE)
###############
    cat(file=fout, "/rightshow {", sep="\n",append=TRUE)
    cat(file=fout, " dup stringwidth pop", sep="\n",append=TRUE)
    cat(file=fout, " -1 mul", sep="\n",append=TRUE)
    cat(file=fout, "gsave 90 rotate", sep="\n",append=TRUE)
    cat(file=fout, " 0 rmoveto", sep="\n",append=TRUE)

    cat(file=fout, "show grestore } def", sep="\n",append=TRUE)

###############
    cat(file=fout, "/xcirc { /xs exch def", sep="\n",append=TRUE)

    cat(file=fout, "    newpath", sep="\n",append=TRUE)

    cat(file=fout, "    xs .5 mul  0 360 arc stroke  } def", sep="\n",append=TRUE)

###############
    cat(file=fout, "/xfillcirc { /xs exch def", sep="\n",append=TRUE)

    cat(file=fout, "    newpath", sep="\n",append=TRUE)
    cat(file=fout, "    xs  .5 mul  0 360 arc fill  } def", sep="\n",append=TRUE)
    
    cat(file=fout, "0 cm 0 cm translate", sep="\n",append=TRUE)
    
###############
    
    cat(file=fout, "/pwidth { 21.6 cm } def", sep="\n",append=TRUE)
    
    
    cat(file=fout, "1 setlinewidth", sep="\n",append=TRUE)
    
###############
    cat(file=fout, "1 cm 1 cm moveto", sep="\n",append=TRUE)
    
    cat(file=fout, "20.55 cm 1 cm lineto", sep="\n",append=TRUE)

    cat(file=fout, "20.55 cm 7.4 cm lineto", sep="\n",append=TRUE)
    
    cat(file=fout, "1  cm 7.4 cm lineto", sep="\n",append=TRUE)
    
    cat(file=fout, "closepath", sep="\n",append=TRUE)
    
    cat(file=fout, "stroke", sep="\n",append=TRUE)
    
###############
    cat(file=fout, "1  cm 7.4 cm moveto", sep="\n",append=TRUE)
    cat(file=fout, "1  cm 13.8  cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "6.6  cm 13.8  cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "6.6  cm 7.4 cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "closepath", sep="\n",append=TRUE)
    cat(file=fout, "stroke", sep="\n",append=TRUE)
###############


    

    cat(file=fout, "gsave", sep="\n",append=TRUE)
    cat(file=fout, " .3 .3 .3 setrgbcolor", sep="\n",append=TRUE)

    cat(file=fout, "19  cm 25.05 cm translate", sep="\n",append=TRUE)

    cat(file=fout, "0  cm 0 cm moveto", sep="\n",append=TRUE)
    cat(file=fout, "0  cm  1.7   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "1.4 cm  1.7   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "1.4  cm  0   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "closepath", sep="\n",append=TRUE)
    cat(file=fout, "fill", sep="\n",append=TRUE)
    cat(file=fout, "grestore", sep="\n",append=TRUE)

    cat(file=fout, "gsave", sep="\n",append=TRUE)
    cat(file=fout, " 1 1 1  setrgbcolor", sep="\n",append=TRUE)

    cat(file=fout, "19 .2 add  cm 25.05 .2 add  cm translate", sep="\n",append=TRUE)

    cat(file=fout, "0  cm 0 cm moveto", sep="\n",append=TRUE)
    cat(file=fout, "0  cm  1.7   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "1.4 cm  1.7   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "1.4  cm  0   cm lineto", sep="\n",append=TRUE)
    cat(file=fout, "closepath", sep="\n",append=TRUE)
    cat(file=fout, "save fill restore", sep="\n",append=TRUE)
    cat(file=fout, " 0 0 0  setrgbcolor", sep="\n",append=TRUE)

    cat(file=fout, "stroke", sep="\n",append=TRUE)


    cat(file=fout, "grestore", sep="\n",append=TRUE)


    cat(file=fout, " 0 0 0  setrgbcolor", sep="\n",append=TRUE)

###############




    
    cat(file=fout, "13.9 cm 1 cm moveto", sep="\n",append=TRUE)
    
    cat(file=fout, "13.9  cm 7.4 cm lineto", sep="\n",append=TRUE)
    
    cat(file=fout, "stroke", sep="\n",append=TRUE)

###############




    cat(file=fout, "/qsiz { .31 cm  } def", sep="\n",append=TRUE)
    
    
    cat(file=fout, ".1 setlinewidth", sep="\n",append=TRUE)

    cat(file=fout, "/Times-Roman findfont 24 scalefont setfont 0 setgray", sep="\n",append=TRUE)

    cat(file=fout, "2 cm 2 cm moveto", sep="\n",append=TRUE)
    itit = paste(sep=" ","(",LAB, ")")

    cat(file=fout, paste(sep=" ", itit, "show"), sep="\n",append=TRUE)

    cat(file=fout, "HEL10", sep="\n",append=TRUE)
    cat(file=fout, "/delbub  { .42 cm } def", sep="\n",append=TRUE)
    cat(file=fout, "/DX { .85  cm } def", sep="\n",append=TRUE)

    cat(file=fout, "/arow { 0 1 4 { delbub mul  y add x exch   qsiz  xcirc   }  for  } def", sep="\n",append=TRUE)

    M = 0
###  21.6 cm
    DX =.85
    pwidth=  21.6
    
    Ax = pwidth - 4.95
    Ay = 8.05
    charshift = "x .1 cm add y .3 cm sub moveto"


    for(i in 0:4)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, " arow", sep="\n",append=TRUE)

        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)

            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            
          }


        
      }

    Ay = 8.05 +  3.4

    for(i in 0:4)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, "   arow", sep="\n",append=TRUE) 
        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)

            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            
          }


      }

    Ay =  8.05  +     3.4 *  2 

    Ax =  1.35 

    for(i in 0:22)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, "   arow", sep="\n",append=TRUE) 
        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)

            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            
            
          }


      }
    
    Ay =  8.05  +     3.4 *  3

    Ax =  1.35 

    for(i in 0:22)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, "   arow", sep="\n",append=TRUE) 
        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)


            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            

            
          }


      }
    
    Ay =  8.05  +     3.4 *  4

    Ax =  1.35 

    for(i in 0:22)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, "   arow", sep="\n",append=TRUE) 
        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)

            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            

            
          }


      }
    Ay =  8.05  +     3.4 *  5

    Ax =  1.35 

    for(i in 0:20)
      {
        x = DX*i +Ax
        y =  Ay
        cat(file=fout, paste(sep=" ",x , "cm /x exch def"),  sep="\n",append=TRUE)
        cat(file=fout, paste(sep=" ",y , "cm /y exch def"),  sep="\n",append=TRUE)
        
        cat(file=fout, "   arow", sep="\n",append=TRUE) 
        M = M+1
        if(M<=N)
          {
            bub = g[M]-1
            cat(file=fout, paste(sep=" ",bub,"delbub mul  y add x exch   qsiz xfillcirc"),  sep="\n",append=TRUE)
            cnum = paste(sep="", "(", M, ".)")
            cat(file=fout,paste(charshift, cnum , "rightshow"  ),  sep="\n",  append=TRUE    )
            

            
          }


      }
    

    
    cat(file=fout, "showpage", sep="\n",append=TRUE)

  }

