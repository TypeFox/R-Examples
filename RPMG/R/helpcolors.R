`helpcolors` <-
function()
  {

   
    cat(file = "", "\n")
    
    pl = c("grey", "lightblue1", "pink", "darkseagreen2", "gold1",
      "chartreuse1", "aquamarine", "plum1", "goldenrod", "maroon1",
      "deepskyblue", "palegreen2", "salmon")


    a = paste("\"", pl, "\"", sep="")

    
    b = as.character(paste(a,  collapse = ","))


    
    
    cat(file = "", paste(collapse = "", 'plite', "=c(", b, ")"), fill = TRUE)

    pd=c('black','darkmagenta','forestgreen','blueviolet','tan3','lightseagreen','deeppink','cyan3','bisque3','magenta1','lightsalmon3','darkcyan','darkslateblue','chocolate4','goldenrod4','mediumseagreen')
    cat(file = "", "\n")
    a = paste("\"", pd, "\"", sep="")

    b = as.character(paste(a,  collapse = ","))


    
    cat(file = "", paste(collapse = "", 'pdark', "=c(", b, ")"), fill = TRUE)
    cat(file = "", "\n")
    
    ##   cat(file = "","Just installed: plite  pdark XMCOL ", sep="\n")

    cat(file = "","### to change the color palette use:  palette(yourpal)  ", sep="\n")
    
    cat(file = "","###  Example: palette(plite)  ", sep="\n")
    
    cat(file = "","###  Example:  palette(pdark)  ", sep="\n")
     cat(file = "", "\n")
    cat(file = "","###  for dark background:  ", sep="\n")

    cat(file = "","par(bg=rgb(0.1,0.1,0.1), fg=rgb(1,1,1),col.axis=rgb(1,1,1),col.lab=rgb(1,1,1),col.main=rgb(1,0,0) ,col.sub=rgb(1,1,1)   )", "\n")
    cat(file = "", "\n")


    
    ##   plot(rnorm(10) , rnorm(10), col=1:10 )
    
    cat(file = "","SHOWPAL(plite)", sep="\n\n")
    cat(file = "","SHOWPAL(plite, NAME=TRUE)", sep="\n\n")
    cat(file = "","SHOWPAL(plite, NAME=TRUE, NUM=TRUE)", sep="\n\n")
    
     cat(file = "", "\n")
    
    cat(file = "","SHOWPAL(pdark)", sep="\n\n")
    
    cat(file = "", "\n")

    
    cat(file = "","XMCOL =  setXMCOL()", sep="\n")

    cat(file = "","SHOWPAL(XMCOL)", sep="\n")
    
  cat(file = "", "\n")
 cat(file = "","SYSCOL =  colors()", sep="\n")
  cat(file = "","SHOWPAL(SYSCOL, ncol=10)", sep="\n")

  cat(file = "", "\n")
 cat(file = "","SYSCOL =  hcl(h=seq(from=0, to=200, by=10) ) ", sep="\n")
  cat(file = "","SHOWPAL(SYSCOL, NAME=TRUE)", sep="\n")

  cat(file = "", "\n")
 cat(file = "","SYSCOL =  hcl(h=seq(from=0, to=360, by=10) ) ", sep="\n")
  cat(file = "","SHOWPAL(SYSCOL, NAME=TRUE)", sep="\n")

  cat(file = "", "\n")
 cat(file = "","SYSCOL = sepia.colors(100)  ", sep="\n")
  cat(file = "","SHOWPAL(SYSCOL, ncol=10) ", sep="\n")

   cat(file = "","SYSCOL =  pastel.colors(100)  ", sep="\n")
  cat(file = ""," ", sep="\n")

    cat(file = "", "\n")




    return()

  }

