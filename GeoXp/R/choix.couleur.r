choix.couleur <- function(graphChoice,listvar=NULL,listnomvar=NULL,varChoice1=NULL,legends,col,pch, spdf=FALSE)
 {
    if((graphChoice == "Barplot")||(graphChoice == "Moran"))
     { if(graphChoice == "Barplot")
        {method <- "Cluster"
         n.col <- length(levels(as.factor(listvar[,which(listnomvar == varChoice1)])))
         labmod <- levels(as.factor(listvar[,which(listnomvar == varChoice1)]))
        }
       else
        {method <-"Quadrant"
         n.col <- 4
         # labmod is not used but it permits to uniform the programm
         labmod <- c("Low-High","High-High","Low-Low","High-Low")
        }

      tt2 <- tktoplevel()

       OnOK <- function()
        {

         tt <- tktoplevel()
         txt <- tktext(tt, height=8)

         tkpack(tklabel(tt,text=paste("Please, could you give a vector of ",n.col," colors like c(\"orange\", \"purple\", etc) or colors()[99:",98+n.col,"]",sep="")) )
         tkpack(txt)

         run <- function()
          {
           code <- tclvalue(tkget(txt,"0.0","end"))
           e <- try(parse(text=code))

            if (inherits(e, "try-error"))
             {
              tkmessageBox(message="Syntax error", icon="error")
              return()
             }


            col2<<-eval(e)
            col3<<-col2

           if (length(col2)!=n.col)
            {
             tkmessageBox(message="Length of vector is not good", icon="error")
             return()
            }
           else
            {
             if(!spdf) {tt2 <- tktoplevel()}

              OnOK3 <- function()
               {

                tt <- tktoplevel()
                txt <- tktext(tt, height=8)

                tkpack(tklabel(tt,text=paste("Please, could you give a vector of ",n.col," symbols like c(3, 4, etc)")))
                tkpack(txt)


                run <- function()
               {
                code <- tclvalue(tkget(txt,"0.0","end"))
                e <- try(parse(text=code))

                 if (inherits(e, "try-error"))
                  {tkmessageBox(message="Syntax error", icon="error")
                    return()
                  }

                pch2<<-eval(e)

                 if (length(pch2)!=n.col)
                  {tkmessageBox(message="Length of vector is not good", icon="error")
                   return()
                  }



                   OnOK5 <- function()
                    {

                      msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
                      tkmessageBox(message=msg)

                      dev.set(2)
                      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')

                      loc <- locator(1)
                      loc$name <- varChoice1

                      legends<<-list(legends[[1]],TRUE,legends[[3]],loc)
                      #print(legends)
                      tkdestroy(tt1)
                     }


                    OnOK6 <- function()
                     {
                       legends<<-list(legends[[1]],FALSE,legends[[3]],"")
                       tkdestroy(tt1)
                     }


                  tt1<-tktoplevel()
                  ifelse(method=="Quadrant",labelText12 <- tclVar("Do you want a legend for quadrants on the map"),
                  labelText12 <- tclVar("Do you want a legend for factors on the map"))
                  label12 <- tklabel(tt1,justify = "center", wraplength = "3i", text=tclvalue(labelText12))
                  tkconfigure(label12, textvariable=labelText12)
                  tkgrid(label12,columnspan=2)

                  point.but <- tkbutton(tt1, text="  Yes  ", command=OnOK5);
                  poly.but <- tkbutton(tt1, text=" No ", command=OnOK6);
                  tkgrid(point.but, poly.but)
                  tkgrid(tklabel(tt1,text="    "))
                  tkfocus(tt1)
                  tkwait.window(tt1)
                  tkdestroy(tt)

                }

              topMenu <- tkmenu(tt)
              tkconfigure(tt, menu=topMenu)
              fileMenu <- tkmenu(topMenu, tearoff=FALSE)
              tkadd(topMenu, "command", label="Run",command=run)
              tkwait.window(tt)
              tkdestroy(tt2)
             }

            OnOK4 <- function()
             {
              tt1<-tktoplevel()
              pch2 <<- pch[1]

                  OnOK5 <- function()
                   {
                    tkdestroy(tt1)
                    msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
                    tkmessageBox(message=msg)

                    dev.set(2)
                    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')

                    loc <- locator(1)
                    loc$name <- varChoice1
                    legends<<-list(legends[[1]],TRUE,legends[[3]],loc)
                   }


                  OnOK6 <- function()
                   {
                    legends<<-list(legends[[1]],FALSE,legends[[3]],"")
                    tkdestroy(tt1)
                   }



              ifelse(method=="Quadrant",labelText12 <- tclVar("Do you want a legend for Quadrant on the map"),
              labelText12 <- tclVar("Do you want a legend for factors on the map"))
              label12 <- tklabel(tt1,justify = "center", wraplength = "3i", text=tclvalue(labelText12))
              tkconfigure(label12, textvariable=labelText12)
              tkgrid(label12,columnspan=2)

              point.but <- tkbutton(tt1, text="  Yes  ", command=OnOK5);
              poly.but <- tkbutton(tt1, text=" No ", command=OnOK6);
              tkgrid(point.but, poly.but)
              tkgrid(tklabel(tt1,text="    "))
              tkfocus(tt1)
              tkwait.window(tt1)
              tkdestroy(tt2)
            }

          if(!spdf)
           {labelText13 <- tclVar("Do you want also symbols on the map ?")
           label13 <- tklabel(tt2,justify = "center", wraplength = "3i", text=tclvalue(labelText13))
           tkconfigure(label13, textvariable=labelText13)
           tkgrid(label13,columnspan=2)

           point2.but <- tkbutton(tt2, text="  Yes  ", command=OnOK3);
           poly2.but <- tkbutton(tt2, text=" No ", command=OnOK4);
           tkgrid(point2.but, poly2.but)
           tkgrid(tklabel(tt2,text="    "))
           tkfocus(tt2)
           tkwait.window(tt2)
           }
           else
           {OnOK4()}
           tkdestroy(tt)
          }

         }

       topMenu <- tkmenu(tt)
       tkconfigure(tt, menu=topMenu)
       fileMenu <- tkmenu(topMenu, tearoff=FALSE)
       tkadd(topMenu, "command", label="Run",
       command=run)
       tkwait.window(tt)
       tkdestroy(tt2)

      }


     OnOK2 <- function()
      {
       tkdestroy(tt2)
       col2 <<- "blue"
       col3 <<- col[1]
       pch2 <<-pch[1]
      }

    ifelse(method=="Quadrant",labelText12 <- tclVar("Do you want different colors on the quadrant on the moran plot and on the map ?"),
    labelText12 <- tclVar("Do you want different colors on the barplot and on the map ?"))
    label12 <- tklabel(tt2,justify = "center", wraplength = "3i", text=tclvalue(labelText12))
    tkconfigure(label12, textvariable=labelText12)
    tkgrid(label12,columnspan=2)

    point.but <- tkbutton(tt2, text="  Yes  ", command=OnOK);
    poly.but <- tkbutton(tt2, text=" No ", command=OnOK2);
    tkgrid(point.but, poly.but)
    tkgrid(tklabel(tt2,text="    "))
    tkfocus(tt2)
    tkwait.window(tt2)
    }
  else
    {
     legends<-list(legends[[1]],FALSE,legends[[3]],"")
     method <- ""
     col2 <- "blue"
     col3 <- col[1]
     pch2 <- pch[1]
     labmod <- ""
     }

 return(list(method=method,col2=col2,col3=col3,pch2=pch2,legends=legends,labmod=labmod))
}