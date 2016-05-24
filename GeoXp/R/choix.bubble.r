choix.bubble<-function(buble,listvar,listnomvar,legends)
{
 if(!buble)
  {
   if ((length(listvar) != 0) && (length(listnomvar) != 0))
    {
    if(listnomvar[1]!="ilocal")
     {varChoix <- choixvarfunc("Choice of variables","Choose a variable",listnomvar)
      bubble <- listvar[,which(listnomvar == varChoix)]
     }
     else
     {varChoix <- "abs(LISA)"
      bubble <- listvar[,which(listnomvar == "ilocal")]
     }
       if ((length(bubble) != 0)&&(min(bubble)>=0))
        {
         if(varChoix != "chi2.quant")
         {buble<-TRUE

         tt2 <- tktoplevel()
         z<-sqrt(abs(bubble)/max(abs(bubble)))*3
         legmap<-NULL

         OnOK <- function()
          {
            tt3 <- tktoplevel()

            bubl<-function()
             {
              tkdestroy(tt3)
               msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
               tkmessageBox(message=msg)
               ifelse(listnomvar=="ilocal",dev.set(3),dev.set(2))
               title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')

               loc <- locator(1)

               legmap<<-c(sqrt(abs(as.numeric(tclvalue(ma)))/max(abs(bubble)))*3,sqrt(abs(as.numeric(tclvalue(mea)))/max(abs(bubble)))*3,
               sqrt(abs(as.numeric(tclvalue(mi)))/max(abs(bubble)))*3,as.numeric(tclvalue(ma)),as.numeric(tclvalue(mea)),
               as.numeric(tclvalue(mi)),varChoix)

               legends<<-list(TRUE,legends[[2]],loc,legends[[4]])
             }

            mi<-tclVar(round(min(bubble),2))
            mea<-tclVar(round(mean(bubble),2))
            ma<-tclVar(round(max(bubble),2))
            entry.Name <-tkentry(tt3,width="5",textvariable=mi)
            entry.Name2 <-tkentry(tt3,width="5",textvariable=mea)
            entry.Name3 <-tkentry(tt3,width="5",textvariable=ma)

            tkgrid(tklabel(tt3,text="Break Points:"))
            tkgrid(tklabel(tt3,text="Small Bubble"),entry.Name)
            tkgrid(tklabel(tt3,text="Middle Bubble"),entry.Name2)
            tkgrid(tklabel(tt3,text="Large Bubble"),entry.Name3)

            autre.but <- tkbutton(tt3, text="     OK     " , command=bubl);
            tkgrid(autre.but,columnspan=2)
            tkgrid(tklabel(tt3,text="    "))
            tkfocus(tt3)
            tkwait.window(tt3)
            tkdestroy(tt2)
          }

        OnOK2 <- function()
          {
            tkdestroy(tt2)
          }


         labelText12 <- tclVar("Do you want a legend for bubbles")
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
          msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
          tkmessageBox(message=msg)
          dev.set(2)
          title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main='red')

          buble<-TRUE
          loc <- locator(1)
           z<-sqrt(abs(bubble+1)/max(abs(bubble)+1))*2.3
           #print(z)
           legmap<-c(3,sqrt(4/5)*2.3,sqrt(3/5)*2.3,sqrt(2/5)*2.3,sqrt(1/5)*2.3,"chi2.quant")
           legends<-list(TRUE,legends[[2]],loc,legends[[4]])
         }
        }
       else
        {
         tkmessageBox(message="Bubbles have not been given or variable is not strictly positive",icon="warning",type="ok")
         buble<-FALSE
         legends<-list(FALSE,legends[[2]],"",legends[[4]])
         legmap<-NULL
         z<-NULL
        }
     }

   else
     {
        tkmessageBox(message="To use Bubbles, the lists wich contain the variables and their names must have been given",icon="warning",type="ok");
        buble<-FALSE
        legends<-list(FALSE,legends[[2]],"",legends[[4]])
        z<-NULL
        legmap<-NULL
     }
  }
  else
  {buble<-FALSE
   legends<-list(FALSE,legends[[2]],"",legends[[4]])
   z<-NULL
   legmap<-NULL
  }

return(list(buble=buble,legends=legends,legmap=legmap,z=z))

}