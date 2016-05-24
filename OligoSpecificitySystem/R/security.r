"security"<-function(){
   if (length(txt4)==1 & length(txt1)==1 & length(txt2)==1) tkmessageBox(message="Error. Please, load before at least 2 databases")
   if (length(txt4)==1 & length(txt1)==1 & length(txt2)==1) stop()
   if (length(txt4)==1 & length(txt1)==1 & length(txt3)==1) tkmessageBox(message="Error. Please, load before at least 2 databases")
   if (length(txt4)==1 & length(txt1)==1 & length(txt3)==1) stop()
   if (length(txt4)==1 & length(txt3)==1 & length(txt2)==1) tkmessageBox(message="Error. Please, load before at least 2 databases")
   if (length(txt4)==1 & length(txt3)==1 & length(txt2)==1) stop()
   }