tz2offset <- function(tz){
      to.ret <- NA

      if(as.character(tz)=="CEST" || as.character(tz)=="UTC" || as.character(tz)=="GMT" || as.character(tz)=="0" || as.character(tz)=="000"  || as.character(tz)=="0000" || as.character(tz)=="+0000"
		|| as.character(tz)=="W. Europe Daylight Time"  || as.character(tz)=="Mitteleurop\u00E4ische Zeit"){
	      to.ret=0
      }
      if(as.character(tz)=="-1" ||  as.character(tz)=="-100" || as.character(tz)=="-0100"){
	      to.ret=-100
      }
      if(as.character(tz)=="1" ||  as.character(tz)=="100" || as.character(tz)=="0100" || as.character(tz)=="+0100" || as.character(tz)=="W. Europe Standard Time"
	  ){
	      to.ret=100
      }
      if(as.character(tz)=="2" ||  as.character(tz)=="200" || as.character(tz)=="0200" || as.character(tz)=="+0200"){
	      to.ret=200
      }
      if(as.character(tz)=="-2" ||  as.character(tz)=="-200" || as.character(tz)=="-0200"){
	      to.ret=-200
      }
      if(as.character(tz)=="3" ||  as.character(tz)=="300" || as.character(tz)=="0300" || as.character(tz)=="+0300"){
	      to.ret=300
      }
      if(as.character(tz)=="-3" ||  as.character(tz)=="-300" || as.character(tz)=="-0300"){
	      to.ret=-300
      }
      if(as.character(tz)=="4" ||  as.character(tz)=="400" || as.character(tz)=="0400" || as.character(tz)=="+0400"){
	      to.ret=400
      }
      if(as.character(tz)=="-4" ||  as.character(tz)=="-400" || as.character(tz)=="-0400"){
	      to.ret=-400
      }
      if(as.character(tz)=="5" ||  as.character(tz)=="500" || as.character(tz)=="0500" || as.character(tz)=="+0500"){
	      to.ret=500
      }
      if(as.character(tz)=="-5" ||  as.character(tz)=="-500" || as.character(tz)=="-0500"){
	      to.ret=-500
      }
      if(as.character(tz)=="6" ||  as.character(tz)=="600" || as.character(tz)=="0600" || as.character(tz)=="+0600"){
	      to.ret=600
      }
      if(as.character(tz)=="-6" ||  as.character(tz)=="-600" || as.character(tz)=="-0600"){
	      to.ret=-600
      }
      if(as.character(tz)=="7" ||  as.character(tz)=="700" || as.character(tz)=="0700" || as.character(tz)=="+0700"){
	      to.ret=700
      }
      if(as.character(tz)=="-7" ||  as.character(tz)=="-700" || as.character(tz)=="-0700"){
	      to.ret=-700
      }
      if(as.character(tz)=="8" ||  as.character(tz)=="800" || as.character(tz)=="0800" || as.character(tz)=="+0800"){
	      to.ret=800
      }
      if(as.character(tz)=="-8" ||  as.character(tz)=="-800" || as.character(tz)=="-0800"){
	      to.ret=-800
      }
      if(as.character(tz)=="9" ||  as.character(tz)=="900" || as.character(tz)=="0900" || as.character(tz)=="+0900"){
	      to.ret=900
      }
      if(as.character(tz)=="-9" ||  as.character(tz)=="-900" || as.character(tz)=="-0900"){
	      to.ret=-900
      }
      if(as.character(tz)=="10" ||  as.character(tz)=="1000" || as.character(tz)=="+1000"){
	      to.ret=1000
      }
      if(as.character(tz)=="-10" ||  as.character(tz)=="-1000"){
	      to.ret=-1000
      }
      if(as.character(tz)=="11" ||  as.character(tz)=="1100" || as.character(tz)=="+1100"){
	      to.ret=1100
      }
      if(as.character(tz)=="-11" ||  as.character(tz)=="-1100"){
	      to.ret=-1100
      }
      if(as.character(tz)=="12" ||  as.character(tz)=="1200" || as.character(tz)=="+1200" || as.character(tz)=="-12" ||  as.character(tz)=="-1200"){
	      to.ret=1200
      }
      if(as.character(tz)=="NULL"){
	      to.ret="+0"
	      warning("There was a TZ value of 'NULL', replaced by 'UTC'")
      }
	  	  #browser()
      if(is.na(to.ret)){
		stop(paste("unknown tz:", tz))
      }
      return(to.ret)

}
