
iClick.VisAssetPrice <- function(DAT,color4="r2b",color5="jet") {
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252")

x=cbind(DAT[,2])
rownames(x)=as.character(DAT[,1])
y=timeSeries::as.timeSeries(x)


colnames(y)=c(names(DAT)[2])
head(y)
YMD=time(y)
yr=unique(lubridate::year(YMD))
charvec <- timeDate::timeCalendar(m = 12, d = 15:31, y = yr[1]-1, FinCenter = "GMT")
fake=timeSeries::as.timeSeries(rnorm(length(charvec)),charvec)
colnames(fake)=colnames(y)
x=rbind(fake,y);head(x,30)
names(x)=names(y)


YMD=time(x)
yr=unique(lubridate::year(YMD))
full.Date=as.Date(seq(from=YMD[1],to=YMD[length(YMD)],by="1 day"))
date=as.POSIXlt(paste(full.Date,"02:00:00"),"GMT")
full.Data=timeSeries::as.timeSeries(data.frame(date,1));

data.Date=as.POSIXlt(paste(YMD,"02:00:00"),"GMT")
real.Data=timeSeries::as.timeSeries(data.frame(data.Date,unclass(x)))

newData=cbind(full.Data,real.Data)[,-1]

dat=data.frame(date,unclass(newData))
colnames(dat)=c("date",names(DAT)[2])
    dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickCalendarPlot(obj.name = "plotType"))
    #Unit = colnames(x)

        # Print Basic Return Statistics:
        if (type == 1) {
        summaryTable=as.matrix(fBasics::basicStats(y)[-c(10:12),])
        rownames(summaryTable)=rownames(fBasics::basicStats(y))[-c(10:12)]
        print(summaryTable)
       }

        #=== Price Series Plot:
        if (type == 2) {
          dev.new();plot.new();
        seriesPlotX(y,ylab="Price", col = "indianred2")
        }

        #=== Cut and Connect
        if (type == 3) {
          dev.new();plot.new();
        print(cutAndStack(y, number=6, overlap = 0.1))
        }

        #=== Calender heatmap:
        if (type == 4) {
          dev.new();plot.new();
        print(calendarHeat(y,color=color4))
        }

        if (type == 5) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[2]),cols =color5,year=yr[2])
        }

        if (type == 6) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[3]),cols=color5,year=yr[3])
        }

        if (type == 7) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[4]),cols = color5,year=yr[4])
        }

        if (type == 8) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[5]),cols = color5,year=yr[4])
        }

        if (type == 9) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[6]),cols = color5,year=yr[6])
        }

        if (type == 10) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[7]),cols = color5,year=yr[7])
        }

        if (type == 11) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[8]),cols = color5,year=yr[8])
        }

        if (type == 12) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[9]),cols = color5,year=yr[9])
        }

        if (type == 13) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[10]),cols = color5,year=yr[10])
        }

        if (type == 14) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[11]),cols = color5,year=yr[11])
        }

        if (type == 15) {
          dev.new();plot.new();
          openair::calendarPlot(dat,pollutant=names(dat)[2],main=paste(names(dat),"in", yr[12]),cols = color5,year=yr[12])
        }

}  #End of dataRefreshCode()

    nAssets = dim(x)[2]
       V0=c("1 Descriptive Statistics Table",
            "2 Price Series Plot",
            "3 Breaking Plot",
            "4 Calender Heatmap, up to 6 years only")

        no=5:((length(yr)-1)+4)
        V1=paste(no," Calender Plot, year ", yr[-1],sep="")


    .oneClickCalendarPlot(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      nAssets),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "14")
                dataRefreshCode()},
        function(...){
                .oneClickCalendarPlot(obj.name = "plotType", obj.value = "15")
                dataRefreshCode()}
        ),



        button.names =c(V0,V1),

        title = "1-Click Visualization: Asset Price"
        )

  .oneClickCalendarPlot(obj.name = "type", obj.value = "1", no = 1)

   # Return Value()
   invisible()
}


.oneClickCalendarPlot.env = new.env()


.oneClickCalendarPlot <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {

    if(!exists(".oneClickCalendarPlot.env")) {
      .oneClickCalendarPlot.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickCalendarPlot.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickCalendarPlot.env)
      }
      return(obj.value)
    }
    if(missing(title)) {
      title = "Control Widget"
    }

    # GUI Settings:
    myPane <- tktoplevel()
    tkwm.title(myPane, title)
    tkwm.geometry(myPane, "+0+0")

    # Buttons:
    framed.button <- ttkframe(myPane,padding=c(3,3,12,12))
    tkpack(framed.button, fill = "x")

    if (missing(button.names)) {
      button.names <- NULL
    }

#looping button names
    for (i in seq(button.names)) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "45")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)

}


#===== Quit Button:
    quitCMD = function() {
    tkdestroy(myPane)

    }

   quitButton<-tkbutton(framed.button, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q", function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="indianred2", font=tkfont.create(weight="bold",size=10))

   tkconfigure(quitButton,underline=0)
   tkpack(quitButton, side = "right",fill = "x",ipady=3)


assign(".oneClickCalendarPlot.values.old", starts, envir = .oneClickCalendarPlot.env)

    # Return Value:
   invisible(myPane)
  }
