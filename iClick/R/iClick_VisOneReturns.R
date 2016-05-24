iClick.VisOneReturns <- function(dat) {

if (class(dat)=="ts"){
  x=timeSeries::as.timeSeries(dat)
}
else if (ncol(dat)==2) {
y=cbind(dat[,2])
rownames(y)=as.character(dat[,1])
x=timeSeries::as.timeSeries(y)
colnames(x)=c(names(dat)[2])
}

    N = ceiling(sqrt(ncol(x)))

    if (max(x)>=1) { x=x*0.01  }

        dataRefreshCode <- function(...)  {
    type = as.integer(.oneClickVisUniReturns(obj.name = "plotType"))
    Unit = colnames(x)

        # Print Basic Return Statistics:
        if (type == 1) {
        summaryTable=as.matrix(fBasics::basicStats(x)[-c(10:12),])
        rownames(summaryTable)=rownames(fBasics::basicStats(x))[-c(10:12)]
        print(summaryTable)
       }

        # Return Series Plot:
        if (type == 2) {
          dev.new();plot.new();
        seriesPlotX(x,ylab="Returns", col = "indianred2") }

        # Cumulate Return Series Plot
        if (type == 3) {
        dev.new()
         seriesPlotX(100*exp(timeSeries::colCumsums(x)),ylab="Cumulated Returns", col = "indianred2")

         abline(h = 100, col = "steelblue")
        }

        # Up-down Plot:
        if (type == 4) {
        dev.new()
plot.new();
        par(mfrow = c(2,1))
        drawupPlotX(x,ylab="Up Returns",col="indianred2")
        drawdownPlotX(x,ylab="Down Returns",col="darkgreen")
        par(mfrow = c(1,1))
        }

        # Histogram Plot:
        if (type == 5) {
          dev.new();
          fBasics::histPlot(x, skip = TRUE,col="indianred2") }

        # Density Plot:
        if (type == 6) {
          dev.new();
          fBasics::densityPlot(x,col="indianred2") }

        # Normal QQ Plot:
        if (type == 7) {
          dev.new();
          qqnormPlotX(x,col="indianred2") }

        # Box-Whisker Plot:
        if (type == 8) {
          dev.new();
          boxPlotX(x,col="indianred2") }

        # ACF/PACF Plot:
        if (type == 9) {
          dev.new();
        par(mfrow = c(2,1))
         fBasics::acfPlot(x)
         fBasics::pacfPlot(x)
        par(mfrow = c(1,1))
        }

        # lagged ACF Plot:
       if (type == 10) {
         dev.new();
         fBasics::lacfPlot(x) }

        # Taylor effect Plot:
       if (type == 11) {
         dev.new();
         fBasics::teffectPlot(x)}

       # Put them together: time series plots
       if (type == 12) {
         dev.new();
        par(mfrow = c(2,2))
        seriesPlotX(x,ylab="Returns", col = "indianred2")
        seriesPlotX(100*exp(timeSeries::colCumsums(x)),ylab="Cumulated Returns", col = "indianred2");abline(h = 100, col = "steelblue")
        drawupPlotX(x,ylab="Up Returns", col = "indianred2")
        drawdownPlotX(x,ylab="Down Returns", col = "darkgreen")
        par(mfrow = c(1,1))

       }

       # Put them together: distribution plots
       if (type == 13) {
         dev.new();
         par(mfrow = c(2,2))
         fBasics::histPlot(x, skip = TRUE)
         fBasics::densityPlot(x)
         qqnormPlotX(x)
         boxPlotX(x)
         par(mfrow = c(1,1))
       }

       # Put them together: Auto-Correlation Plots
       if (type == 14) {
         dev.new();
        par(mfrow = c(2,2))
         fBasics::acfPlot(x)
         fBasics::pacfPlot(x)
         fBasics::lacfPlot(x)
         fBasics::teffectPlot(x)
        par(mfrow = c(1,1))
       }


}  #End of dataRefreshCode()

    nAssets = dim(x)[2]

    .oneClickVisUniReturns(
        dataRefreshCode,
        names       = c("Selected Asset"),
        minima      = c(      0),
        maxima      = c(      nAssets),
        resolutions = c(      1),
        starts      = c(      0),

        button.functions = list(
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "1")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "2")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "3")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "4")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "5")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "6")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "7")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "8")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "9")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "10")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "11")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "12")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "13")
                dataRefreshCode()},
        function(...){
                .oneClickVisUniReturns(obj.name = "plotType", obj.value = "14")
                dataRefreshCode()}
        ),

        button.names = c(
            " 1 Descriptive Statistics Table",
            " 2 Return Series Plot",
            " 3 Cumulated Return Series Plot",
            " 4 Draw Up/Down Plots",
            " 5 Histogram Plot",
            " 6 Density Plot",
            " 7 Q-Q Plot for Normality",
            " 8 Box-Whisker Plot",
            " 9 ACF/PACF Plots",
            "10 Lagged ACF Plot",
            "11 Taylor Effects Plot",
            "12 Put them together: Time series plots",
            "13 Put them together: Distribution plots",
            "14 Put them together: Auto-Correlation plots"),

        title = "1-Click Visualization: Single Asset Returns"
        )

  .oneClickVisUniReturns(obj.name = "type", obj.value = "1", no = 1)

   # Return Value()
   invisible()
}


.oneClickVisUniReturns.env = new.env()


.oneClickVisUniReturns <-
  function(names, minima, maxima, resolutions, starts,button.functions, button.names, no, set.no.value, obj.name, obj.value,reset.function, title)
  {

#    if (!require(tcltk, quietly = TRUE))
#      stop("\n -- Package tcltk not available -- \n\n")

    if(!exists(".oneClickVisUniReturns.env")) {
      .oneClickVisUniReturns.env <<- new.env()
    }
    if(!missing(obj.name)){
      if(!missing(obj.value)) {
        assign(obj.name, obj.value, envir = .oneClickVisUniReturns.env)
      } else {
        obj.value <- get(obj.name, envir = .oneClickVisUniReturns.env)
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

#    framed.button <- scrollable_frame(framed.button0,350,450)
#    tkpack(framed.button0, fill = "x")
    tkpack(framed.button, fill = "x")

#loop through button names
    for (i in seq(button.names)) {
      button.fun <-button.functions[[i]]
      plotButtons<-tkbutton(framed.button, text = button.names[i], command = button.fun, anchor = "nw",relief="ridge",width = "45")
      tkconfigure(plotButtons,foreground="blue",font=tkfont.create(size=10,weight="bold"))
      tkpack(plotButtons,fill = "x", pady=1)
}

  #===== Quit Button:
    quitCMD = function() {tkdestroy(myPane)}

   quitButton<-tkbutton(framed.button, text = "Quit", command = quitCMD, anchor = "center",relief="ridge",width = "8")
   tkbind(myPane,"Q",function() tcl(quitButton,"invoke"))
   tkfocus(quitButton)
   tkconfigure(quitButton,foreground="indianred2", font=tkfont.create(weight="bold",size=10))

   tkconfigure(quitButton,underline=0)
   tkpack(quitButton, side = "right",fill = "x",ipady=3)

#    assign(".oneClickVisUniReturns.values.old", starts, envir = .oneClickVisUniReturns.env)

    # Return Value:
    invisible(myPane)
  }

