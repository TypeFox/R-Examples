ARSTAN <-
function (){
  n.sink<-sink.number();  while(n.sink>0) {sink(); n.sink=n.sink-1}
  f.apply (FUN=arstan)
}


arstan <-
function (a=NULL){
if(is.null(a)) {
a<-tk_choose.files(multi=FALSE)
if(length(a)!=1) return(tk_messageBox(type="ok", "You should select at least one file!", caption="Problems")->a)
}

file.name<-(getFileName(a))
folder.work<-substr(a, 1,nchar(a)- nchar(basename(a)))
Output<-paste(folder.work,  sep="\\")
outputFolder<-paste(file.name,  sep="\\")
output=paste(file.name, "out", sep=".")
setwd(folder.work)
dir.create(file.name, showWarnings = FALSE)
setwd(file.name)
sink(output, append=FALSE)
Detrender.Options(Stc=stc)
sink()
sink(output, append=TRUE)

cat("Input: ", a, "\nOutput: ", Output, outputFolder, "/", output, "\n",sep ="")

#rwl<-readRwl(a,n.header=n.lines.header, info=FALSE)      #read the file [i]
rwl<-readRwl(a,n.header=NULL, info=FALSE)      #read the file [i]
TreeIds(rwl, stc=stc)                                                           # Give the number of trees and the no of cores per tree

TrwLessThan(rwl, TRW=0)                                                         # MISSING RINGS
RwlInfo(rwl, print=TRUE)                                                        # Rwl information

if (makeSegPlot){                                                               # Make graph with series length
jpeg(paste(file.name, "SegPlot.jpg", sep="-"), width = 700, height = 700, quality = 85 )
segPlot(rwl, nS=45,  main=file.name)                                           #main -> graph title #
dev.off()}

if (remove.shorter.series){                                                     # Delete series shorter than
rwl = DeleteSeriesShorterThan(rwl, filename=file.name, YEAR=delete.shorter.series)
}

# Plot RAW RWL
jpeg(paste(file.name, "RAW.jpg", sep="-"), width=900, height=700, quality=85)   #graph to JPEG file#
plotRwl(rwl,file.name=file.name, save.csv=T)                                   #main -> graph title #
dev.off()

################################################################################
#  DETRENDING  #################################################################
################################################################################
if (makeFirstDetrending) {


SaveDetrendJPG = function (rwl, detrend, folderName = "Detrend", work.dir = NULL, 
    detrend.method = "", select.series = 1:(ncol(rwl))) 
{


if (length(detrend.method)==1) detrend.method<-rep(detrend.method,ncol(rwl))
    if (is.null(work.dir)) 
        work.dir = getwd()
    dir.create(folderName, showWarnings = FALSE)
    setwd(folderName)
            seriesnames = colnames(rwl)
        yr.vec = as.numeric(rownames(rwl))
    for (i in select.series) {

        jpeg(paste(seriesnames[i], ".jpg", sep = ""), width = 1200, 
            height = 600, quality = 100)
        plot(yr.vec, rwl[, i], type = "l", xlab = "Years", ylab = "", 
            main = seriesnames[i], las = 1, col = "blue")
        mtext(paste(detrend.method[i], sep = ""), line = 0.5, side = 3, 
            adj = 1, cex = 0.9, col = "blue", font = 1)
        mtext("Detrender", line = 0.5, side = 3, adj = 0, cex = 0.9, 
            col = "blue", font = 1)
        lines(yr.vec, detrend[, i], col = 2)
        dev.off()
    }
    setwd(work.dir)
}





        cat("\nDETRENDING\n\nFirst detrending [", first.detrending.method, 
            "]\n", sep = "")
            
            
       # eval(parse(text = paste("assign('detrend1', apply(rwl, 2, detrendeR:::RemoveTrend, method=method1, BandwidthPerc=nPerc1, Bandwidth = n1, P=p1))", 
       #    sep = "")), envir = .GlobalEnv)
			
        detrend1 = apply(rwl, 2, detrendeR:::RemoveTrend, method = method1, 
            BandwidthPerc = nPerc1, Bandwidth = n1, P = p1)
      #  if (save.detrend1) 
      #      saveDetrendJPG(rwl, detrend1, folderName = "FirstDetrend", 
      #         detrend.method = first.detrending.method)
        if (interactive.detrend) {
           # eval(parse(text = paste("detrendeR:::InteractiveDetrending(rwl, detrend1, , folderName = 'FirstDetrend', method=method1,  n=n1, nPerc=nPerc1, p=p1)", 
            #    sep = "")), envir = .GlobalEnv)
				
          detrendeR:::InteractiveDetrending(rwl, detrend1, folderName = "FirstDetrend", 
                method = method1, n = n1, nPerc = nPerc1, p = p1)
         
            if (detrendeR:::.get("DETRENDING_INTERACTIVE_FLAG")) {
            #   eval(parse(text = paste("detrend1<<-detrendeR:::.get('DETRENDING_INTERACTIVE_OUTPUT')", 
            #      sep = "")), envir = .GlobalEnv)
        
                detrend1 = detrendeR:::.get("DETRENDING_INTERACTIVE_OUTPUT")
  								
                first.detrending.method <<- as.vector(detrendeR:::.get("DETRENDING_INTERACTIVE_OUTPUT_CHANGE")[-1,3])
                first.detrending.method
                
            }
           # else {
           #     if (save.detrend1) 
           #       saveDetrendJPG(rwl, detrend1, folderName = "FirstDetrend", 
           #         detrend.method = first.detrending.method)
           # }
        }

		 if (save.detrend1) {
		 
                  SaveDetrendJPG(rwl, detrend1, folderName = "FirstDetrend", 
                    detrend.method = first.detrending.method)
		}
        if (min(detrend1, na.rm = TRUE) < 0) {
            cat("\n")
            TrwLessThan(detrend1, TRW = 0.01)
        }
        rw1 = rwl/detrend1
        write.rwl(as.data.frame(detrend1), fname = paste(file.name, "cv1", sep = "."), 
            long.names = TRUE)
        write.rwl(rw1, fname = paste(file.name, "in1", sep = "."), 
            long.names = TRUE)
        RwlInfo(rw1, print = TRUE)
		
        if (makeSecondDetrending) {
            cat("\nSecond detrending [", second.detrending.method, 
                "]\n", sep = "")
            #eval(parse(text = paste("assign('detrend2', apply(rw1, 2, detrendeR:::RemoveTrend, method=method2, BandwidthPerc=nPerc2, Bandwidth = n2, P=p2))", 
            #sep = "")), envir = .GlobalEnv)
            
			detrend2 = apply(rw1, 2, detrendeR:::RemoveTrend, method = method2, 
                BandwidthPerc = nPerc2, Bandwidth = n2, P = p2)
          #  if (save.detrend2) 
          #      saveDetrendJPG(rw1, detrend2, folderName = "SecondDetrend", 
          #        detrend.method = second.detrending.method)
            if (interactive.detrend) {
 #eval(parse(text = paste("detrendeR:::InteractiveDetrending(rw1, detrend2, , folderName = 'SeconDetrend', method=method2,  n=n2, nPerc=nPerc2, p=p2)", 
 #               sep = "")), envir = .GlobalEnv)
				
      
             detrendeR:::InteractiveDetrending(rw1, detrend2, folderName = "SecondDetrend", 
                  method = method2, n = n2, nPerc = nPerc2, p = p1)
                if (detrendeR:::.get("DETRENDING_INTERACTIVE_FLAG")) {
                  detrend2 = detrendeR:::.get("DETRENDING_INTERACTIVE_OUTPUT")
                  second.detrending.method = as.vector(detrendeR:::.get("DETRENDING_INTERACTIVE_OUTPUT_CHANGE")[-1,3])
                }
          #      else {
          #        if (save.detrend2) 
          #          saveDetrendJPG(rw1, detrend2, folderName = "SecondDetrend", 
          #            detrend.method = second.detrending.method)
          #      }
            }
        if (save.detrend2) 
                     SaveDetrendJPG(rw1, detrend2, folderName = "SecondDetrend", 
                      detrend.method = second.detrending.method)
            if (min(detrend2, na.rm = TRUE) < 0) {
                cat("\n")
                TrwLessThan(detrend2, TRW = 0.01)
            }
            rw2 = rw1/detrend2
            write.rwl(as.data.frame(detrend2), fname = paste(file.name, "cv2", 
                sep = "."), long.names = TRUE)
            write.rwl(rw2, fname = paste(file.name, "in2", sep = "."), 
                long.names = TRUE)
            RwlInfo(rw2, print = TRUE)
        } else {
            rw2 <- rw1
        }
    } else {
        rw2 <- rwl
    }

# Save standard chronology
crnStd <- Chron(rw2, prefix = paste(file.name, "STD", sep="-"), biweight = .get("biweightMean"), prewhiten = FALSE,
                  stc=stc)
PrintPlotChrono(crnStd, rwl=rw2, file.name = file.name, crono.type="STD" )

if (run.win.analysis){cat("\nRunning analysis of detrended series\n\n")
runWinRw2<-Run.Win(rw2, winLength=winLength, stc=stc, step = stepWin)
for (i in c(0.75, 0.80, 0.85, 0.90)) {EPS.resume(runWinRw2, EPS=i)}
}

if (make.common.EPS){EPS.common.interval(rw2, stc=stc, out=FALSE)}
if (make.select.period.EPS) { if (!(is.null(first.year.common)||is.null(last.year.common))){
EPS.common.interval(rw2, first.year.common=first.year.common,
                    last.year.common=last.year.common, stc=stc, out=FALSE) } }

if (makeAr){
            ArFunction(rw2, order.max = arMAX)                                  #AUTOREGRESSIVE MODELING
            res = apply(rw2, 2, ar.func, order.max = arMAX)
            crnRes <- Chron(res, prefix = paste(file.name, "RES", sep="-"), biweight = .get("biweightMean"), prewhiten = FALSE, stc=stc)
            PrintPlotChrono(crnRes, rwl=res, file.name = file.name, crono.type="RES" )
            cat("\nResidual series\n")
            RwlInfo(res, print=TRUE)
            write.rwl(as.data.frame(res), fname= paste(file.name, "res", sep="."),long.names=TRUE)
  if (run.win.analysis){cat("\nRunning analysis of residual series\n\n")
            runWinRes<-Run.Win(res, winLength=winLength, stc=stc, step = stepWin)
            for (i in c(0.75, 0.80, 0.85, 0.90)) {EPS.resume(runWinRes, EPS=i)}}
          
if (make.common.EPS){  EPS.common.interval(res, stc=stc, out=FALSE)  }
if (make.select.period.EPS) { if (!(is.null(first.year.common)||is.null(last.year.common))) {
              EPS.common.interval(res, first.year.common=first.year.common,
              last.year.common=last.year.common, stc=stc, out=FALSE) }  }
}
cat("[", date(),"]",sep="")
sink()
}

