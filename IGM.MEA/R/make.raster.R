# Diana Hall
# purpose: make raster plot with Robject complete with burst + ns data made from IGM.main
#  prompts user for input or the file path may be passed
make.raster<-function(RobjectFile=NULL,
                      outputdir=NULL ,
                      well.for.raster=NULL, 
                      interval.for.raster=NULL,
                      show.bursts=F, 
                      show.burst.number=F, 
                      show.networkspikes=F,
                      show.ns.number=F ){
  
  #data("parameters")
  # development
  # well.for.raster=NULL; RobjectFile=NULL; outputdir=NULL
  # show.ns.number=F; show.networkspikes=F; show.burst.number=F;
  # show.bursts=F; interval.for.raster=NULL; 
  
  
  # user input
  if ( is.null(RobjectFile) ){
    have.data=F; counter=1
    while(have.data==F){
      RobjectFile<-tk_choose.files(caption="Select R-object (.RData) for raster plot")
      data.type<-strsplit( basename(RobjectFile), split="[.]")[[1]][2] 
      
      #laod and check
      if ( !is.element( data.type,  c("rdata","RData","rData","Rdata") ) ){
        counter=counter+1
        tkmessageBox(message = ".RData file needed!", icon = "error", type = "ok")
      } else {
        t=load( RobjectFile, verbose=T )
        s=list(); s[[1]]<-get(t)
        have.data=T
        analysis<-list()
        analysis$output.dir<-dirname(RobjectFile) 
        analysis$Routput.dir<-dirname(RobjectFile)
      }
      if (counter==4){
        stop()
      } 
    }
  } else{
    analysis<-list()
    analysis$output.dir<-dirname(RobjectFile) 
    analysis$Routput.dir<-dirname(RobjectFile)
    t=load( RobjectFile, verbose=T )
    s=list(); s[[1]]<-get(t)
    have.data=T
  }
  
  if( is.null(outputdir) ){
    outputdir= analysis$output.dir 
  } else{
    dir.create(outputdir, showWarnings = FALSE)
  }
  
  setwd( outputdir )
  
  
  
  # error checks
  # check for correct data
  if ( !is.element("allb",names(s[[1]]) ) ){
    tkmessageBox(message = "No burst data in Robject!\n")
    stop("No burst data in Robject!\n")
  } 
  if ( !is.element("ns.all",names(s[[1]]) ) ){ 
    
    tkmessageBox(message = "No network spike data in Robject!\n")
    stop("No network spike data in Robject!\n")
  }
  
  ## ready data for plot
  # correct spike times for recording start offset
  
  
  # initialize loop and keep raster up
  # want.new.raster<-T;show.bursts<-T;show.burst.number<-T;show.networkspikes<-T; show.ns.number<-T;

      # +++++++++++++++++++++well.for.raster
      if ( !exists("well.for.raster") || is.null( well.for.raster ) ){
        well.for.raster = unique(s[[1]]$cw)[1]
      }
      well.for.raster<-toupper(well.for.raster) #ensure user entered upper case well names
      # check if its among available channels
      if ( !is.element( well.for.raster, unique(s[[1]]$cw) ) ){
        
        print(c("Improper 'well.for.raster' parameter specification.",
                "Specify a well begining with A-F and ending with 1-8",
                "for example well.for.raster<-'A4' ",
                "Specified well may not appear in recording") )
        well.for.raster = unique(s[[1]]$cw)[1]
        print(paste( "Using well ", well.for.raster,sep="")  )
        
      }
      
      
      
      # +++++++++++++++++++++ show bursts
      if ( !exists("show.bursts") || is.na(show.bursts) ){ show.bursts=F }
      if( !is.element(show.bursts, c(T,F) ) ){
        show.bursts=T
      }
      
      # +++++++++++++++++++++ show network spikes
      if ( !exists("show.networkspikes") || is.na(show.networkspikes) ){ 
        show.networkspikes=F }
      if(!is.element(show.networkspikes, c(T,F) ) ){
        show.networkspikes=F
      }
      
      # +++++++++++++++++++++ show bursts number
      if ( !exists("show.burst.number") || is.na(show.burst.number) ){ 
        show.burst.number=F }
      if(!is.element(show.burst.number, c(T,F) ) ){
        show.burst.number=F
      }
      
      
      # ++++++++++++++++++++++++++interval check
      if (is.null(interval.for.raster)||length(interval.for.raster)<2 ) {
        if ( is.element('rec.time', names(s[[1]])) ){
          interval.for.raster<-c(0, floor(s[[1]]$rec.time[2]-s[[1]]$rec.time[1]) )
        } else {
          interval.for.raster<-c(min(unlist( lapply(s[[1]]$spikes, min) )),
                                 max(unlist( lapply(s[[1]]$spikes, max) )) )
        }
         
      }
      
      if ( any( is.na(interval.for.raster) ) ){
        if( is.na(interval.for.raster[1]) ){
          interval.for.raster[1]<-0
        }
        if( is.na(interval.for.raster[2]) ){
          interval.for.raster[2]<-floor(s[[1]]$rec.time[2]-s[[1]]$rec.time[1])
        }
      }
      # error check on times chosen
      if ( 0>interval.for.raster[1]  ){
        interval.for.raster[1]<- 0
        print( paste("Beginning of raster interval preceeds recording start",
                     "resetting start of raster interval to start of recording", 
                     sep="\n" )  )
        
      }
      
      if ( s[[1]]$rec.time[2]<interval.for.raster[2] ){
        interval.for.raster[2]<-floor(s[[1]]$rec.time[2])
        
        print( paste("End of raster interval exceeds recording end",
                     "resetting end of raster interval to end of recoding", 
                     sep="\n" ) )
        
      }
      if ( interval.for.raster[2]<interval.for.raster[1] ){
        interval.for.raster<-c( ceiling(s[[1]]$rec.time[1]), floor( s[[1]]$rec.time[2])  )
        
      }
      
      # +++++++++++++++++++++ show.ns.number
      if ( !exists("show.ns.number") || is.null(show.ns.number) ){ 
        show.ns.number=F }
      if(!is.element(show.ns.number, c(T,F) ) ){
        show.ns.number=F
      }
      
      
      # each spike has a name, e.g. "A2_321", "A2_322" that's 1st and 2nd spike of channel A2_32
      if ( show.networkspikes ){
        if ( !( s[[1]]$ns.all[[well.for.raster]]$brief[1]==0 || is.na( s[[1]]$ns.all[[well.for.raster]]$brief[1]) ) ){
          raster.ns.t<-s[[1]]$ns.all[[well.for.raster]]$measures[ ,c("time","durn","peak.val") ]
          if (is.matrix(raster.ns.t) ){
            raster.ns.t[,"time"]=raster.ns.t[,"time"]
            index.want<-which(raster.ns.t[,"time"]<interval.for.raster[2]& raster.ns.t[,"time"]>interval.for.raster[1])
            if( length(index.want)>0){
              raster.ns<-raster.ns.t[index.want, ]
            } else{
              raster.ns=NULL
            }
            
          } else{
            raster.ns.t["time"]=raster.ns.t["time"]
            index.want<-which(raster.ns.t["time"]<interval.for.raster[2]& raster.ns.t["time"]>interval.for.raster[1])
            if ( length(index.want)>0  ){
              raster.ns<-raster.ns.t
            } else{
              raster.ns=NULL
            }
            
          }
          
        } else{
          raster.ns=NULL
        }
      } else{
        raster.ns=NULL
      }
      
      #######################################################
      
      
      
      ss<-summary(s[[1]])
      
      #make plot title
      plot.title.first.line<-paste( unlist(strsplit(basename(s[[1]]$file),split="_"))[1:4], 
                                    collapse="_" )
      #ns plot title line
      if (!is.null(raster.ns) && show.networkspikes ){
        plot.title.ns.line<-paste("green = network spike, # Elect at peak " ,sep="" ) 
      } else if ( is.null(raster.ns)  && show.networkspikes ){
        plot.title.ns.line<-paste(  "no network spikes" ) 
      } else {
        plot.title.ns.line<-paste(  " " ) 
      }
      # burst title line
      if ( show.bursts && show.burst.number ){
        plot.title.b.line<-paste("red (horz)= burst,  ", 
                                 "blue=# spikes/burst",sep="" ) 
      } else if ( show.bursts && !show.burst.number ){
        plot.title.b.line<-paste("red= bursts ",sep="" )  
      } else {
        plot.title.b.line<-paste(  " " ) 
      }
      plot.title<-paste(plot.title.first.line,
                        plot.title.b.line,
                        plot.title.ns.line,
                        sep="\n" )
      
      
      ####plot interesting responses for all files
      RasterPlotPath = paste(analysis$output.dir, "/rasterPlots.pdf", sep="" )
      pdf(file= RasterPlotPath ) 
      
      
      
      if (show.networkspikes ){
        if ( is.matrix(raster.ns) ){
          
          .plot.mm.s( s[[1]], beg = interval.for.raster[1], 
                     label.cells = T,
                     end = interval.for.raster[2] , 
                     show.bursts = show.bursts, 
                     whichcells = well.for.raster ,
                     show.burst.number=show.burst.number,
                     main= plot.title,
                     panel.first=for(cur.l in 1:length(raster.ns[,"time"]) ){
                       lines(xy.coords( c(raster.ns[cur.l,"time"], raster.ns[cur.l,"time"] ), 
                                        c(0,.98) ) , 
                             col="green") }
                     # panel.first=abline(v=raster.ns[,"time"], col="green")
          ) 
          # fix #Electrodes so no overlap exists
          if(show.ns.number){
            n.digits=sum(floor(raster.ns[,3]/10))+length(raster.ns[,3])
            close.neighbors=which(c(raster.ns[-1,1],tail(raster.ns[-1,1]+1.1, n=1))-raster.ns[,1]<1)+1
            for (i in 1:length(raster.ns[,1]) ){
              x.coord=raster.ns[i,"time"]
              if( is.element(i, close.neighbors)){
                deltay_t=(interval.for.raster[2]-interval.for.raster[1])/75
                x.coord=x.coord+deltay_t
              }
              
              print(x.coord)
              cur.peak.val=raster.ns[i,"peak.val"]
              print( cur.peak.val )
              text( x=x.coord, y=1.02  , 
                    labels =cur.peak.val , pos=NULL, col="green" )
            }#end of label for loop
          } # end of if(show.ns.number)
          mtext(text=summary(s[[1]]), outer=T, side=3 )  
        } else if( is.vector(raster.ns) ) {
          .plot.mm.s( s[[1]], beg = interval.for.raster[1], 
                     label.cells = T,
                     end = interval.for.raster[2] , 
                     show.bursts = show.bursts, 
                     whichcells = well.for.raster ,
                     show.burst.number=show.burst.number,
                     main= plot.title,
                     
                     panel.first=lines(xy.coords(c(raster.ns["time"], 
                                                   raster.ns["time"] ), c(0,.98) ) , 
                                       col="green")
          ) 
          if(show.ns.number){
            #use vector indexing      
            x.coord=raster.ns["time"]
            cur.peak.val=raster.ns["peak.val"]
            text( x=x.coord, y=1.02 , 
                  labels =cur.peak.val , pos=NULL, col="green" )
          }# end of show.ns.number
          mtext(text=summary(s[[1]]), outer=T, side=3 )  
        } else {
          # raster.ns is null
          .plot.mm.s( s[[1]], beg = interval.for.raster[1], 
                     label.cells = T,
                     end = interval.for.raster[2] , 
                     show.bursts = show.bursts, 
                     whichcells = well.for.raster ,
                     show.burst.number=show.burst.number,
                     main= plot.title
          )
          mtext(text=summary(s[[1]]), outer=T, side=3 )  
        }
        #end of if show.networkspikes
      } else{
        # if don't show ns
        .plot.mm.s( s[[1]], beg = interval.for.raster[1], 
                   label.cells = T,
                   end = interval.for.raster[2] , 
                   show.bursts = show.bursts, 
                   whichcells = well.for.raster ,
                   show.burst.number=show.burst.number,
                   main= plot.title
        )
        mtext(text=summary(s[[1]]), outer=T, side=3 )  
      }
      
      
      
      
      dev.off()
      # pop open file
      system(paste("open ", RasterPlotPath, sep="") )
      

  
  
} # end of function

