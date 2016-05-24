read1segy<-function(fname, Iendian=1 , HEADONLY=FALSE, BIGLONG=FALSE)
  {

    #####   given a file name with one PASSCAL SEGY, read it in and return a list
    
 if(missing(Iendian)) { Iendian=1 }
  if(missing(HEADONLY)) {HEADONLY=FALSE }
  if(missing(BIGLONG)) { BIGLONG=FALSE}

 DATAendian =c("little", "big", "swap")
  Kendian =c(1,2,3)

 

if(is.character(Iendian))
  {
    endianVEC =Iendian
  }
  else
    {
      endianVEC = DATAendian[match(Iendian , Kendian )]
    }


theENDIAN =  endianVEC
 ## 
SEGYhead.names = c("lineSeq", "reelSeq", "event_number", "channel_number",
        "energySourcePt", "cdpEns", "traceInEnsemble", "traceID",
        "vertSum", "horSum", "dataUse", "sourceToRecDist",
        "recElevation", "sourceSurfaceElevation", "sourceDepth",
        "datumElevRec", "datumElevSource",
        "sourceWaterDepth","recWaterDepth",
        
        "elevationScale", "coordScale",

        "sourceLongOrX", "sourceLatOrY","recLongOrX", "recLatOrY",
        
        "coordUnits", "weatheringVelocity", "subWeatheringVelocity",
        "sourceUpholeTime", "recUpholeTime", "sourceStaticCor",
        "recStaticCor", "totalStatic", "lagTimeA", "lagTimeB",
        "delay", "muteStart", "muteEnd",
        
        "sampleLength", "deltaSample", "gainType", "gainConst",
        "initialGain",
        
        "correlated", "sweepStart", "sweepEnd", "sweepLength",
        "sweepType", "sweepTaperAtStart", "sweepTaperAtEnd",
        "taperType", "aliasFreq", "aliasSlope",
        "notchFreq","notchSlope", "lowCutFreq", "hiCutFreq",
        "lowCutSlope", "hiCutSlope",
        
        "year", "day", "hour", "minute", "second", "timeBasisCode",
        
        "traceWeightingFactor", "phoneRollPos1", "phoneFirstTrace",
        "phoneLastTrace", "gapSize", "taperOvertravel",
        
        "station_name", "sensor_serial", "channel_name",
        "totalStaticHi", "samp_rate", "data_form", "m_secs",
        "trigyear", "trigday", "trighour", "trigminute", "trigsecond",
        "trigmills", "scale_fac", "inst_no", "not_to_be_used",
        "num_samps", "max", "min")

formsegy=c("long", "long", "long", "long", "long", "long", "long", "short", "short", "short", "short",
  "long", "long", "long", "long", "long", "long", "long", "long", "short", "short", "long", "long",
  "long", "long", "short", "short", "short", "short", "short", "short", "short", "short", "short",
  "short", "short", "short", "short", "short", "short", "short", "short", "short", "short", "short",
  "short", "short", "short", "short", "short", "short", "short", "short", "short", "short", "short",
  "short", "short", "short", "short", "short", "short", "short", "short", "short", "short", "short",
  "short", "short", "short", "short", "char", "char", "char", "short", "long", "short", "short",
  "short", "short", "short", "short", "short", "short", "float", "short", "short", "long", "long", "long")

charlen = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  6, 8, 4, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

segynotes=c("0 Sequence numbers within line",
" 4 Sequence numbers within reel ",
"   8 Original field record number or trigger number ",
"  12 Trace channel number within the original field record ",
"  16 X ",
"  20 X ",
"  24 X ",
"  28 Trace identification code: seismic data = 1 ",
"  30 X ",
"  32 X ",
"  34 X ",
"  36 X ",
"  40 X ",
"  44 X ",
"  48 X ",
"  52 X ",
"  56 X ",
"  60 X ",
"  64 X ",
"  68 Elevation Scaler: scale = 1 ",
"  70 Coordinate Scaler: scale = 1 ",
"  72 X ",
"  76 X ",
"  80 X ",
"  84 X ",
"  88 Coordinate Units:  = 2 (Lat/Long) ",
"  90 X ",
"  92 X ",
"  94 X ",
"  96 X ",
"  98 X ",
" 100 X ",
" 102 Total Static in MILLISECS added to Trace Start Time (lower 2 bytes)",
" 104 X ",
" 106 X ",
" 108 X ",
" 110 X ",
" 112 X ",
" 114 Number of samples in this trace (unless == 32767) ",
" 116 Sampling interval in MICROSECONDS (unless == 1) ",
" 118 Gain Type: 1 = Fixed Gain ",
" 120 Gain of amplifier ",
" 122 X ",
" 124 X ",
" 126 X ",
" 128 X ",
" 130 X ",
" 132 X ",
" 134 X ",
" 136 X ",
" 138 X ",
" 140 X ",
" 142 X ",
" 144 X ",
" 146 X ",
" 148 X ",
" 150 X ",
" 152 X ",
" 154 X ",
" 156 year of Start of trace ",
" 158 day of year at Start of trace ",
" 160 hour of day at Start of trace ",
" 162 minute of hour at Start of trace ",
" 164 second of minute at Start of trace ",
" 166 Time basis code: 2 = GMT ",
" 168 X ",
" 170 X ",
" 172 X ",
" 174 X ",
" 176 X ",
" 178 X ",
" 180 Station Name code (5 chars + \\0) ",
" 186 Sensor Serial code (7 chars + \\0) ",
" 194 Channel Name code (3 chars + \\0) ",
" 198 Total Static in MILLISECS added to Trace Start Time (high 2 bytes)",
" 200 Sample interval in MICROSECS as a 32 bit integer ",
" 204 Data Format flag: 0=16 bit, 1=32 bit integer ",
" 206 MILLISECONDS of seconds of Start of trace ",
" 208 year of Trigger time ",
" 210 day of year at Trigger time ",
" 212 hour of day at Trigger time ",
" 214 minute of hour at Trigger time ",
" 216 second of minute at Trigger time ",
" 218 MILLISECONDS of seconds of Trigger time ",
" 220 Scale Factor (IEEE 32 bit float) ",
" 224 Instrument Serial Number ",
" 226 X ",
" 228 Number of Samples as a 32 bit integer",
" 232 Maximum value in Counts ",
" 236 Minimum value in Counts ")

HEAD = vector(mode="list")
 
 if(BIGLONG)
    {
      
    ishort = 2
  iint  = 4
  ilong = 8
  ifloat = 4
  idouble = 8

  }
  else
    {

    ishort = 2
  iint  = 4
  ilong = 4
  ifloat = 4
  idouble = 8
  }

sizes = c(ilong, ishort , iint, ifloat)

m1 = match( formsegy, c("long",  "short", "char",  "float") )
isize = rep(iint, length(formsegy))
isize = sizes[m1]

##  cat(paste("In read1segy:", fname), sep="\n")
 
zz <- file(fname , "rb")

for(i in 1:length(formsegy) )
  {
##    print(paste(i,SEGYhead.names[i], formsegy[i]) )
       if(formsegy[i]=="char")
         {
           ## fchar =  paste(format(SEGYhead.vals[i] , width=charlen[i]-1), "\\0" , sep="" )
           
           fchar=readChar(zz, charlen[i], useBytes = FALSE)
           HEAD[[i]] = fchar
         }
       else
           {


               if(formsegy[i]!='float')
                   {
                       A1 =  readBin(zz, integer() , n = 1, size = isize[i] , signed = TRUE,
                           endian = theENDIAN)
                       
                       
                   }
               else
                   {
                       
                      ##   print(i, isize[i] )
                       A1 =  readBin(zz, numeric() , n = 1, size = isize[i] , signed = TRUE,
                           endian = theENDIAN)
                   }
               
               
               HEAD[[i]] = A1
               
           }
   }
 
 names(HEAD) = SEGYhead.names
 
 segyLIST=list(HEAD=HEAD, amp=NULL)
 
 if(HEADONLY==TRUE)
   {
     close(zz)
     invisible(segyLIST)
   }
 else
   {

 
N = HEAD[[ which(SEGYhead.names== "num_samps")  ]]
D1 = readBin(zz, integer() , n = N , size =iint ,  endian = theENDIAN, signed = TRUE)

 segyLIST$amp = D1

    close(zz)

invisible(segyLIST)
}
}


