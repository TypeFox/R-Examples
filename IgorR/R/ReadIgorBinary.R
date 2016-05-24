# ReadIgorBinary.R
# 
# Copyright (C) 2007 Gregory Jefferis <gsxej2@cam.ac.uk>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#################

# 2005-11-14
# File to read the igor binary file format into R
# There are various issues to deal with
# 1) Formats are either v2 or v5
# 2) Can be either Mac or PC

# Tested with the following from Igor Tech Note # 3 package
# zeroPointWave.ibw
# version5.ibw
# version2.ibw
# double.ibw
# along with my own testuint and testuword

# 2006-06-15
# Now reads the Igor packed experiment format into a list
# So far just works on the open/close data folders and reads 
# v2 and v5 multidimensional data.  Still need to check if matrices and
# higher dimensional arrays are laid out in the same order.
# Can supply a regex to read in only a selection and can also request that 
# only headers rather than data are loaded. 
# Have also tested it with one of Rachel's .pxp files and that worked fine 
# as well

# 2006-11-23
# Now reads text waves 
# tested with SamplesMac/textWave.ibw and SampleWin/textWave.ibw
# and some neuromatic pxp sweeps

# 2006-11-25
# Fixed problems reading waves with notes / extended units etc
# Now able to load a large Neuromatic pxp file (about 6 or 7Mb, 25s)
# Improved error checking for offsets of text waves
# Made return of v2 or v5 read more consistent

# 2009-07-09
# Fixed handling of dates in ReadNclampLogTable

#' Read binary files in the Igor Binary Wave format (IBW)
#' 
#' @param wavefile either a character vector containing the path to a file or an
#'   R \link{connection}.
#' @param Verbose if \code{TRUE}, print status information while reading the
#'   file.
#' @param ReturnTimeSeries if \code{TRUE}, return as an R time series (package
#'   \code{\link{ts}}).
#' @param MakeWave if \code{TRUE}, assign wave to a list in the global user
#'   environment.
#' @param HeaderOnly if \code{TRUE}, only return the header of the Igor Wave.
#' @return A vector containing the wave data or, if \code{MakeWave == TRUE},
#'   returns the name of a new R vector containing the data which has been made
#'   in the user environment
#' @author jefferis
#' @export
#' @examples
#' # return a list containing the wave
#' wavedata=read.ibw(system.file("igor","version5.ibw",package="IgorR")) 
#' sum(wavedata)
#' 
#' # make a list containing the wave's data in the users's environment
#' wavename=read.ibw(system.file("igor","version5.ibw",package="IgorR"),MakeWave=TRUE) 
#' sum(get(wavename))
read.ibw<-function(wavefile,Verbose=FALSE,ReturnTimeSeries=FALSE,
    MakeWave=FALSE,HeaderOnly=FALSE){
  if (is.character(wavefile)) {
    # NB setting the encoding to "macintosh" resolves some problems with utf-8 incompatible chars
    # in the mac or windows-1252 encodings
    wavefile <- file(wavefile, "rb",encoding="macintosh")
    on.exit(close(wavefile))
  }
  
  # read one byte
  bytes=readBin(wavefile,"integer",2,size=1)
  if(bytes[1]==0){
    endian="big"; version=bytes[2]
  } else {
    endian="little"; version=bytes[1]
  }
  if(Verbose) cat("version = ",version,"endian = ",endian,"\n")
  
  if(version==5) {
    rval=.ReadIgorBinary.V5(wavefile,Verbose=Verbose,endian=endian,HeaderOnly=HeaderOnly)
  } else if(version==2){
    rval=.ReadIgorBinary.V2(wavefile,Verbose=Verbose,endian=endian,HeaderOnly=HeaderOnly)
  }
  else stop(paste("Unable to read from Igor Binary File:",summary(wavefile)$description,"with version:",version))
  # Store Igor wave version number
  attr(rval,'BinHeader')$version=version
  
  if(ReturnTimeSeries){
    rval=WaveToTimeSeries(rval)
  }
  
  # makes a wave with a specified name in the user environment
  if(MakeWave){
    OriginalWaveName=attr(rval,"WaveHeader")$WaveName
    WaveName=make.names(OriginalWaveName)
    if(WaveName!=OriginalWaveName)
      warning("Had to change wave name from ",OriginalWaveName," to ",WaveName,"to keep R happy")
    assign(WaveName,rval,envir = .GlobalEnv)
    invisible(WaveName)
  } else invisible(rval)
}

#' Reads an Igor Pro Packed Experiment (.pxp) file
#' 
#' Note that PXP files are only partially documented so some contents cannot be 
#' parsed (e.g. image data). This function currently reads data records (Igor 
#' waves and variables), history, procedures, recreation macros and plain text 
#' notebooks. Formatted notebooks cannot be read.
#' 
#' \code{IgorPlatform} will determine in which encoding text is read 
#' (WINDOWS-1252 for windows and macintosh for macintosh). Unique abbreviations 
#' are acceptable. Defaults to \code{"windows"} on Windows, \code{"macintosh"}
#' otherwise. Note that Igor Pro 5.5 added a PlatformRecord to the PXP file
#' format which is used to determine the file's platform of origin when
#' available. Since this is information straight from the horse's mouth it will
#' override the \code{IgorPlatform} argument.
#' @param pxpfile character vector naming a PXP file or an R \link{connection}.
#' @param regex if \code{TRUE}, only read records (e.g. waves) in the PXP file
#'   whose names match a \link{regex}.
#' @param ReturnTimeSeries if \code{TRUE}, Igor waves are returned as a
#'   \code{link{ts}} object with sensible x scaling (\code{FALSE} by default).
#' @param Verbose whether to print information to console during loading 
#'   (numeric values are also allowed 0=none, 1=basic, 2=all).
#' @param StructureOnly (TODO) if \code{TRUE}, only the structure of the PXP
#'   file for inspection.
#' @param ExtractText whether to extract procedures, recreation macros, history 
#'   and plain text notebooks (\code{FALSE} by default).
#' @param IgorPlatform OS on which Igor file was saved (windows or macintosh).
#' @param ... optional parameters passed to \link{read.ibw}.
#' @return A list containing all the individual waves or variables in the PXP 
#'   file.
#' @author jefferis
#' @export
#' @import bitops
#' @examples 
#' r=read.pxp(system.file("igor","testexpt.pxp",package="IgorR"))
read.pxp<-function(pxpfile,regex,ReturnTimeSeries=FALSE,Verbose=FALSE,
    StructureOnly=FALSE,ExtractText=FALSE,IgorPlatform=NULL,...){
  if (is.character(pxpfile)) {
    # NB setting the encoding to "MAC" resolves some problems with utf-8 incompatible chars
    # in the mac or windows-1252 encodings
    pxpfile <- file(pxpfile, "rb",encoding="LATIN1")
    on.exit(close(pxpfile))
  }
  filename=summary(pxpfile)$description
  
  if(is.logical(Verbose))
    Verbose=ifelse(Verbose,1,0)

  # Check if this file needs to be byte swapped
  firstShort=readBin(pxpfile,"integer",1,size=2)
  firstShort=bitAnd(firstShort,0x7FFF)
  if (bitAnd(firstShort,0xFF00) != 0) endian="swap"
  else endian = .Platform$endian

  # default to windows on windows, mac otherwise
  if(is.null(IgorPlatform))
    IgorPlatform=ifelse(.Platform$OS.type=="windows",'windows','macintosh')
  else
    IgorPlatform=match.arg(IgorPlatform,choices=c('windows','macintosh'))

  encoding = ifelse(IgorPlatform=='windows','WINDOWS-1252','macintosh')

  root=list() # we will store data here
  currentNames="root"
  recordStartPos=0
  fileSize=file.info(filename)$size
  while(recordStartPos<fileSize){
    seek(pxpfile,recordStartPos)
    ph=.ReadPackedHeader(pxpfile,endian)
    if(Verbose>1) cat("recordStartPos =",recordStartPos,",Record type =",ph$recordType,
          "of length =",ph$numDataBytes,"\n")
    
    recordStartPos=seek(pxpfile)+ph$numDataBytes
    #if(!is.na(recTypeFuns[ph$recordType])) x=match.fun(recTypeFuns[ph$recordType])(pxpfile,endian)
    if (ph$recordType==1){
      # variables
      vh=.ReadVarHeader(pxpfile,endian)
      
      vars=list()
      if(vh$numSysVars>0) vars$sysVars=readBin(pxpfile,n=vh$numSysVars,size=4,what=numeric(),signed=TRUE,endian=endian)
      if(vh$numUserVars>0) vars=c(vars,.ReadUserVar(pxpfile,endian,n=vh$numUserVars))
      if(vh$numUserStrs>0) vars=c(vars,
            .ReadUserStr(pxpfile,endian,n=vh$numUserStrs,IgorPlatform=IgorPlatform),
            attr(vh,"version"))
      if(Verbose>1) print(vh)
      if(Verbose>1) print(vars)
      el=paste(paste(currentNames,collapse="$"),sep="$","vars")
      eval(.myparse(text=paste(el,"<-vars")))
    } else if (ph$recordType==3){
      # wave record; verbose made for wave reading if we are passed
      # Verbose = 2
      x=.ReadWaveRecord(pxpfile,endian,Verbose=ifelse(Verbose==2,TRUE,FALSE),
          HeaderOnly=StructureOnly,ReturnTimeSeries=ReturnTimeSeries,...)
      if(!is.null(x)){
        wavename=attr(x,"WaveHeader")$WaveName
        if(is.null(wavename)) {
          # assume this is a wave name
          wavename=x
        }
        el=paste(paste(currentNames,collapse="$"),sep="$",wavename)
        # store the record if required
        if(missing(regex) || any( grep(regex,el) )){
          eval(.myparse(text=paste(el,"<-x")))
        }
        if (Verbose>0) cat("el:",el,"\n")
      }
    } else if (ph$recordType == 2 && ExtractText){
        # Experiment History
        history = .readCharsWithEnc(pxpfile, ph$numDataBytes, encoding)
        el = paste(paste(currentNames, collapse="$"), sep="$", "history")
        eval(.myparse(text=paste(el,"<-history")))
        if(Verbose > 1)
          cat("history: ", substr(history,0,20), " ...\n")
    } else if (ph$recordType == 4 && ExtractText){
        # Recreation Macro
        recmacro = .readCharsWithEnc(pxpfile, ph$numDataBytes, encoding)
        el = paste(paste(currentNames, collapse="$"), sep="$", "recmacro")
        eval(.myparse(text=paste(el,"<-recmacro")))
        if(Verbose > 1)
          cat("recreation macro: ", substr(recmacro,0,20), " ...\n")
    } else if (ph$recordType == 5 && ExtractText){
        # Procedure Text
        mainproc = .readCharsWithEnc(pxpfile, ph$numDataBytes, encoding)
        el = paste(paste(currentNames, collapse="$"), sep="$", "mainproc")
        eval(.myparse(text=paste(el,"<-mainproc")))
        if(Verbose > 1)
          cat("main procedure: ", substr(mainproc,0,20), " ...\n")
    } else if (ph$recordType == 8 && ExtractText){
        # packed procedures and notebooks
        file = .ReadPackedFile(pxpfile, ph$numDataBytes, encoding,  Verbose)
        if(length(file) > 0){
          el = paste(paste(currentNames, collapse="$"), sep="$", file$name)
          eval(.myparse(text=paste(el,"<-file$data")))
          if(Verbose > 1)
            cat("packed file ", file$name, ": ", substr(file$data,1,20), " ...\n")
        }
    } else if (ph$recordType==9){
      # Open Data Folder
      foldername=.ReadDataFolderStartRecord(pxpfile,endian)
      # backquote protects funny folder names with no effect on valid names
      currentNames=c(currentNames, sprintf("`%s`", foldername))
    } else if (ph$recordType==10){
      # Close Data Folder
      currentNames=currentNames[-length(currentNames)]
    } else if (ph$recordType==20){
      # Platform info
      ipi=.ReadPlatformInfo(pxpfile,endian,ph$numDataBytes)
      if(ipi$platform==1) IgorPlatform='macintosh'
      else if(ipi$platform==2) IgorPlatform='windows'
    }
  }
  invisible(root)
}

.DetermineIgorWaveType<-function(WaveHeader,endian=NULL){
  if(endian=="big") WaveHeader$type=rev(WaveHeader$type)
  WaveHeader$typeBits=as.integer(rawToBits(WaveHeader$type[1]))
  
  # // From IgorMath.h
  #define NT_CMPLX 1      // Complex numbers.
  #define NT_FP32 2     // 32 bit fp numbers.
  #define NT_FP64 4     // 64 bit fp numbers.
  #define NT_I8 8       // 8 bit signed integer. Requires Igor Pro 2.0 or later.
  #define NT_I16  0x10    // 16 bit integer numbers. Requires Igor Pro 2.0 or later.
  #define NT_I32  0x20    // 32 bit integer numbers. Requires Igor Pro 2.0 or later.
  #define NT_UNSIGNED 0x40  // Makes above signed integers unsigned. Requires Igor Pro 3.0 or later.
  
  # default is character ie text wave, which I assume to be ascii
  what="character";size=1
  if(WaveHeader$typeBits[6]) { what="integer"; size=4 }
  if(WaveHeader$typeBits[5]) { what="integer"; size=2 }
  if(WaveHeader$typeBits[4]) { what="integer"; size=1 }
  if(WaveHeader$typeBits[3]) { what="double"; size=8 }
  if(WaveHeader$typeBits[2]) { what="double"; size=4 }
  if(WaveHeader$typeBits[1]) what="complex"
  WaveHeader$signAdjustment=0 # note that readBin handles unsigned data poorly,
  # so have to do it myself
  if(WaveHeader$typeBits[7]) { what="integer"; WaveHeader$signAdjustment=2^(size*8)}
  
  WaveHeader$what=what;
  WaveHeader$size=size;
  return(WaveHeader)
}

.ConvertIntToUInt<-function(x,adjustment=2^32){
  x=as.numeric(x)
  signs=sign(x)
  x[signs<0]=x[signs<0]+adjustment
  x
}
.readNullTermString<-function(con,totalLength,encoding,...){
  if(totalLength==0) return ("")
  # first read in exactly the number of raw bytes we've been told
  # note totalLength is the length including the null terminator (if present)
  rawc=readBin(con, what='raw', n=totalLength)
  # then read a string from that - readBin will drop any null terminator
  s=readBin(rawc, what="character")
  if(missing(encoding)) return(s)
  else return(iconv(s,from=encoding,to=""))
}

.readCharsWithEnc<-function(con,nchars,encoding,targetenc='UTF-8'){
  s=readChar(con,nchars,useBytes=TRUE)
  if(missing(encoding)) return(s)
  else return(iconv(s,from=encoding,to=targetenc))
}

.readIgorDate<-function(con,endian){
  # Igor uses date as seconds since midnight 1 Jan 1904
  # Unfortunately does not seem to record time zone, so have to use current
  # readBin can only read signed shorts or bytes,
  # so must add 2^32 to read as unsigned long 
  i=as.numeric(readBin(con,endian=endian,what=integer(),size=4))+2^32
  i=.convertIgorDate(i)
  i
}

igor_date_origin<-as.numeric(ISOdate(1904,1,1,hour=0,tz=""))

.convertIgorDate<-function(dateval){
  dateval=dateval+igor_date_origin
  class(dateval)<-"POSIXct"
  dateval
}

# enum PackedFileRecordType {
#   kUnusedRecord = 0,
#   kVariablesRecord,   //  1: Contains system numeric variables (e.g., K0) and user numeric and string variables.
#   kHistoryRecord,     //  2: Contains the experiment's history as plain text.
#   kWaveRecord,      //  3: Contains the data for a wave
#   kRecreationRecord,    //  4: Contains the experiment's recreation procedures as plain text.
#   kProcedureRecord,   //  5: Contains the experiment's main procedure window text as plain text.
#   kUnused2Record, //6
#   kGetHistoryRecord,    //  7x6: Not a real record but rather, a message to go back and read the history text.
#   kPackedFileRecord,    //  8x7: Contains the data for a procedure file or notebook in a packed form.
#   kDataFolderStartRecord, //  9x8: Marks the start of a new data folder.
#   kDataFolderEndRecord  // 10: Marks the end of a data folder.
#   kPlatformRecord=20    // 20: Info about Igor that wrote experiment. Added for Igor Pro 5.5D.
#   
#   /*  Igor writes other kinds of records in a packed experiment file, for storing
#     things like pictures, page setup records, and miscellaneous settings. The
#     format for these records is quite complex and is not described in PTN003.
#     If you are writing a program to read packed files, you must skip any record
#     with a record type that is not listed above.
#   */
# };

#' Private functions in IgorR Package
#' @name IgorR-private
NULL

#' @title Read the a short record header from the current location in a PXP file
#' 
#' Note that the recordType will be one of the constants from Igor's
#' enum PackedFileRecordType
#' @param con an R connection to the file we are reading
#' @param endian either little or big
#' @return a list containing information about the current record
#' @rdname IgorR-private
#' @author jefferis
.ReadPackedHeader<-function(con,endian){
  recordType=readBin(con,size=2,what=integer(),signed=FALSE,endian=endian)
  version=readBin(con,size=2,what=integer(),endian=endian)
  numDataBytes=readBin(con,size=4,what=integer(),endian=endian) # is 4 correct?
  return(list(recordType=recordType,version=version,numDataBytes=numDataBytes))
}

.ReadDataFolderStartRecord<-function(con,endian){
  folderName=.readNullTermString(con,32)
  #cat("Start of Data Folder",folderName,"\n")
  folderName
}
.ReadDataFolderEndRecord<-function(con,endian){
  #cat("End of Data Folder")
}
.ReadWaveRecord<-function(con,endian,Verbose=FALSE,...){  
  x=read.ibw(con,Verbose=Verbose,...)
  #cat("Wave Record", attr(x,"WaveHeader")$WaveName,"\n")
  x
}

.ReadPackedFile<-function(con, recordSize, encoding, Verbose){

  # discard first header part
  numBytes   = 32
  readChar(con, numBytes)
  recordSize = recordSize - numBytes

  # read the filename
  file       = list()
  file$name  = readChar(con, numBytes)
  recordSize = recordSize - numBytes

  # discard last header part
  numBytes   = 90
  readChar(con, numBytes)
  recordSize = recordSize - numBytes

  file$data  = .readCharsWithEnc(con, recordSize, encoding)

  if(nchar(file$data) <= 1) { # assume it is a formatted notebook with custom binary header

    if(Verbose > 1)
      cat("Ignoring formatted notebook ", file$name, "\n")

    return(list())
  }

  return(file)
}

#struct PlatformInfo {      // Data written for a record of type kPlatformRecord.
# short platform;       // 0=unspecified, 1=Macintosh, 2=Windows.
# short architecture;     // 0=invalid, 1=PowerPC, 2=Intel.
# double igorVersion;     // e.g., 5.00 for Igor Pro 5.00.
# char reserved[256 - 12];  // Reserved. Write as zero.
#};

.ReadPlatformInfo<-function(con,endian,size){
  # I didn't find this documented anywhere, but this record can be preceded
  # by some additional information besides the PlatformInfo struct
  readBin(con,n=size-256,size=1,what=raw())
  l=list()
  l$platform=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian)
  l$architecture=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian)
  l$igorVersion=readBin(con,size=8,what=numeric(),endian=endian)
#  print(l)
  # skip remaining bytes
  readBin(con,n=256-12,size=1,what=raw())
  l
}

.ReadVarHeader<-function(con,endian){
  l=list()
  vhversion=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian) # Version number is 2 for this header.
  l$numSysVars=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian)  # Number of system variables (K0, K1 . . .).
  l$numUserVars=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian) # Number of user numeric variables -- may be zero.
  l$numUserStrs=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian) # Number of user string variables -- may be zero.
  if(vhversion==2){
    l$numDependentVars==readBin(con,size=2,what=integer(),signed=TRUE,endian=endian) # Number of dependent numeric variables -- may be zero.
    l$numDependentStrs==readBin(con,size=2,what=integer(),signed=TRUE,endian=endian) # Number of dependent string variables -- may be zero.
  }
  attr(l,"version")=vhversion
  l
}
.ReadUserVar<-function(con,endian,n=1){
  # typedef struct UserNumVarRec {
  #   char name[31+1];          // Name of the variable as a C string including null terminator.
  #   short type;             // Must be 1, signifying numeric.
  #   VarNumRec num;            // Type and value of the variable if it is numeric.
  # } UserNumVarRec;
  l=list()
  if(n<0) {return(NULL)}
  for(i in 1:n){
    varname=.readNullTermString(con,32)
    vtype=readBin(con,size=2,what=integer(),signed=TRUE,endian=endian)
    if(vtype!=1) warning(paste("inappropriate vartype",vtype," for varname",varname))
    # struct VarNumRec {
    #   short numType;            // Type from IgorMath.h.
    #   double realPart;          // The real part of the number.
    #   double imagPart;          // The imag part if the number is complex.
    #   long reserved;            // Reserved - set to zero.
    # };
    # make a fake Igor Binary Wave Header and then use that to figure out numeric type
    fakewh=list(type=readBin(con,n=2,size=1,what="raw"))
    fakewh=.DetermineIgorWaveType(fakewh,endian)
    x=readBin(con,n=2,size=8,what=numeric(),endian=endian)
    if(fakewh$what=='complex'){
      x=complex(real=x[1],imaginary=x[2])
    } else{
      x=x[1]
    }
    readBin(con,n=1,size=4,what=integer(),endian=endian)
    l[[varname]]=x
  }
  l
}

.ReadUserStr<-function(con,endian,n=1,version=1,IgorPlatform){
  # typedef struct UserStrVarRec1 {
  #   char name[31+1];          // Name of the string variable.
  #   short strLen;           // The real size of the following array.
  #   char data[1];
  # } UserStrVarRec1;
  
  # typedef struct UserStrVarRec2 {
  #   char name[31+1];          // Name of the string variable.
  #   long strLen;            // Number of bytes in the string.
  #   char data[1];           // Start of string data. This is not a C string - no null terminator.
  # } UserStrVarRec2;         
  l=list()
  if(n<0) {return(NULL)}
  for(i in 1:n){
    #cat("seek=",seek(con),"\n")
    varname=.readNullTermString(con,32)
    #cat("seek=",seek(con),"\n")
    #print(varname)
    #print(version)
    strLen=readBin(con,n=1,size=version*2,what=integer(),endian=endian)
    #cat("strLen=",strLen,"\n")
    x=.readCharsWithEnc(con,strLen,encoding=ifelse(IgorPlatform=='windows','WINDOWS-1252','macintosh'))
    if(varname!=""){
      l[[varname]]=x
    }
  }
  l
}

# keep.source FALSE is a recent addition
if(R.version$major>2) {
  .myparse<-function(text) parse(text=text, keep.source=FALSE) 
} else {
  .myparse<-function(text) parse(text=text)
}

.ReadIgorBinary.V2<-function(con,Verbose=FALSE,ReturnTimeSeries=NULL,endian=NULL,HeaderOnly=FALSE){
  myread=function(what="integer",size=4,...) readBin(con,endian=endian,what=what,size=size,...)
  # File pointer should be positioned at the start of the header
  # this should be 2 bytes into the wave data 
  # (ie after endian and version bytes).  Store that position:
  startPos=seek(con)

  # Read binary header 
  BinHeader2=list()
  BinHeader2$wfmSize=myread()
  BinHeader2$noteSize=myread(); myread()
  BinHeader2$checksum=myread(size=2)

  if(Verbose) print("Hello!")
  if(Verbose) print(BinHeader2)
  if(Verbose) print(BinHeader2)
  # Read WaveHeader
  seek(con,where=startPos+16-2) # is this necessary?
  WaveHeader2=list()
  WaveHeader2$type=myread(what="raw",size=1,n=2)
  myread()
  # The platform default is utf-8, but some mac or windows-1252 chars are invalid in this
  # encoding (e.g. mu) so read as mac and convert later to windows=1252
  defaultEncoding="macintosh"
  WaveHeader2$WaveName=.readNullTermString(con,18+2,encoding=defaultEncoding)
  myread(n=2) # Skip 8 bytes
  WaveHeader2$dataUnits=.readNullTermString(con,3+1,encoding=defaultEncoding)
  WaveHeader2$xUnits=.readNullTermString(con,3+1,encoding=defaultEncoding)
  WaveHeader2$npts=myread()
  myread(size=2) # skip aModified
  WaveHeader2$hsA=myread(what="double",size=8)
  WaveHeader2$hsB=myread(what="double",size=8)
  myread()
  WaveHeader2$fsValid=myread(size=2)!=0
  WaveHeader2$topFullScale=myread(what="double",size=8)
  WaveHeader2$botFullScale=myread(what="double",size=8)
  myread(size=1,n=10)
  WaveHeader2$CreationDate=.readIgorDate(con,endian)
  WaveHeader2$platform=myread(size=1,signed=FALSE)  
  myread(size=1)
  WaveHeader2$modDate=.readIgorDate(con,endian)
  myread()
  
  # reencode character fields
  # if(WaveHeader2$platform==1 || (WaveHeader2$platform==0 && endian=="big")){
  if(WaveHeader2$platform==1 || (WaveHeader2$platform==0)){
    encoding="macintosh"
  } else {
    encoding="WINDOWS-1252"
  }
  if(Verbose) print(paste("encoding is:", encoding,"\n"))
  charfields=sapply(WaveHeader2,is.character)
  WaveHeader2[charfields]=lapply(WaveHeader2[charfields],iconv,from=encoding,to="")
  
  WaveTypeInfo=.DetermineIgorWaveType(WaveHeader2,endian=endian)
  if(Verbose) cat("signed=",(WaveTypeInfo$signAdjustment==0),"\n")
  
  if(HeaderOnly) {
    x=NA; attr(x,"WaveHeader")=WaveHeader2
    return(x)
  }
  # Note that only unsigned 16 bit data can be read by readBin
  WaveData=myread( what=WaveTypeInfo$what,size=WaveTypeInfo$size,n=WaveHeader2$npts,
    signed=(WaveTypeInfo$signAdjustment==0))
  # Handle adjustment for unsigned integer data
  if(WaveTypeInfo$signAdjustment!=0  && WaveTypeInfo$size>2){
    WaveData=.ConvertIntToUInt(WaveData,WaveTypeInfo$signAdjustment)
  }
  
  # Finish up
  attr(WaveData,"dataUnits")=WaveHeader2$dataUnits
  attr(WaveData,"dimUnits")=WaveHeader2$dimUnits[WaveHeader2$nDim>0]

  if(BinHeader2$noteSize>0) {
    #attr(WaveData,"Note")=.readCharsWithEnc(con,BinHeader2$noteSize,encoding)
  }
  attr(WaveData,"WaveHeader")=WaveHeader2
  attr(WaveData,"BinHeader")=BinHeader2
  attr(WaveData,"start")=attr(WaveData,"WaveHeader")$hsB
  attr(WaveData,"deltat")=attr(WaveData,"WaveHeader")$hsA
  WaveData
}

.ReadIgorBinary.V5<-function(con,Verbose=FALSE,ReturnTimeSeries=NULL,endian=NULL,HeaderOnly=FALSE){
  # File pointer should be positioned at the start of the header
  # this should be 2 bytes into the wave data 
  # (ie after endian and version bytes).  Store that position:
  startPos=seek(con)
  
  # Read the BinHeader5 (64 bytes)
  # 1 byte  char
  # 2 bytes short
  # 4 bytes int, long, float, Handle, any kind of pointer
  # 8 bytes double
  MAXDIMS=4

# typedef struct BinHeader5 {
#   short version;            // Version number for backwards compatibility.
#   short checksum;           // Checksum over this header and the wave header.
#   long wfmSize;           // The size of the WaveHeader5 data structure plus the wave data.
#   long formulaSize;         // The size of the dependency formula, if any.
#   long noteSize;            // The size of the note text.
#   long dataEUnitsSize;        // The size of optional extended data units.
#   long dimEUnitsSize[MAXDIMS];    // The size of optional extended dimension units.
#   long dimLabelsSize[MAXDIMS];    // The size of optional dimension labels.
#   long sIndicesSize;          // The size of string indicies if this is a text wave.
#   long optionsSize1;          // Reserved. Write zero. Ignore on read.
#   long optionsSize2;          // Reserved. Write zero. Ignore on read.
# } BinHeader5;

  BinHeader5Def=data.frame(name=I(c("checksum","wfmSize","formulaSize","noteSize","dataEUnitsSize",
      "dimEUnitsSize","dimLabelsSize","sIndicesSize","optionsSize1","optionsSize2")),
    size=c(2,rep(4,9)),what="integer",n=c(rep(1,5),rep(4,2),rep(1,3)))
  BinHeader5=list()
  for(i in seq(nrow(BinHeader5Def)))
    BinHeader5[[BinHeader5Def$name[i]]]=readBin(con,endian=endian,what=BinHeader5Def$what[i],size=BinHeader5Def$size[i],n=BinHeader5Def$n[i])
  if(Verbose) print(BinHeader5)
  
  # Read the WaveHeader5 (320 bytes)
  seek(con,where=startPos+64-2)
  WaveHeader5=list()
  myread=function(what="integer",size=4,...) readBin(con,endian=endian,what=what,size=size,...)
  
  # next WaveHeader5 pointer
  myread()
  WaveHeader5$creationDate=.readIgorDate(con,endian)
  WaveHeader5$modDate=.readIgorDate(con,endian)
  WaveHeader5$npts=myread()
  WaveHeader5$type=myread(what="raw",size=1,n=2)
  
  myread(size=1,n=10)
  WaveHeader5$WaveName=.readNullTermString(con,32)
  myread(n=2,size=4)
  WaveHeader5$nDim=myread(n=MAXDIMS)
  #highestDim=max(which(WaveHeader5$nDim>0))
  dims=WaveHeader5$nDim[WaveHeader5$nDim>0]
  nDims=length(dims)
  #if(any(WaveHeader5$nDim[-1]>0)) warning(paste("Ignoring additional dimensions for wave with dims",WaveHeader5$nDim))

  WaveHeader5$sfA=myread(n=MAXDIMS,what="double",size=8)
  WaveHeader5$sfB=myread(n=MAXDIMS,what="double",size=8)  
  # Units
  WaveHeader5$dataUnits=.readNullTermString(con,3+1)
  WaveHeader5$dimUnits=replicate(MAXDIMS,.readNullTermString(con,3+1))
  
  # Platform
  # skip to correct location
  seek(con,2+2+8+8+4*(MAXDIMS*2+2),origin="current")
  # unsigned char platform;
  # // 0=unspecified, 1=kMacPlatform, 2=kWinPlatform

  WaveHeader5$platform=myread(size=1,signed=FALSE)
  if(WaveHeader5$platform==1 || (WaveHeader5$platform==0 && endian=="big")){
    encoding="macintosh"
  } else {
    encoding="WINDOWS-1252"
  }
  # reencode character fields
  charfields=sapply(WaveHeader5,is.character)
  WaveHeader5[charfields]=lapply(WaveHeader5[charfields],iconv,from=encoding,to="")

  # Read the wave data 
  seek(con,startPos+64+320-2)
  WaveTypeInfo=.DetermineIgorWaveType(WaveHeader5,endian=endian)
  if(Verbose) print(WaveHeader5)
  if(Verbose) cat("data position = ",seek(con),"\n")

  # Have the option to return the header info alone - as a sort of
  # preview
  if(HeaderOnly) {
    x=NA; attr(x,"WaveHeader")=WaveHeader5
    return(x)
  }

  if(WaveTypeInfo$what=="character"){
    # gist is we are going to have a char vector containing npts items;
    # each item is going to be nchars longs, where nchars is a vector
    # generated by reading npts longs from right after the end
    # of the txt wave.
    
    dataStartPos=seek(con) #store current location
    # skip to the end of the wave data
    seek(con,sum(unlist(BinHeader5[2:7]))-320,"current")
    if(Verbose) cat("Reading offsets at:",seek(con),"\n")
    offsets=myread(n=WaveHeader5$npts)
    if(Verbose) cat("offets:",offsets,"\n")
    nchars=diff(c(0,offsets))
    if(sum(nchars)>BinHeader5$wfmSize  || any(nchars<0) ) {
      warning(paste("Unable to read text wave",WaveHeader5$WaveName))
      return (NULL)
    }
    # return to start of data
    seek(con,dataStartPos)
#     if(T) WaveData="test"
#     else{
    WaveData=readChar(con,nchars)
    # convert from appropriate text encoding to R default
    if(WaveHeader5$platform==1 || (WaveHeader5$platform==0 && endian=="big")){
      if(Verbose) cat("Converting text wave from macintosh encoding to R default\n")
      WaveData=iconv(WaveData,from="macintosh",to="")
    } else {
      if(Verbose) cat("Converting text wave from WINDOWS-1252 encoding to R default\n")
      WaveData=iconv(WaveData,from="WINDOWS-1252",to="")
    }
#     }
  } else {
    # Note that only unsigned 16 bit data can be read by readBin
    WaveData=myread( what=WaveTypeInfo$what,size=WaveTypeInfo$size,n=WaveHeader5$npts,
      signed=(WaveTypeInfo$signAdjustment==0))
    # Handle adjustment for unsigned integer data
    if(WaveTypeInfo$signAdjustment!=0  && WaveTypeInfo$size>2){
      # nb must convert int to numeric because R does not handle ints >
      # 2^31-1
      WaveData=.ConvertIntToUInt(WaveData,WaveTypeInfo$signAdjustment)
    }
  }
  # Set units 
  attr(WaveData,"dataUnits")=WaveHeader5$dataUnits
  attr(WaveData,"dimUnits")=WaveHeader5$dimUnits
  # Read anything else
  if(BinHeader5$formulaSize>0) {
    attr(WaveData,"Dependency Formula")= .readNullTermString(con,BinHeader5$formulaSize,encoding)
  }
  if(BinHeader5$noteSize>0) {
    attr(WaveData,"Note")=.readCharsWithEnc(con,BinHeader5$noteSize,encoding)
  }
  
  if(BinHeader5$dataEUnitsSize>0) {
    attr(WaveData,"dataUnits")=.readNullTermString(con,BinHeader5$dataEUnitsSize,encoding)
  }
  
  if(any(BinHeader5$dimEUnitsSize>0)) {
    x=.readCharsWithEnc(con,BinHeader5$dimEUnitsSize,encoding)
    attr(WaveData,"dimUnits")[BinHeader5$dimEUnitsSize>0]=x[BinHeader5$dimEUnitsSize>0]
  }
  # Trim units down to active dimensions
  attr(WaveData,"dimUnits")=attr(WaveData,"dimUnits")[WaveHeader5$nDim>0]

  if(any(BinHeader5$dimLabelsSize>0)) {
    x=.readCharsWithEnc(con,BinHeader5$dimLabelsSize,encoding)
    attr(WaveData,"dimLabels")=x[WaveHeader5$nDim>0]
  }
  
  # Finish up
  # Re-dimension
  if(nDims>1) dim(WaveData)=dims
  attr(WaveData,"WaveHeader")=WaveHeader5
  attr(WaveData,"BinHeader")=BinHeader5
  return (WaveData)
}

#' Convert an Igor wave (or list of waves) loaded by read.ibw into an R time series 
#' 
#' Where there are multiple waves, they are assumed to be of compatible lengths 
#' so that they can be joined together by cbind.
#' 
#' @param WaveData, a wave or list of waves
#' @param ReturnOriginalDataOnError If we can't make a time series, return
#'  return original data (default TRUE)
#' @return a time series or multi time series (ts, mts)
#' @author jefferis
#' @export
WaveToTimeSeries<-function(WaveData,ReturnOriginalDataOnError=TRUE){
  if(is.list(WaveData)) {
    # process separate waves into multi wave time series
    l=lapply(WaveData,WaveToTimeSeries)
    return(do.call(cbind,l))
  }
  tspw = tsp.igorwave(WaveData)
  res=try(ts(WaveData,start=tspw[1],frequency=tspw[3]),silent=TRUE)
  if(inherits(res,'try-error')){
    warning(res)
    WaveData
  } else res
}

#' Return tsp attribute of igor wave (start, end, frequency)
#' 
#' Note that end = (npts-1) * deltat
#' @param wave Igor wave loaded by read.ibw or read.pxp
#' @return numeric vector with elements start, end, frequency
#' @author jefferis
#' @seealso \code{tsp}
#' @export
tsp.igorwave<-function(wave){
  bh=attr(wave,"BinHeader")
  if(is.null(bh)) stop("Can't find BinHeader")
  wh=attr(wave,"WaveHeader")
  if(is.null(wh)) stop("Can't find WaveHeader")
  if(bh$version==2){
    # WaveHeader2
    start=wh$hsB
    deltat=wh$hsA
    npts=wh$npts
  } else if(bh$version==5){
    # WaveHeader5, can be multi-dimensional
    start=wh$sfB[1]
    deltat=wh$sfA[1]
    npts=wh$nDim[1]
  } else {
    stop(paste("Unsupported Igor Wave Version Number,",bh$version))
  }
  return(c(start=start,end=deltat*(npts-1),frequency=1/deltat))
}
