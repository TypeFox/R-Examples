#####################################################################
## Thomas Hoffmann                                                  #
## EXPORTED:    07/08/2005                                          #
## MODIFIED:    01/25/2006                                          #
##                                                                  #
## DESCRIPTION:                                                     #
##  Creating the command file that is passed to pbat.               #
##  Also contains the routine that gets pbatdata.zip from the       #
##   internet, gets a time stamp / tests if a string has been       #
##   stamped, and some of the vector testing/pasting routines.      #
#####################################################################

PBATDATAURL <- "http://www.biostat.harvard.edu/~clange/pbatdata.zip";

#####################################################################
## getPbatdata()                                                    #
## Gets the 'pbatdata.txt' file, from the internet.                 #
#####################################################################
getPbatdata <- function() {
  return(); ## sigh... this file was removed!!!

#   ## Give the user a chance to say yes or no:
#   msgStr <- paste("Can I attempt to download 'pbatdata.txt' from '",
#                   PBATDATAURL,
#                   "'? This file is needed.", sep="");
#   if( isPackageLoaded( 'tcltk' ) ) {
#     if( "yes" != tclvalue(tkmessageBox(title="pbatdata.txt",message=msgStr,icon="question",type="yesno")) )
#       return();
#   }else{
#     if( 'yes' != textMessageBox( msgStr, c('yes','no') ) )
#       return();
#   }
#
#   ## Carry on with downloading
#   pbatpath <- str.getpath(pbat.get());
#   ######zipfile <- paste( pbatpath, "/pbatdata.zip", sep="" );  ## not used
#   ######if( pbatpath=="" ) zipfile <- "./pbatdata.zip";         ## ?codetools.
#   ###download.file( PBATDATAURL, zipfile );
#   download.file( PBATDATAURL, "./pbatdata.zip" );
#   Sys.sleep(2); ## Just in case this is why it was getting corrupted
#   pbatdatafile <- zip.file.extract( file="pbatdata.txt", zipname="pbatdata.zip" );
#   destfile <- paste( pbatpath, "/pbatdata.txt", sep="" );
#   if( pbatpath=="" ) destfile <- "pbatdata.txt";
#   if( file.exists(pbatdatafile) ){
#     file.copy(pbatdatafile, destfile);
#   }
}

#####################################################################
## checkAndGetPbatdata()                                            #
## Gets the 'pbatdata.txt' into the cwd, trying _everywhere_!       #
## Exported out of it's original function 7/25/2006 for use with    #
##  the power routines.                                             #
## Note: this will error if it cannot be found (it's not going to   #
##  compute the prob table).                                        #
#####################################################################
checkAndGetPbatdata <- function() {
  ## First try to copy from the pbat directory.
  pbatdatafile <- paste( str.getpath(pbat.get()), "/pbatdata.txt", sep="" );
  if( file.exists( pbatdatafile ) )
    file.copy( from=pbatdatafile, to=paste(getwd(),"/pbatdata.txt",sep="") );

  ## If we can't find it, then try other things
  if( !file.exists( paste(getwd(),"/pbatdata.txt",sep="") ) ) {

    ## See if it's anywhere in the path
    newLoc <- pathFindFile("pbatdata.txt");
    if( newLoc != "" ){
      ## found it! copy it over!
      file.copy( from=newLoc, to=paste(getwd(),"/pbatdata.txt",sep="") ); ## seems to have to be in cwd - more than just the path somewhere
    }else{
      ## Last thing to try is to download from the internet
      getPbatdata(); ## but puts it in pbat dir first
      if( file.exists( pbatdatafile ) ) {
        file.copy( from=pbatdatafile, to=paste(getwd(),"/pbatdata.txt",sep="") );
      }
    }

    ## So make sure that it finally got copied in
    if( !file.exists( paste(getwd(),"/pbatdata.txt",sep="") ) ) {
      stop( paste("'pbatdata.txt was not found in the current",
                  " working directory '", getwd(),
                  "', or in the pbat directory '", pbatdatafile, "'",
                  ", or anywhere in your current path,",
                  " and it could not be downloaded online. ",
                  " Please see",
                  " http://www.biostat.harvard.edu/~clange/Downloading%20PBAT.htm",
                  " for more details.", sep="" ) );
    }

    ## Old coding below (changed order to look)

    ## Now, see if we can find it in the pbat directory location, and go from there...
    ##getPbatdata();
    ##if( file.exists( pbatdatafile ) )
    ##  file.copy( from=pbatdatafile, to=paste(getwd(),"/pbatdata.txt",sep="") );
    ##if( !file.exists( paste(getwd(),"/pbatdata.txt",sep="") ) ) {
    ##  if( pathFindFile( "pbatdata.txt" ) == "" ){
    ##    stop( paste("'pbatdata.txt was not found in the current",
    ##                " working directory '", getwd(),
    ##                "', or in the pbat directory '", pbatdatafile, "'",
    ##                ", or anywhere in your current path,",
    ##                " and it could not be downloaded online. ",
    ##                " Please see",
    ##                " http://www.biostat.harvard.edu/~clange/Downloading%20PBAT.htm",
    ##                " for more details.", sep="" ) );
    ##  }
    ##}
  }

}

#####################################################################
## getTimeStamp(...)                                                #
## Gets a unique time-stamp string for the output!                  #
#####################################################################
getTimeStamp <- function() {
  zpad <- function(n, pad=2) {
    if( nchar(n)<pad )
      return( paste( rep("0",pad-nchar(n)), n, sep="" ) );
    return(n);
  }

  options(digits.secs=6) ## millisecond timing out to the millionth of a second... try to eliminate that race condition
  d <- as.POSIXlt( Sys.time() );
  return( paste( 1900+d$year, zpad(d$mon), zpad(d$mday), zpad(d$hour), zpad(d$min), zpad(floor(d$sec*1000000)), sep="" ) );  ## R2.3 change... seconds decide to have bloody decimal points... why can't we just be consistent between releases???? WHY??????
}

## Returns if all characters in a string are numbers (decimals not allowed since we're checking for time-stamps ultimately)
isStringAllNumeric <- function( str ){
  if( nchar(str)==0 ) return(TRUE); ## I guess so

  for( i in 1:nchar(str) ) {
    ch <- substring(str,i,i);
    #if( !(  ( '0'<=ch & ch<='9') | ch=='.'  ) )
    if( !( '0'<=ch & ch<='9' ) )
      return(FALSE);
  }
  return(TRUE);
}

isTimeStamped <- function( str, extLen=3 ) {
  timeLen <- strlen(getTimeStamp());
  strLen <- strlen(str);
  strPbatLen <- strlen("pbat")+1;  # 10/07/2005
  if( strLen < 5+timeLen ) return( FALSE ); # too short
  possTimeStr <- substring( str, strPbatLen, strPbatLen+timeLen-1 );
  #print( possTimeStr );
  ###if( !is.na( as.numeric( possTimeStr ) ) )
  ###  return(TRUE);
  ###return(FALSE);

  return( isStringAllNumeric(possTimeStr) );  ## No more warnings!
}

# Assumes isTimeStamped( str ) _already_ returned TRUE
extractTimeStamp <- function( str, extLen=3 ) {
  timeLen <- strlen(getTimeStamp());
  strLen <- strlen(str);
  strPbatLen <- strlen("pbat")+1;  # 10/07/2005
  if( strLen < 5+timeLen ) return( FALSE ); # too short
  possTimeStr <- substring( str, strPbatLen, strPbatLen+timeLen-1 );
  return( possTimeStr );
}


## Moved outside of pbat.create.commandfile(...) 9/?/05
############################################
## Whether 'subcol' is contained in 'col'. #
############################################
isVecContained <- function( subcol, col ) {
  for( i in 1:length(subcol) ) {
    if( sum( subcol[i]==col ) != 1 )
      return( FALSE );
  }
  return( TRUE );
}

# Moved outside of pbat.create.commandfile(...) 9/20/05
#################################################################
## paste a vector of strings together                           #
## SQUOTE    if TRUE, surrounds each string with a single quote #
## COMMASEP  if TRUE, comma seperates the values                #
#################################################################
pasteVector <- function( vector, SQUOTE=FALSE, COMMASEP=FALSE ) {
  if( length(vector) < 1 ) return("");

  squote <- function(str) {return(str)};
  if( SQUOTE==TRUE )
    squote <- function(str){return(paste("'",str,"'",sep=""));};

  if( length(vector) == 1 ) return( squote(vector[1]) );
  strRet = squote(vector[1]);
  for( i in 2:length(vector) ) {
    if( COMMASEP==FALSE ) {
      strRet <- paste(strRet, squote(vector[i]));
    }else{
      strRet <- paste(strRet, ", ", squote(vector[i]), sep="");
    }
  }
  return(strRet);
} # Status: debugged
##### And just a wrapper to make commands simpler later #######
csPasteVector <- function( vector ){
  return( pasteVector( vector, SQUOTE=TRUE, COMMASEP=TRUE ) );
} # Status: debugged

pasteVector2 <- function( vector, sep=" " ){
  if( length(vector) < 1 ) return("");

  if(length(vector) == 1 ) return(vector);
  strRet <- vector[1];
  for( i in 2:length(vector) )
    strRet <- paste(strRet, sep, vector[i], sep="");
  return(strRet);
}

## Moved outside pbat.create.commandfile(...) 5/31/06
#####################################
## write the command (unless empty) #
#####################################
writeCommand <- function( commandStr, vals, end=FALSE, outfile=NULL ) {
  if(  !( length(vals)==1 && (is.na(vals)||vals=="") )  ) {
    if( end==FALSE ) {
      writeLines( paste( commandStr, pasteVector(vals) ), con=outfile );
    }else{
      writeLines( paste( commandStr, pasteVector(vals), "end" ), con=outfile );
    }
  }
} # Status: debugged

####################################################################
## Prints helpful error msg if 'subcol' is not contained in 'col'. #
## if AT.MOST.SINGLETON=TRUE, then 'subcol' can be of at           #
##  most length 1                                                  #
###################################################################
errorVecNotContained <- function( command, subcol, col,
                                 AT.MOST.SINGLETON=FALSE ){
  if( length(subcol)==1 && subcol=="" ) return(); # empty set is contained!

  if( !isVecContained(subcol,col) ) {
    stop( paste( "For the option '", command,
                "', the values that you specified {",
                csPasteVector(subcol),
                "} did not match the possible values {",
                csPasteVector(col),
                "}.",
                sep=""
                ) );
  }
  if( length(subcol)>1 && AT.MOST.SINGLETON==TRUE ) {
    stop( paste("For the option '", command,
                "', the values that you specified {",
                csPasteVector(subcol),
                "} was of too long a length.  Keep in mind it must take at most a single string value from the following collection {",
                csPasteVector(col), "}.",
                sep=""
                ) );
  }
}

##########################################################################
## Prints out a useful error message if their is a nonempty intersection #
##  between 'col1' and 'col2'                                            #
##########################################################################
errorIfAnyMatch <- function( col1, col2, nameCol1, nameCol2 ) {
  if( length(col1)==1 && col1=="" ) return(); # collection is empty
  if( length(col2)==1 && col2=="" ) return();

  for( i in 1:length(col1) )
    if( sum(col1[i]==col2)>0 )
      stop( paste("There should not be any overlap in the following two options collections: ",
                  nameCol1, "{", csPasteVector(col1), "}, ",
                  nameCol2, "{", csPasteVector(col2), "}.",
                  sep="" ) );
}

################################################
## write the command from human-readable input #
## The beauty of this function is it includes  #
##   the error handling routines within.       #
## stop() is like throwing an exception, which #
##   if not caught goes all the way to the user#
##   and halts the program.                    #
###############################################
writeCommandStrMatch <- function( commandStr, str, strVec, vals=c(0:(length(strVec)-1)), outfile=NULL ){
  if( sum(str==strVec)!=1 ) {
    ## Error! New: print out a message here for ease!
    stop( paste( "'", commandStr, "' can only take on the following values: ",
                csPasteVector( strVec ),
                ".  You passed the invalid value '", str, "'.",
                sep="" ) );
  }

  writeCommand( commandStr, vals[which(str==strVec)], outfile=outfile );
  return( TRUE ); # success!
}
###### Again just another wrapper to make life a little easier. #######
writeCommandStrMatch1 <- function( commandStr, str, strVec, outfile ) {
  return(  writeCommandStrMatch( commandStr, str, strVec, vals=c(1:length(strVec)), outfile=outfile )  );
} # NOT DEBUGGED

################################################################
## Prints error messages if outside a range (Closed interval). #
##  If the value is NULL or "", then there is no max/min       #
################################################################
errorRangeCheck <- function( commandStr, value, min=1, max=NULL, IS.INTEGER=TRUE ) {
  hasMin <- !is.null(min) && min!="";
  hasMax <- !is.null(max) && max!="";

  minStr <- "-INF"; ## 12/28/2008
  maxStr <- "INF";  ## 12/28/2008

  if( (hasMin && value<min) || (hasMax && value>max) ) {
    minStr <- as.character(min); if(!hasMin) minStr <- "INF";
    maxStr <- as.character(max); if(!hasMax) maxStr <- "INF";
    stop( paste("For the option '", commandStr, "', the value must be in the range [",
                minStr, ",", maxStr, "]. The value you supplied was ", value, ".",
                sep="") );
  }

  ##if( IS.INTEGER && value!=floor(value) )
  if( IS.INTEGER && as.numeric(value)!=floor(as.numeric(value)) )
    stop( paste("For the option '", commandStr,
                "', the value must be _an integer_ in the range [",
                minStr, ",", maxStr, "]. The value you supplied was ", value, ".",
                sep="") );
}

## 08/31/2010
## All possible subsets
all.possible.subsets = function(snps){
  sublist = list()
  for(i in 1:length(snps))
    sublist[[i]] = c(FALSE, TRUE) #c(TRUE,FALSE)
  snpgrid = expand.grid(sublist)

  haplolist = list()
  for(g in 1:nrow(snpgrid))
    if(!all(snpgrid[g,] == FALSE))
      haplolist[[length(haplolist) + 1]] = snps[as.logical(snpgrid[g,])]
  return(haplolist)
}
all.possible.subsets(paste("m", 1:4, sep=""))

#####################################################################
## pbat.create.commandfile(...)                                     #
## LOTS of options and info! See .Rd file!                          #
##                                                                  #
## This includes all of the debugging code and everything.          #
## Hopefully some good debugging here :).                           #
####################################################################
pbat.create.commandfile <- function(
       pedfile, phefile="",
       snps="",
       phenos="", time="",
       preds="", preds.order="",
       inters="",
       groups.var="", groups="",
       fbat="gee",
       censor="",
       max.pheno=1, min.pheno=1,
       null="no linkage, no association", alpha=0.05,
       trans.pheno="none", trans.pred="none", trans.inter="none",
       scan.pred="all", scan.inter="all",
       scan.genetic="additive",
       offset="gee",
       screening="conditional power", distribution="default",
       logfile="",
       max.gee=1,
       max.ped=7, min.info=0,
       haplos=NULL, incl.ambhaplos=TRUE, infer.mis.snp=FALSE,
       sub.haplos=FALSE, length.haplos=2, adj.snps=TRUE,
       overall.haplo=FALSE, cutoff.haplo=FALSE,
       output="normal",
       max.mating.types=10000,
       commandfile="",
       future.expansion=NULL,
       LOGFILE.OVERRIDE=TRUE, ## But FALSE for the GUI!?
       monte=0,
       mminsnps=NULL, mmaxsnps=NULL,
       mminphenos=NULL, mmaxphenos=NULL,
       env.cor.adjust=FALSE,
       gwa=FALSE,
       snppedfile=FALSE,
       extended.pedigree.snp.fix=FALSE,
       new.ped.algo=FALSE,
       cnv.intensity=2, cnv.intensity.num=3
                                    )
{
  ##cat( "cnv.intensity", cnv.intensity, "cnv.intensity.num", cnv.intensity.num, "\n" )
  ##cat( "createCommandfile pedfile", pedfile, "\n" );
  ## Fix up a couple of variables passed in
  gwa <- (gwa==TRUE); ## in case it was a string
  snppedfile <- (snppedfile==TRUE);

  ## NEW 08/31/2010, if extended.pedigree.snp.fix=TRUE, turn off new.ped.algo...
  if(extended.pedigree.snp.fix)
    new.ped.algo=FALSE

  if(overall.haplo==TRUE && screening=="conditional power")
    stop("To use overall.haplo=TRUE, you must set screening='wald'.")

  ##-----------------------------
  ## fix up extensions / naming -
  ##-----------------------------
  ####pedfile <- str.file.extension( pedfile, ".ped" );
  pedfile.ext <- file.extension( pedfile );
  #print( "pedfile" );
  #print( pedfile );
  #print( "pedfile.ext" );
  #print( pedfile.ext );
  if( pedfile.ext!="ped" & pedfile.ext!="pped" & pedfile.ext!="cped" ) {
    pedfile <- paste( pedfile, ".ped", sep="" );
    pedfile.ext <- "ped";
  }


  #print( "pedfile" );
  #print( pedfile );

  if( phefile=="" ) {
    phefile <- paste( substring(pedfile,1,strlen(pedfile)-3), "phe", sep="" );
    ## But we don't need a phe file for AffectionStatus
    if( !file.exists( phefile ) ) phefile <- "";
  }else{
    phefile <- str.file.extension( phefile, ".phe" );
  }
  if( logfile=="" ) {
    if( isTimeStamped( pedfile ) ) {
      logfile <- paste( substring(pedfile,1,strlen(pedfile)-4),
                        "", sep="" );
    }else{
      logfile <- paste( substring(pedfile,1,strlen(pedfile)-4),
                       getTimeStamp(), "", sep="" );
    }
  }
  ## New overriding
  if( LOGFILE.OVERRIDE ) {
    potlogfile <- str.extract.afterb(logfile, '/' );
    if( potlogfile=="" ) {
      potlogfile <- str.extract.afterb(logfile, '\\' );
    }else{
      #addi 09/05/2006
      potlogfile2 <- str.extract.afterb(logfile, '\\' );
      if( potlogfile2!="" )
        potlogfile <- potlogfile2;
      #idda
    }

    if( potlogfile != "" )
      logfile <- potlogfile;
  }


  if( commandfile=="" )
    commandfile <- paste( substring(logfile,1,strlen(logfile)-3), "txt", sep="" );
  #print( "COMMANDFILE" );
  #print( commandfile );
  #print( "LOGFILE" );
  #print( logfile );

  #warning( "COMMented out this part of logfile stuff too..." );
#  if( file.exists(logfile) )
#    stop( paste("Logfile '",logfile,"' already exists!  (Note, the default name for the 'logfile' option is the prefix of the pedigree file with a .txt suffix.", sep="" ) );

  # Take certain strings to lowercase (but NOT the phenos stuff!)
  fbat <- tolower(fbat);
  null <- tolower(null);
  trans.pheno <- tolower(trans.pheno);
  trans.pred <- tolower(trans.pred);
  trans.inter <- tolower(trans.inter);
  scan.pred <- tolower(scan.pred);
  scan.inter <- tolower(scan.inter);
  scan.genetic <- tolower(scan.genetic);
  offset <- tolower(offset);
  screening <- tolower(screening);
  distribution <- tolower(distribution);
  output <- tolower(output);

  # Note: these functions are included within this function
  #        because they are operating on variables in this
  #        function.  Slightly confusing, but I think it makes
  #        the most sense this way.


  #-----------------
  # Some debugging -
  #-----------------

  # first make sure certain files exist

  ## The infamous pbatdata.txt file...
  checkAndGetPbatdata();
#   print("haplos")
#   print(haplos)

  if( !file.exists(pedfile) )
    stop( paste("The pedigree file '",
                pedfile, "' does not exist.  Current working directory is '",
                getwd(), "'.", sep="") );
  if( phefile!="" && !file.exists(phefile) )
    stop( paste("The phenotype file '",
                phefile, "' does not exist.  Current working directory is '",
                getwd(), "'.", sep="") );

  # other debugging
  if( phenos[1]!="" & time[1]!="" )
    stop( "Both 'phenos' and 'time' cannot have values set to them.  See the help file for more details." );

  # New debugging -- Christoph explains...
  if( scan.genetic=="all" &&
     !is.null(mminsnps) && !is.null(mmaxsnps) && !is.null(mminphenos) && !is.null(mmaxphenos) ) {
    stop( "Multimarker mode (mminsnps, mmaxsnps, mminphenos, mmaxphenos) is not supported for scan.genetic='all'. Try doing each one in turn." );
  }
  if( gwa==TRUE ) {
    if( output!="short" ){
      output <- "short";
      warning( "gwa mode only supported with short output format; short output format enforced." );
    }
    if( !snppedfile ){
      msg <- "Ensure that you do not really have just snps. This will go much faster (storage enhancment) if it is and you specify 'snppedfile=TRUE'. Stop this and try again if that is the case.";
      print( msg ); ## need to try to get it to the user as fast as possible
      warning( msg );
    }
  }

  # much more advanced debugging!

  # pedigree file information
  #ped.b <- read.badheader( pedfile );
  #posSnps <- ped.b$header[3:length(ped.b$header)];
  ##posSnps <- read.badheader( pedfile )$header;
  posSnps <- NULL;
  if( pedfile.ext != "pped" )
    posSnps <- read.badheader( pedfile, onlyHeader=TRUE )$header;
  ####print( "got past here" );

  # phenotype file information
  posPhenos <- c();
  if( !is.null(phefile) && phefile!="" ){
    phe <- read.phe( phefile, sym=TRUE );
    posPhenos <- names(phe);
  }
  posPhenos <- c( posPhenos, "AffectionStatus" ); ## When, where, and _why_ did this get lost??

  # check containment of various options...
  if( !is.null(posSnps) )
    errorVecNotContained( "snps", snps, posSnps );

  errorVecNotContained( "phenos", phenos, posPhenos );
  errorVecNotContained( "time", time, posPhenos, AT.MOST.SINGLETON=TRUE );
  ##errorVecNotContained( "inters", inters, phenos );
  errorVecNotContained( "inters", inters, preds );  ## 01/18/2006
  #############errorVecNotContained( "groups", groups, posPhenos );
  ##print( "groups" );
  ##print( groups ); print( groups.var );

  errorVecNotContained( "preds", preds, posPhenos ) ########################################################################################################
  errorVecNotContained( "censor", censor, posPhenos )
  errorVecNotContained( "groups", groups.var, posPhenos );

  errorIfAnyMatch( groups.var, phenos, "groups", "phenos" );
  errorIfAnyMatch( groups.var, time, "groups", "time");
  errorIfAnyMatch( groups.var, censor, "groups", "censor" );
  errorIfAnyMatch( groups.var, preds, "groups", "preds" );
  errorIfAnyMatch( preds, censor, "preds", "censor" );
  errorIfAnyMatch( preds, time, "preds", "time" );
  errorIfAnyMatch( preds, phenos, "preds", "phenos" );
  errorIfAnyMatch( censor, time, "censor", "time" );
  errorIfAnyMatch( censor, phenos, "censor", "phenos" );
  errorIfAnyMatch( time, phenos, "time", "phenos" );

  ## Enforce haplotype mode for multiprocessing
  ## 01/18/2006 rewrite - this should _always_ be done!
  ## 01/29/2007 son of a ... apparently not always...
  ##if( pbat.getNumProcesses() > 1 && is.null(haplos) ) {
  if( is.null(haplos) && extended.pedigree.snp.fix==FALSE ) {
    haplos <- list();
    if( is.null(snps) || snps[1]=="" ) {
      # simple hack - we need to insert the names into the haplotypes...
      # - it's really not pretty but it solves our problem so simply
      #   now that ew have the addition of multiple processing built
      #   in to pbat.
      # - on a second note, it's really the only way to be able to fix
      #   this on such a low level to guarantee that the command-line
      #   stuff will also work without a massive changes?

      ####### 09/08/2006 - what is going on here?
      #junk <- read.ped( pedfile ); ##, lowercase=FALSE ); ## that lowercase...
      #allSnps <- names( as.pedlist( junk ) );
      #haplos[[1]] <- NULL;
      #if( !is.null(allSnps) )
      #  haplos[[1]] <- allSnps[7:length(allSnps)];  ## 01/18/06 fix - list
      ####### 09/08/2006 - what is going on here?

      # it now suffices just to do 'haplos 1 end' for them all; this is pure symbolic!
      haplos[[1]] <- snps; ## additien 09/08/2006ish
    }else{
      haplos[[1]] <- snps;
    }
    snps <- "";
    sub.haplos <- TRUE;
    length.haplos <- 1;
    adj.snps <- TRUE;
  }

  ## caveat -- extended pedigree snp fix can only be done in single mode...
  if( extended.pedigree.snp.fix=="TRUE" && pbat.getmode()$mode!="single" )
    stop( "The extended pedigree SNP fix can only be done in single mode currently." );

  ####print( "got past here 3" );

  # Verification of the haplos structure
  if( !is.null(haplos) && !is.list(haplos) )
    stop( "Haplos must be a list of string vectors (or objects that can coerced into strings." );
  if( !is.null(haplos) && is.list(haplos) && length(haplos)>0 ) {
    # check to make sure the snps are in, and no overlap.
    for( i in 1:length(haplos) ){
      # make sure the snps are in the list
      errorVecNotContained( "haplos dataframe", haplos[[i]], posSnps );

#       # then make sure there is no overlap (can't have a snp in more than one block!)
#       ##print( length(haplos) );
#       if( i<length(haplos) ) {
#         for( j in (i+1):length(haplos) ) {
#           errorIfAnyMatch(haplos[[i]], haplos[[j]],
#                           paste("Haplotype block ",i," (",names(haplos)[i],") ",sep="" ),
#                           paste("Haplotype block ",j," (",names(haplos)[j],") ",sep="" ) );
#         }
#       }
    }

    ## Now, what if sub.haplos?
    ## NEW, NEW, NEW
    if(sub.haplos==TRUE && adj.snps==FALSE){
      sub.haplos = FALSE
      haplos2 = list()
      for(h in 1:length(haplos)){
        haplosaddi = all.possible.subsets(haplos[[h]])
        for(hh in 1:length(haplosaddi))
          haplos2[[length(haplos2) + 1]] = haplosaddi[[hh]]
      }
      haplos = haplos2
    }
  }

  # Some simple range checking...
  errorRangeCheck( "max.pheno", max.pheno, min=min.pheno );
  errorRangeCheck( "min.pheno", min.pheno, max=max.pheno );
  errorRangeCheck( "alpha", alpha, min=0, max=1, IS.INTEGER=FALSE );
  errorRangeCheck( "max.gee", max.gee );
  errorRangeCheck( "max.ped", max.ped );
  errorRangeCheck( "min.info", min.info, min=0 );
  errorRangeCheck( "length.haplos", length.haplos );
  errorRangeCheck( "max.mating.types", max.mating.types );
  ##warning( "Range checking for 'max.mating.types' is _NOT_ realistic!" );

  ## Newer range checking
  monte <- as.integer( monte ); ## force to be an integer
  errorRangeCheck( "monte", monte, min=0 );

  errorRangeCheck( "cutoff.haplo", cutoff.haplo, min=0, max=1, IS.INTEGER=FALSE );  ##

  ####print( "got past here 5" );

  ## Oy-vay! pbat needs a bloody phe file!!! arghhhhhhhhhhhhhhhhhhhhhhHHHhhhhhhhhhhhh!
  if( phefile == "" ) {
    phefile <- "p2bat_empty.phe"
    ped <- read.ped( pedfile, sym=FALSE )
    write.phe( phefile, as.phe( data.frame(pid=ped$pid, id=ped$id, p2batGenerated=rep(1,nrow(ped) ) ) ) );
  }


  #-------------------------
  # now do the actual work -
  #-------------------------

  outfile <- file( commandfile, "w" );
  on.exit( close(outfile) );

  if( logfile!="" )
    writeCommand( "logfile", logfile, outfile=outfile );

  ## apparently this can only be written when specified!
  ## EDIT: No, the filthy beast needs to be before the pedfile??? What the hell???
  ##if( snppedfile )
  writeCommand( "snppedfile", as.integer(snppedfile), outfile=outfile );

  ####print( "TESTING 1" );
  ##writeCommand( "pedfile", pedfile, outfile=outfile);     # (1)
  if( pedfile.ext=="ped" ) {  # (1)
    writeCommand( "pedfile", pedfile, outfile=outfile );
  }else if( pedfile.ext=="pped" ) {
    writeCommand( "ppedfile", pedfile, outfile=outfile );
  }else if( pedfile.ext=="cped" ) {
    writeCommand( "cnv", paste(cnv.intensity,cnv.intensity.num), outfile=outfile ); ## NEW, ARGHHHHHHH...
    writeCommand( "cpedfile", pedfile, outfile=outfile );
    writeCommand( "shortoutput", "1", outfile=outfile ) ## cnv suffering
  }else{
    stop( paste( "The pedigree file '", pedfile, "' has the extension '", pedfile.ext, "', which is not supported (must be 'ped' or 'pped'." ) );
  }
  ####print( "TESTING 2" );
  if( phefile != "" )
    writeCommand( "phenofile", phefile, outfile=outfile );  # (3)

  #if( snps!="" )
  writeCommand( "snps", c(snps), outfile=outfile, end=TRUE ); # (2)

  writeCommand( "censor", c(censor), outfile=outfile, end=TRUE ); # (4)
  if( time=="" && fbat=="logrank" ) { # (5)
    stop( "time-to-onset variable is required for pbat-logrank." );
  }else{
    writeCommand( "phenos", c(time), outfile=outfile, end=TRUE );
  }
  writeCommand( "phenos", phenos, outfile=outfile, end=TRUE );

  if( !is.null(preds) && preds[1]!="" ) { # (6)   ## 01/27/2006
    ##if( length(preds)!=length(preds.order) ) {
    ##  warning("'preds' and 'preds.order' must be of the same length. This information will be ignored.");
    ##}else{
    ##  writeCommand( "preds", pasteVector(c(preds,preds.order,"end")) );
    ##}

    ## 01/18/2006 bugfix - alteration in the pbat syntax?
    if( (length(preds)!=length(preds.order)) || preds.order[1]=="" ) {
      warning("'preds' and 'preds.order' are not of the same length; all order's will be forced to 1.'");
      preds.order = rep(1,length(preds));
    }
    writeCommand( "preds", pasteVector( c( preds, "end", preds.order, "end" ) ), outfile=outfile );
  }

  ## 01/18/2006 bugfix - didn't have end
  writeCommand( "inters", inters, outfile=outfile, end=TRUE );      # (7)

  if( groups.var!="" ){                    # (8)
    if( is.null(groups) || groups[1]=="" ) {
      ## Then we have to extract the groups values from the file!
      junk.phe <- read.phe( phefile, sym=FALSE ); ## WOW - difficult bug to find after adding the sym option
      groups <- unique( junk.phe[[groups.var]] );
    }

    writeCommand( "groups", c(groups.var, "end", groups, "end"), outfile=outfile );
  }

  writeCommandStrMatch1( "fbat", fbat, c("gee","pc","logrank"), outfile=outfile );

  if( censor=="" && fbat=="logrank" ) stop( "Need a censoring variable." );

  # (10-11)
  if( fbat!="logrank" ) {
    ## 01/25/2005 - amazing - this bug only shows up in multiprocessing mode
    writeCommand( "max", max.pheno, outfile=outfile );
    writeCommand( "min", min.pheno, outfile=outfile );
  }

  writeCommandStrMatch( "null", null, c("no linkage, no association", "linkage, no association"),
                        outfile=outfile, vals=c(1,2) );  ## 01/25/2006
  writeCommand( "alpha", alpha, outfile=outfile );        # (13)

  writeCommandStrMatch( "transpheno", trans.pheno, c("none","ranks","normal score"), outfile=outfile );
  writeCommandStrMatch( "transpred",  trans.pred , c("none","ranks","normal score"), outfile=outfile );
  writeCommandStrMatch( "transinter", trans.inter, c("none","ranks","normal score"), outfile=outfile );

  if( fbat!="logrank" ) {
    writeCommandStrMatch( "scanpred", scan.pred, c("all","subsets"), outfile=outfile );
    writeCommandStrMatch( "scaninter", scan.inter, c("all","subsets"), outfile=outfile );
  }

  writeCommandStrMatch( "scangenetic", scan.genetic,
                        c("additive","dominant","recessive",
                          "heterozygous advantage","all"), outfile=outfile );

  if( offset!="default" & offset!="" )
    writeCommandStrMatch("offset", offset, c("none","max power","gee + marker score","gee"), outfile=outfile );

  writeCommandStrMatch1( "screening", screening, c("conditional power","wald"), outfile=outfile );
  #warning( "replace with is.factor??" );
  ## AHHHHH... he changed this!!! Wow, it's backwards?
  ##writeCommandStrMatch( "distribution", distribution, c("continuous","categorical"), outfile=outfile );

  #warning( "I commented out logfile." );
  #writeCommand( "logfile", logfile );  # (23)

  # (24) NA
  if( fbat=="gee" )
    writeCommand( "maxgee", max.gee, outfile=outfile );

  writeCommand( "maxped", max.ped, outfile=outfile );  # (25)
  writeCommand( "mininfo", min.info, outfile=outfile );  # (26)

  if( is.list(haplos) ) {
    if( length(haplos)==0 ){
      ## Pure symbolic
      writeCommand( 'haplos 1 end', "", outfile=outfile );
    }else{
      ## Not pure symbolic
      vec <- c(  as.character(length(haplos))  );
      for( colnum in 1:length(haplos) ) {
        for( rownum in 1:length(haplos[[colnum]]) ) {
          vec <- c(vec, as.character(haplos[[colnum]][rownum]));
        }
        vec <- c(vec,"end");
      }
      writeCommand( "haplos", vec, outfile=outfile );
    }
  }

  writeCommandStrMatch( "ambhaplos", incl.ambhaplos, c(TRUE,FALSE), outfile=outfile ); #backwards
  writeCommandStrMatch( "infermissnp", infer.mis.snp, c(FALSE,TRUE), outfile=outfile );

  if( sub.haplos==TRUE ) {
    writeCommandStrMatch( "subhaplos", sub.haplos, c(FALSE,TRUE), outfile=outfile );
    writeCommand( "lengthhaplos", length.haplos, outfile=outfile );  # (31)
    writeCommandStrMatch( "adjsnps", adj.snps, c(FALSE,TRUE), outfile=outfile );
  }

  writeCommandStrMatch( "overallhaplo", overall.haplo, c(FALSE,TRUE), outfile=outfile );
  #writeCommandStrMatch( "cutoffhaplo", cutoff.haplo, c(FALSE,TRUE), outfile=outfile );
  writeCommand( "cutoffhaplo", as.numeric(cutoff.haplo), outfile=outfile ); # 12/29/06

  # (35-36)
  if( output=="short" )
    writeCommand( "shortoutput", "1", outfile=outfile );
  if( output=="detailed" )
    writeCommand( "detailedoutput", "1", outfile=outfile );

  writeCommand( "maxmatingtypes", max.mating.types, outfile=outfile );  # (37)

  if( fbat=="logrank" )
    writeCommand( "splus", "1", outfile=outfile ); # (38)

  if( !is.null(future.expansion) ) {
    for( i in 1:length(future.expansion) )
      writeLines( future.expansion[i], con=outfile );
  }

  ## 12/29/2006 additions

  ## monte carlo method -- 0 indicates not to use it
  writeCommand( "montecarloiteration", monte, outfile=outfile );

  ## MFBAT -- multi-phenotype multi-marker tests
  if( !is.null(mmaxsnps) || !is.null(mminsnps) || !is.null(mmaxphenos) || !is.null(mminphenos) ) {
    if( is.null(mminsnps) ) mminsnps <- 1;
    if( is.null(mmaxsnps) ) mmaxsnps <- mminsnps;
    if( is.null(mminphenos) ) mminphenos <- 1;
    if( is.null(mmaxphenos) ) mmaxphenos <- mminphenos;

    writeCommand( "MFBAT", "1", outfile=outfile );
    writeCommand( "mminsnps", mminsnps, outfile=outfile );
    writeCommand( "mmaxsnps", mmaxsnps, outfile=outfile );
    writeCommand( "mminphenos", mminphenos, outfile=outfile );
    writeCommand( "mmaxphenos", mmaxphenos, outfile=outfile );

    ## and debugging
    if( pbat.getmode()$mode != "single" )
      stop( "Multi-marker / multi-phenotype tests are not supported under any mode but 'single' at this time. Please use pbat.setmode('single') and rerun." );
  }

  ## environmental correlation adjust (GFBAT)
  if( !is.na(env.cor.adjust) && is.null(env.cor.adjust) )
    writeCommand( "GFBAT", as.integer(env.cor.adjust==TRUE), outfile=outfile );

  ## genome-wide acceleration
  if( gwa==TRUE )
    writeCommand( "gwa", 1, outfile=outfile );

  ## 04/22/2007
  if( new.ped.algo==TRUE )
    writeCommand( "newpedalgo", 1, outfile=outfile );

  ## hmm... this was changed a bit...
  writeCommandStrMatch( "distribution", distribution, c("default","jiang","murphy","naive","observed"), outfile=outfile );

  ##writeCommand( "channing", "1", outfile=outfile ) ## ahhhhhhhhhh!!!

  return( logfile );  # for future processing!
}
