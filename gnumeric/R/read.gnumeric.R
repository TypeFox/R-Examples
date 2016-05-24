### install.packages("XML", repos = "http://www.omegahat.org/R")
library(XML);

## sheet.name defaults to 'Sheet1' instead of NA, which would get the
## 'current sheet' from the .gnumeric file.
##
## Requirements:
##
## ssconvert (spreadsheet converter program, comes with gnumeric)
##
##            We use its undocumented flag:
##            --export-range='Sheet1!A1:Z100' or
##            --export-range='A1:Z100' (read whichever sheet was 'current'
##            in gnumeric when the file was saved)
##
##            The sheet name may require extra quoting like: --export-range='"Sheet1"!A1:Z100'
##
##            Wed Feb 9 22:40:12 2011: ssconvert for gnumeric 1.10.13
##            does not accept --export-range='"Sheet1"!A1:Z100', but
##            with double quotes outside and single quotes inside like
##            --export-range="'Sheet1'!A1:Z100" appears to work
##
##            But: (on windows) quotes around the filename
##            ("file.gnumeric") have to be double quotes. Maybe:
##            windows does not understand single quotes.
##
##            Maybe we should use the text exporter: it has a sheet option.
##
## quiet=TRUE uses '2>/dev/null' redirection (when .Platform$OS.type=="unix" )
##
##
## ??? Character encoding issues? gnumeric appears to use utf-8 (at
## least under an utf-8 locale)
##
read.gnumeric.sheet <-
  function(file,       ## 'filename.gnumeric'
           head=FALSE, ## as in read.csv, but first row seen by
                       ## read.csv is decided by row in top.left
           sheet.name='Sheet1', ## name of the sheet as appears in
                                ## gnumeric
           top.left='A1',       ## top left cell to request from
                                ## ssconvert (gnumeric utility)
           bottom.right=NA,    ## bottom right cell. The default
                               ## reads all rows, which is
                               ## slow. Speed up by giving a more
                               ## accurate upper bound, or
                               ## accurate value.
           drop.empty.rows="bottom", ## Drop rows containing only NAs
                                     ## and empty strings from result
           drop.empty.columns="right",## Drop columns containing only
                                      ## NAs and empty strings from
                                      ## result
           colnames.as.sheet=FALSE,   ## When TRUE, names(result) are
                                      ## set to those of
                                      ## gnumeric. (Overwrites names
                                      ## from read.csv)
           rownames.as.sheet=colnames.as.sheet, ## When TRUE, rownames
                                                ## are set to gnumeric
                                                ## row
                                                ## indices. (head=TRUE
                                                ## is taken into
                                                ## account, but e.g. skip!=0 not)
           quiet=TRUE,                ## Redirect stderr of ssconvert
                                      ## to /dev/null and do not print
                                      ## command executed
           LANG='C',                  ## Environment for ssconvert under unix
           locale='C',                ## value for: ssconvert -O "locale=C"
           import.encoding=NA,        ## --import-encoding for ssconvert
           field.format='automatic',  ## For: ssconvert -O "format=value".
                                      ## Options: raw automatic preserve
           ... ## passed to read.csv
           )
{

  ## 'all': drop even those between row (or columns) contining data
  if ( ! drop.empty.rows %in% c('none','top','bottom','both', 'all' ) ){
    stop( "drop.empty.rows is not in c('none','top','bottom','both', 'all' )" )
  }
  if ( ! drop.empty.columns %in% c('none','left','right','both', 'all' ) ){
    stop( "drop.empty.columns is not in c('none','left','right','both', 'all' )" )
  }
  if ( ! field.format %in% c('raw', 'automatic', 'preserve') ){
    stop( "field.format is not in c('raw', 'automatic', 'preserve')" );
  }


  ### build command
  SHEET='';
  if ( !is.na(sheet.name) ){
    ### bug: Should check for "'" inside sheet.name
    sheet.name.with.quotes = paste(sep='',  "'", sheet.name, "'" );
    SHEET=paste(sheet.name.with.quotes, sep='' );
  }


  ssconvert = "ssconvert";
  ssconvert.full.path=Sys.which( ssconvert );
  ## Wed Sep  1 08:11:00 2010 : On Solaris Sys.which returns a message instead of "",
  ## thus we need to check the value returned does exist.
  if ( ( ssconvert.full.path == "") || ! file.exists( ssconvert.full.path ) ){
    stop("Required program '",ssconvert,"' not found." );
  }

  IMPORT.ENCODING='';
  if ( ! is.na( import.encoding ) ){
    IMPORT.ENCODING = paste(" --import-encoding='", import.encoding, "' ", sep='' );
  }
  
  ## --export-range needed because I know of no other way to select the
  ## sheet. This in turn forces to also provide top.left and
  ## bottom.right, even when we just want 'all the sheet'
  #  cmd <- paste(ssconvert.full.path,
  #               " --export-type=Gnumeric_stf:stf_csv ",
  #               " --export-range=", '"', SHEET , top.left ,":", bottom.right, '" ',
  #               IMPORT.ENCODING,
  #               ' "', file, '"',
  #               " fd://1 ", sep='');

  ## With --export-type=Gnumeric_stf:stf_assistant we can select the
  ## sheet without using --export-range. 
  
  if (is.na(bottom.right) && (is.na(top.left)|| top.left=="A1" ) ){
     range="";
  } else {
    if ( is.na(top.left) ){
      top.left = 'A1';
    }
    if ( is.na( bottom.right ) ){
      bottom.right = 'IV65536';
    }
    range=paste(sep='', " --export-range=", '"', SHEET , "!",top.left ,":", bottom.right, '" ');
  }

  cmd <- paste(ssconvert.full.path,
               " --export-type=Gnumeric_stf:stf_assistant ",
               ' -O "locale=',locale,
                   ' format=',field.format,
                   ' separator=, eol=unix sheet=',
               sheet.name.with.quotes ,'"',
               range , 
               IMPORT.ENCODING,
               ' "', file, '"',
               " fd://1 ", sep='');

  if ( .Platform$OS.type == "unix" ){
    ## force decimal point in ssconvert output under e.g. hungarian locale
    cmd=paste( "LANG=",LANG," ", cmd, sep='' );
  }
  
  if ( ! quiet ){
    cat(cmd,"\n");
  } else {
    ## redirect stderr of cmd to /dev/null
    if ( .Platform$OS.type == "unix" ){
      ## the 'grep ,' filter is a temporary workaround for ods
      ## to ignore diagnostic messages (may be removed in 2010)
      cmd = paste( cmd, " 2> /dev/null " ); ## unix
    } else {
      ## ( .Platform$OS.type == "windows" )
      ## ??? 2>NIL:  ??? or similar? 
    }
  }

  ### read data
  x=read.csv( pipe( cmd  ) , header=head, encoding="UTF-8", ... )

  ### optionally rename columns and rows to correspond to gnumeric
  ### cell names.
  if ( colnames.as.sheet || rownames.as.sheet ){
    ABC=LETTERS;
    ## COLNAMES: A .. AA, AB, .. IV
    COLNAMES= as.vector( t(outer(c('',ABC[1:9]), ABC, paste, sep='')))[1:256]; 

    left=''
    i=1;
    while ( substr( top.left,i,i) %in% ABC ){
      i=i+1;
    }
    left = substr( top.left, 1, (i-1) )
    top.str = substr( top.left, i, nchar(top.left) );
    top=as.integer( top.str );

    if ( colnames.as.sheet ){
      left.index = match( left, COLNAMES, nomatch=0 )
      stopifnot( left.index > 0 )
      names(x)=COLNAMES[left.index:(left.index+length(x)-1)];
    }

    if ( rownames.as.sheet ){

      if ( head ){
        top = top+1;
      }
      
      rownames(x) = top:(top+length(x[[1]])-1)
    }
  }

  ### optionally drop empty columns and rows 
  if ( drop.empty.columns != "none" || drop.empty.rows!="none" ){
    ## drop empty columns and rows
    last.col=length(x);
    last.row=length(x[[1]]);

    is.empty <- function(x){ all( is.na(x) | x=='' )  }
    m=as.matrix(x);

    i=1:last.row; ## indices of rows to keep

    if ( drop.empty.rows != "none" ){
      bi=!apply(m,1,is.empty);
      i=i[bi]; ## indices of non-empty rows
      
      if ( drop.empty.rows == "bottom" ){
        ## extend i to include all rows from row 1 to last non-empty
        ## (includes empty rows between non-empty ones as well)
        i = 1:i[length(i)]; 
      } else if ( drop.empty.rows == "top" ){
        i = i[1]:last.row
      } else if ( drop.empty.rows == "both" ){
        i = i[1]:i[length(i)]
      }
    }


    j=1:last.col;
    if ( drop.empty.columns  != "none"){
      bj=!apply(m,2,is.empty)
      j=j[bj]
      if ( drop.empty.columns == "right" ){
        j = 1:j[length(j)]
      } else if ( drop.empty.columns == "left" ){
        j = j[1]:last.col
      } else if ( drop.empty.columns == "both" ){
        j = j[1]:j[length(j)]
      }
    }
    x=x[i,j]
  }
  
  x
}



## Like read.gnumeric.sheet, but bottom.right is mandatory and we drop
## no rows or columns by default.
read.gnumeric.range <-
  function(file,       
           head=FALSE, 
           sheet.name='Sheet1', 
           top.left='A1',
           bottom.right, 
           drop.empty.rows="none", 
           drop.empty.columns="none",
           colnames.as.sheet=FALSE,
           rownames.as.sheet=colnames.as.sheet,
           quiet=TRUE,
           LANG='C',
           locale='C',
           import.encoding=NA,
           field.format='automatic',
           ... ## passed to read.csv
           )
{
  read.gnumeric.sheet(file,
                      head,
                      sheet.name,
                      top.left,
                      bottom.right,
                      drop.empty.rows,
                      drop.empty.columns,
                      colnames.as.sheet, 
                      rownames.as.sheet,
                      quiet,
                      LANG=LANG,
                      locale=locale,
                      import.encoding=import.encoding,
                      field.format=field.format,
                      ... );
}

## Note: read.gnumeric.sheet.names works, but is superfluous since
## we have read.gnumeric.sheet.info
##
## return sheet names from a .gnumeric file
## read.gnumeric.sheet.names <- function(filename){
##  doc = xmlTreeParse(filename);
##  x1 <- doc$doc$children$Workbook['SheetNameIndex'][[1]];
##  x2=xmlSApply( x1, xmlValue );
##  names(x2)=NULL;
##  x2
##}


## return a data.frame with columns sheet.name (string), width
## (integer, may be zero), height (integer, may be zero) bottom
## (string or NA) to be passed to read.gnumeric.sheet if not
## NA. bottom is NA if width or height is zero, i.e. when the sheet is
## empty.
read.gnumeric.sheet.info <- function(file){
  doc = xmlTreeParse(file);
  x1 <- doc$doc$children$Workbook['Sheets'][[1]];
  names=c();
  widths=c();
  heights=c();

  bottoms=c();
  ABC=LETTERS;
  ## COLNAMES: A .. AA, AB, .. IV
  COLNAMES= as.vector( t(outer(c('',ABC[1:9]), ABC, paste, sep='')))[1:256]; 

  for ( i in 1:length(x1) ){
    xx <- x1[i];
    width  <-  1+as.integer(xmlValue( xx$Sheet['MaxCol']$MaxCol[[1]] )) ;
    height <-  1+as.integer(xmlValue( xx$Sheet['MaxRow']$MaxRow[[1]] )) ;
    name   <-  xmlValue( xx$Sheet['Name']$Name[[1]]     ) ;
    bottom <- if ( width==0 || height==0 ){
      NA
    } else {
      paste( COLNAMES[width], height , sep='' );
    }
    
    names=c(names,name);
    widths=c(widths,width);
    heights=c(heights,height);
    bottoms=c(bottoms,bottom);
  }
  data.frame(sheet.name=names,
             width=widths,
             height=heights,
             bottom.right=bottoms,
             stringsAsFactors=FALSE );
}

## Read nonempty sheets from a .gnumeric file. Returns a list of
## data.frames. Names in the list are the sheet names, values are read
## using read.gnumeric.sheet
read.gnumeric.sheets <- function(file,
                                 head=FALSE,
                                 drop.empty.rows="none",
                                 drop.empty.columns="none",
                                 colnames.as.sheet=FALSE,
                                 rownames.as.sheet=colnames.as.sheet,
                                 quiet=TRUE,
                                 LANG='C',                  
                                 locale='C',                
                                 import.encoding=NA,        
                                 field.format='automatic',  
                                 ...  ## passed to read.csv
                                 ){
  si <- read.gnumeric.sheet.info( file );
  res=list();
  for ( i in 1:nrow(si) ){
    if ( !is.na(si[i,'bottom.right']) ){
      x=read.gnumeric.sheet(file=file,head=head,
        sheet.name=si[i,'sheet.name'],
        bottom.right=NA, #si[i,'bottom.right'],
        drop.empty.rows=drop.empty.rows,
        drop.empty.columns=drop.empty.columns,
        colnames.as.sheet=colnames.as.sheet,
        rownames.as.sheet=rownames.as.sheet,
        quiet=quiet,
        LANG=LANG,
        locale=locale,
        import.encoding=import.encoding,
        field.format=field.format,
        ...
        );
      xx=list(x); names(xx)=si[i,'sheet.name'];
      res=c(res,xx);
    }
  }
  res;
}

### ### gnumeric.raw.date.epoch may vary depending on a code in .gnumeric
### ### files. Possibly also on the input file type for non .gnumeric files.
### ### For variations on good values for 'epoch' see:
### ### help(as.Date)  (search for 'Excel')
### gnumeric.raw.date.as.date <- function(x,epoch='1899-12-30'){
###   v=as.numeric(as.character(x));
###   v+as.Date(epoch);
### }
### 
### gnumeric.raw.datetime.as.datetime <- function(x, epoch="1899-12-29 23:59:59", tz='UTC'){
###   as.POSIXct(as.numeric(as.character(x))*(60*60*24),  origin=epoch, tz=tz)
### }


