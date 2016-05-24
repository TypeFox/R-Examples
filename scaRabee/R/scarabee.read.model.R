
#Copyright (c) 2009-2014 Sebastien Bihorel
#All rights reserved.
#
#This file is part of scaRabee.
#
#    scaRabee is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    scaRabee is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with scaRabee.  If not, see <http://www.gnu.org/licenses/>.
#

scarabee.read.model <- function(files=NULL, runtype=NULL){
  
  ## Checks inputs
  if (is.null(files))
    stop('files argument is NULL.')
  
  if (is.null(files$model))
    stop('files argument does not have any model level.')
  
  if (!file.exists(files$model))
    stop('model file does not exist.')
  
  if (is.null(runtype))
    stop('runtype argument is NULL.')
  
  # Read model file
  cat('\nProcessing model file:\n')
  text <- scan(file=files$model,what='character',sep='\n',quiet=TRUE)
  
  # Strip comments
  text <- gsub('#.*','',text)
  
  # Check for tags
  if (!any(grepl('\\$',text)))
    stop('model file could not be read because of the lack of tags.')
  
  # Check for tags followed by text
  if (!any(grepl('\\$',text)))
    stop('model file could not be read because of the lack of tags.')
  
  # Add tags at the beginning and end of the files
  text <- c('$START', text, '$END')
  
  # Process text to put each $ block code on one element of text
  text <- gsub('^\\$','&^?*$',text)
  text <- paste(text,collapse='%#!@')
  text <- unlist(strsplit(text,'&^?*',fixed=TRUE))
  text <- gsub('%#!@',';',text)[-c(1,2,length(text))]
  
  # Check for $ANALYSIS tag(s)
  if (!any(grepl('^\\$ANA',text)))
    stop('model file must contain an $ANALYSIS block.')
  if (sum(grepl('^\\$ANA',text))>1)
    stop('model file must contain a single $ANALYSIS tag.')
  
  # Extract $ANALYSIS tag
  anatext <- text[grepl('^\\$ANA',text)]
  text <- text[-grepl('^\\$ANA',text)]
  
  # Check tags
  if (any(grepl('^\\$\\w+ [ ]*\\w',text)))
    stop('model file must not contain lines with code following the block tag.')
  
  if (!any(grepl('^\\$\\w+ [ ]*\\w',anatext)))
    stop('$ANALYSIS tag must be followed by a name for the analysis.')
  
  # Get analysis name
  analysis.name <- unlist(strsplit(anatext,'[ ;]'))[2]
  
  # Identify type of run
  hasoutput <- explicit <- ode <- dde <- FALSE
  if (any(grepl('^\\$OUT',text))) hasoutput <- TRUE
  if (any(grepl('^\\$ODE',text))) ode <- TRUE
  if (any(grepl('^\\$DDE',text))) dde <- TRUE
  
  # Check combinations of hasoutput, ode, and dde
  if (all(c(hasoutput,ode,dde)==FALSE))
    stop('model file does not contain any $OUTPUT, $ODE, or $DDE tags.')
  
  if (ode & dde)
    stop('model file must contain either a $ODE or a $DDE block.')
  
  if ((ode & !hasoutput) | (dde & !hasoutput))
    stop('model file must contain a $OUTPUT block.')
  
  if (sum(grepl('^\\$OUT',text))>1 | sum(grepl('^\\$ODE',text))>1 | 
      sum(grepl('^\\$DDE',text))>1)
    stop('model file must contain a single $OUTPUT, $ODE, or $DDE tags.')
  
  if (!(ode | dde)) explicit <- TRUE
  
  if (explicit){
    cat('  Model with explicit solution detected.\n\n')
  } else if (ode){
    cat('  Model using ordinary differential equations detected.\n\n')
  } else {
    cat('  Model using delay differential equations detected.\n\n')
  }
  
  # Check block requirements for ode or dde model
  if (ode | dde){
    if (!any(grepl('^\\$DER',text)))
      stop('model file must contain a $DERIVED block.')
    if (!any(grepl('^\\$IC',text)))
      stop('model file must contain a $IC block.')
    if (!any(grepl('^\\$SCA',text)))
      stop('model file must contain a $SCALE block.')
    if (sum(grepl('^\\$DER',text))>1 | sum(grepl('\\$IC',text))>1 | 
        sum(grepl('^\\$SCA',text))>1)
      stop(paste('model file must contain a single $DERIVED, $IC, and ',
                 '$SCALE tags.',sep=''))
  }
  if (dde){
    if (!any(grepl('^\\$LAG',text)))
      stop('model file must contain a $LAGS block.')
    if (sum(grepl('^\\$LAG',text))>1)
      stop('model file must contain a single $LAGS tags.')
  }
  
  # Check block requirements for any model
  if (!runtype=='simulation'){
    if (!any(grepl('^\\$VAR',text)))
      stop('model file must contain a $VARIANCE block.')
  }
  
  # Process code blocks
  code <- list(template='explicit.model')
  if (ode | dde){
    # Parameter block
    code$deriv <- text[grep('^\\$DER',text)]
    code$deriv <- sub('\\$DER[A-Za-z]*;','',code$deriv)
    
    # Initial condition block
    code$ic <- text[grep('^\\$IC',text)]
    code$ic <- sub('\\$IC[A-Za-z]*;','',code$ic)
    
    # Scale block
    code$scale <- text[grep('^\\$SCA',text)]
    code$scale <- sub('\\$SCA[A-Za-z]*;','',code$scale)
    
    # Differential equation block
    if (ode){
      code$template <- 'ode.model'
      code$de <- text[grep('^\\$ODE',text)]
      code$de <- sub('\\$ODE[A-Za-z]*;','',code$de)
    } else {
      code$template <- 'dde.model'
      code$de <- text[grep('^\\$DDE',text)]
      code$de <- sub('\\$DDE[A-Za-z]*;','',code$de)
    }
    
    # Lag block
    if (dde){
      code$lags <- text[grep('^\\$LAG',text)]
      code$lags <- sub('\\$LAG[A-Za-z]*;','',code$lags)
    }
  }
  
  # output block
  code$output <- text[grep('^\\$OUT',text)]
  code$output <- sub('\\$OUT[A-Za-z]*;','',code$output)
  
  # variance block
  if (!runtype=='simulation'){
    code$var <- text[grep('^\\$VAR',text)]
    code$var <- sub('\\$VAR[A-Za-z]*;','',code$var)
  }
  
  # secondary parameters block
  if (any(grepl('^\\$SEC',text))){
    if (sum(grepl('^\\$SEC',text))>1)
      stop('model file must contain a single $SECONDARY tag.')
    
    code$sec <- text[grep('^\\$SEC',text)]
    code$sec <- sub('\\$SEC[A-Za-z]*;','',code$sec)
  }
  
  # Remove ; in elements of code
  code <- lapply(code,function(x)gsub(';','\n',x))
  
  # check for init, scale, dadt, and f objects
  if (ode | dde){
    if (!any(grepl('^[ ]*init[ ]*<-[A-Za-z0-9_({[ .]',
                   unlist(strsplit(code$ic,'\n')))))
      stop('init object should be created in the $IC block.')
    
    if (!any(grepl('^[ ]*scale[ ]*<-[A-Za-z0-9_({[ .]',
                   unlist(strsplit(code$scale,'\n')))))
      stop('scale object should be created in the $SCALE block.')
    
    if (!any(grepl('^[ ]*dadt[ ]*<-[A-Za-z0-9_({[ .]',
                   unlist(strsplit(code$de,'\n')))))
      stop(sprintf('dadt object should be created in the $%s block.',
                   ifelse(ode,'ODE','DDE')))
  }
  if (!any(grepl('^[ ]*y[ ]*<-[A-Za-z0-9_({[ .]',
                 unlist(strsplit(code$output,'\n')))))
    stop('y object should be created in the $OUTPUT block.')
  
  if (!runtype=='simulation'){
    if (!any(grepl('^[ ]*v[ ]*<-[A-Za-z0-9_({[ .]',
                   unlist(strsplit(code$var,'\n')))))
      stop('v object should be created in the $VARIANCE block.')
  }
  
  return(c(name=analysis.name,code))
  
}
