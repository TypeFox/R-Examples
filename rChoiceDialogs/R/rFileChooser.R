
#' Get default set of filters
#'
#' @name getDefaultFilters
#' @return Default set of filters for jchoose.files
#' @seealso {\code{\link{jchoose.files}}, \code{\link{tkchoose.files}}, \code{\link{rchoose.files}}}
#' @export getDefaultFilters
#' @author  Alex Lisovich, Roger Day

getDefaultFilters<-function(){
  res<-rbind(
	c("R or S files (*.R,*.q,*.ssc,*.S)", "*.R;*.q;*.ssc;*.S"),
	c("Enhanced metafiles (*.emf)",       "*.emf"),            
	c("Postscript files (*.ps)",          "*.ps"),             
	c("PDF files (*.pdf)",                "*.pdf"),            
	c("Png files (*.png)",                "*.png"),            
	c("Windows bitmap files (*.bmp)",     "*.bmp"),            
	c("Jpeg files (*.jpeg,*.jpg)",        "*.jpeg;*.jpg"),     
	c("Text files (*.txt)",              "*.txt"),            
	c("R images (*.RData,*.rda)",         "*.RData;*.rda"),    
	c("Zip files (*.zip)",                "*.zip"),           
	c("All files (*.*)",                 "*.*")
  );             
  return(as.matrix(res));
}


#' Choose a list of files interactively using tcltk
#'
#' Provides the same functionality as choose.files from utils package for Windows,
#' but relies on tcltk package and therefore is system independent provided tcltk is installed.
#'
#'
#' @note tkchoose.files() is called internally by rchoose.files() if it's appropriate for a given platform/graphics combination.
#' Calling tkchoose.files() directly forces the package to use tcl tk based dialog regardless of system capabilities and therefore may fail.
#' Use the direct call to tkchoose.files() only if it seems beneficial to bypass the rchoose.files() decision logic.
#'
#' @name tkchoose.files
#' @param default Which filename to show initially
#' @param caption The caption on the file selection dialog
#' @param multi Whether to allow multiple files to be selected
#' @param filters A matrix of filename filters. If NULL, all files are shown.
#' Default is filters=getDefaultFilters().
#' @param index Which row of filters to use by default.
#' @return A character vector giving zero or more file paths.  If user cancels operation, character(0) is returned.
#'
#' @seealso {\code{\link{getDefaultFilters}}, \code{\link{rchoose.files}}}
#' @examples
#' \dontrun{
#' tkchoose.files();
#' }
#' @export tkchoose.files
#' @author  Alex Lisovich, Roger Day


tkchoose.files<-function(default = "", caption = "Select files",
             multi = TRUE, filters = getDefaultFilters(),
             index = nrow(filters)){
	if (is.null(filters))
		filters=c("All","*.*");

    	if (inherits(filters, "character")) 
		filters <- matrix(filters, ncol = 2);
    	if (!inherits(filters, "matrix") || mode(filters) != "character" || 
		ncol(filters) != 2) 
		stop("'filters' must be a n*2 matrix of characters!")

	filters[,2]<-gsub(";",",",filters[,2],fixed=TRUE);
	for (i in 1:nrow(filters))		
		filters[i,1]<-gsub(paste("(",filters[i,2],")",sep=""),"",filters[i,1],fixed=TRUE);

    	filters <- paste("{\"", filters[, 1], "\" {\"", gsub(",", 
		"\" \"", filters[, 2]), "\"}}", sep = "", collapse = " ");


	tclres<-tcltk::tclvalue(tcltk::tkgetOpenFile(title=caption,initialdir = default,initialfile = "",
										multiple=multi,filetypes=filters));


	res<-character(0);	
	if(nzchar(tclres)){
		if (multi){
			starts<-gregexpr("{",tclres,fixed=TRUE)[[1]]+1;
			stops<-gregexpr("}",tclres,fixed=TRUE)[[1]]-1;
			if(starts[1]<1)
				return(tclres);
			res<-NULL;
			for (i in 1:length(starts)){
				res<-c(res,as.character(substr(tclres,starts[i],stops[i])));
			}
			
		} else
			 res<-tclres;
	}
	return(res);
}

#' Choose a list of files interactively using rJava
#'
#' Provides the same functionality as choose.files from utils package,
#' but relies on Java and rJava package and therefore is system independent provided Java 1.5 and higher is installed.
#'
#' @note jchoose.files() is called internally by rchoose.files() if it's appropriate for a given platform/graphics combination.
#' Calling jchoose.files() directly forces the package to use Java based dialog regardless of system capabilities and therefore may fail.
#' Use the direct call to jchoose.files() only if it seems beneficial to bypass the rchoose.files() decision logic.
#'
#' @name jchoose.files
#' @param default Which filename or directory to show initially. Default is current work directory.
#' @param caption The caption on the file selection dialog
#' @param multi Whether to allow multiple files to be selected
#' @param filters A matrix of filename filters. If NULL, all files are shown.
#' Default is filters=getDefaultFilters().
#' @param index Which row of filters to use by default.
#' @param modal Indicates how the modality of the dialog is implemented.
#' If TRUE, the modal dialog is used and if FALSE, R repeatedly checks for dialog status (active or not).
#' The latter is used to refresh R Gui window  on Windows. Default is canUseJavaModal().
#' @return A character vector giving zero or more file paths. If user cancels operation, character(0) is returned.
#'
#' @seealso {\code{\link{getDefaultFilters}}, \code{\link{rchoose.files}}, \code{\link{canUseJavaModal}}}
#' @examples
#' \dontrun{
#' jchoose.files();
#' }
#' @export jchoose.files
#' @author  Alex Lisovich, Roger Day

jchoose.files<-function(default = getwd(), caption = "Select files",
             multi = TRUE, filters = getDefaultFilters(),
             index = nrow(filters),modal=canUseJavaModal()){



    if (inherits(filters, "character")) {
        filters <- matrix(filters, ncol = 2)
	  index=nrow(filters);
    }
    if (!inherits(filters, "matrix") || mode(filters) != "character" || 
        ncol(filters) != 2) 
        stop("'filters' must be a n*2 matrix of characters!")

	chooser<-new(J("rjavautils.rJavaFileChooser"));

	chooser$setDialogTitle(caption);
	chooser$setMultiSelectionEnabled(multi);
	chooser$setCurrentDirectory(.jnew("java.io.File",default));

	if(!is.null(filters)){
		active_filter<-NULL;
		chooser$setAcceptAllFileFilterUsed(FALSE);
		for (i in 1:nrow(filters)){
			extensions<-unlist(strsplit(filters[i,2],split=";"));
			extensions<-gsub("*.","",extensions,fixed=TRUE);
			filter<-new(J("rjavautils.ExtensionFileFilter"),filters[i,1],.jarray(extensions));
			chooser$addChoosableFileFilter(filter);
			if(index==i)
				active_filter<-filter;
		}
		if(!is.null(active_filter))
			chooser$setFileFilter(active_filter);
	}



	resVal<-chooser$showOpenDialog(NULL,modal);
	if(!modal){
		while (chooser$isRunning()){};
		chooser$dispose();
	}

	if(chooser$selectionApproved()){
		if(multi) {
			files<-as.list(chooser$getSelectedFiles());
			fnames<-unlist(lapply(files, function(file) {file$getCanonicalPath()}));
		} else {
			fnames<-chooser$getSelectedFile()$getCanonicalPath();
		}
	} else {
		fnames<-character(0);
	}



	return(fnames);
}

#' Choose a list of files interactively using the command line
#'
#' Allows to choose files or directories using the using the command line based interaction
#' providing the same functionality as GUI counterparts without using any graphical framework.
#'
#' @note cmdchoose.files() is called internally by rchoose.files() if neither Java nor TclTk are available
#' Calling cmdchoose.files() directly forces the package to use command line interaction regardless of system capabilities and therefore may fail.
#' Use the direct call to cmdchoose.files() only if it seems beneficial to bypass the rchoose.files() decision logic.
#'
#' @name cmdchoose.files
#' @param default Which filename or directory to show initially. Default is current work directory.
#' @param caption The caption on the file selection dialog
#' @param multi Whether to allow multiple files to be selected
#' @param dir.only If TRUE (default is FALSE) works as directory chooser.
#' @param filters A matrix of filename filters. If NULL, all files are shown.
#' Default is filters=getDefaultFilters().
#' @param index Which row of filters to use by default.
#' @return A character vector giving zero or more file paths. If user cancels operation, character(0) is returned.
#'
#' @seealso {\code{\link{getDefaultFilters}}, \code{\link{rchoose.files}}, \code{\link{canUseJavaModal}}}
#' @examples
#' \dontrun{
#' cmdchoose.files();
#' }
#' @export cmdchoose.files
#' @author  Alex Lisovich, Roger Day

cmdchoose.files = function(default = getwd(), caption = "Select files", multi = TRUE, dir.only=FALSE,
                           filters = getDefaultFilters(), index = nrow(filters))
{
  
  getChoice = function(prompt="",multi=TRUE){
    cat(prompt);
    
    done<-FALSE
    while(!done){
      s<-readline(": ");
      s<-unlist(strsplit(s," ",fixed=TRUE));
      suppressWarnings(s<-as.integer(s));
      if(length(s)>0 && sum(is.na(s))==0){
        if(!multi && length(s)>1){
          cat(gettext("!!!invalid input: multiple selections (should use multi=TRUE)\n\n"));          
        } else {
          done=TRUE;
        }
      } else {
        cat(gettext("!!!invalid input: not a numeric choice\n\n"));
      }
    }
    return(s);  
  }
  
  
  display.selection = function(path=getwd(),pattern=NULL,dir.only=FALSE){
    all.items<-dir(path=path,pattern=pattern);
    all.items.fullnames<-dir(path=path,pattern=pattern,full.names=TRUE);
    are.dirs<-file.info(all.items.fullnames)$isdir;
    
    dirs<-all.items[are.dirs];
    if(dir.only){
      files<-dirs;
    }else{
      files<-all.items[!are.dirs];
    }
    
    cat("===Look in:",path,"===\n");
    cat("1:\t../\n");
    if(length(dirs)>0){
      cat(paste(c(2:(length(dirs)+1)),":\t./",dirs,"\n",sep=""),sep="");
    }
    if(length(files)>0){
      cat(paste(c(2:(length(files)+1))+length(dirs),":\t",files,"\n",sep=""),sep="");
    }
    
    return(list(ndirs=length(dirs)+1,items=c("../",dirs,files)));
  } 
  
  if(is.matrix(filters)) {
    pattern = filters[index, 2]
    ### Convert from e.g. "*.txt;*.csv;*.tsv"  to "*.txt$|*.csv$|*.tsv$"
    pattern = paste(gsub(';', '$|', pattern), '$', sep="")
  }
  else if(is.vector(filters) & length(filters)==2 & mode(filters)=='character') {
    ### filters = c("CSV files (*.txt,*.csv,*.tsv)","*.txt;*.csv;*.tsv"));
    pattern = filters[2]
    ### Convert from e.g. "*.txt;*.csv;*.tsv"  to "*.txt$|*.csv$|*.tsv$"
    pattern = paste(gsub(';', '$|', pattern), '$', sep="")
  }
  else if(is.vector(filters) & length(filters)==1 & mode(filters)=='character') {
    pattern = filters  ###
    if(substring(pattern, nchar(pattern)) != '$')
      pattern = paste(gsub(';', '$|', pattern), '$', sep="")
  }
  else {
    warning("file extension pattern not recognized.")
    pattern = ".*"
  }
  
  cur.path<-default;
  prev.path<-"";
  res<-character(0);
  done<-FALSE;
  
  
  while(!done){
    if(cur.path!=prev.path){
      prev.path<-cur.path;
      cat("===",caption,"===\n");
      items<-display.selection(cur.path,pattern,dir.only);
    }  
    choices<-getChoice("\nEnter number(s) separated by spaces, or 0 to cancel",multi=multi);
    if(min(choices)<0 || max(choices)>length(items$items)){
      cat(gettext("!!!invalid input: selection out of range\n\n"));
      next;
    }
    #browser();
    numdirs<-sum(choices<=items$ndirs);
    if(numdirs>0){ #at least one directoriy selected
      if(sum(choices<=items$ndirs)>1){
        cat(gettext("!!!invalid input: more than one directory was selected\n\n"));
      } else if (length(choices)>1){
        cat(gettext("!!!invalid input: selection contains both directories and files\n\n"));
      } else if (choices[1]==1){ #go up one level
        prev.path<-cur.path;
        cur.path<-unlist(strsplit(cur.path,"/",fixed=TRUE));
        cur.path<-paste(cur.path[-length(cur.path)],collapse="/");        
      } else if (choices[1]>0) {#drill down
        prev.path<-cur.path;
        cur.path<-paste(cur.path,items$items[choices[1]],sep="/");
      } else {
        done=TRUE; #0 was entered
      }
    }else{ #files selected
      return(normalizePath(paste(cur.path,items$items[choices],sep="/")));
    }
  }#while
  return(res);
}



#' Choose a directory interactively using rJava
#'
#' Provides the same functionality as choose.dir from utils package,
#' but relies on Java and rJava package and therefore is system independent provided Java 1.5 and higher is installed.
#'
#'
#' @note jchoose.dir() is called internally by rchoose.dir() if it's appropriate for a given platform/graphics combination.
#' Calling jchoose.dir() directly forces the package to use Java based dialog regardless of system capabilities and therefore may fail.
#' Use the direct call to jchoose.dir() only if it seems beneficial to bypass the rchoose.dir() decision logic.
#'
#' @name jchoose.dir
#' @param default Which filename or directory to show initially. Default is current work directory.
#' @param caption The caption on the file selection dialog
#' @param modal Indicates how the modality of the dialog is implemented.
#' If TRUE, the modal dialog is used and if FALSE, R repeatedly checks for dialog status (active or not).
#' The latter is used to refresh R Gui window  on Windows. Default is canUseJavaModal().
#' @return A character vector giving zero or more file paths. If user cancels operation, character(0) is returned.
#'
#' @seealso {\code{\link{rchoose.dir}}, \code{\link{canUseJavaModal}}}
#' @examples
#' \dontrun{
#' jchoose.dir();
#' }
#' @export jchoose.dir
#' @author  Alex Lisovich, Roger Day


jchoose.dir<-function(default = getwd(), caption = "Select Directory",modal=canUseJavaModal()){

	chooser<-new(J("rjavautils.rJavaFileChooser"));

	chooser$setFileSelectionMode(J("javax.swing.JFileChooser")$DIRECTORIES_ONLY);
	chooser$setDialogTitle(caption);
	chooser$setCurrentDirectory(.jnew("java.io.File",default));


	resVal<-chooser$showOpenDialog(NULL,modal);
	if(!modal){
		while (chooser$isRunning()){};
		chooser$dispose();
	}

	if(chooser$selectionApproved()){
		dirname<-chooser$getSelectedFile()$getCanonicalPath();
	} else {
		dirname<-character(0);
	}



	return(dirname);
}


#' Choose a list of files interactively
#'
#' Provides the same functionality as choose.files from the utils package for Windows,
#' but is intended to be system independent.
#'
#' @name rchoose.files
#' @param default Which filename or directory to show initially. Default is current work directory.
#' @param caption The caption on the file selection dialog
#' @param multi Whether to allow multiple files to be selected
#' @param filters A matrix of filename filters. If NULL, all files are shown.
#' Default is filters=getDefaultFilters().
#' @param index Which row of filters to use by default.
#' @return A character vector giving zero or more file paths. If user cancels operation, character(0) is returned.
#'
#' @seealso {\code{\link{getDefaultFilters}}, \code{\link{jchoose.files}}, 
#' \code{\link{tkchoose.files}}, \code{\link{canUseJava}}, \code{\link{canUseTclTk}}}
#' @examples
#' \dontrun{
#' rchoose.files();
#' }
#' @export rchoose.files
#' @author  Alex Lisovich, Roger Day

rchoose.files<-function(default = getwd(), caption = "Select files",
             multi = TRUE, filters = getDefaultFilters(),
             index = nrow(filters)){

  if (!interactive()) 
    stop("rchoose.files() cannot be used non-interactively")
  
	if(canUseJava()){
		return(jchoose.files(default,caption,multi,filters,index));
	} else if(canUseTclTk()) {
		return(tkchoose.files(default,caption,multi,filters,index));
	} else {
	  return(cmdchoose.files(default,caption,multi,FALSE,filters,index));
	}

	return(character(0));
}



#' Choose directory interactively
#'
#' Provides the same functionality as choose.dir from the utils package for Windows,
#' but is intended to be system independent.
#'
#' @name rchoose.dir
#' @param default Which filename or directory to show initially. Default is current work directory.
#' @param caption The caption on the file selection dialog
#' @return A character vector giving zero or more file paths. If user cancels operation, character(0) is returned.
#' @seealso {\code{\link{getDefaultFilters}}, \code{\link{jchoose.files}}, 
#' \code{\link{tkchoose.files}}, \code{\link{canUseJava}}, \code{\link{canUseTclTk}}}
#' @examples
#' \dontrun{
#' rchoose.dir();
#' }
#' @export rchoose.dir
#' @author  Alex Lisovich, Roger Day


rchoose.dir<-function(default=getwd(), caption="Select Directory"){
  if (!interactive()) 
    stop("rchoose.files() cannot be used non-interactively")
  
  if(canUseJava()){
		return(jchoose.dir(default,caption));
	} else if(canUseTclTk()) {
		return(tcltk::tk_choose.dir(default,caption));
	} else {
	  return(cmdchoose.files(default,caption,TRUE,TRUE));
	}

	return(character(0));

}


