find_npp <- 
function(){
    exe <- 
    file.path( utils::readRegistry("SOFTWARE\\Notepad++", "HLM", view="32-bit")[[1]]
             , "Notepad++.exe"
             )
    if(file.exists(exe)) return(normalizePath(exe))
    exe <- file.path("C:", "Program Files (x86)", "Notepad++", "Notepad++.exe", fsep="\\")
    if(file.exists(exe)) return(normalizePath(exe))
    exe <- file.path("C:", "Program Files", "Notepad++", "Notepad++.exe", fsep="\\")
    if(file.exists(exe)) return(normalizePath(exe))
    return(invisible(NULL))
}
find_npptor <- 
function(){
    exe <- file.path(Sys.getenv("PROGRAMFILES"), "NppToR", "NppToR.exe", fsep="\\")
    if(file.exists(exe)) return(normalizePath(exe))
    exe <- file.path(Sys.getenv("PROGRAMFILES(X86)"), "NppToR", "NppToR.exe", fsep="\\")
    if(file.exists(exe)) return(normalizePath(exe))
    exe <- file.path(Sys.getenv("APPDATA"), "NppToR", "NppToR.exe", fsep="\\")
    if(file.exists(exe)) return(normalizePath(exe))
    return(invisible(NULL))
}

#' Launch NppToR
#' 
#' @param ... passed on as arguments to npptor
#' @param exe path to the NppToR exacutable
#' @param startup should the '-startup' parameter be passed to npptor?
#' 
#' @export
npptor <- 
function( ...
        , exe = getOption("wingui::npptor", find_npptor())
        , startup = getOption("wingui::startup", is_r_startup())
        ){
    c( "-rhome", R.home()
     , "-npp", getOption("wingui::Notepad++", find_npp())
     , if(startup) "-startup"
     )

    if(!is.null(exe))
        system2(exe, list(...), wait=F, invisible=F)
}

#' Launch Notepad++
#' 
#' @param file  file to open in Notepad++
#' @param new   open in a new instance?
#' @param exe   Path to Notepad++ executable.
#' 
#' @export
npp <- 
function( file  = NULL
        , new   = FALSE
        , exe   = getOption("wingui::Notepad++", find_npp())
        ){
	args <- NULL
	if(is.character(file))
		args <- shQuote(normalizePath(file))
	if(new){
		args <- c("-multiInst", args)
	}
    if(!is.null(exe))
        system2(exe, args, wait=F, invisible=F)
}
