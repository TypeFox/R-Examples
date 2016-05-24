#' Function to load the package built-in RData
#'
#' \code{xRDataLoader} is supposed to load the package built-in RData.
#'
#' @param RData which built-in RData to load. It can be one of "GWAS2EF", "GWAS_LD", "IlluminaHumanHT", "IlluminaOmniExpress", "ig.DO", "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.string", "org.Hs.PCommons_DN", "org.Hs.PCommons_UN"
#' @param RData.customised a file name for RData-formatted file. By default, it is NULL. It is designed when the user wants to import customised RData that are not listed in the above argument 'RData'. However, this argument can be always used even for those RData that are listed in the argument 'RData' 
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @param RData.location the characters to tell the location of built-in RData files. By default, it remotely locates at \url{https://github.com/hfang-bristol/RDataCentre/blob/master/XGR}. For the user equipped with fast internet connection, this option can be just left as default. But it is always advisable to download these files locally. Especially when the user needs to run this function many times, there is no need to ask the function to remotely download every time (also it will unnecessarily increase the runtime). For examples, these files (as a whole or part of them) can be first downloaded into your current working directory, and then set this option as: \eqn{RData.location="."}. Surely, the location can be anywhere as long as the user provides the correct path pointing to (otherwise, the script will have to remotely download each time)
#' @return 
#' any use-specified variable that is given on the right side of the assigement sign '<-', which contains the loaded RData.
#' @note If there are no use-specified variable that is given on the right side of the assigement sign '<-', then no RData will be loaded onto the working environment.
#' @export
#' @import dnet
#' @import igraph
#' @importFrom GenomicRanges findOverlaps distance mcols seqnames as.data.frame
#' @importFrom grDevices colorRampPalette dev.cur rgb
#' @seealso \code{\link{xRDataLoader}}
#' @include xRDataLoader.r
#' @examples
#' ImmunoBase <- xRDataLoader(RData.customised='ImmunoBase')
#' \dontrun{
#' org.Hs.eg <- xRDataLoader(RData='org.Hs.eg')
#' ig.HPPA <- xRDataLoader(RData='ig.HPPA')
#' org.Hs.egHPPA <- xRDataLoader(RData='org.Hs.egHPPA')
#' org.Hs.egHPPA <- xRDataLoader(RData.customised='org.Hs.egHPPA')
#' }

xRDataLoader <- function(RData=c(NA,"GWAS2EF", "GWAS_LD", "IlluminaHumanHT", "IlluminaOmniExpress", "ig.DO", "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA", "ig.HPMI", "ig.HPPA", "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP", "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA", "org.Hs.egHPMI", "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1", "org.Hs.egMsigdbC2BIOCARTA", "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall", "org.Hs.egMsigdbC2CP", "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME", "org.Hs.egMsigdbC3MIR", "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM", "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF", "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH", "org.Hs.egPS", "org.Hs.egSF", "org.Hs.string", "org.Hs.PCommons_DN", "org.Hs.PCommons_UN"), RData.customised=NULL, verbose=T, RData.location="https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0")
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    RData <- match.arg(RData)

	######################################################################################
	# RData now is primarily hosted in github
	######################################################################################
	my_https_downloader <- function (url, method=c("auto","internal","wininet","libcurl","wget","curl"), quiet=T, mode=c("w","wb","a","ab"), cacheOK=T, extra=getOption("download.file.extra")){
	
		## https://stat.ethz.ch/R-manual/R-devel/library/utils/html/download.file.html
		method <- match.arg(method)
		mode <- match.arg(mode)
	
		## specify the temporary image files
		tdir <- tempdir()
		destfile <- file.path(tdir, "temp.RData")
		## remove the existing temporary RData file
		unlink(destfile, recursive=T, force=T)
	
		if(base::grepl("^https?://", url)){
			isR32 <- base::getRversion() >= "3.2"
			if(.Platform$OS.type == "windows"){
				if(isR32){
					method <- "wininet"
				}else{
					seti2 <- utils::"setInternet2"
					internet2_start <- seti2(NA)
					if(!internet2_start){
						on.exit(suppressWarnings(seti2(internet2_start)))
						suppressWarnings(seti2(TRUE))
					}
					method <- "internal"
				}
				suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
			}else{
				if(isR32 && capabilities("libcurl")){
					method <- "libcurl"
				}else if(nzchar(Sys.which("wget")[1])){
					method <- "wget"
				}else if(nzchar(Sys.which("curl")[1])){
					method <- "curl"
					orig_extra_options <- getOption("download.file.extra")
					on.exit(options(download.file.extra = orig_extra_options))
					options(download.file.extra = paste("-L", orig_extra_options))
				}else if(nzchar(Sys.which("lynx")[1])) {
					method <- "lynx"
				}else{
					stop("no download method found")
				}
				suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
			}
		}else{
			suppressWarnings(utils::download.file(url, destfile=destfile, method=method, quiet=quiet, mode=mode, cacheOK=cacheOK, extra=extra))
		}
	
		if(file.exists(destfile) & file.info(destfile)$size!=0){
			res_RData <- get(load(destfile))
			res_flag <- T
		}else{
			res_RData <- NULL
			res_flag <- F
		}
		
		res <- list(RData = res_RData,
					flag = res_flag)
		
		invisible(res)
	}
	######################################################################################
	######################################################################################


    if(is.na(RData) & !is.null(RData.customised)){
		RData <- RData.customised
	}else if(is.na(RData) & is.null(RData.customised)){
		stop("There is no input! Please input one of two parameters ('RData' or 'RData.customised').\n")
	}
	
	RData <- gsub('.RData$', "", RData, ignore.case=T, perl=T)
	RData <- gsub(".RDa$", "", RData, ignore.case=T, perl=T)
	
    ###############################
    ## make sure there is no "/" at the end
    path_host <- gsub("/$", "", RData.location)
    if(path_host=="" || length(path_host)==0 || is.na(path_host)){
        path_host <- "https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0"
    }
    
    ## load 
    load_remote <- paste(path_host, "/", RData, ".RData", sep="")
    load_local1 <- file.path(path_host, paste("data/", RData, ".RData", sep=""))
    load_local2 <- file.path(path_host, paste(RData, ".RData", sep=""))
    load_package <- RData
    
    ## first, load data from the package itself
    if(length(suppressWarnings(tryCatch(eval(parse(text=paste("data(",load_package,", package='XGR')",sep=""))), error=function(e) e, warning=function(w) w)))==2){
        ## second, load local R files
        RData_local <- c(load_local1, load_local2)
        load_flag <- sapply(RData_local, function(x){
            if(.Platform$OS.type=="windows") x <- gsub("/", "\\\\", x)
            ifelse(file.exists(x), TRUE, FALSE)
        })
        ## otherwise, load remote R files
        if(sum(load_flag)==0){
        	
        	flag_failed <- F
        	if(length(grep('^https',load_remote,perl=T))){
        		if(length(grep('github',load_remote,perl=T))){
        			load_remote <- paste(load_remote, "?raw=true", sep="")
        		}
        		res <- my_https_downloader(load_remote, mode="wb")
        		if(res$flag==F){
        			flag_failed <- T
        		}else{
        			eval(parse(text=paste(RData, " <- res$RData", sep="")))
        		}
        	}else{
        		if(class(suppressWarnings(try(load(url(load_remote)), T)))=="try-error"){
        			flag_failed <- T
        		}
        	}
        	
			if(flag_failed){
				load_remote <- paste("https://github.com/hfang-bristol/RDataCentre/blob/master/XGR/1.0.0/", RData, ".RData?raw=true", sep="")
				
				res <- my_https_downloader(load_remote, mode="wb")
				if(res$flag==F){
					stop("Built-in Rdata files cannot be loaded. Please check your internet connection or their location in your local machine.\n")
				}else{
					eval(parse(text=paste(RData, " <- res$RData", sep="")))
				}
			}
		
			load_RData <- load_remote
        }else{
            load_RData <- RData_local[load_flag]
            load(load_RData)
        }
    }else{
        load_RData <- sprintf("package 'XGR' version %s", utils::packageVersion("XGR"))
    }
    
	out <- base::get(RData)
	
    if(verbose){
        now <- Sys.time()
        message(sprintf("'%s' (from %s) has been loaded into the working environment (at %s)", RData, load_RData, as.character(now)), appendLF=T)
    }
    
    invisible(out)
}
