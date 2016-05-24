# compress2matrix is a function to read 'compressed format'
#  text files giving the species found at particular sites
#  and rewrite them into a sites by species matrix

# Compressed format data has a line of data for each site
#  consisting of tab-separated records, the first two records
#  give the site name; the remaining records are pairs
#  separated by a single space giving species names and counts
#  (species names can be repeated)

# path gives the directory within which to process all the files
#  that end in '.BB' (which must be capitalized)
# skip is the number of lines to skip at the beginning of each
#  file, defaulting to 3
# verbose = TRUE (the default) prints diagnostic messages
# lump.queries = TRUE treats taxon names/codes that are followed
#  by question marks as if they were the same as the same
#  name/code not followed by a question mark. FALSE (the default)
#  treats queries as separate taxa.

# The return value of the function is a list including a vector
#  of site names, a vector of taxon names/codes, and a sites-by-
#  species matrix. As a side-product, the function saves (in the 
#  directory specified by path) a matrix (tab-delimited text
#  file ending in .MTX) for each locality (.BB file) and the
#  overall sites-by-species matrix as a file named 'SiteBySpp.MTX'

compress2matrix <- function (path, skip = 3, verbose = TRUE,
                             lump.queries = FALSE,
                             exclude.taxa = NULL){
  if(is.null(path)) path <- getwd()
  filelist <- dir(path)
  filelist <- filelist[grep('.BB$', filelist)]
  if(length(filelist) == 0)
    stop('there are no .BB files to process')
  
  pa.all <- vector(mode = 'list', length = length(filelist))
  cover.all <- vector(mode = 'list', length = length(filelist))
      
  # loop through list of .BB files to process  
  for(i in 1:length(filelist)){
    raw <- scan(paste(path, filelist[i], sep=""),
                what = '', sep = '\n', skip = skip)
    raw <- toupper(raw)
    if(lump.queries) raw <- gsub('?', '', raw, fixed = TRUE)
    if(verbose) cat(paste('successfully scanned', filelist[i]))
    if(verbose) cat('\n')
    
    # remove extra tab characters and then split each line on tabs
    raw <- gsub('\t{2,}', '', raw)
    raw <- gsub('\t$', '', raw)
    splits <- strsplit(raw, '\t')
    
    # row lables for new matrix (site names) are made up
    #  from the first pair of items in each data line
    #  of the compressed format file
    rowlbls <- lapply(splits, function(x){x[1:2]})
    rowlbls <- unlist(lapply(rowlbls,
                      function(x){paste(x, collapse='')}))
    
    data <- lapply(splits, function(x){x[3:length(x)]})
    # 'data' is a list each element of which is the raw data from
    #  a line in the compressed format file, i.e. one census
   
    # using subroutine 'adding()'... to convert each element
    #  of the list 'data' to a table of counts

#if(filelist[i] == '14_1.BB') browser()
# example of a leading space throwing things off in prior version

    addeddata <- lapply(data, adding)
    addeddata <- lapply(addeddata, as.data.frame)
    names(addeddata) <- rowlbls
    data.out <- rbind.all(addeddata)
    if(verbose) {cat(paste(nrow(data.out),
                          'rows/census localities and',
                          ncol(data.out),
                          'columns/species\n'))
    }

    newname <- gsub('BB', 'MTX', filelist[i])
    write.table(data.out, paste(path, newname, sep=""))

   
    #exclude.taxa <- c("IF", "ICO", "ICY", "IM", "ID", "BL",
    #                  "WA", "X50PH", "X100PH", "AMB", "CHAR",
    #                  "IND1", "IREP", "PREP")
    data.excl <- data.out[,!colnames(data.out) %in% exclude.taxa]
    pa <- apply(data.excl > 0, 2, function(x){sum(x,
    	                                          na.rm = TRUE)})
    cover <- apply(data.excl, 2, function(x){sum(x,
    	                                         na.rm = TRUE)})
    cover.all[[i]] <- cover
    pa.all[[i]] <- pa
    
    pa.pc <- 100 * (pa / sum(pa))
    cover.pc <- 100 * (cover / sum(cover))

    write.table(rbind(data.excl, rbind(pa, pa.pc,
                                       cover, cover.pc)),
                paste(path, newname, 'excl.dat.plus.sums.txt',
                      sep = ''))
    write.table(rbind(pa, pa.pc, cover, cover.pc),
                paste(path, newname, 'excl.dat.sums.txt',
                      sep = ''))
    
    if(i == 1){
      alltaxa <- colnames(data.out)
      site.sums <- list()
      ncensuses <- nrow(data.out)
    }
    else{
      alltaxa <- c(alltaxa, colnames(data.out))
      ncensuses <- ncensuses + nrow(data.out)
    }
    site.sums[[i]] <- colSums(data.out, na.rm = TRUE)
    names(site.sums[[i]]) <- colnames(data.out)
    if(verbose) print(site.sums[[i]])
    if(verbose) cat('__________\n')

  } # end loop through files
  
  sites <- gsub('.BB', '', filelist)
  names(site.sums) <- sites
  taxa <- unique(alltaxa)
  site.sp.mtx <- mergelist(site.sums)
  rownames(site.sp.mtx) <- sites
  site.sp.mtx <- site.sp.mtx[,sort(colnames(site.sp.mtx))]
  write.table(site.sp.mtx, paste(path, 'SiteBySpp.MTX', sep=""))
  pa.all <- mergelist(pa.all)
  cover.all <- mergelist(cover.all)
  rownames(pa.all) <- gsub('.BB', '', filelist)
  rownames(cover.all) <- gsub('.BB', '', filelist)
  write.table(pa.all, paste(path, 'PA.MTX', sep = ''))
  write.table(cover.all, paste(path, 'Cover.MTX', sep = ''))
   
  cat('Processed a total of', ncensuses,
      'censuses,\ngiving counts of',
      length(taxa), 'unique taxa at',
      length(sites),'sites\n')
      
  return(list(sites = sites, taxa = taxa,
              site.sp.mtx = site.sp.mtx,
              pa = pa.all, cover = cover.all))
} # end function

############################################################
# Subroutines

# adding takes a string of taxon-count pairs, splits it on
#  spaces, and then adds up the counts for each taxon 
adding <- function(x){
  x <- gsub('^[[:blank:]]+', '', x, perl = TRUE)
  x <- gsub('[[:blank:]]+$', '', x, perl = TRUE)
  spcountmatrix <- matrix(unlist(strsplit(x,
                                          split = '[[:blank:]]+',
                                          perl = TRUE)),
                          ncol = 2, byrow = TRUE)
  spcountlist <- split(spcountmatrix[,2], spcountmatrix[,1])
  spsums <- lapply(spcountlist,
                   function(x){sum(as.numeric(x))})
  return(spsums)
}

#rbind.all code modified from:
# From: Sundar Dorai-Raj <sundar.dorai-raj <at> pdf.com># Subject: Re: merging# Newsgroups: gmane.comp.lang.r.general# Date: 2006-05-30 20:09:27 GMT
rbind.all <- function(x) {
  cn <- unique(unlist(lapply(x, colnames)))  for(i in seq(along = x)) {    if(any(m <- !cn %in% colnames(x[[i]]))) {
      na <- matrix(NA, nrow(x[[i]]), sum(m))      dimnames(na) <- list(rownames(x[[i]]), cn[m])      x[[i]] <- cbind(x[[i]], na)    }
  }  do.call(rbind, x)}

mergelist <- function(x){
  frame.out <- t(x[[1]])
  for(i in 2:length(x))
    frame.out <- merge(frame.out, t(x[[i]]),
                       all = TRUE, sort = FALSE)
  return(frame.out)
}