
globalVariables(c("DEBUG", "entries", "this"))

library(digest) # digest
library(tools) # md5sum


repo_open <- function(root="~/.R_repo", force=F)
{
    DEBUG <- F

    handleErr <- function(err, ...)
        {
            pars <- list(...)
            if(length(pars)>0)
                lpars <- paste0(paste(pars[[1]], collapse=", "), ".")
            switch(err,
                   "DEBUG" = {
                       message(pars[[1]])
                   },
                   "ID_NOT_FOUND" = {
                       stop(paste0("Item not found: ", lpars))
                   },
                   "ID_EXISTING" = {
                       stop(paste0("There is already an item with this name: ", lpars))
                   },
                   "ID_RESERVED" = {
                       stop(paste0("Id not valid (reserved): ", lpars))
                   },
                   "TAG_RESERVED" = {
                       warning(paste0("Reserved TAG used: ", lpars))
                   },
                   "EMPTY_REPO" = {
                       stop("Repo is empty.")
                   },
                   "MISS_OBJ_HAS_URL" = {
                       stop(paste0("The file object could not be found. ",
                                   "However, it can be downloaded using pull."))
                   },
                   "NO_URL" = {
                       stop("The object has no associated URL.")
                   },
                   "LAZY_FOUND" = {
                       message("lazydo found precomputed resource.")
                   },
                   "LAZY_NOT_FOUND" = {
                       message("lazydo is building resource from code.")
                   },
                   "LAZY_NOT_EXPR" = {
                       stop("expr must be of class expression.")
                   },
                   "LAZY_NAME" = {
                       message(paste0("Cached item name is: ", lpars))
                   },
                   "DATA_ALREADY_THERE" = {
                       stop(paste0("There is existing content for ", lpars, ". ",
                                   "Use replace=T to overwrite."))
                   },
                   "ATTACHMENT_FILE_NOT_FOUND" = {
                       stop(paste0("Attachment file not found"))
                   }
                   )

        }

    
    getFile <- function(name)
        {
            entry <- getEntry(name)

            if(substr(normalizePath(entry$dump, mustWork=F), 1, nchar(root)) == root) {
                fpath <- normalizePath(entry$dump, mustWork=F)
                } else fpath <- normalizePath(file.path(root, entry$dump), mustWork=F)

            return(fpath)
        }

    setData <- function(name, obj, asattach)
    {
        w <- findEntryIndex(name)
        newdata <- list()

        if(!asattach) {
            if(!is.null(dim(obj)))
                dims <- dim(obj) else dims <- length(obj)
        } else dims <- NULL

        rmData(name, "temp")            
        
        tryCatch({
            fdata <- storeData(name, obj, asattach)
        }, error = function(e) {
            print(e)
            rmData(name, "undo")
        }, finally = {
            rmData(name, "finalize")
        }
        )                
        
        newdata[["dump"]] <- relativePath(fdata[["path"]])
        newdata[["size"]] <- fdata[["size"]]
        newdata[["checksum"]] <- md5sum(path.expand(fdata[["path"]]))
        newdata[["dims"]] <- dims

        return(newdata)
    }
    
    checkTags <- function(tags, name=NULL)
    {

        dups <- which(duplicated(tolower(tags)))
        if(length(dups)>0)
            warning(paste0("The following tags are duplicated (not case sensitive): ",
                        paste0(tags[dups], collapse=", ")))
        
        if(!is.null(name)) {
            e <- getEntry(name)
            comm <- intersect(tolower(e$tags), tolower(tags))
            if(length(comm)>0)
                warning(paste0("The following tags are already present (not case sensitive): ",
                            paste0(comm, collapse=", ")))
        }
        reservedTags <- c("stash", "attachment")
        if(any(tolower(tags) %in% reservedTags))
            handleErr("TAG_RESERVED", tags[tolower(tags) %in% reservedTags])                 

        return(unique(tags))
    }
    
    checkVersions <- function(name)
        {
            names <- sapply(entries, get, x="name")
            ## searching for names ending with # and a number
            w <- regexpr(paste0(name, "#[[:digit:]]+$"), names)
            
            ## extract the numbers
            v <- as.numeric(gsub(paste0(name,"#"),"",names[w!=-1]))
            if(length(v)>0)
                newname <- paste0(name, "#", max(v)+1)
            else
                newname <- paste0(name, "#1")
            
            return(list(w=which(w!=-1), v=v, new=newname))
        }


    relativePath <- function(path)
    {
        sep <- .Platform$file.sep
        root <- get("root",thisEnv)                
        relpath <- gsub(paste0(root, sep), "", path ,fixed=T)
        return(relpath)
    }
    
    stopOnEmpty <- function(doreturn=F) {
        if(length(entries)<1) {
            if(doreturn)
                return(1) else handleErr("EMPTY_REPO")
        }
        return(0)
    }

    setEntry <- function(name, newEntry)
    {
        stopOnNotFound(name)
        e <- findEntryIndex(name)
        entries[[e]] <- newEntry
        assign("entries", entries, thisEnv)
    }
    
    stopOnNotFound <- function(names=NULL, tags=NULL)
        {
            stopOnEmpty()
            allnames <- sapply(entries, get, x="name")
            w <- match(names, allnames)
            if(all(is.na(w)))
                handleErr("ID_NOT_FOUND", names[is.na(w)])
        }
    
    getEntry <- function(name) {
            e <- findEntryIndex(name)
            if(is.null(e))
                return(invisible(NULL))
            
            return(entries[[e]])
        }

    runWithTags <- function(f, tags, names, askconfirm, tagfun="OR", ...) {
        if(!is.null(tags))
            e <- findEntries(tags, tagfun) else {
                dbnames <- sapply(entries, get, x="name")
                e <- match(names, dbnames)
                if(any(is.na(e))) {
                    errmsg <- paste0("The following names could not be found: ",
                           paste(names[which(is.na(e))], collapse=", "))
                    stop(errmsg)
                }
            }
                            
            if(length(e)<1)
                stop("Tag or name list does not match any entry.")
            entr <- entries
            names <- sapply(entr[e], get, x="name")
            if(askconfirm) {
                cat(paste0("Matched entries:\n",
                               paste(names, collapse="\n"), "\n"))
                n <- readline("Type \"yes\" to confirm: ")
            } else n <- "yes"
            if(n == "yes") {
                for(i in 1:length(e))
                    get("this", thisEnv)[[f]](name=names[[i]], ...)
                return(invisible())
            } else
                {
                    message("Nothing done.")
                    return(invisible())
                }
        }

    
    hmnRead <- function(bytes) {

        values <- c(
            bytes,
            bytes/2^10,
            bytes/2^20,
            bytes/2^30,
            bytes/2^40,
            bytes/2^50
            )
        names(values) <- c("B","kB", "MB", "GB", "TB", "PB")
        okvals <- values[values>1]
        m <- okvals[which.min(okvals)] ## preserves name
        final <- format(round(m,2), scientific=F)
        return(paste0(final, " ", names(m)))
        }
    
    checkIndexUnchanged <- function() {
        if(DEBUG) {
            message(paste0("checkIndexUnchanged: cur MD5 is ", md5sum(repofile)))
            message(paste0("checkIndexUnchanged: stored MD5 is ", indexMD5))
        }
        
      if(indexMD5 != md5sum(repofile))
          stop(format(paste0("Repo index has been modified outside this session. ",
                             "For security reasons I will stop here. Please open the ",
                             "repo again to sync with the latest changes ",
                             "(which may include deletions). You may want to run \"check\" ",
                             "on this session first.")))
    }
    
    storeIndex <- function() {
        saveRDS(entries, repofile)
        if(DEBUG) {
            message(paste0("storeIndex: stored MD5 is ", indexMD5))
            message(paste0("StoreIndex: new MD5 is ", md5sum(repofile)))
        }
        ## NOTE: do not indexMD5 <- md5sum(repofile)... doesn't work
        assign("indexMD5", md5sum(repofile), thisEnv)
        }

    buildpath <- function(resname)
        {
            resname <- paste0(sample(c(0:9,letters), 32, T),collapse="")
            return(list(root,
                        substr(resname, 1, 2),
                        substr(resname, 3, 4),
                        substr(resname, 5, 6),
                        resname))
        }

    checkName <- function(name)
        {
            if(length(entries)<1)
                return(T)
            names <- sapply(entries, get, x="name")
            return(!(name %in% names))
        }

    rmData <- function(name, phase)
    {
        fpath <- getFile(name)
        fpath_temp <- paste0(fpath, ".remove_me")

        if(!file.exists(fpath) && !file.exists(fpath_temp)){
            ## this should never happen, unless file was removed from
            ## something else:
            warning(paste("File to be removed was not found:", fpath))
            return(invisible(0))
        }
        
        if(phase=="temp") {
            file.rename(fpath, fpath_temp)
        } 
        if(phase=="finalize") {
            file.remove(fpath_temp)
        }
        if(phase=="undo") {
            file.rename(fpath_temp, fpath)
        }
        
        return(invisible(NULL))
    }

    storeData <- function(name, obj, attach=F)
        {
            opath <- buildpath(name)
            if(!file.exists(do.call(file.path, opath[1:2])))
                dir.create(do.call(file.path, opath[1:2]))
            if(!file.exists(do.call(file.path, opath[1:3])))
                dir.create(do.call(file.path, opath[1:3]))
            if(!file.exists(do.call(file.path, opath[1:4])))
                dir.create(do.call(file.path, opath[1:4]))

            fpath <- do.call(file.path, opath)

            if(!attach) {
                saveRDS(obj, fpath)
                } else {
                    didwork <- file.copy(obj, fpath)
                    if(!didwork)
                        stop(paste0("There was an error while trying to write file: ",
                                    fpath))
                }
                
            return(list(path=fpath, size=file.info(fpath)$size))
        }

    depgraph <- function(depends=T, attached=T, generated=T)
        {
            stopOnEmpty()
            if(!any(c(depends, attached, generated)))
              stop("One of depends, attached or generated must be true.")
            
            nodes <- unique(unlist(sapply(entries, get, x="name")))
            ## if(generated) {
            ##   srcs <- unique(unlist(sapply(entries, get, x="source")))
            ##   nodes <- c(nodes, srcs)
            ## }
            n <- length(nodes)
            depgraph <- matrix(0,n,n)
            for(i in 1:length(entries))
                {
                  e <- entries[[i]]
                  if(depends) {
                    w <- match(e$depends, nodes)
                    depgraph[i, w] <- 1
                  }
                  if(attached) {
                    w <- match(e$attachedto, nodes)
                    depgraph[i, w] <- 2
                  }
                  if(generated) {
                      w <- match(e$source, nodes)
                      depgraph[i, w] <- 3
                    }
                }
            rownames(depgraph) <- colnames(depgraph) <- nodes
            return(depgraph)
        }

    checkFoundEntry <- function(e)
        {
            if(is.null(e))
                cat("Entry not found.\n")
            if(e==-1)
                cat("Repo is empty.\n")
        }

    findEntryIndex <- function(name)
        {
            if(is.null(entries) | length(entries)<1) {
                return(-1)
            }
            names <- sapply(entries, get, x="name")
            w <- match(name, names)
            if(length(w)<1){
                return(NULL)
            }
            return(w)
        }

    findEntries <- function(tags=NULL, tagfun="OR", find=NULL)
        {
            if(!is.null(tags)) {
                   tagsets <- lapply(entries, get, x="tags")

                   if(is.character(tagfun) && tagfun =="AND")
                       tagfun <- function(x, tags=tags)all(tags %in% x)
                   if(is.character(tagfun) && tagfun=="NOT")
                       tagfun <- function(x, tags=tags)all(!(tags %in% x))
                   if(is.character(tagfun) && tagfun=="OR")
                       tagfun <- function(x, tags=tags)any(tags %in% x)

                   if(class(tagfun)!="function")
                       stop("tagfun must be either a function or one of OR, AND, NOT")

                   w <- sapply(tagsets, tagfun, tags)
               } else {
                   strmat <- entriesToMat(1:length(entries))
                   w <- apply(strmat, 1, function(l)
                       length(grep(find, l, ignore.case=T))>0)
               }
               
            return(which(w))
        }

    isAttachment <- function(name)
        {
            w <- findEntryIndex(name)
            return("attachment" %in% entries[[w]]$tags)
        }

    attachments <- function(name)
        {
            r <- match(name,  sapply(entries, get, x="attachedto"))
            if(is.na(r))
                return(NULL)            
            return(r)
        }

    dependants <- function(name)
        {
            r <- sapply(sapply(entries, get, x="depends"), match, x=name)
            w <- which(!is.na(r))
            if(length(w)<1)
                return(NULL)            
            return(w)
        }
    
    compressPath <- function(path)
        {
            hp <- path.expand("~")
            return(gsub(paste0("^",hp), "~", path))
        }

    entriesToMat <- function(w)
        {
            entr <- entries[w]

            labels <- c("ID", "a@><", "Dims", "Tags", "Size")
            names <- sapply(entr, get, x="name")

            a <- matrix(NA, length(names), length(labels))
            colnames(a) <- labels

            attachs <- depends <- hasattach <- allows <- rep(" ", length(entr))
                            
            tagsets <- lapply(entr, get, x="tags")
            attachs[sapply(tagsets, is.element, el="attachment")] <- "x"
            depends[sapply(lapply(entr, get, x="depends"), length)>0] <- "x"
            allows[!sapply(lapply(names, dependants), length)>0] <- "x"
            hasattach[!sapply((sapply(names, attachments)), is.null)] <- "x"            

            flags <- paste0(attachs, hasattach, depends, allows)
            
            descriptions <- sapply(entr, get, x="description")
            prefixes <- rep("", length(names))
            prefixes[attachs == "x"] <- "@"                        

            a[,"ID"] <- paste0(prefixes, names)
            a[,2] <- flags
            a[,"Dims"] <- sapply(lapply(entr, get, x="dims"), paste, collapse="x");
            a[a[,"Dims"]=="", "Dims"] <- "-"            
            a[,"Tags"] <- sapply(tagsets, paste, collapse=", ")
            a[,"Size"] <- sapply(lapply(entr, get, x="size"), hmnRead)

            return(a)
        }

    
    cutString <- function(text, len, dotsafter=T)
        {
            if(nchar(text) == len)
                return(text)
            if(nchar(text) < len)
                return(format(text, width = len))
            if(dotsafter) {
                text <- substr(text, 1, len-3)
                return(paste0(text ,"..."))
            } else {
                text <- substr(text, nchar(text)-(len-4), nchar(text))
                return(paste0("...", text))
            }
            
        }

    
    me <- list(
        dependencies = function(depends=T, attached=T, generated=T, plot=T)
        {
          deps <- depgraph(depends, attached, generated)
          if(plot) {
              if (requireNamespace("igraph", quietly = TRUE)) {
                  deps2 <- deps
                  rownames(deps2) <- colnames(deps2) <- basename(rownames(deps))
                  g <- igraph::graph.adjacency(deps2, weighted=c("type"))
                  igraph::plot.igraph(g, edge.label=c("depends", "attached", "generated")
                               [igraph::get.edge.attribute(g,"type")])
              } else {
                  stop("The suggested package igraph is not installed.")
              }              
          }
          invisible(depgraph())
        },

        check = function()
            {
                stopOnEmpty()
                entr <- entries

                warn <- 0
                for(i in 1:length(entr))
                {
                    cat(paste0("Checking ", entr[[i]]$name, "..."))
                    if(file.exists(getFile(entr[[i]]$name))){
                        md5s <- md5sum(getFile(entr[[i]]$name))
                        if(md5s != entr[[i]]$checksum) {
                            cat(" changed!")
                            warning("File has changed!")
                        } else cat(" ok.")
                    } else {
                        cat(" not found!")
                        warn <- warn + 1
                    }
                    cat("\n")
                }
                if(warn > 0)
                    warning(paste0("There were ", warn, " missing files!"))

        cat("\nChecking for extraneous files in repo root... ")
                allfiles <- file.path(root, list.files(root, recursive=T))
                dumps <- sapply(sapply(entries, get, x="name"), getFile)
                junk <- setdiff(path.expand(allfiles), path.expand(dumps))
                junk <- setdiff(junk, repofile)
                if(length(junk)>0){
                    cat("found some:\n")
                    cat(paste(junk, collapse="\n"))
                } else cat("ok.")
                cat("\n")
                invisible()
            },

        pies = function() {
            sizes = sapply(entries, get, x="size")
            names(sizes) <- sapply(entries, get, x="name")
            if(length(sizes)>10) {
                sizes <- sort(sizes, decreasing=T)
                sizes <- c(sizes[1:9], sum(sizes[10:length(sizes)]))
                names(sizes)[10] <- "Other"
            }
            pie(sizes)
        },

        copy = function(destrepo, name, tags=NULL)
        {            
            if(!("repo" %in% class(destrepo)))
                stop("destrepo must be an object of class repo.")
            if(!xor(missing(name), is.null(tags)))
                stop("You must specify either names or tags.")

            if(length(name) > 1 | !is.null(tags))
                runWithTags("copy", tags, name, T, destrepo) else {
                    e <- findEntryIndex(name)
                    checkFoundEntry(e)
                    entr <- entries[[e]]
                    obj <- get("this", thisEnv)$get(name)
                    destrepo$put(obj, name, entr$description, entr$tags, entr$source, entr$depends)
                }
        },

        handlers = function()
        {
            h <- list()
            for(i in 1:length(entries))
                {
                    fbody <- paste0('function(f="get", ...)',
                                    'get("this", thisEnv)[[f]](name ="', entries[[i]]$name, '",...)')
                    h[[i]] <- eval(parse(text=fbody))
                }
            h[[length(h)+1]] <- get("this", thisEnv)
            names(h) <- c(sapply(entries, get, x="name"), "repo")
            return(h)
        },
                
        tags = function()
        {
            entr <- entries
            tagset <- unique(unlist(lapply(entr, get, x="tags")))
            return(tagset)
        },

        sys = function(name, command)
        {
            stopOnNotFound(name)
            e <- getEntry(name)
            syscomm <-paste0(command, " ", file.path(root, e[["dump"]]))
            message(paste("Running system command:", syscomm))
            system(syscomm)
        },

        find = function(what, all=F, show="ds")
        {
        get("this", thisEnv)$print(find=what, all=all, show=show)
        },
        
        print = function(tags=NULL, tagfun="OR", find=NULL, all=F, show="ds")
        {
            ## TODO: Part of the code is now in function entriesToMat,
            ## should be removed from here.

            if(!is.null(tags) & !is.null(find))
                stop("Please provide either tags or find.")
            
            stopOnEmpty()
            
            entr <- entries
            if(!is.null(tags) | !is.null(find)) {
                w <- findEntries(tags=tags, tagfun=tagfun, find=find)
                if(length(w)<1)
                    {
                        message("No matches.")
                        return(invisible())
                    } else {
                        entr <- entr[w]
                    }
            }

            
            labels <- c("ID", "a@><", "Dims", "Tags", "Size")
            names <- sapply(entr, get, x="name")

            a <- matrix(NA, length(names), length(labels))
            colnames(a) <- labels

            attachs <- depends <- hasattach <- allows <- rep(" ", length(entr))
                            
            tagsets <- lapply(entr, get, x="tags")
            attachs[sapply(tagsets, is.element, el="attachment")] <- "x"
            depends[sapply(lapply(entr, get, x="depends"), length)>0] <- "x"
            allows[!sapply(lapply(names, dependants), length)>0] <- "x"
            hasattach[!sapply((sapply(names, attachments)), is.null)] <- "x"

            flags <- paste0(attachs, hasattach, depends, allows)

            descriptions <- sapply(entr, get, x="description")
            
            #tagsets <- lapply(tagsets, setdiff, y="attachment")
            #tagsets <- lapply(tagsets, setdiff, y="hide")
            #tagsets <- lapply(tagsets, setdiff, y="stash")

            prefixes <- rep("", length(names))
            prefixes[attachs == "x"] <- "@"                        
            
            a[,"ID"] <- paste0(prefixes, names)
            a[,2] <- flags
            a[,"Dims"] <- sapply(lapply(entr, get, x="dims"), paste, collapse="x"); a[a[,"Dims"]=="", "Dims"] <- "-"            
            a[,"Tags"] <- sapply(tagsets, paste, collapse=", ")
            a[,"Size"] <- sapply(lapply(entr, get, x="size"), hmnRead)
            ##a[,"URL"] <- sapply(entr, get, x="URL")

            h <- rep(F,length(entr))
            hidden <- sapply(tagsets, is.element, el="hide")

            if(!all)
                h[hidden] <- T

            if(length(entr)>1 & all(h))
            {
                message("All matched entries are hidden, use all=T.")
                return(invisible(NULL))
            }

                       
            cols <- c(T, sapply(c("f","d","t","s"), grepl, show))
            m <- as.data.frame(a[!h,cols], nm="")

            if(sum(!h)>1)
                print(m, quote=F, row.names=F) else print(t(m), quote=F, row.names=F)

            invisible(m)
        },

        export = function(name, where=".", tags=NULL)
        {
            if(!xor(missing(name), is.null(tags)))
                stop("You must specify either names or tags.")

            if(!is.null(tags) | length(name)>1){
                runWithTags("export", tags, name, T, where)
            } else {
                ipath = do.call(file.path, buildpath(name))
                if(isAttachment(name))
                    fname <- name else fname <- paste0(name, ".RDS")
                file.copy(ipath, file.path(where, fname))
            }
        },
        
        info = function(name = NULL, tags = NULL)
        {
            stopOnEmpty()
            
            if(!is.null(name))
                stopOnNotFound(name)
            
            if(!xor(is.null(name), is.null(tags))) {
                labels <- c("Root:", "Number of items:", "Total size:")
                maxw <- max(sapply(labels, nchar))
                vals <- c(compressPath(root), length(entries),
                          hmnRead(sum(sapply(entries, get, x="size"))))
                lines <- paste(format(labels, width=maxw), vals, sep=" ")
                for(i in 1:length(lines))
                    cat(lines[[i]], "\n")
                return(invisible(NULL))
            }            

            
            if(!is.null(tags) | length(name)>1){
                runWithTags("info", tags, name, askconfirm=F)
            } else {            
                                        #e <- get("findEntryIndex",thisEnv)(name)
                e <- findEntryIndex(name)
                if(is.null(e))
                    stop("Identifier not found.")

                labels <- c("ID:", "Description:", "Tags:",
                            "Dimensions:", "Timestamp:",
                            "Size on disk:", "Provenance:",
                            "Attached to:", "Stored in:", 
                            "MD5 checksum:", "URL:")
                maxlen <- max(sapply(labels, nchar))

                if(is.null(entries[[e]]$attachedto))
                    att <- "-" else att <- paste(entries[[e]]$attachedto, collapse=", ")
                if(is.null(entries[[e]]$URL))
                    url <- "-" else url <- entries[[e]]$URL

                vals <- c(entries[[e]]$name, entries[[e]]$description,
                          paste0(entries[[e]]$tags, collapse=", "),
                          paste(entries[[e]]$dims, collapse="x"),
                          as.character(entries[[e]]$timestamp),
                          hmnRead(entries[[e]]$size),
                          paste(entries[[e]]$source, collapse=", "),
                          att, entries[[e]]$dump, entries[[e]]
                          $checksum, url)
                cat(paste0(format(labels, width=maxlen+1), vals, "\n"), sep="")
                cat("\n")
            }
        },

        rm = function(name = NULL, tags = NULL, force = F)
        {
          checkIndexUnchanged()                   
            
            if(!xor(missing(name),missing(tags)))
                stop("You must specify either a name or a set of tags.")

            if(!is.null(tags) | length(name)>1){
                runWithTags("rm", tags, name, !force)
            } else {            
                e <- get("findEntryIndex",thisEnv)(name)
                if(is.null(e))
                    return(invisible(NULL))

                rmData(name, "temp")
                rmData(name, "finalize")

                assign("entries", entries[-e], thisEnv)                
                storeIndex()
            }
      },

        bulkedit = function(outfile=NULL, infile=NULL)
            {
                if(!xor(is.null(infile), is.null(outfile)))
                    stop("Please provide exactly one of infile or outfile.")

                if(!is.null(outfile)) {
                    if(file.exists(outfile))
                        stop("File already exists.")
                
                    outf <- file(outfile,"at")
                    writeLines(paste0(digest(entries),
                                    " EDIT NEXT LINES ONLY. FIELDS MUST BE TAB-SEPARATED. ",
                                    "TAGS MUST BE COMMA-SEPARATED."), outfile)
                    for(i in 1:length(entries)){
                        src <- entries[[i]]$source
                        if(is.null(src)) src <- "NULL"
                        line <- paste0(c(
                            entries[[i]]$name,
                            entries[[i]]$description,
                            paste0(entries[[i]]$tags, collapse=", "),
                            src,
                            entries[[i]]$attachedto,
                            entries[[i]]$depednds),
                                       collapse="\t"
                                       )
                        writeLines(line, outf)
                    }
                    close(outf)
                } else {
                    checkIndexUnchanged()
                                
                    if(!file.exists(infile))
                        stop("Can't find input file.")
                
                    indata <- strsplit(readLines(infile), "\t")

                    csum <- substr(indata[[1]],1,32) 
                    if(csum != digest(entries))
                        stop(paste0("Checksum mismatch: it seems that input data were ",
                            "made for entries that have changed in the meantime. \n",
                            ## "Current checksum is: ", digest(entries), "\n",
                            ## "Checksum in input file is: ", indata[[1]], "\n",
                            "Overwriting could be dangerous, so I will stop here. ",
                            "Call bulkedit again to create a new input file.")
                            )
                    indata <- indata[-1]

                    rset <- get("this", thisEnv)$set

                    for(i in 1:length(indata)) {
                        src <- indata[[i]][[4]]
                        if(src=="NULL")
                            src <- NULL
                        entries[[i]]$name <- indata[[i]][[1]]
                        entries[[i]]$description <- indata[[i]][[2]]
                        entries[[i]]$tags <- strsplit(gsub(" ", "", entries[[i]]$tags), ",")
                        entries[[i]]$source <- src
                    }

                    assign("entries", entries, thisEnv)                
                    storeIndex()
                    message("Entries updated.")
                }
            },

        get = function(name)
        {          
            if(checkName(name)){                
                enames <- sapply(entries, get, x="name")
                x <- agrep(name, enames)
                if(length(x)>0) {
                    x <- x[abs(sapply(enames[x],nchar) - nchar(name))<=3]
                    message(paste0(
                        "Maybe you were looking for: ",
                        paste0(enames[x], collapse=", ")
                    ))
                }
                handleErr("ID_NOT_FOUND", name)
                return(invisible())
            }
            entry <- getEntry(name)
            root <- get("root",thisEnv)
            if(substr(normalizePath(entry$dump, mustWork=F), 1, nchar(root)) == root) {
                newpath <- relativePath(normalizePath(entry$dump, mustWork=F))
                message(paste0("This resource was indexed in a deprecated format. ",
                               "Now updating position from:\n", entry$dump, "\nto:\n",
                               newpath))
                entry$dump <- newpath
                setEntry(name, entry)
                storeIndex()
                }

            if(isAttachment(name))
              stop("Get is not valid for attachments.")

            f <- getFile(name)
            if(!file.exists(f) && !is.null(entry$URL))
                handleErr("MISS_OBJ_HAS_URL")
            data <- readRDS(f)
            
            return(data)
        },
                
        entries = function()
        {
            return(get("entries",thisEnv))
        },

        tag = function(name = NULL, newtags, tags = NULL)
        {
            if(!xor(is.null(name), is.null(tags)))
                stop("You must provide either name or tags.")
            if(!is.null(tags) | length(name)>1)
                runWithTags("tag", tags, name, F, newtags) else
                    get("this", thisEnv)$set(name, addtags=newtags)                   
        },

      lazydo = function(expr, force=F, env=parent.frame())
      {
          if(!is.expression(expr))
              handleErr("LAZY_NOT_EXPR")
          
          src <- as.character(expr)
          resname <- digest(src)

          if(checkName(resname) || force)
          {
              handleErr("LAZY_NOT_FOUND")
              res <- eval(expr, envir=env)
              get("this", thisEnv)$stash(res, resname)
              get("this", thisEnv)$set(resname,
                                       description=quote(expr),
                                       addtags="lazydo")
              handleErr("LAZY_NAME", resname)
              return(res)
          } else {
              handleErr("LAZY_FOUND")
              return(get("this", thisEnv)$get(resname))
           }
      },


        untag = function(name = NULL, rmtags, tags = NULL)
        {
            if(!xor(is.null(name), is.null(tags)))
                stop("You must provide either name or tags.")
            if(!is.null(tags) | length(name)>1)
                runWithTags("untag", tags, name, F, rmtags) else {
                    currtags <- getEntry(name)$tags
                    w <- rmtags %in% currtags
                    if(!any(w))
                        warning(paste0("Tag/s ", paste0(rmtags[!w], collapse=", "),
                                       " not present in entry ", name))
                    currtags <- setdiff(currtags, rmtags)
                    get("this", thisEnv)$set(name, tags = currtags)
                }
        },
        
        set = function(name, obj=NULL, newname=NULL, description=NULL,
            tags=NULL, src=NULL, depends=NULL, addtags=NULL, URL=NULL)
        {
            checkIndexUnchanged()                    
            
            if(missing(name) | (missing(newname) & missing(obj) & missing(description) &
                                 missing(tags) & missing(addtags)  & missing(src) & missing(URL)))
                stop("You must provide name and one of: obj, description, tags or addtags, src.")
            if(!missing(tags) & !missing(addtags))
                stop("You can not specify both tags and addtags.")
            
            if(checkName(name))
                handleErr("ID_NOT_FOUND", name)

            w <- findEntryIndex(name)
            entr <- entries[[w]]

            entr$timestamp <- Sys.time()
            if(!is.null(newname))
                entr$name <- newname
            if(!is.null(description))
                entr$description <- description
            if(!is.null(tags)) {                
                entr$tags <- checkTags(tags)
            }
            if(!is.null(addtags)) {
                entr$tags <- unique(c(entr$tags, checkTags(addtags, name)))
            }
            if(!is.null(src))
                entr$src <- src

            if(!missing(URL))
                entr$URL <- URL

            if(!is.null(obj)) {
                newinfo <- setData(entr$name, obj, isAttachment(entr$name))
                entr$dump <- newinfo[["dump"]]
                entr$size <- newinfo[["size"]]
                entr$checksum <- newinfo[["checksum"]]
                entr$dims <- newinfo[["dims"]]
            }
            
      entries[[w]] <- entr
      assign("entries", entries, thisEnv)
      storeIndex()
    },


        
        attach = function(filepath, description, tags, src=NULL, replace=F, to=NULL)
        {
            get("this", thisEnv)$put(filepath, basename(filepath),
                                     description, tags, src, replace=replace, asattach=T, to=to)
        },


        stash = function(object, rename = deparse(substitute(object)))
        {
            name <- deparse(substitute(object))
                if(!stopOnEmpty(T)){
                    e <- getEntry(rename)
                    if(!is.null(e))
                        if(!"stash" %in% e$tags)
                            stop(paste("A non-stash entry by the same name already exists,",
                                       "try setting the rename parameter."))
                }
            
            get("this", thisEnv)$put(object, rename, "Stashed object",
                                     c("stash", "hide"), replace=T)
        },

        stashclear = function(force=F)
        {
            get("this", thisEnv)$rm(tags=c("stash", "hide"), force=force)
        },

        pull = function(name, replace=F) {
            e <- getEntry(name)
            if(is.null(e$URL))
                handleErr("NO_URL", name)
            if(file.exists(e$dump) && !replace)
                handleErr("DATA_ALREADY_THERE", name)
            tf <- tempfile()
            download.file(e$URL, tf)

            if(isAttachment(name)) {
                get("this", thisEnv)$set(name, obj=tf)
            } else get("this", thisEnv)$set(name, obj=readRDS(tf))            
        },
        
      put = function(obj, name, description, tags, src=NULL,                       
                     depends=NULL, replace=F, notes=NULL, asattach=F,
                     to=NULL, addversion=F, URL=NULL)
        {
            checkIndexUnchanged()

            if(addversion)
                stop("addversion is deprecated, use replace=\"addversion\"")
            
            if(replace == "addversion") {
                ## This code is to cope with new interface after
                ## removing addversion parameter
                addversion = T 
                replace = F
            }
            
            if(missing(obj) | missing(name) | missing(description) | missing(tags))
                stop("You must provide all of: obj, name, description, tags.")
            
            
            if(!is.null(to))
                asattach <- T

            if(asattach)
                if(!file.exists(obj))
                    handleErr("ATTACHMENT_FILE_NOT_FOUND")
            
            if(name == "repo")
                handleErr("ID_RESERVED")

            notexist <- checkName(name)
            if(!notexist & !replace & !addversion)
                handleErr("ID_EXISTING", name)
            
            if(!asattach) {
                if(!is.null(dim(obj)))
                    dims <- dim(obj) else
                dims <- length(obj)
            } else {
                dims <- NULL
                tags <- unique(c(tags, "attachment"))
            }

            if(!is.null(src)) 
                stopOnNotFound(src)
            
            if(!is.null(to))
                stopOnNotFound(to)

            ## if(is.null(src))
            ##     src <- NA

            storedfrom <- src
            
            repoE <- list(name = name,
                          description = description,
                          tags = tags,
                          class = class(obj),
                          dims = dims,
                          timestamp = Sys.time(),
                          dump = NULL,
                          size = NULL,
                          checksum = NULL,
                          source = src,
                          depends = depends,
                          attachedto = to,
                          URL = URL)            
            
            if(!notexist & addversion) {
                newname <- checkVersions(name)$new
                get("this", thisEnv)$set(name, newname=newname)
                get("this", thisEnv)$tag(newname, "hide")
            }            

            entr <- get("entries", thisEnv)

            if(!notexist & replace) {
                ei <- findEntryIndex(name)
                rmData(name, "temp")
                oldEntr <- entr[[ei]]
            } else ei <- length(entries)+1
            
            tryCatch({
                fdata <- get("storeData", thisEnv)(name, obj, asattach)
            }, error = function(e) {
                print(e)
                if(!notexist & replace) 
                    rmData(name, "undo")
                stop("Error writing data.")
            }, finally = {
                repoE["size"] <- fdata[["size"]]
                repoE["checksum"] <- md5sum(path.expand(fdata[["path"]]))
                repoE["dump"] <- relativePath(fdata[["path"]])

                ## rmData must be called before overwriting the old
                ## entry (particularly the dump field)
                if(!notexist & replace)
                    rmData(name, "finalize")

                entr[[ei]] <- repoE
                assign("entries", entr, thisEnv)
                get("storeIndex", thisEnv)()
                
            }
            )

            if(asattach)
                get("this", thisEnv)$tag(name, "hide")
        },

      
      cpanel=function()
      {
          repo_cpanel(get("this", thisEnv)$root())
      },
        
        test=function()
        {
            print(ls(envir=thisEnv))
        },

        append = function(id, txtorfunc)
        {
          checkIndexUnchanged()
                                
          notexist <- checkName(id)
          if(notexist)
              stop("Identifier not found.")

          if(class(txtorfunc)=="function")
              txtorfunc <- paste0("\n",
                                  paste(deparse(txtorfunc), collapse="\n"),
                                  "\n")

          if(class(txtorfunc)!="character")
              stop("txtorfunc must be an object of class function or character")
          
          #e <- findEntryIndex(id)
          curobj <- this$get(id)
          this$set(id, obj=paste0(curobj, txtorfunc))
        },        
        
        root = function()
        {
            return(get("root",thisEnv))
        }

        )


    root <- normalizePath(root, mustWork=F)
    REPOFNAME <- "R_repo.RDS"
    repofile <- file.path(root, REPOFNAME)
    thisEnv <- environment()
    class(me) <- append(class(me),"repo")
    assign('this', me, envir=thisEnv)
    assign('entries', list(), envir=thisEnv)
    
    if(file.exists(repofile))
        {
            message(paste0("Found repo index in \"",
                           repofile, "\"."))
            assign("entries", readRDS(repofile), thisEnv)
        } else {

            if(!file.exists(root))
                {
                    if(force)
                        n <- "yes" else {
                            cat(paste0(
                                "Repo root \"", get("root",thisEnv),
                                "\" does not exist. Create it? "))
                            n <- readline("Type \"yes\" to proceed: ")
                        }
                    if(tolower(n) == "yes") {
                        dir.create(root)
                        message("Repo root created.")
                    } else message("Nothing done.")
                }

            storeIndex()
            message("Repo created.")
        }

    indexMD5 <- md5sum(repofile)
    
    return(me)
}

