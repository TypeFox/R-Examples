rattleInfo <- function(all.dependencies=FALSE,
                       include.not.installed=FALSE,
                       include.not.available=FALSE,
                       include.libpath=FALSE)
{

  # TODO: Add in support for BIOC

  cran.repos <- "http://cran.rstudio.org"
  bioc.repos <- ""

  # Using installed.packages() can be a "very slow way to find
  # information on one or a small number of packages" (Brian Riply
  # 2012). This is stated in the man page and I am very aware of
  # it. Brian also note: "In addition, many of you are using it to
  # find out if a package is installed, when actually you want to know
  # if it is usable (it might for example be installed for a different
  # architecture or require a later version of R), for which you need
  # to use require()." This was particularly relevant within
  # packageIsAvailable() and there I use a better way of checking for
  # an installed package.  Here I think it might still remain
  # appropriate to use installed.packages().
  
  iv <- utils::installed.packages()
  av <- available.packages(contriburl=contrib.url(cran.repos))
  have.av <- nrow(av) != 0
  # not a cran repos bv <- available.packages(contriburl=contrib.url(cran.repos))

  riv <- iv["rattle", "Version"]
  if (have.av) rav <- av["rattle", "Version"]
  
  cat(sprintf("Rattle: version %s", riv))
  if (have.av) cat(sprintf(" CRAN %s", rav))
  cat("\n")

  # Record the packages that can be upgraded

  up <- if (have.av && compareVersion(rav, riv) == 1) "rattle" else NULL
    
  cat(sprintf("%s\n", sub(" version", ": version", version$version.string)))

  cat("\n")
  si <- Sys.info()
  for (i in seq_along(si))
    cat(sprintf("%s%s: %s\n", toupper(substr(names(si)[i], 1, 1)),
                substring(names(si)[i], 2), si[i]))

  cat("\nInstalled Dependencies\n")

  deps2vec <- function(deps)
  {
    if (is.na(deps)) return(NULL)
    strsplit(gsub("\\n", " ", gsub(' ?\\([^\\)]+\\)', '', deps)), ", ?")[[1]]
  }
    
  if (all.dependencies)
  {
    if (! "pkgDepTools" %in% rownames(iv))
    {
      source("http://bioconductor.org/biocLite.R")
      pkg <- "pkgDepTools"
      biocLite("pkgDepTools")
    }
    if (! "Rgraphviz" %in% rownames(iv))
    {
      source("http://bioconductor.org/biocLite.R")
      biocLite("Rgraphviz")
    }

    # 150711 There does not seem to be a way to get both suggest and
    # depend links using pkgDepTools::makeDepGraph which I used to
    # deploy here. It's either one or the other. Rattle only has
    # suggests links. So I want to get what Rattle suggests and then
    # find all the depends in cran.deps as the packages that are
    # reported on. Instead of going through the repository and build a
    # dependency graph, we've already dounloaded the available package
    # information so use it here instead.

    pkg.deps <- function(pkg, pkgs, av)
    {
      if (pkg %in% pkgs) return(pkgs)

      if (! pkg %in% rownames(av)) return(c(pkg, pkgs))

      for (p in union(deps2vec(av[pkg, "Suggests"]), deps2vec(av[pkg, "Depends"])))
      {
        pkgs <- pkg.deps(p, union(pkg, pkgs), av)
      }
      return(union(pkg, pkgs))
    }
    
    if (have.av)
      deps <- pkg.deps("rattle", NULL, av)
    else
      deps <- pkg.deps("rattle", NULL, iv)
  }    
  else
    deps <- union(deps2vec(iv["rattle", "Depends"]), deps2vec(iv["rattle", "Suggests"]))

  for (p in sort(setdiff(deps, 'rattle')))
  {
    if (have.av && ! p %in% rownames(av))
    {
      if (include.not.available) cat(sprintf("%s: not available\n", p))
    }
    else if (! p %in% rownames(iv))
    {
      if (include.not.installed) cat(sprintf("%s: not installed\n", p))
    }
    else
      cat(sprintf("%s: version %s%s%s%s", p, iv[p,"Version"],
                  ifelse(have.av && compareVersion(av[p,"Version"], iv[p,"Version"]) == 1,
                         {
                           up <- c(up, p);
                           sprintf(" upgrade available %s", av[p,"Version"])
                         },
                         ""),
                  ifelse(include.libpath, paste("\t", iv[p,"LibPath"]), ""),
                  "\n"))
  }

  cat("\nThat was",
      if (include.not.available)
        length(deps)
      else
        sum(sapply(deps, function(p) p %in%
                     if (have.av && include.not.installed) rownames(av) else rownames(iv))),
      "packages.\n")
  
  if (! is.null(up))
  {
    cat(sprintf(paste('\nUpdate the packages with either',
                      'of the following commands:\n\n ',
                      '> install.packages(c("%s"))\n\n ',
                      '> install.packages(rattleInfo(%s%s%s%s%s%s%s))\n\n'),
                paste(strwrap(paste(up, collapse='", "'),
                              width=60, exdent=23), collapse="\n"),
                ifelse(all.dependencies, "all.dependencies=TRUE", ""),
                ifelse(all.dependencies &&
                       (include.not.installed ||
                        include.not.available ||
                        include.libpath), ", ", ""),
                ifelse(include.not.installed, "include.not.installed=TRUE", ""),
                ifelse(include.not.installed &&
                       (include.not.available ||
                        include.libpath), ", ", ""),
                ifelse(include.not.available, "include.not.available=TRUE", ""),
                ifelse(include.not.available &&
                       include.libpath, ", ", ""),
                ifelse(include.libpath, "include.libpath=TRUE", "")))
    if (isWindows() && "rattle" %in% up)
      cat("Detach rattle (and other attached packages) before updating:\n\n ",
          '> detach("rattle")\n\n')
    cat("Alternatively update all installed packages:\n\n ",
        '> update.packages()\n\n')

  }

  invisible(up)

}
