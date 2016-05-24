#' Import DHS data.
#' 
#' This step by step function guides users to import data from a Demographic and Health Survey (DHS) 
#' and create an object of class \code{\link[=prevR-class]{prevR}}.
#' 
#' @param file.sav DHS data (one individual per line) in SPSS format (.sav), 
#' downloaded from \url{http://www.dhsprogram.com/}. Could also be directly a data.frame.
#' @param file.dbf GPS position of clusters in DATABASE format (.dbf), downloaded from 
#' \url{http://www.dhsprogram.com/}.  Could also be directly a data.frame.
#' 
#' @note If you don't provide the precise path of files, \R will check the working directory 
#' (see \code{\link[base]{setwd}}). To specify the file path, see \code{\link[base]{file.path}}.
#' 
#' This function was developed specifically for importing DHS. 
#' For a generic function for creating an object of class \code{\link[=prevR-class]{prevR}}, 
#' see \code{\link{as.prevR}}.
#' 
#' @seealso \code{\link{as.prevR}}, \code{\link{prevR-class}}.
#' 
#' @examples 
#' \dontrun{
#' import.dhs("data.sav", "gps.dbf")
#' }
#' 
#' @keywords manip
#' @export

import.dhs  = function(file.sav, file.dbf){
  make.clust.dbf <-function (file, ind) {
    if(is.data.frame(file))
      temp.clust <- file
    else {
      temp.clust <- foreign::read.dbf(file)
      if(is.data.frame(temp.clust$dbf)) temp.clust  = temp.clust$dbf
    }
    temp.var <- attr(temp.clust, "names")
    message("A window will open presenting the data contained in the file. Thank you to identify the following variables: \n- Cluster number (needed) \n- Longitude (decimal format in degrees, needed) \n- Latitude (decimal format in degrees, needed) \n- Clusters with missing coordinates (optional) \n- Clusters type (optional)\n Once the names of these variables identified, close the window so that the program can continue. \n\n Are you ready?",domain="R-prevR")
    menu(gettext("Yes",domain="R-prevR"))
    edit(temp.clust)
    ok <- 0
    while (ok != 1) {
      message("Please indicate the following variables:\n",domain="R-prevR")
      message("* Cluster number (usually called DHSCLUST):",domain="R-prevR")
      nb.clust <- menu(temp.var)
      message("* Longitude (decimal value, usually called LONGNUM):",domain="R-prevR")
      long <- menu(temp.var)
      message("* Latitude (decimal value, usually called LATNUM):",domain="R-prevR")
      lat <- menu(temp.var)
      message("* Clusters with missing coordinates (optionnal, usually called SOURCE, type 0 if none):",domain="R-prevR")
      c.source <- menu(temp.var)
      message("* Type of cluster (optionnal, type 0 if none):",domain="R-prevR")
      c.type <- menu(temp.var)
      message("\n----------------------------------------------------\n")
      message("Please check the following informations:\n",domain="R-prevR")
      message("* Cluster number:",domain="R-prevR")
      message(temp.var[nb.clust], if (nb.clust == 0)  gettext("Not available - WARNING: this variable must be specified!",domain="R-prevR"))
      message("* Longitude (decimal value):",domain="R-prevR")
      message(temp.var[long], if (long == 0) gettext("Not available - WARNING: this variable must be specified!",domain="R-prevR"))
      message("* Latitude (decimal value):",domain="R-prevR")
      message(temp.var[lat], if (lat == 0) gettext("Not available - WARNING: this variable must be specified!",domain="R-prevR"))
      message("* Clusters with missing coordinates (optionnal):",domain="R-prevR")
      message(temp.var[c.source], if (c.source == 0) gettext("Not available",domain="R-prevR"))
      message("* Type of cluster (optionnal):",domain="R-prevR")
      message(temp.var[c.type], if (c.type == 0) gettext("Not available",domain="R-prevR"))
      message("\n----------------------------------------------------\n")
      if (nb.clust == 0) 
        message("WARNING: cluster number must be specified!",domain="R-prevR")
      if (long == 0) 
        message("WARNING: longitude must be specified!",domain="R-prevR")
      if (lat == 0) 
        message("WARNING: latitude number must be specified!",domain="R-prevR")
      if (nb.clust == 0 | long == 0 | lat == 0) {
        alarm()
        message("\n----------------------------------------------------\n")
        message("WARNING: some problems were found (see above). You have to start again. Are you ready?",domain="R-prevR")
        menu(gettext("Yes",domain="R-prevR"))
      }
      else {
        message("Are these data correct?",domain="R-prevR")
        ok <- menu(gettext(c("Yes","No"),domain="R-prevR"))
      }
    }
    clust <- data.frame(id = temp.clust[nb.clust], x = temp.clust[long],  y = temp.clust[lat], c.type = temp.clust[c.type], c.source = temp.clust[c.source])
    colNames   = c("id","x","y","c.type","c.source")[c(T, T , T, c.type!=0, c.source!=0)]
    names(clust) = colNames
    
    if (is.numeric(clust$id))
      clust$id <- as.integer(clust$id)
    if (!is.integer(clust$id))
      clust$id <- as.integer(as.character(clust$id))
    if (!is.numeric(clust$x)) 
      clust$x <- as.numeric(as.character(clust$x))
    if (!is.numeric(clust$y)) 
      clust$y <- as.numeric(as.character(clust$y))
    if (c.type != 0)
      if (!is.factor(clust$c.type))
        clust$c.type <- factor(clust$c.type)
    if (c.source != 0)
      if (!is.factor(clust$c.source))
        clust$c.type <- factor(clust$c.source)
    
    #Deleting clusters with missing coordinates
    if (c.source != 0) {
      message("\n----------------------------------------------------\n")
      message("Please select the modality corresponding to missing coordinates (usually coded MIS).",domain="R-prevR")
      modalites <- attr(clust$c.source, "levels")
      c.mis <- select.list(modalites, multiple = TRUE, title = gettext("Missing coordinates",domain="R-prevR"))
      for (i in 1:length(c.mis)) {
        clust$x[clust$c.source == c.mis[i]] <- NA
        clust$y[clust$c.source == c.mis[i]] <- NA
      }
      n.missing <- length(clust[is.na(clust$x)|is.na(clust$y),'id'])
      sprintf(ngettext(n.missing,"%s cluster was deleted due to missing values.","%s clusters were deleted due to missing values.",domain="R-prevR"),n.missing) -> warning.mess
      warning(warning.mess, call.=F)
      clust <- clust[!is.na(clust$x)&!is.na(clust$y),]
      # Deleting c.source, not needed anymore
      if (c.type != 0) clust <- clust[c("id","x","y","c.type")] else clust <- clust[c("id","x","y")]
    }
    
    # Checking clusters with x==0 and y==0
    if (c.source==0) {
      n.missing <- length(clust[clust$x==0 & clust$y==0,'id'])
      sprintf(ngettext(n.missing,"%s cluster has latitude and longitude equal to 0. You should check if there is a variable in the dataset indicating clusters with missing coordinates.","%s clusters have latitude and longitude equal to 0. You should check if there is a variable in the dataset indicating clusters with missing coordinates.",domain="R-prevR"),n.missing) -> warning.mess
      warning(warning.mess, call.=F)
    }
    
    # Deleting individuals from a cluster who was deleted
    ind = ind[ind$id %in% clust$id,]
    # Partie merge en clusters et donnees GPS     
    clust.ind = merge(clust,ind,by="id")
    clust.ind = na.omit(clust.ind)
    sp        = split(clust.ind,clust.ind$id)
    id        = names(sp)
    n         = sapply(sp,nrow)
    pos       = sapply(sp, function(df) length(df$case[df$case == "Positive"]))
    x         = sapply(sp, function(df) unique(df$x))
    y         = sapply(sp, function(df) unique(df$y))
    clusters = data.frame(id = id, x =x, y = y, n = n, pos = pos)
    if(is.element("c.type",names(clust.ind))){
      c.type    = sapply(sp, function(df) unique(df$c.type))
      clusters  = cbind(clusters, c.type = c.type)
    }
    if(is.element("wcase",names(clust.ind))){
      wn        = sapply(sp, function(df) sum(df$wcase))
      wpos      = sapply(sp, function(df) sum(df$wcase[df$case == "Positive"]))
      clusters  = cbind(clusters, wn = wn, wpos = wpos)
    }
    return(clusters)
  }
  
  ################################################################################
  ################################################################################
  
  make.ind.spss <- function (file) {
    if (is.data.frame(file))
      temp.ind <- file
    else
      temp.ind <- foreign::read.spss(file, use.value.labels = TRUE, to.data.frame = TRUE)
    temp.var <- paste(attr(temp.ind, "names"), attr(temp.ind, "variable.labels"), sep = " - ")
    ok <- 0
    while (ok != 1) {
      message("1.1 VARIABLES SELECTION\n",domain="R-prevR")
      
      message("* Cluster number (0 if not, will be calculated from individual identification number):",domain="R-prevR")
      clust <- menu(temp.var)
      
      id = 0
      if(clust == 0){
        message("* Individual identification number (0 if not):",domain="R-prevR")
        id <- menu(temp.var)
        if(id == 0){
          message("* Individual identification number or Cluster number must be specified.",domain="R-prevR")
          next
        }
      }
      
      message("* Analyzed variable:",domain="R-prevR")
      result <- menu(temp.var)
      
      message("* Statistical weight (0 if not. All persons will have the same weight of 1.):",domain="R-prevR")
      weight <- menu(temp.var)
      
      message("\n----------------------------------------------------\n")
      message("Please check the following informations:\n",domain="R-prevR")
      message("* Individual identification number:",domain="R-prevR")
      message(temp.var[id], if (id == 0) gettext("Not available",domain="R-prevR"))
      message("* Cluster number:",domain="R-prevR")
      message(temp.var[clust], if (clust == 0) gettext("Not available - It will be calculated with individual identification number.",domain="R-prevR"))
      message("* Analyzed variable:",domain="R-prevR")
      message(temp.var[result], if (result == 0) gettext("Not available - WARNING: this variable must be specified!",domain="R-prevR"))
      message("* Statistical weight:",domain="R-prevR")
      message(temp.var[weight], if (weight == 0) gettext("Not available - All persons will have a weight of 1.",domain="R-prevR"))
      message("\n----------------------------------------------------\n")
      
      if (id == 0 & clust == 0) 
        message("WARNING: if cluster number not specified, you need to specify individual identification number!!",domain="R-prevR")
      
      if (result == 0) 
        message("WARNING: you have to specify the analyzed variable!",domain="R-prevR")
      
      
      if ((id == 0 & clust == 0) | result == 0) {
        alarm()
        message("\n----------------------------------------------------\n")
        message("WARNING: some problems were found (see above). You have to start again. Are you ready?",domain="R-prevR")
        menu(gettext("Yes",domain="R-prevR"))
      }
      else {
        message("Are these data correct?",domain="R-prevR")
        ok <- menu(c(gettext("Yes",domain="R-prevR"), gettext("No",domain="R-prevR")))
      }
    }
    
    ind        = data.frame(temp.ind[id], temp.ind[clust], temp.ind[result], temp.ind[weight])  
    colNames   = c("id","cluster","original.result","weight")[c(id!=0, clust!= 0, result!=0, weight != 0)]
    names(ind) = colNames
    
    ok <- 0
    message("\n----------------------------------------------------\n")
    message("Please specify the modalities corresponding to a positive result (the analysed phenomenon occured), a negative result (not occured) and an undetermined result (considered as a missing value).",domain="R-prevR")
    modalites <- attr(ind$original.result, "levels")
    
    
    while (ok != 1) {
      pos <- select.list(modalites, multiple = TRUE, title = gettext("Positive result",domain="R-prevR"))
      neg <- select.list(modalites, multiple = TRUE, title = gettext("Negative result",domain="R-prevR"))
      und <- select.list(modalites, multiple = TRUE, title = gettext("Undetermined result",domain="R-prevR"))
      message("\n----------------------------------------------------\n")
      if (length(neg) + length(pos) + length(und) != length(modalites)) {
        alarm()
        message("WARNING: You specified the same modality two times or you forgot one. Please start again.\nAre you ready?",domain="R-prevR")
        menu(gettext("Yes",domain="R-prevR"))
      }
      else {
        message("Please check the following informations:\n",domain="R-prevR")
        message("\n* Positive result:",domain="R-prevR")
        message("- ", paste(pos, collapse = "\n- "))
        message("\n* Negative result:",domain="R-prevR")
        message("- ", paste(neg, collapse = "\n- "))
        message("\n* Undetermined result:",domain="R-prevR")
        message("- ", paste(und, collapse = "\n- "))
        message("\nAre these data correct?",domain="R-prevR")
        ok <- menu(gettext(c("Yes", "No"),domain="R-prevR"))
      }
    }
    
    
    ind$result <- NA
    for (i in 1:length(neg)) ind$result[ind$original.result == neg[i]] <- 1
    for (i in 1:length(pos)) ind$result[ind$original.result == pos[i]] <- 2
    ind$result <- factor(ind$result,levels = c(1,2), labels = c("Negative","Positive"))
    
    
    if (weight!=0) {
      message("\n----------------------------------------------------\n")
      message("Often, in DHS, the weight variable have to be divided by a factor, usually 1'000'000.",domain="R-prevR")
      message(gettextf("The mean value of the weight variable is %.2f.",mean(ind$weight[ind$weight > 0]),domain="R-prevR"))
      message("\nIf this value is close to 1, the variable does not have to be modified.\nIf it is close to 1'000'000, then it must be divided by this factor.\nElse consult the survey documentation.",domain="R-prevR")
      message("\nDoes the weight variable have to be divided by a factor?",domain="R-prevR")
      choix <- menu(gettext(c("No modification", "Divided by 1'000'000","Divided by an other factor"),domain="R-prevR"))
      if (choix > 1) {
        if (choix == 2) 
          division.factor <- 1e+06
        else {
          message("Wich factor?",domain="R-prevR")
          division.factor <- 1
          division.factor <- as.integer(de(division.factor))
        }
        ind$weight <- ind$weight/division.factor
      }
    }
    
    
    if (clust != 0) {
      if (is.character(ind$cluster)) 
        ind$cluster <- as.integer(ind$cluster)
      if (is.factor(ind$cluster)) 
        ind$cluster <- as.integer(as.character(ind$cluster))
    } else {
      if (!is.character(ind$id)) 
        ind$id <- as.character(ind$id)
      exemple1 <- ind$id[1]
      exemple2 <- ind$id[length(ind$id)%/%2]
      exemple3 <- ind$id[length(ind$id)]
      message("\n----------------------------------------------------\n")
      message("The cluster variable is not specified. It must be then calculated from the individual identification number. Usually in DHS, the cluster number corresponds to the first three digits of the individual identification number. You can see here three individual identification numbers selected at the beginning, the medium and the end of the file:",domain="R-prevR")
      message("- Individual identification number 1: ", exemple1,domain="R-prevR")
      message("- Individual identification number 2: ", exemple2,domain="R-prevR")
      message("- Individual identification number 3: ", exemple3,domain="R-prevR")
      message("\nLocate for each of them the cluster number.\nAmong the various proposals below, which extracts the good cluster numbers?\n",domain="R-prevR")
      choix <- c()
      for (i in 1:(nchar(exemple1) - 2)) choix <- c(choix, paste(substr(exemple1, i, i + 2),substr(exemple2, i, i + 2),substr(exemple3, i, i + 2),sep=" / "))
      lim = 0
      while(lim==0)
        lim <- menu(choix, graphics = TRUE, title=gettext("Select the correct cluster numbers:",domain="R-prevR"))
      ind$cluster <- as.integer(substr(ind$id, lim, lim + 2))
    }
    df = data.frame(id = ind$cluster, case = ind$result)
    if (weight!=0) df = cbind(df, wcase = ind$weight) 
    return(df)
  }
  
  ################################################################################
  ################################################################################
  
  message("\n\n\nSTEP 1/3: DEFINE ANALYZED VARIABLE\n\n\n",domain="R-prevR")
  ind = make.ind.spss(file.sav)
  message("\n\n\nSTEP 2/3: CLUSTERS INFORMATION\n\n\n",domain="R-prevR")
  clusters = make.clust.dbf(file.dbf,ind)
  message("\n\n\nSTEP 3/3: BOUNDARY OF THE COUNTRY\n\n\n",domain="R-prevR")
  message("Do you want to define a boundary?",domain="R-prevR")
  bound = menu(gettext(c("Yes","No"),domain="R-prevR"))
  boundary = NULL
  if(bound==1) boundary = create.boundary(multiple = F)
  col = names(clusters)
  names(col) = col
  as.prevR(clusters,col,boundary = boundary,proj = "+proj=longlat +ellps=WGS84")
}