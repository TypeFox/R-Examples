#### 1- Construction functions ####
constLatex <- function(directory, filename = NULL, identifier = NULL, param = NULL, tabular = NULL, extra_text = NULL, 
                       subsection = NULL, label = NULL, 
                       width = 0.9, trim = c(0,0,0,0), plotPerPage = 3, 
                       width.legend = 0.35, trim.legend = c(0,0,0,0), 
                       title = "", date = "", author = "", verbose = TRUE){
  
  ls.text <- list()
  
  #### preparation ####
  
  ## dossiers ou sont present les graphiques a afficher
  res_init <- initDirPat_constLatex(directory = directory, param = param, identifier = identifier, verbose = verbose)
  
  directory.plot <- res_init$directory.plot
  n.directory.plot <- res_init$n.directory.plot
  names_dirs <- res_init$names_dirs
  param <- res_init$param
  identifier <- res_init$identifier
  n.identifier <- res_init$n.identifier
  
  ## gestion des sections, subsections, subsubsection, label
  res_init <- initSection_constLatex(directory.plot = directory.plot, names_dirs = names_dirs, param = param, 
                                     subsection = subsection, label = label)
  
  subsection <- res_init$subsection
  index_subsection <- res_init$index_subsection
  label <- res_init$label
  param <- res_init$param
  
  ## initialisation des parametres latex
  if(!is.list(trim)){
    trim <- lapply(1:100, function(x){trim})
  }
  if(length(width) == 1){
    width <- rep(width, 100)
  }
  
  ## initialisation de donnes cliniques
  if(!is.null(tabular)){
    
    if(class(tabular) != "list"){tabular <- list(tabular)}
    n.tabular <- length(tabular)
    
    ## conversion de vector en data.frame (si necessaire)
    tabular <- lapply(tabular, function(x){  
      if(is.null(dim(x))){ 
        res <- data.frame(matrix(NA, ncol = length(x), nrow = 1), stringsAsFactors = FALSE)
        res[1,] <- x ; names(res) <- names(x) ; res
      }else{x}}
    )
    
    ## conversion des factors en character    
    tabular <- lapply(tabular, function(x){ 
      df.tempo <- data.frame(rapply(x, as.character, classes = "factor", how = "replace"), stringsAsFactors = FALSE)
      names(df.tempo) <- names(x)
      return(df.tempo)
    })
    
  }
  
  ## initialisation du texte a afficher  
  if(!is.null(extra_text)){
    validDim_vector(value1 = unlist(extra_text), value2 = identifier, name1 = "extra_text", name2 = "identifier", type = "length", method = "constLatex")
  }
  
  ## display 
  if(verbose == TRUE){
    space.index_subsection <- initSpace(c("num", index_subsection))
    space.subsection <- initSpace(c("subsection", subsection))
    
    cat("%%%% num", space.index_subsection[1], " : subsection", space.subsection[1], " | directory \n", sep = "")
    cat("% (label) \n\n", sep = "")
    
    for(iter_subsection in 1:length(subsection)){
      cat("%%%% ", iter_subsection, space.index_subsection[iter_subsection + 1], " : ", subsection[iter_subsection], sep = "")
      
      cat(paste(space.subsection[iter_subsection + 1], " | /", names_dirs[index_subsection == iter_subsection], "\n", 
                "% label : \"", gsub("\n", "", label[[iter_subsection]], fixed = TRUE), "\" \n\n", 
                sep = ""), sep = "")      
    }
    
    cat("% graph of patient ")       
  }
  
  #### preambule ####
  text.preamble <- c("%endTrace \n \n \\documentclass[a4paper]{article} \n \n \n", 
                     "%%%%%%%%%%%%%%%%% Preamble %%%%%%%%%%%%%%%%% \n", 
                     "\\usepackage[utf8]{inputenc} \n", 
                     "\\usepackage{amssymb} \n", 
                     "\\usepackage{amsmath} \n", 
                     "\\usepackage{titlesec} \n", 
                     "\\usepackage{geometry} \n", 
                     "\\usepackage{enumitem} \n", 
                     "\\usepackage{graphicx} \n", 
                     "\\usepackage{color} \n", 
                     "\\usepackage{xspace} \n", 
                     "\\usepackage{hyperref} \n", 
                     "\\usepackage{caption} \n \n \n", 
                     "%%%% Margin %%%% \n", 
                     "\\geometry{\n", 
                     "a4paper, \n", 
                     "total = {210mm, 297mm}, \n", 
                     "left = 10mm, \n", 
                     "right = 10mm, \n", 
                     "top = 5mm, \n", 
                     "bottom = 12.5mm, \n", 
                     "}\n", 
                     "\\setlength{\\textfloatsep}{5pt} % space between last top float or first bottom float and the text \n", 
                     "\\setlength{\\intextsep}{5pt} % space left on top and bottom of an in-text float \n", 
                     "\\setlength{\\abovecaptionskip}{10.0pt} % space above caption \n", 
                     "\\setlength{\\belowcaptionskip}{0pt} % space below caption \n", 
                     "\\titlespacing\\section{0pt}{0ex}{5ex} % {left space}{above space}{below space} \n", 
                     "\\titlespacing\\subsection{25pt}{3ex}{1ex} % {left space}{above space}{below space} \n \n \n",                     
                     "%%%% Margin %%%% \n", 
                     paste("\\graphicspath{{./", utils::tail(strsplit(directory, split = "/", fixed = TRUE)[[1]], 1), "/}} \n \n \n", sep = ""), 
                     "%%%% Title %%%% \n", 
                     paste("\\title{", title, "} \n", sep = ""), 
                     paste("\\date{", date, "} \n", sep = ""), 
                     paste("\\author{", author, "} \n \n \n", sep = ""), 
                     "%%%%%%%%%%%%%%%%% End Preamble %%%%%%%%%%%%%%%%% \n \n \n"
  )
  
  #### begin doc ####
  text.begin <- c("%", "\n", 
                  "\\begin{document} \n \n ", 
                  "\\maketitle \n \n ", 
                  "\\setcounter{tocdepth}{1} \n", 
                  "\\tableofcontents \n", 
                  "\\setcounter{tocdepth}{3} \n", 
                  "\\clearpage \n"
  ) 
  
  #### boucle sur les patients ####
  for(iter_pat in 1:n.identifier){
    
    if(verbose == TRUE){cat(iter_pat, "(", identifier[iter_pat], ") ", sep = "")}
    
    #### definition des graphiques a afficher pour chaque patient 
    res_init <- initPlot_constLatex(directory = directory, names_dirs = names_dirs, directory.plot = directory.plot, tabular = tabular, 
                                    identifier = identifier[iter_pat], plotPerPage = plotPerPage)
    
    ls.plot <- res_init$ls.plot
    ls.newplot <- res_init$ls.newplot
    ls.endplot <- res_init$ls.endplot
    ls.legend <- res_init$ls.legend
    ls.slices <- res_init$ls.slices
    ls.tabular <- res_init$ls.tabular
    
    #### affichage clinique ####
    ls.text[[iter_pat]] <- paste("%", "\n \n \\section", paste("{patient ", identifier[iter_pat], "} \n \n ", sep = ""), sep = "")
    
    if(length(ls.tabular) > 0){ # donnees cliniques disponibles pour le patient
      
      for(iter_table in 1:length(ls.tabular)){
        
        if(ncol(ls.tabular[[iter_table]]) == 0){next}
        
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                 "\\begin{table}[!h] \n", 
                                 "\\centering \n", 
                                 paste("\\begin{tabular}{|", paste(rep("c", ncol(ls.tabular[[iter_table]])), collapse = ""), "|} \n", sep = ""), 
                                 paste(paste(names(ls.tabular[[iter_table]]), collapse = " & "), " \\\\ \\hline \n", sep = ""))
        
        for(iter_rowtable in 1:nrow(ls.tabular[[iter_table]])){
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                   paste(paste(ls.tabular[[iter_table]][iter_rowtable,], collapse = " & "), " \\\\  \n", sep = ""))
        }
        
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                 "\\end{tabular} \n", 
                                 "\\end{table} \n \n", 
                                 "\\smallskip  \n \n")
        
      }
    }
    
    #### extra text ####
    if(!is.null(extra_text)){
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], extra_text[[iter_pat]]) 
    }
    
    if(length(ls.tabular) > 0 || !is.null(extra_text)){ # si donnees cliniques 
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], '\n \n \\clearpage \n \n')    
    }else{ # sinon si la page est pleine
      if(unlist(lapply(ls.plot, function(x){ if(all(is.na(x))){NULL}else{length(x)}} ))[1] >= plotPerPage){
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], '\n \n \\vspace{-5ex} \n \n')
      }
    }
    
    #### affichage des images ####
    count.subsection <- 0
    
    for(iter_subsection in 1:length(subsection)){
      
      ## pas de graphique on saute la subsection
      if(all(is.na(ls.plot[[iter_subsection]]))){next}else{count.subsection <- count.subsection + 1} 
      
      ## debut de la subsection
      ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                               paste("\\subsection", paste("{", subsection[iter_subsection], "} \n  ", sep = ""), sep = "")
      )
      
      ## initialisation
      iter_subsubsection <- 0
      n.plot_tempo <- length(ls.plot[[iter_subsection]])
      count_plot <- 1
      count_figure <- 0
      
      #       if(length(ls.tabular) == 0 && is.null(extra_text) && count.subsection == 1 && n.plot_tempo >= plotPerPage){ # si section sur la meme page que subsection et premiere subsection et la figure arrive en bas
      #         ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], '\n \n \\vspace{-3ex} \n \n')
      #       }
      
      for(iter_plot in 1:n.plot_tempo){ #### gestion des graphiques
        
        ## debut d une nouvelle page
        if(ls.newplot[[iter_subsection]][iter_plot] == 1){
          
          ## debut d une nouvelle figure
          count_figure <- count_figure + 1 
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                   "\\begin{figure}[!h] \n", 
                                   "\\centering \n"
          ) 
          
        } 
        
        ## parametres latex entre legend et image
        if(ls.legend[[iter_subsection]][iter_plot] == TRUE){
          trim_tempo <- trim.legend
          width_tempo <- width.legend
        }else{
          trim_tempo <- trim[[count_plot]]
          width_tempo <- width[[count_plot]]
        }
        
        ## affichage de l image 
        ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                 paste("\\includegraphics[trim =  ", trim_tempo[1], "mm ", trim_tempo[2], "mm ", trim_tempo[1], "mm ", trim_tempo[1], "mm, clip, width = ", width_tempo, "\\textwidth]", sep = ""), 
                                 paste("{", ls.plot[[iter_subsection]][iter_plot], "} \n", sep = "")
        )
        
        ## affichage de la legende
        if(ls.endplot[[iter_subsection]][iter_plot] == 1){
          
          ls.text[[iter_pat]] <- c(ls.text[[iter_pat]], 
                                   paste("\\caption[foo bar]{Patient ", identifier[iter_pat], " ", ls.slices[[iter_subsection]][count_figure], " - ", label[[iter_subsection]], "} \n", sep = ""), 
                                   #paste("\\label{fig:Sweave_", identifier[iter_pat], "_", param[iter_subsection], "} \n", sep = ""), 
                                   "\\vspace{-1cm} \n \n", 
                                   "\\end{figure} \n \n",                                   
                                   '\\clearpage \n \n'
          )          
        }
        count_plot <- count_plot + 1         
        
        
      } # iteration plot
    } # iteration subsection
  } # iteration patient
  if(verbose == TRUE){cat("\n")}
  
  #### end doc ####
  text.end <- c("%", "\n \n \\end{document}  ")    
  
  #### create Latex file ####
  if(!is.null(filename)){
    
    file.create(paste(filename, "tex", sep = "."))
    
    sink(paste(filename, "tex", sep = "."))
    cat(text.preamble, sep = "")
    cat(text.begin, sep = "") 
    for(iter_list in 1:length(ls.text)){
      cat(ls.text[[iter_list]], sep = "")
    }
    cat(text.end, sep = "")
    sink()
    
  }
  
  #### export ####
  
  return(list(text.preamble = text.preamble, 
              text.begin = text.begin, 
              ls.text = ls.text, 
              text.end = text.end)
  )
}


#### 2- Calculation functions ####

calcAUPRC <- function(x, y, subdivisions = 10000, performance = NULL, ci = TRUE, alpha = 0.05, 
                      method = "Kronrod", reltol = .Machine$double.eps^0.25){
  
  
  #### tests
  ## packages
  # initPackage("ROCR", method =  "calcAUPRC")	
  
  ## arguments
  if (optionsMRIaggr("checkArguments")) {
  
  validInteger(value = subdivisions, validLength = 1, min = 0, method = "calcAUPRC")
  validLogical(value = ci, validLength = 1, method = "calcAUPRC")
  if(!is.null(performance)){
    validClass(value = performance, validClass = "performance", superClasses = FALSE, method = "calcAUPRC")
    
    if( any(c("Precision", "Recall") %in% c(performance@x.name, performance@y.name) == FALSE) ){
      stop("calcAUPRC : wrong specification of \'performance\' \n", 
           "\'performance\' must contains \"Precision\" and \"Recall\" performance measures \n", 
           "measures in \'performance\' : \"", performance@x.name, "\" \"", performance@y.name, "\" \n")   
    }
  }
  validNumeric(value = alpha, validLength = 1, min = 0, max = 1, method = "calcAUPRC")
  validCharacter(value = method, validLength = 1, 
                 validValues = c("integrate", "Kronrod", "Richardson", "Clenshaw", "Simpson", "Romberg"),
                 refuse.NULL = TRUE, method = "calcAUPRC")
  if(method != "integrate"){
    initPackage("pracma", argument =  "method != \"integrate\"", method =  "calcAUPRC")	
  }
  validNumeric(value = reltol, validLength = 1, min = 0, method = "calcAUPRC")
  
  }
  
  #### initialisation
  if(is.null(performance)){    
    performance <- ROCR::performance(ROCR::prediction(x, y), x.measure =  "rec", measure =  "prec")
  }
  
  M.PR <- cbind(performance@x.values[[1]], 
                performance@y.values[[1]])
  n <- nrow(M.PR)
  colnames(M.PR) <- c(performance@x.name,performance@y.name)
  M.PR <- M.PR[c(-1,-nrow(M.PR)),, drop = FALSE]
  
  # deal with multiple values at the same threshold
  M.PR <- cbind(Recall = unique(M.PR[,"Recall"]), 
                Precision = tapply(M.PR[,"Precision"], M.PR[,"Recall"], function(x){max(x)}))
  if(nrow(M.PR) == 1){return(M.PR[,"Precision"])} # cas avec un seuil 
  
  # conversion to a function 
  f  <- stats::approxfun(x = M.PR[,"Recall"], y = M.PR[,"Precision"])
  
  #### computation : integration
  if(method == "integrate"){
    f.int <- stats::integrate(f, min(M.PR[,"Recall"]), max(M.PR[,"Recall"]), subdivisions = subdivisions, rel.tol = reltol)$value
  }else{
    f.int <- pracma::integral(fun = f, xmin = min(M.PR[,"Recall"]), xmax = max(M.PR[,"Recall"]), 
                              method = method, reltol = reltol)
  }
  
  AUPRC <- f.int / (max(M.PR[,"Recall"]) - min(M.PR[,"Recall"]))
  
  if(ci == TRUE){
    eta <- log(AUPRC / (1 - AUPRC))
    tau <- 1 / sqrt(n * AUPRC * (1 - AUPRC))    
    
    inf_tempo <- exp(eta - stats::qnorm(p = 1 - alpha / 2) * tau)
    sup_tempo <- exp(eta + stats::qnorm(p = 1 - alpha / 2) * tau)
  }else{
    inf_tempo <- NA
    sup_tempo <- NA
  }
    
  #### export
  return(c(AUPRC = AUPRC, 
           IC_inf = inf_tempo / ( 1 + inf_tempo), 
           IC_sup = sup_tempo / ( 1 + sup_tempo)
          ))
  
}

EDK <-function(x, bandwidth, power = 2){ 
  1 / (2 * pi * bandwidth^2)^(1 / power) * exp( - (x / bandwidth)^power) 
}

logit<-function (x){
  return(log(x / (1 - x)))
}

inv.logit<-function (x){
  return(1 / (1 + exp(-x)))
}

# clone of the rtnorm function of the msm package
rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf) 
{
  if (length(n) > 1) 
    n <- length(n)
  mean <- rep(mean, length = n)
  sd <- rep(sd, length = n)
  lower <- rep(lower, length = n)
  upper <- rep(upper, length = n)
  lower <- (lower - mean) / sd
  upper <- (upper - mean) / sd
  ind <- seq(length = n)
  ret <- numeric(n)
  alg <- ifelse(lower > upper, -1, ifelse(((lower < 0 & upper == Inf) | (lower == -Inf & upper > 0) | (is.finite(lower) & 
                                                                                                         is.finite(upper) & (lower < 0) & (upper > 0) & (upper - 
                                                                                                                                                           lower > sqrt(2 * pi)))), 0, ifelse((lower >= 0 & (upper > 
                                                                                                                                                                                                               lower + 2 * sqrt(exp(1)) / (lower + sqrt(lower^2 + 4)) * 
                                                                                                                                                                                                               exp((lower * 2 - lower * sqrt(lower^2 + 4)) / 4))), 
                                                                                                                                                                                              1, ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1)) / (-upper + 
                                                                                                                                                                                                                                                              sqrt(upper^2 + 4)) * exp((upper * 2 - -upper * sqrt(upper^2 + 
                                                                                                                                                                                                                                                                                                                    4)) / 4)), 2, 3))))  
  ind.nan <- ind[alg == -1]
  ind.no <- ind[alg == 0]
  ind.expl <- ind[alg == 1]
  ind.expu <- ind[alg == 2]
  ind.u <- ind[alg == 3]
  ret[ind.nan] <- NaN
  while (length(ind.no) > 0) {
    y <- stats::rnorm(length(ind.no))
    done <- which(y >= lower[ind.no] & y <= upper[ind.no])
    ret[ind.no[done]] <- y[done]
    ind.no <- setdiff(ind.no, ind.no[done])
  }
  stopifnot(length(ind.no) == 0)
  while (length(ind.expl) > 0) {
    a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expl), a) + lower[ind.expl]
    u <- stats::runif(length(ind.expl))
    done <- which((u <= exp( - (z - a)^2 / 2)) & (z <= upper[ind.expl]))
    ret[ind.expl[done]] <- z[done]
    ind.expl <- setdiff(ind.expl, ind.expl[done])
  }
  stopifnot(length(ind.expl) == 0)
  while (length(ind.expu) > 0) {
    a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expu), a) - upper[ind.expu]
    u <- stats::runif(length(ind.expu))
    done <- which((u <= exp( - (z - a)^2 / 2)) & (z <= -lower[ind.expu]))
    ret[ind.expu[done]] <- -z[done]
    ind.expu <- setdiff(ind.expu, ind.expu[done])
  }
  stopifnot(length(ind.expu) == 0)
  while (length(ind.u) > 0) {
    z <- stats::runif(length(ind.u), lower[ind.u], upper[ind.u])
    rho <- ifelse(lower[ind.u] > 0, exp((lower[ind.u]^2 - 
                                           z^2) / 2), ifelse(upper[ind.u] < 0, exp((upper[ind.u]^2 - 
                                                                                      z^2) / 2), exp(-z^2 / 2)))
    u <- stats::runif(length(ind.u))
    done <- which(u <= rho)
    ret[ind.u[done]] <- z[done]
    ind.u <- setdiff(ind.u, ind.u[done])
  }
  stopifnot(length(ind.u) == 0)
  ret * sd + mean
}


# clone of the dtnorm function of the msm package
dtnorm <- function (x, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = FALSE) 
{
  ret <- numeric(length(x))
  ret[x < lower | x > upper] <- if (log) 
    -Inf
  else 0
  ind <- x >= lower & x <= upper
  if (any(ind)) {
    denom <- stats::pnorm(upper, mean, sd) - stats::pnorm(lower, mean, 
                                            sd)
    xtmp <- stats::dnorm(x, mean, sd, log)
    if (log) 
      xtmp <- xtmp - log(denom)
    else xtmp <- xtmp / denom
    ret[x >= lower & x <= upper] <- xtmp[ind]
  }
  ret
}

#### 3- initialisation functions ####

initDirPat_constLatex <- function(directory, param, identifier, verbose){ 
  
  ## definition de directory
  validPath(directory, name = "directory", method = "initDirPat_constLatex")
  directory.plot <- list.dirs(directory)[-1]
  
  n.directory.plot <- length(directory.plot)
  if(n.directory.plot == 0){
    stop("constLatex : no directory inside the root directory \n", 
         "\'directory\' must contain one directory per parameter \n", 
         "proposed directory : ", directory, "\n") 
  }  
  
  names_dirs <- unlist(lapply(strsplit(directory.plot, split = "/"), function(x){x[length(x)]}))
  
  if(verbose == TRUE){
    cat("%% ", n.directory.plot, " directories found in the root directory \n", sep = "")
  }
  
  ## dir restrected to param
  if(!is.null(param)){
    
    validCharacter(value = param, validLength = NULL, validValues = names_dirs, refuse.NULL = TRUE, method = "constLatex")
    
    directory.plot <- directory.plot[sapply(param, function(x){which(names_dirs == x)})]
    n.directory.plot <- length(directory.plot)
    names_dirs <- unlist(lapply(strsplit(directory.plot, split = "/"), function(x){x[length(x)]}))    
    if(verbose == TRUE){
      cat("% ", n.directory.plot, " directories will be considered according to \'param\' argument \n", sep = "")
    }
  }
  
  ## param
  if(is.null(param)){
    param <- names_dirs
  }
  
  if(length(param) != length(names_dirs)){
    stop("constLatex : \'param\' and \'names_dirs\' must have same length \n", 
         " length(param)  : ", length(param), " \n", 
         " length(names_dirs) : ", length(names_dirs), " \n", 
         " param : ", paste(param, collapse = " "), " \n", 
         " names_dirs : ", paste(names_dirs, collapse = " "), " \n")
  }
  
  
  ## identifiants
  Allidentifier <- unique(unlist(lapply(directory.plot, function(x){res <- strsplit(list.files(x), split = "_", fixed = TRUE);
  unique(unlist(lapply(res, "[", 1)))})))
  
  if(is.null(identifier)){
    identifier <- Allidentifier  
  }
  n.identifier <- length(identifier)
  
  if(verbose == TRUE){
    cat("%% ", length(Allidentifier), " identifiers found in the subdirectories \n", sep = "")
    cat("% ", n.identifier, " identifiers will be considered according to \'param\' argument \n \n", sep = "")
  }
  
  if(n.identifier == 0){
    stop("constLatex : no graphics inside the root directory \n", 
         "proposed directory : ", directory, "\n", 
         "proposed directory.plot : ", paste(directory.plot, collapse = "\n                     "), "\n")
  }  
  
  #### export ####
  return(list(directory.plot = directory.plot, 
              n.directory.plot = n.directory.plot, 
              names_dirs = names_dirs, 
              param = param, 
              identifier = identifier, 
              n.identifier = n.identifier))
  
}

initSection_constLatex <- function(directory.plot, names_dirs, param, subsection, label){
  
  #### initialisation de subsection avec les noms des dossiers figures
  if(is.null(subsection)){
    subsection <- gsub(pattern = "_", replacement = " ", x = names_dirs)
  }
  n.subsection <- length(subsection)
  
  index_subsection <- 1:n.subsection
  
  #### initialisation des legendes des graphiques avec les noms des dossiers figures
  
  if(is.null(label)){
    label <- gsub(pattern = "_", replacement = " ", x = names_dirs)
  }
  
  if(is.list(label) == TRUE){
    stop("constLatex : wrong specification of \'label\'\n", 
         "\'label\' must be a vector and not a list \n", 
         "is(label) : ", paste(is(label), collapse = " "), " \n") 
  }
  
  validDim_vector(value1 = param, value2 = label, name1 = "param", name2 = "label", type = "length", method = "constLatex")
  
  #### export ####
  return(list(subsection = subsection, 
              index_subsection = index_subsection, 
              label = label, 
              param = param))
}


initPlot_constLatex <- function(directory, names_dirs, directory.plot, identifier, tabular, 
                                plotPerPage){
  
  n.directory.plot <- length(directory.plot)
  
  ls.plot <- list() # list of the plot to display for each patient
  ls.newplot <- list() # list of the indicator for the creation of a new figure in latex
  ls.endplot <- list() # list of the indicator for the end of a figure in latex
  ls.legend <- list() # list of the legend files
  ls.slices <- list() # list of the first-last slices for each figure
  
  #### iteration sur les sous dossiers
  for(iter_dir in 1:n.directory.plot){ 
    
    split_tempo <- strsplit(list.files(directory.plot[iter_dir]), split = "_", fixed = TRUE)
    
    if(length(split_tempo) == 0){ # pas de fichier dans le repertoire 
      ls.legend[[iter_dir]] <- NA
      ls.plot[[iter_dir]] <- NA
      ls.newplot[[iter_dir]] <- NA 
      ls.endplot[[iter_dir]] <- NA 
      ls.slices[[iter_dir]] <- NA 
      next      
    } 
    
    index_plot <- which(unlist(lapply(split_tempo, 
                                      function(x){x[1] == identifier})))      
    n.index_plot <- length(index_plot) # nb de fichiers images du patient
    
    if(n.index_plot == 0){ # pas de fichier dans le repertoire correspondants au patient      
      ls.legend[[iter_dir]] <- NA
      ls.plot[[iter_dir]] <- NA
      ls.newplot[[iter_dir]] <- NA 
      ls.endplot[[iter_dir]] <- NA 
      ls.slices[[iter_dir]] <- NA 
      next      
    } 
    
    legend_tempo <- unlist(lapply(split_tempo[index_plot], 
                                  function(x){length(grep(pattern = "legend", x = x)) > 0}))
    ls.legend[[iter_dir]] <- legend_tempo
    
    if(length(index_plot) == 0){
      stop("constLatex : an legend file was found but with no corresponding image \n", 
           "directory : ", directory.plot[iter_dir], " \n", 
           "legend file : \"", list.files(directory.plot[iter_dir])[legend_tempo > 0], "\" \n") 
    }
    
    # test fichier multiplot
    position_slices <- unlist(lapply(split_tempo[index_plot], 
                                     function(x){grep(x, pattern = "slice", fixed = TRUE)}))
    
    if(length(position_slices) == 0){ # si ce n est pas un fichier multiplot
      ls.legend[[iter_dir]] <- rep(FALSE, n.index_plot)
      ls.plot[[iter_dir]] <- paste(utils::tail(strsplit(directory, split = "/", fixed = TRUE)[[1]], 1), "/", 
                                   names_dirs[iter_dir], "/", 
                                   list.files(directory.plot[iter_dir])[index_plot], 
                                   sep = "")
      ls.newplot[[iter_dir]] <- seq(0, n.index_plot - 1) %% plotPerPage == 0
      ls.endplot[[iter_dir]] <- rep(FALSE, n.index_plot)
      ls.endplot[[iter_dir]][intersect(which(ls.newplot[[iter_dir]]) - 1, 1:n.index_plot)] <- TRUE
      ls.endplot[[iter_dir]][n.index_plot] <- TRUE
      ls.slices[[iter_dir]] <- rep("", n.index_plot)        
      next
    }
    
    # extraction des coupes
    slice_tempo <- lapply(1:n.index_plot, 
                          function(x){strsplit(split_tempo[[index_plot[x]]][position_slices[x]], 
                                               split = "slice", fixed = TRUE)[[1]][[2]]})
    
    # verifie qu il s agit du meme plot
    if(length(index_plot) > 1){
      test.plot <- sapply(1:n.index_plot, function(x){
        x_tempo <- gsub(pattern = "-legend", replacement = "", fixed = TRUE, list.files(directory.plot[iter_dir])[index_plot[x]])
        gsub(pattern = slice_tempo[[x]], replacement = "", fixed = TRUE, x_tempo)
      })
      
      if(any(test.plot != test.plot[1])){
        stop("constLatex : different types of plot in the same directory \n", 
             "directory : ", directory.plot[iter_dir], " \n", 
             "files : \"", list.files(directory.plot[iter_dir])[index_plot[1]], "\" \n", 
             "        \"", paste(list.files(directory.plot[iter_dir])[index_plot[which(test.plot != test.plot[1])]], collapse = "\"\n        \""), "\"\n")
      }
    }      
    
    # mettre en ordre les figures suivant les coupes (avec les figures en premier et la legend en dernier)
    order_tempo <- order(as.numeric(unlist(lapply(slice_tempo, 
                                                  function(x){strsplit(x, split = "-", fixed = TRUE)[[1]][1]}))) + 10^6 * legend_tempo)
    index_plot <- index_plot[order_tempo]
    
    
    # store results
    ls.legend[[iter_dir]] <- c(ls.legend[[iter_dir]][legend_tempo == 0], ls.legend[[iter_dir]][legend_tempo > 0])          
    ls.plot[[iter_dir]] <- paste(utils::tail(strsplit(directory, split = "/", fixed = TRUE)[[1]], 1), "/", 
                                 names_dirs[iter_dir], "/", 
                                 list.files(directory.plot[iter_dir])[index_plot], 
                                 sep = "")
    ls.newplot[[iter_dir]] <- seq(0, n.index_plot - 1) %% plotPerPage == 0
    ls.endplot[[iter_dir]] <- rep(FALSE, n.index_plot)
    ls.endplot[[iter_dir]][intersect(which(ls.newplot[[iter_dir]]) - 1, 1:n.index_plot)] <- TRUE
    ls.endplot[[iter_dir]][n.index_plot] <- TRUE
    
    slice_tempo.first <- unlist(lapply(slice_tempo[order_tempo], function(x){strsplit(x, split = "-", fixed = TRUE)[[1]][[1]]}))
    slice_tempo.last <- unlist(lapply(slice_tempo[order_tempo], function(x){strsplit(x, split = "-", fixed = TRUE)[[1]][[2]]}))
    
    ls.slices[[iter_dir]] <- c(paste("( slices ", slice_tempo.first[ls.newplot[[iter_dir]] == TRUE], "-", 
                                     slice_tempo.last[ls.endplot[[iter_dir]] == TRUE], " ", 
                                     if(any(legend_tempo)){"with legend "}, ")", sep = "")
    )
  }
  
  names(ls.legend) <- names_dirs
  names(ls.newplot) <- names_dirs    
  names(ls.endplot) <- names_dirs    
  names(ls.plot) <- names_dirs
  names(ls.slices) <- names_dirs
  
  #### data clinique
  if(!is.null(tabular)){
    ## identification du patient
    ls.tabular <- lapply(tabular, function(x){
      table_tempo <- x[x[,1] == identifier,-1, drop = FALSE];
      names(table_tempo) <- gsub(pattern = "_", replacement = "\\_", x = names(table_tempo))
      names(table_tempo) <- gsub(pattern = "%", replacement = "\\%", x = names(table_tempo))
      return(table_tempo)
    }
    )
    
  }else{
    ls.tabular <- NULL
  }
  
  return(list(ls.legend = ls.legend, 
              ls.newplot = ls.newplot, 
              ls.endplot = ls.endplot, 
              ls.plot = ls.plot, 
              ls.slices = ls.slices, 
              ls.tabular = ls.tabular))
}
