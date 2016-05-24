hasLaTeX <- function(pdfLatexPath) {
  return(file.exists(paste0(pdfLatexPath, "/pdflatex.exe")) |
         file.exists(paste0(pdfLatexPath, "/pdflatex")));
}

### Function to escape special latex characters, based on
### http://stackoverflow.com/questions/5406071/r-sweave-latex-escape-variables-to-be-printed-in-latex
sanitizeLatexString <- function(str) {
  return(gsub('([#$%&~_\\^\\\\{}])', '\\\\\\1', str, perl = TRUE));
}

rnwString.initiate <- function(studyName, authorName,
                               docClassArgs = 'a4paper,portrait,11pt',
                               newPage = TRUE, pageMargins=15) {
  res <- paste0("\\documentclass[", docClassArgs, "]{article}
                    
% For adjusting margins
\\usepackage[margin=", pageMargins, "mm]{geometry}
% For printing on rotated page
\\usepackage{pdflscape}
% For resizing correlationtables that become too large
\\usepackage{adjustbox}

% !Rnw weave = knitr
                    
\\title{", studyName, "}
\\author{", authorName, "}
\\begin{document}
\\raggedright
\\noindent
\\maketitle
GENERATED ON ", date(), "
\\tableofcontents\n");
  if (newPage) {
    res <- paste0(res, "\\newpage\n");
  }
  return(res);
}

rnwString.terminate <- function(rnwString) {
  return(paste0(rnwString, "\n\\end{document}"));
}

rnwString.generate <- function(rnwString, rnwPath, fileName, pdfLatexPath,
                               envir = parent.frame()) {
    
  ### Write the knitr text to a file
  fileConn<-file(paste0(rnwPath, '/', fileName, '.rnw'));
  writeLines(rnwString, fileConn);
  close(fileConn);
  
  ### Knit the .tex file
  knit(paste0(rnwPath, '/', fileName, '.rnw'),
       paste0(rnwPath, '/', fileName, '.tex'), envir=envir);

  if (!hasLaTeX(pdfLatexPath)) {
    stop('In path "', pdfLatexPath, '", the file pdflatex.exe (Windows) or ',
         'pdflatex (MacOS or Ubuntu (Linux)) does not exist! Please ',
         'locate the file and provide its path (without the last ',
         'slash) See ?rnwString for more information.');
  }
  
  ### First run to generate TOC
  tryCatch(
    texOutput <- system(paste0('"', pdfLatexPath, '/pdflatex" "',
                               paste0(rnwPath, '/', fileName, '.tex'), '" ',
                               '-output-directory "', rnwPath, '"'),
                        intern=TRUE)
    , error = function(e) {
      cat(paste("Error returned by pdflatex: ", e));
    }
  );
  
  ### Then generate the final PDF
  tryCatch(
    texOutput <- system(paste0('"', pdfLatexPath, '/pdflatex" "',
                               paste0(rnwPath, '/', fileName, '.tex'), '" ',
                               '-output-directory "', rnwPath, '"'),
                        intern=TRUE)
    , error = function(e) {
      cat(paste("Error returned by pdflatex: ", e));
    }
  );  
}
