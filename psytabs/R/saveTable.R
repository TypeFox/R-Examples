saveTable <-
function (x, file, HTML=FALSE, post.editing = TRUE) {
  if(HTML) {
    target <- R2HTML::HTMLInitFile(getwd(),filename=file)
    R2HTML::HTML("<br>Table X.<br>",file=target)
    R2HTML::HTML(x, file=target)
    R2HTML::HTMLEndFile()
  } else {
    if(length(grep(".rtf", file)) == 0 & length(grep(".doc", file)) == 0) {
      file <- paste(file, ".rtf", sep="")
    } else if (length(grep(".docx", file)) > 0) {
      file <- paste(file, ".rtf", sep="")
    }
    output <- file
    rtf <- rtf::RTF(output, font.size=12)
    rtf::addParagraph(rtf,"Table X.")
    rtf::addTable(rtf, x, row.names=TRUE, NA.string="")
    rtf::done(rtf)
  }
  
  if(post.editing == TRUE) {
    
    if(rownames(x)[1] == c("Configural")) {
      #Post-editing
      rtf.code <- readLines(file)
      
      rtf.code <- gsub("paperw12240", "paperw16840", rtf.code)
      rtf.code <- gsub("paperh15840", "paperh11907", rtf.code)
      rtf.code <- gsub("clwWidth424", "clwWidth600", rtf.code)
      rtf.code <- gsub("clwWidth761", "clwWidth1000", rtf.code)
      rtf.code <- gsub("clwWidth536", "clwWidth700", rtf.code)
      rtf.code <- gsub("clwWidth874", "clwWidth1100", rtf.code)
      
      writeLines(rtf.code, file)
      
    }
    
    if(length(class(x)) > 1 && class(x)[2] == "efaTable.m.sd") {
      if(class(x)[3] == "2factors") {
        rtf.code <- readLines(file)
        rtf.code <- rtf.code[-(14:25)]
        rtf.code[13] <- 

"\\trowd\\trgaph100\\trleft0 \\trhdr\\clvertalc\\clshdrawnil\\clwWidth1104\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx1
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx2
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx3
\\clmgf\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx4
\\clmrg\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx5
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx6
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b  \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\qc\\u0955?\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\widctlpar\\intbl\\row

\\trowd\\trgaph100\\trleft0 \\clvertalc\\clshdrawnil\\clwWidth1104\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrb\\brdrs\\brdrw15\\cellx1
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx2
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrb\\brdrs\\brdrw15\\cellx3
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx4
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx5
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrb\\brdrs\\brdrw15\\cellx6
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b  \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b M\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b SD\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b ML1\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b ML2\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b h2\\b0\\cell
\\widctlpar\\intbl\\row"

        writeLines(rtf.code, file)
      }
      
    }
    
    if(length(class(x)) > 1 && class(x)[2] == "efaTable") {
      if(class(x)[3] == "2factors") {
        rtf.code <- readLines(file)
        rtf.code <- rtf.code[-(14:21)]
        rtf.code[13] <- 
          
"\\trowd\\trgaph100\\trleft0 \\trhdr\\clvertalc\\clshdrawnil\\clwWidth1104\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx1
\\clmgf\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx4
\\clmrg\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx5
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\cellx6
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b  \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\qc\\u0955?\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b \\b0\\cell
\\widctlpar\\intbl\\row

\\trowd\\trgaph100\\trleft0 \\clvertalc\\clshdrawnil\\clwWidth1104\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrb\\brdrs\\brdrw15\\cellx1
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx4
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrt\\brdrs\\brdrw15\\clbrdrb\\brdrs\\brdrw15\\cellx5
\\clvertalc\\clshdrawnil\\clwWidth941\\clftsWidth3\\clheight260\\clpadl100\\clpadr100\\gaph\\clbrdrb\\brdrs\\brdrw15\\cellx6
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b  \\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b ML1\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b ML2\\b0\\cell
\\pard\\ql\\widctlpar\\intbl\\fi0\\f2\\fs24\\b h2\\b0\\cell
\\widctlpar\\intbl\\row"

        writeLines(rtf.code, file)
      }
      
    }
  
    #For every file
    rtf.code <- readLines(file)
    #Remove bold formatting of column headers
    rtf.code[12:(12+(ncol(x)+1)*2)] <- gsub("\\b ", " ", rtf.code[12:(12+(ncol(x)+1)*2)], fixed=TRUE)
    rtf.code[12:(12+(ncol(x)+1)*2)] <- gsub("\\b0\\", "\\", rtf.code[12:(12+(ncol(x)+1)*2)], fixed=TRUE)

    writeLines(rtf.code, file)
  }
  
}
