permPATH2HTML = function(dat, dir, fname, title=NULL, bgcolor="#BBBBEE"){
  target = HTMLInitFile(dir, filename=fname, Title=title, BackGroundColor=bgcolor)
  HTML(xtable(dat), file=target)
  HTMLEndFile()
}
        
