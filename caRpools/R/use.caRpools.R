use.caRpools = function(type=NULL, file="CaRpools-extended-PDF.Rmd", miaccs = "MIACCS.xls", check=TRUE, work.dir=NULL)
{
  # Main function which is used to start the RMD conversion to either HTML or PDF output. 
  # start Markdown

  
  if(is.null(work.dir))
  {
    setwd(getwd())
  }
  else
  {
    setwd(work.dir)
  }
  
  # first check if everything is wokring?
  if(identical(check,TRUE))
  {
    test.all = check.caRpools(template = file, miaccs=miaccs)
    if(identical(test.all,FALSE))
    {
      stop("CaRpools package verification failed.")
    }
  }
  
  
  
  if(is.null(type))
  {
    miaccs.file = miaccs
    rmarkdown::render(input = file)
  }
  else if(type=="PDF")
  {
    miaccs.file = miaccs
    rmarkdown::render(input = file, output_format = rmarkdown::pdf_document(latex_engine="xelatex",
                                                                 fig_height = 6,
                                                                 fig_width= 11,
                                                                 highlight = NULL,
                                                                 keep_tex = TRUE,
                                                                 toc = TRUE,
                                                                 toc_depth = 3)
                      )
  }
  else if(type=="HTML")
  {
    miaccs.file = miaccs
    rmarkdown::render(input = file, output_format = rmarkdown::html_document(
                                                                 fig_height = 8,
                                                                 fig_width= 11,
                                                                 keep_md = TRUE,
                                                                 toc = TRUE,
                                                                 toc_depth = 3)
                      )
  }

  
}