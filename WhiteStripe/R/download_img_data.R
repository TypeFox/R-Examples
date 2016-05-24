#' @title Download T1 and T2 data
#' @description Download T1 and T2 data for Examples
#' @return Logical indicator if the files were downloaded.
#' @export
download_img_data = function(){
  stubs = c("T1Strip.nii.gz",
            "T2Strip.nii.gz")
  img_files = system.file(stubs,
                          package="WhiteStripe")
  
  if (!all(file.exists(img_files))){
    for (istub in stubs){
      url = paste0("http://muschellij2.github.io/WhiteStripe/", istub)
      urlfile <- file.path(system.file(package="WhiteStripe"), 
                         istub)
      download.file(url, urlfile, quiet=TRUE)
    }
  }
  img_files = system.file(stubs,
                          package="WhiteStripe")
  return(all(file.exists(img_files)))
}