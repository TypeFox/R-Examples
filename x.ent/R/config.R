options("encoding" = "UTF-8")
xconfig <- function()
{
    opencpu$browse("/library/x.ent/www/config.html")  
}
save_config <- function(data = "")
{
  if(data == ""){
    stop("Data is empty!")
  }
  {
    #message = paste("Save data successful!", R.Version()$version.string)  
    write(data,paste(.libPaths()[1],"/x.ent/www/config/ini.json",sep="")) 
  }
}