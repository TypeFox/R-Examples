DAISIE_tutorial = function()
{
   filename = system.file("DAISIE_tutorial.pdf",package = "DAISIE")
   os = .Platform$OS.type
   if(os == "windows")
   {
       shell.exec(filename)
   }
   if(os == "unix")
   {
       system(paste("open",filename,sep = " "))
   }
}