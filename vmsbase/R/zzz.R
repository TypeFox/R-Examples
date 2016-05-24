.onAttach <- function(lib, pkg)
{
  packageStartupMessage("

    ------------------------------------------------------------------

                          -- VMSbase 2.0 --

               VMSbase comes with ABSOLUTELY NO WARRANTY

        This is free software, and you are welcome to redistribute    
             it under certain conditions, as outlined at

               http://www.gnu.org/licenses/gpl-2.0.html

        ----------------------------------------------------------    

     Type in the console ' gui_main() ' to start the main VMSbase GUI 

    ------------------------------------------------------------------

")
}
# 
# 
# .onUnload <- function(libpath)
# {
#     library.dynam.unload("vmsbase", libpath)
# }
# 
# ## no S4 methodology here; speedup :
# .noGenerics <- TRUE
