"boa.chain.reset" <-
function()
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   boa.chain(work = boa.chain("master"),
              work.support = boa.chain("master.support"), work.sync = TRUE)
   invisible()
}
