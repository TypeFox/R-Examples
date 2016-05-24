
.onAttach <- function(lib, pkg)  {

  packageStartupMessage(sprintf("Loaded TSdist v%s. See ?TSdist for help, citation(\"TSdist\") for use in publication.\n",
            utils::packageDescription("TSdist")$Version ) );
      
  ## Register the distance measures into package proxy
  pr_DB$set_entry(FUN=TSDistances, 
                  names="TSDistances", 
                  loop=TRUE, type="metric", distance=TRUE,
               description="A function that computes different 
                 distance measures for time series.")

}
