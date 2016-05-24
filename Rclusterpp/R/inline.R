inlineCxxPlugin <-
    Rcpp:::Rcpp.plugin.maker(
			include.before = "#include <Rclusterpp.h>", 
      libs           = Rclusterpp::RclusterppLdFlags(FALSE),
			package        = "Rclusterpp",
			LinkingTo      = c("Rclusterpp", "RcppEigen", "Rcpp"),
			Makevars       = NULL, 
			Makevars.win   = NULL
    )
