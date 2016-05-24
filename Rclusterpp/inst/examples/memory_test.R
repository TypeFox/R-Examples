#!/usr/bin/env Rscript

library(Rclusterpp)
suppressMessages(require(inline))

allocate_matrix <- function(size) {
	fx <- cxxfunction( signature(n = "integer"), '
		typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrix;

		int size = as<int>(n);
		
		matrix m(size, size);
		m.setZero();
		
		return wrap( m.rows() );
	',
	plugin = "Rclusterpp")
	return(fx(size))
}

sizes = c(10, 100, 1000, 10000, 50000)
for (s in sizes) {
	message("Creating ",s," x ",s," matrix...",appendLF=FALSE)
	allocate_matrix(as.integer(s))
	message(" successfully")
}
