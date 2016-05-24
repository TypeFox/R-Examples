get_num_cores <-
function()    .Call("get_num_cores", PACKAGE="RhpcBLASctl")
get_num_procs <-
function()    .Call("get_num_procs", PACKAGE="RhpcBLASctl")

blas_get_num_procs <-
function()    .Call("blas_get_num_procs", PACKAGE="RhpcBLASctl")

blas_set_num_threads <-
function(threads) invisible(.Call("blas_set_num_threads",as.integer(threads),PACKAGE="RhpcBLASctl"))

omp_get_num_procs <-
function()    .Call("Rhpc_omp_get_num_procs", PACKAGE="RhpcBLASctl")

omp_get_max_threads <-
function()    .Call("Rhpc_omp_get_max_threads", PACKAGE="RhpcBLASctl")

omp_set_num_threads <-
function(threads) invisible(.Call("Rhpc_omp_set_num_threads",as.integer(threads),PACKAGE="RhpcBLASctl"))

