### SHELL> mpiexec -np 4 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))

### Check.
if(comm.size() != 4){
  comm.stop("4 processors are requried.")
}
fn <- "iris.txt"

### Manually distributed iris.
da <- iris[get.jid(nrow(iris)),]

### Dump data.
comm.write.table(da, file = fn, quote = FALSE, sep = "\t",
                 row.names = FALSE)
# comm.write.csv(da, file = fn)

### Read back in.
# .pbd_env$SPMD.IO$max.read.size <- 50
# .pbd_env$SPMD.IO$balance.method <- "block.cyclic"
da.gbd <- comm.read.table(fn, header = TRUE, sep = "\t", quote = "")
# da.gbd <- comm.read.csv(file = fn)
comm.print(c(nrow(da), nrow(da.gbd)), all.rank = TRUE)

### Read in common.
# .pbd_env$SPMD.IO$max.read.size <- 50
da.common <- comm.read.table(fn, header = TRUE, sep = "\t",
                             quote = "", read.method = "common")
# da.common <- comm.read.csv(fn, read.method = "common")
comm.print(c(nrow(da.common), sum(da.common != iris)))

### Finish.
finalize()
