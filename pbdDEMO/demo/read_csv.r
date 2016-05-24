library(pbdDEMO, quietly = TRUE)
init.grid()


if(comm.size() != 4){
  comm.stop("This example requries 4 processors.")
}


# Read in example csv file
x <- system.file("extra/data/x.csv", package = "pbdDEMO")

### Same as calling:
# blacs_gridinit(ICTXT=3, NPROW=2, NPCOL=1L)
# dx <- ddmatrix(1:100, 10, 10, bldim=c(6, 10), ICTXT=3)
dx <- read.csv.ddmatrix(
  file=x, sep=",", nrow=10, ncol=10, header=TRUE, bldim=4, num.rdrs=2, ICTXT=0
)

dx

# Recombine on process 0 and print to show that everything worked.
x <- as.matrix(dx, proc.dest=0)

comm.print(x)


finalize()
