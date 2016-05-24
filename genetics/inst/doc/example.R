# narrow window, few few digits, quiet
options(width=50)
options(digits=2)
options(verbose=FALSE)


library(genetics)

# Load the data from a CSV file
data <- read.csv("example_data.csv")

# Convert genotype columns to genotype variables
data <- makeGenotypes(data)

## Annotate the genes 
marker(data$a1691g) <-
  marker(name="A1691G",
         type="SNP",
         locus.name="MBP2",
         chromosome=9, 
         arm="q", 
         index.start=35,
         bp.start=1691,
         relative.to="intron 1")


marker(data$c104t) <-
  marker(name="C-104T",
         type="SNP",
         locus.name="MBP2",
         chromosome=9, 
         arm="q", 
         index.start=35,
         bp.start=-104,
         relative.to="intron 1")

marker(data$c2249t) <-
  marker(name="C2249T",
         type="SNP",
         locus.name="MBP2",
         chromosome=9, 
         arm="q", 
         index.start=35,
         bp.start=2249,
         relative.to="intron 1")

# Look at some of the data
data[1:5,]

# Get allele information for c104t
summary(data$c104t)


# Check Hardy-Weinberg Equilibrium
HWE.test(data$c104t)

# Check Linkage Disequilibrium
ld <- LD(data)
ld # text display

pdf(file="LD.pdf")
LDtable(ld) # graphics display
dev.off()

summary(lm( DELTA.BMI ~ homozygote(c104t,'C') +
                        allele.count(a1691g, 'G') +
                        c2249t, data=data))

