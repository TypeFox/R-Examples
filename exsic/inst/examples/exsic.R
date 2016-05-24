# Example
load(system.file("data/config.rda", package="exsic"))
###########################################################
# This runs the example file


# Read input file
df = system.file("samples/exsic.csv", package="exsic")
# read only first 10 records
data = read.exsic(df)[1:10,]

# Prepare output file
td = tempdir()
of = file.path(td,"out.html")


# Example 1: mostly default parameters
# Prepare exsiccatae indices
exsic(data, html = of) 


# Example 2: using another format
of = file.path(td,"out_PK.html")
exsic(data, html = of, format = format.PK) 

