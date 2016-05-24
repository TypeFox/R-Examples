## Example code for calling Doptimize function from R

library(oapackage)

# Set variables

N = 32   # number of runs
k=7      # number of factors
nrestarts=16  # number of restarts 

# Create a D-optimal design
alpha1=1; alpha2=0; alpha3=0      # parameters of the optimization function

# run the function
p = Doptimize(N, k, nrestarts, alpha1, alpha2, alpha3) 

print('resulting design:')
print(p)

dd = Defficiencies(p)
print(sprintf('design: D-efficiency %f, Ds-effciency %f', dd[1], dd[2]))

# Create an design that is suitable under effect hieracry
alpha2=2
p = Doptimize(N, k, nrestarts, alpha1, alpha2, alpha3) 

dd = Defficiencies(p)
print(sprintf('design: D-efficiency %f, Ds-effciency %f', dd[1], dd[2]))


