library(kitagawa)
#
# Set a few parameters
Rc <- 0.0508   # m, radius of water-sensing (2in)
Lc <- 146.9    # m, length of grouted region (482ft)
Rs <- 3*Rc     # m, radius of screened region (6in)
Ls <- 9.14     # m, length of screened region (30ft)
#
# calculate the sensing volume for the well parameters
Volw <- sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
#
# Set some frequencies
Frqs <- 10**seq.int(from=-4,to=0,by=0.1) # log10-space
#
# Calculate the response
Rsp <- well_response(omega=Frqs, T.=1e-6, S.=1e-5, Vw.=Volw, Rs.=Rs, Ku.=40e9, B.=0.2, freq.units="Hz")
#
# Plot it
kitplot(Rsp)
