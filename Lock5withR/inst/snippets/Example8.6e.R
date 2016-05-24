SST <- sum( ~((Ants - grandMean)^2), data = SA); SST
SSM <- sum( ~M2, data = SA ); SSM    # also called SSG
SSE <- sum( ~E2, data = SA ); SSE

