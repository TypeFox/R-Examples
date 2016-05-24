times <- faithful$eruptions
kdeFaithfulRect <- densityplot(~times,kernel="rectangular",
    main="Rectangular kernel")
kdeFaithfulTri <- densityplot(~times,kernel="triangular",
    main="Triangular kernel")
kdeFaithfulNormal <- densityplot(~times,
    main="Normal kernel")
kdeFaithfulNormal2 <- densityplot(~times,adjust=0.25,
    main="Normal kernel; adjust=0.25")
density(times)       # display some information about the kde
