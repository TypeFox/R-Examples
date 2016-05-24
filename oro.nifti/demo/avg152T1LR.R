mniLR <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "mniLR"))
image(mniLR)
orthographic(mniLR)
