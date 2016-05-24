mniRL <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                             "mniRL"))
image(mniRL)
orthographic(mniRL)
