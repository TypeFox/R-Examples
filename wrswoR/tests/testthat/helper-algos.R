funcnames <- ls("package:wrswoR", pattern = "^sample_int_")

funcs <- lapply(setNames(nm = funcnames), get, pos = "package:wrswoR")
