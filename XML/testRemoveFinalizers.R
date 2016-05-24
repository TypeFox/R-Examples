x = try(dyn.load(sprintf("testRemoveFinalizers%s", .Platform$dynlib.ext)))
q(status = if(inherits(x, "try-error")) 1L else 0L)



