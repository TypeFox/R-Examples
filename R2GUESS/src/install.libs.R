guess.binary <- if (WINDOWS) "GUESS.exe" else "GUESS"
if (!file.exists(guess.binary)) {
    message(sprintf("%s not found", guess.binary))
} else {
    bin.arch <- if (nzchar(R_ARCH)) paste0("bin", R_ARCH) else "bin"
    dest <- file.path(R_PACKAGE_DIR, bin.arch)
    message(sprintf("Installing %s to %s", guess.binary, dest))
    dir.create(dest, showWarnings=FALSE, recursive=TRUE)
    file.copy(guess.binary, dest, overwrite=TRUE)
}

