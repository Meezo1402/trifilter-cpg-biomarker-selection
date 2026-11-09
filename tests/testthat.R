#!/usr/bin/env Rscript
# Run the unit tests:  Rscript tests/testthat.R
library(testthat)
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]), ".."))
if (is.na(root) || !file.exists(file.path(root, "R", "utils.R"))) root <- getwd()
for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
