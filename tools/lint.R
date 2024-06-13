#!/usr/bin/env Rscript

Sys.setenv(CI = "") # https://github.com/jimhester/lintr/issues/166
library(lintr)

# Disable camel case linter, most of the Shiny functions are camelCase...
linters <- with_defaults(camel_case_linter = NULL, line_length_linter(120))

lints <- 0
for (f in list.files(pattern = "\\.R(md)?$", ignore.case = TRUE, recursive = TRUE)) {
    l <- lint(f, linters = linters)
    lints <- lints + length(l)
    print(l)
}

# TODO disable until all errors are fixed
# # Non-zero exist status for CI
# quit(status = lints)
