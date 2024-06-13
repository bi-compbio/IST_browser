#!/usr/bin/env Rscript
library(styler)

style_file(
  c(
    "R/run.R",
    "inst/app/ui.R",
    "inst/app/server.R",
    "inst/app/global.R"
  )
)
