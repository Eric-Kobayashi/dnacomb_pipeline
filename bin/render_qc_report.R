#!/usr/bin/env Rscript
# Produce QC results for DNAComb pipeline by rendering an Rmd document
library(magrittr)
library(argparser)
library(rmarkdown)

p <- arg_parser("Render DNAComb QC report", name = "render_qc_report.R") %>%
  add_argument("--roots", "Root paths to DNAComb output", nargs = Inf) %>%
  add_argument("--libraries", "Library TSV file(s)", nargs = Inf) %>%
  add_argument("--counts", "Pipeline record-count TSV") %>%
  add_argument("--rmd", "Rmd template", default = "dnacomb_qc_report.Rmd") %>%
  add_argument("--output", "Output base path", default = "dnacomb_qc_report")

args <- parse_args(p)

render(
  input = args$rmd,
  output_file = basename(args$output),
  output_dir = dirname(args$output),
  output_format = "all",
  clean = TRUE,
  params = list(
    root_paths = args$roots,
    library_files = args$libraries,
    record_counts_file = args$counts
  ),
  envir = new.env(parent = globalenv())
)
