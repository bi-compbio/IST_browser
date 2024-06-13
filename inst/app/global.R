# capture.output(sessionInfo(), file=stderr())
# capture.output(.libPaths(), file=stderr())

library(IST)
library(data.table)
library(tibble)

# IST results to display
#
# first check if ist.result.path was specified
ist.results.path <- getShinyOption("ist.result.path")
if (!is.null(ist.results.path)) {
  message("Reading IST results file in ", ist.results.path)
  ist.result <- readRDS(getShinyOption("ist.result.path"))
  # if not, just use the example IST data
} else if (exists(data("sample.results.ist", package = "IST"))) {
  message("Reading sample IST results as bundled in the IST package")
  data("sample.results.ist", package = "IST")
  ist.result <- sample.results.ist
  # if this example data is not available, the IST install is corrupt or missing
} else {
  stop("No IST results object to read: please install the IST package")
}

# Precalculate Gene/Pathway/Signature metadata
genesByPathway <- as.data.table(ist.result@ist.pathways@pathways.table.bin)
genesByPathway[, c("label") := .(IST::vec.ensembl2symbol[gene.id])]
genesVsSignatures <- getTabSignatures(ist.result) %>%
  copy() %>%
  setnames(old = "ortholog", new = "gene.id")
genesVsSignatures[, c("label") := .(IST::vec.ensembl2symbol[gene.id])]
genesVsSignatures <- genesVsSignatures[, c("label", "sig.id", "logFC")]
genesVsSignaturesAndPathways <- genesVsSignatures[genesByPathway, on = "label", nomatch = 0, allow.cartesian = TRUE]
genesVsSignaturesAndPathways <- genesVsSignaturesAndPathways[, c("label", "path.id", "sig.id", "logFC")]

defaultPathways <- function() {
  getMetaPathways(ist.result)[path.class == "default"]$path.id
}

defaultSignatures <- function() {
  getMetaSig(ist.result)[sig.class %in% c("default", "positivecontrol")]$sig.id
}

defaultSignatureMetadata <- function() {
  n <- names(getMetaSig(ist.result))
  n[n %in% c(
    "sig.id",
    "contrast.name",
    "contrast.type",
    "treatment",
    "organism.name",
    "animal.model",
    "model.timepoint",
    "model.meta",
    "tissue",
    "disease",
    "study.id",
    "n.signif"
  )]
}
