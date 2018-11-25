# Install the package ####
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("transcriptogramer", version = BiocManager::version())

# Load the package ####
library(transcriptogramer)

# Load the required data ####
load("data.RData")

# Run the analysis ####
t <- transcriptogramPreprocess(association = edges, ordering = Hs800,
                               radius = 125)
groups <- c(rep("control", 3), rep("aza5", 3), rep("aza10", 3))
idx <- groups %in% (c("control", "aza10"))
t <- transcriptogramStep1(object = t, expression = exp[, idx],
                          dictionary = GPL570)
t <- transcriptogramStep2(object = t)
levels <- groups[idx] %in% "control"
t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.02,
                             title = "aza10 x control Microarray",
                             species = "Homo sapiens",
                             boundaryConditions = FALSE)
rdp <- clusterVisualization(t, onlyGenesInDE = TRUE)
t <- clusterEnrichment(t, species = "Homo sapiens", nCores = TRUE,
                           algorithm = "parentchild", pValue = 0.05,
                           onlyGenesInDE = TRUE)
terms <- Terms(t)
top5terms <- lapply(1:4, function(i){
  head(terms[terms$ClusterNumber == i, c(1, 2, 6)], 5)
})
top5terms
