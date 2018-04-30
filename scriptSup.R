# Install package
source("https://bioconductor.org/biocLite.R")
biocLite("transcriptogramer")

# Load package
library(transcriptogramer)

# Load data
load("data.RData")

# Run analysis
t <- transcriptogramPreprocess(association = edges, ordering = Hs800,
                               radius = 125)
t <- transcriptogramStep1(object = t, expression = exp[, c(1:3, 7:9)],
                          dictionary = GPL570)
t <- transcriptogramStep2(object = t)
levels <- c(rep(TRUE, 3), rep(FALSE, 3))
t <- differentiallyExpressed(object = t, levels = levels, pValue = 0.02,
                             title = "aza10 x aza0 Microarray",
                             species = "Homo sapiens")
rdp <- clusterVisualization(t)
terms <- clusterEnrichment(t, species = "Homo sapiens", nCores = T,
                           algorithm = "parentchild", pValue = 0.05)
terms <- terms[order(terms$pValue),]
top5terms <- lapply(1:4, function(i){
  head(terms[terms$ClusterNumber == i, c(1, 2, 6)], 5)
})
top5terms