library(here)

source("hap2prob.R")
g <- hap2geno(haptab, upos)
saveRDS(g, here("Hapfiles", "called_geno.rds"))
