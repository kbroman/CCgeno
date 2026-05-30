# read in hapfiles.txt with URLs for haplotype files
# split into three sets: seq, mrca, gigamuga

library(here)

hapfiles <- readLines(here("Hapfiles/hapfiles_orig.txt"))
mrca <- grep("MRCA", hapfiles)
gm <- grep("Gigamuga", hapfiles)
muga <- grep("MUGA", hapfiles)

sq <- seq_along(hapfiles) %wnin% c(mrca, gm, muga)

strains_sq <- sapply(strsplit(hapfiles[sq], "[\\/\\-]"), "[", 7)
table(table(strains_sq))

strains_mrca <- sapply(strsplit(hapfiles[mrca], "[\\/\\_]"), "[", 8)
table(table(strains_mrca))

strains_gm <- sapply(strsplit(hapfiles[gm], "[\\/\\-]"), "[", 8)
table(table(strains_gm))

# save the sequenced hapfiles
cat(hapfiles[sq], file=here("Hapfiles/hapfiles_seq.txt"), sep="\n")
cat(hapfiles[mrca], file=here("Hapfiles/hapfiles_mrca.txt"), sep="\n")
cat(hapfiles[gm], file=here("Hapfiles/hapfiles_gm.txt"), sep="\n")
