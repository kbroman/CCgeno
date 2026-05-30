# download all haplotype files

library(here)

urls <- readLines(here("Hapfiles/hapfiles_seq.txt"))

for(file in urls) {
    download.file(file, file.path(here("Hapfiles"), basename(file)))
}




urls <- readLines(here("Hapfiles/hapfiles_mrca.txt"))

dir <- here("Hapfiles", "MRCA")
if(!dir.exists(dir)) dir.create(dir)

for(file in urls) {
    download.file(file, file.path(dir, basename(file)))
}


urls <- readLines(here("Hapfiles/hapfiles_gm.txt"))

dir <- here("Hapfiles", "GM")
if(!dir.exists(dir)) dir.create(dir)

for(file in urls) {
    download.file(file, file.path(dir, basename(file)))
}

# the MRCA hap files don't contain M and Y, but don't have out-of-order problems
# the GigaMUGA hap files don't contain M and Y, and only the last 6 have chr X
#      and the GM files have a ton of out-of-order / overlap problems
