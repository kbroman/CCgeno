# download all haplotype files

library(here)

urls <- readLines(here("Hapfiles/hapfiles_seq.txt"))

for(file in urls) {
    download.file(file, file.path(here("Hapfiles"), basename(file)))
}
