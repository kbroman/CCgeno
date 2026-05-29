## Collaborative Cross genotype information

Genotype information on the [Collaborative Cross]() mouse strains
is available in several forms, but not immediately useful for analysis
with [R/qtl2](https://kbroman.org/qtl2).

GigaMUGA data on many strains are provided in standard R/qtl2 input at
the [qtl2data repository on github](), but not all strains are
represented.

Here, we converted the `.hap` files at
<https://csbio.unc.edu/CCstatus/CCGenomes/> to a genotype probability
array suitable for use with R/qtl2.
