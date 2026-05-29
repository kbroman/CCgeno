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

**Note**: There were a bunch of cases where
tha haplotype files were internally inconsistent, with sets of
overlapping haplotypes. If the end of one segment overlapped the start
of the next, we shifted it back so that they didn't overlap.


### Usage

The genotype probabilities as well as the map of locations is in
[`cc_genoprob.RData`](cc_genoprob.RData).





### License

Code licensed under the [MIT license](License.md).
