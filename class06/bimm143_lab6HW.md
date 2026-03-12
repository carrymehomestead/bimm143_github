# bimm143_class5R
Malibu Slattery (A18488012)

``` r
library(bio3d)
something <- function(s1){
  readprotein<- read.pdb(s1)
  #can use the function return(protein) to see what it gives me
  trimmprotein<-trim.pdb(readprotein, chain="A", elety="CA")
  tempfactor<-trimmprotein$atom$b
  plotgraph<-plotb3(tempfactor, sse=trimmprotein, typ="l", ylab="Bfactor")
}
object<- something("4AKE")
```

      Note: Accessing on-line PDB file

![](bimm143_lab6HW_files/figure-commonmark/unnamed-chunk-1-1.png)

**Documentation**

- What are the inputs to the function?
  - *The inputs to the function include reading the pdb file, trimming
    the file to only include the chain A and electrolyte CA, calling out
    the “b” column, which is the temperature factor (John Hopkins
    University), and finally plotting abundance of the residues, which
    are individual amino acids (Lee, RCSB), against the temperature
    factor.*
- What the function does and how to use it?
  - *The function combines a series of commands into a repeatable strip
    that can have any input and you simply create the string of
    commands, each output within the body the input to the next command
    and treat the object as its own replaceable function.*
- What is the output of the function?
  - *The output is the plot of the temperature factor as it relates to
    the abundance of different residues in a protein.*

## References

*John Hopkins University. (n.d.). 1 1.1. introduction to protein data
bank format. Biostat.
https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf*

*Lee, C. (n.d.).
Https://education.seattlepi.com/residue-mean-biology-6172.html. Seattle
PI. Retrieved January 24, 2026, from
https://education.seattlepi.com/residue-mean-biology-6172.html.*

*U.S. National Science Foundation (DBI-2321666), the US Department of
Energy (DE-SC0019749), the National Cancer Institute, National Institute
of Allergy and Infectious Diseases, and National Institute of General
Medical Sciences of the National Institutes of Health . (n.d.). PDB
statistics: PDB data distribution by residue count. RCSB PDB.
https://www.rcsb.org/stats/distribution-residue-count*

Below is the information I’m referring to and ctrl+p from. Ignore.

``` r
#s1 <- read.pdb("4AKE") # kinase with drug
#s2 <- read.pdb("1AKE") # kinase no drug
#s3 <- read.pdb("1E4Y") # kinase with drug
#s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
#s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
#s1.b <- s1.chainA$atom$b
#s2.b <- s2.chainA$atom$b
#s3.b <- s3.chainA$atom$b
#plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
#plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
#plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```
