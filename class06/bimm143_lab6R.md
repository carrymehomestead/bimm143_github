# bimm143_lab6R
Malibu Slattery (A18488012)

- [Background](#background)
- [Our first function](#our-first-function)
- [A second function](#a-second-function)
- [A Protein generator function](#a-protein-generator-function)

## Background

All functions in R have at least 3 things:

- A **name** that we use to call the function
- One or more input **arguments**.
- The **body** the lines of R code that do work.

## Our first function

Let’s write a silly wee function to add numbers called `add()` to add
some numbers (the input arguments).

``` r
add <- function(x,y) {
  x + y
}
```

Now we can use this function

``` r
add(100, 1)
```

    [1] 101

``` r
add(x=c(100,1,100), y=1)
```

    [1] 101   2 101

> Q. What if I give a multiple element vector to `x` and `y`?

``` r
add(x=c(100,1), y= c(100,1))
```

    [1] 200   2

``` r
#add(x=c(100,1), y=1, z=1)
```

> Q. What happens when I only give one input to the add function?

``` r
addnew <- function(x,y=1){
  x + y
}
```

``` r
addnew(x=100)
```

    [1] 101

``` r
addnew(c(100,1), 100)
```

    [1] 200 101

If we write our function with input arguments having no default value,
then the user will be required to set them up when they use the
functipm. We can give our input arguments “default” values by setting
them equal to some sensible value. e.g. y= 1 in the `addnew()` function.

## A second function

Let’s try something more interesting…make a sequence generating tool.

The `sample()` function can be a useful starting point here.

``` r
sample(1:10, size = 9)
```

    [1]  3  8  6  4  5  1 10  2  9

``` r
sample(1:10, size =12, replace = TRUE)
```

     [1]  7  1  1  6  9  9  5  4 10  7  4  1

> Q. Write code for the `sample()` function that generates nucleotide
> sequences of length 6.

``` r
sample(x=c("A","T","G","C"), 6, replace = TRUE)
```

    [1] "T" "C" "G" "T" "G" "G"

> Q. Write a first function `generate_dna()` that returns a *user
> specified* length DNA sequence.

``` r
generate_dna <- function(len){
          sample(x=c("A","T","G","C"), size = len, replace = TRUE)
}
```

``` r
generate_dna(12)
```

     [1] "C" "G" "C" "T" "G" "C" "A" "C" "C" "C" "T" "G"

``` r
generate_dna(len=100)
```

      [1] "A" "A" "T" "T" "T" "A" "T" "T" "C" "G" "G" "T" "T" "T" "C" "A" "T" "C"
     [19] "G" "G" "G" "G" "C" "G" "C" "T" "G" "C" "G" "T" "T" "A" "G" "T" "G" "G"
     [37] "A" "A" "A" "C" "A" "G" "T" "G" "G" "T" "A" "G" "A" "G" "T" "G" "T" "A"
     [55] "C" "C" "C" "A" "T" "T" "C" "C" "G" "A" "G" "G" "T" "G" "T" "T" "A" "G"
     [73] "A" "C" "G" "T" "C" "G" "G" "T" "C" "C" "C" "C" "T" "C" "A" "T" "A" "A"
     [91] "A" "T" "T" "A" "G" "A" "T" "G" "T" "T"

> **Key Points** Every function in R looks fundamentally the same in
> terms of structure

    name <-function(input) {
     body
    }

> Functions can have multiple inputs. These can be **required**
> arguments or so called **optional** arguments. With optional arguments
> having a default value.

> Q. Modify and improve our function `generate_dna()` to return its
> generated sequence in a more standard form like “AGATTC” rather than
> in a vector “A”, “G”….etc…

``` r
generate_dna <- function(len=6, fasta = FALSE){
         
  ans<- sample(x=c("A","T","G","C"), 
               size = len, replace = TRUE)
  if(fasta) {
  ans<- paste(ans, collapse = "")
}
  return(ans)
}
  
generate_dna()
```

    [1] "T" "T" "A" "C" "A" "C"

The `paste()` function’s job is to join up or stick together (a.k.a.
paste) strings together.

``` r
paste(c("alice", "barry"), "loves R")
```

    [1] "alice loves R" "barry loves R"

``` r
paste("alice", "love R", sep=" DOES NOT ")
```

    [1] "alice DOES NOT love R"

Flow control means where the R brain goes in your code.

``` r
good_mood <- FALSE

if(good_mood) {
  cat("Great!")
} else {
  cat("Bummer!")
}
```

    Bummer!

``` r
generate_dna <- function(len=6, fasta = FALSE){
         
  ans<- sample(x=c("A","T","G","C"), 
               size = len, replace = TRUE)
  if(fasta) {
  cat("Single-element vector output")
  ans<- paste(ans, collapse = "")
  } else {
  cat("Multi-element vector output")
}
  return(ans)
}
  
generate_dna(fasta = TRUE)
```

    Single-element vector output

    [1] "TTCTAC"

## A Protein generator function

> Q. Write a function that generates a user-specified length protein
> sequence.

> Q. Then use this function to generate different sequences between the
> lengths of 6 and 12.

> Q. Are any of your sequences unique? (i.e.: not found anywhere in
> nature.)

``` r
aa<- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
```

``` r
generate_protein <- function(len){
  
  # The amino-acids to sample from:
  aa<- c("A", "R", "N", "D", "C", "Q", "E",
         "G", "H", "I", "L", "K", "M", "F", 
         "P", "S", "T", "W", "Y", "V")
 
   #Draw n=len amino acids to make our sequence
  ans<- sample(aa, size =len, replace = T)
  ans<- paste(ans, collapse="")
  return(ans)
}
```

``` r
myseq<-generate_protein(10)
myseq
```

    [1] "YDRPYASHEP"

``` r
generate_protein(6)
```

    [1] "TCKDCT"

``` r
generate_protein(7)
```

    [1] "HLVHQTW"

``` r
generate_protein(8)
```

    [1] "WRTIHNHY"

``` r
generate_protein(9)
```

    [1] "YMSQNIIED"

``` r
generate_protein(10)
```

    [1] "YSCNILMSFD"

``` r
generate_protein(11)
```

    [1] "YHGEYHFDIFT"

``` r
generate_protein(12)
```

    [1] "QGPAQKEVCKPF"

``` r
for(i in 6:12) {
  #**Fasta formatting ID line**
  cat(">", i, sep="", "\n")
  #Protein sequence line
  cat(generate_protein(i), "\n")
}
```

    >6
    PGFVFQ 
    >7
    SQDKLTR 
    >8
    PGQTDKTN 
    >9
    VIAGTETPS 
    >10
    HKNPSNNAFF 
    >11
    YEEHTHQAMNN 
    >12
    KPYRQSWWELKA 

9 is the “lucky number”, so that it’s safe enough to be unique without
expending too much energy.
