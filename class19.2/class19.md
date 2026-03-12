# class19


Useful functions for this lab… bio3d read.fasta conserv()

Useful programs for this lab… blastp hmmer

``` r
library(bio3d)

fast<-read.fasta("A18488012_mutant_seq.fa")
fast
```

                   1        .         .         .         .         .         60 
    wt_healthy     MKAPAVLAPGILVLLFTLVQRSNGECKEALAKSEMNVNMKYQLPNFTAETPIQNVILHEH
    mutant_tumor   MKAPAVLAPGILVLLFTLVQRSNGECKEALAKSEMNVNMKYQLPNFTAETPIQNVILHEH
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     HIFLGATNYIYVLNEEDLQKVAEYKTGPVLEHPDCFPCQDCSSKANLSGGVWKDNINMAL
    mutant_tumor   HIFLGATNYIYVLNEEDLQKVAEYKTGPVLEHPDCFPCQDCSSKANLSGGVWKDNINMAL
                   ************************************************************ 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     VVDTYYDDQLISCGSVNRGTCQRHVFPHNHTADIQSEVHCIFSPQIEEPSQCPDCVVSAL
    mutant_tumor   VVDTYYDDQLISCGSVNRGTCQRHVFPHNHTADIQSEVHCIFSPQIEEPSQCPDCVVSAL
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     GAKVLSSVKDRFINFFVGNTINSSYFPDHPLHSISVRRLKETKDGFMFLTDQSYIDVLPE
    mutant_tumor   GAKVLSSVKDRFINFFVGNTINSSYFPDHPLHSISVRRLKETKDGFMFLTDQSYIDVLPE
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     FRDSYPIKYVHAFESNNFIYFLTVQRETLDAQTFHTRIIRFCSINSGLHSYMEMPLECIL
    mutant_tumor   FRDSYPIKYVHAFESNNFIYFLTVQRETLDAQTFHTRIIRFCSINSGLHSYMEMPLECIL
                   ************************************************************ 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     TEKRKKRSTKKEVFNILQAAYVSKPGAQLARQIGASLNDDILFGVFAQSKPDSAEPMDRS
    mutant_tumor   TEKRKKRSTKKEVFNILQAAYVSKPGAQLARQIGASLNDDILFGVFAQSKPDSAEPMDRS
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     AMCAFPIKYVNDFFNKIVNKNNVRCLQHFYGPNHEHCFNRTLLRNSSGCEARRDEYRTEF
    mutant_tumor   AMCAFPIKYVNDFFNKIVNKNNVRCLQHFYGPNHEHCFNRTLLRNSSGCEARRDEYRTEF
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .         .         .         480 
    wt_healthy     TTALQRVDLFMGQFSEVLLTSISTFIKGDLTIANLGTSEGRFMQVVVSRSGPSTPHVNFL
    mutant_tumor   TTALQRVDLFMGQFSEVLLTSISTFIKGDLTIANLGTSEGRFMQVVVSRSGPSTPHVNFL
                   ************************************************************ 
                 421        .         .         .         .         .         480 

                 481        .         .         .         .         .         540 
    wt_healthy     LDSHPVSPEVIVEHTLNQNGYTLVITGKKITKIPLNGLGCRHFQSCSQCLSAPPFVQCGW
    mutant_tumor   LDSHPVSPEVIVEHTLNQNGYTLVITGKKITKIPLNGLGCRHFQSCSQCLSAPPFVQCGW
                   ************************************************************ 
                 481        .         .         .         .         .         540 

                 541        .         .         .         .         .         600 
    wt_healthy     CHDKCVRSEECLSGTWTQQICLPAIYKVFPNSAPLEGGTRLTICGWDFGFRRNNKFDLKK
    mutant_tumor   CHDKCVRSEECLSGTWTQQICLPAIYKVFPNSAPLEGGTRLTICGWDFGFRRNNKFDLKK
                   ************************************************************ 
                 541        .         .         .         .         .         600 

                 601        .         .         .         .         .         660 
    wt_healthy     TRVLLGNESCTLTLSESTMNTLKCTVGPAMNKHFNMSIIISNGHGTTQYSTFSYVDPVIT
    mutant_tumor   TRVLLGNESCTLTLSESTMNTLKCTVGPAMNKHFNMSIIISNGHGTTQYSTFSYVDPVIT
                   ************************************************************ 
                 601        .         .         .         .         .         660 

                 661        .         .         .         .         .         720 
    wt_healthy     SISPKYGPMAGGTLLTLTGNYLNSGNSRHISIGGKTCTLKSVSNSILECYTPAQTISTEF
    mutant_tumor   SISPKYGPMAGGTLLTLTGNYLNSGNSRHISIGGKTCTLKSVSNSILECYTPAQTISTEF
                   ************************************************************ 
                 661        .         .         .         .         .         720 

                 721        .         .         .         .         .         780 
    wt_healthy     AVKLKIDLANRETSIFSYREDPIVYEIHPTKSFISGGSTITGVGKNLNSVSVPRMVINVH
    mutant_tumor   AVKLKIDLANRETSIFSYREDPIVYEIHPTKSFISGGSTITGVGKNLNSVSVPRMVINVH
                   ************************************************************ 
                 721        .         .         .         .         .         780 

                 781        .         .         .         .         .         840 
    wt_healthy     EAGRNFTVACQHRSNSEIICCTTPSLQQLNLQLPLKTKAFFMLDGILSKYFDLIYVHNPV
    mutant_tumor   EAGRNFTVACQHRSNSEIICCTTPSLQQLNLQLPLKTKAFFMLDGILSKYFDLIYVHNPV
                   ************************************************************ 
                 781        .         .         .         .         .         840 

                 841        .         .         .         .         .         900 
    wt_healthy     FKPFEKPVMISMGNENVLEIKGNDIDPEAVKGEVLKVGNKSCENIHLHSEAVLCTVPNDL
    mutant_tumor   FKPFEKPVMISMGNENVLEIKGNDIDPEAVKGEVLKVGNKSCENIHLHSEAVLCTVPNDL
                   ************************************************************ 
                 841        .         .         .         .         .         900 

                 901        .         .         .         .         .         960 
    wt_healthy     LKLNSELNIEWKQAISSTVLGKVIVQPDQNFTGLIAGVVSISTALLLLLGFFLWLKKRKQ
    mutant_tumor   LKLNSELNIEWKQAISSTVLGKVIVQPDQNFTGLIAGVVSISTALLLLLGFFLWLKKRKQ
                   ************************************************************ 
                 901        .         .         .         .         .         960 

                 961        .         .         .         .         .         1020 
    wt_healthy     IKDLGSELVRYDARVHTPHLDRLVSARSVSPTTEMVSNESVDYRATFPEDQFPNSSQNGS
    mutant_tumor   IKDLGSELVRYDARVHTPHLDRLVSARSVSPTTEMVSNESVDYRATFPEDQFPNSSQNGS
                   ************************************************************ 
                 961        .         .         .         .         .         1020 

                1021        .         .         .         .         .         1080 
    wt_healthy     CRQVQYPLTDMSPILTSGDSDISSPLLQNTVHIDLSALNPELVQAVQHVVIGPSSLIVHF
    mutant_tumor   CRQVQYPLTDMSPILTSGDSDISSPLLQNTVHIDLSALNPELVQAVQHVVIGPSSLIVHF
                   ************************************************************ 
                1021        .         .         .         .         .         1080 

                1081        .         .         .         .         .         1140 
    wt_healthy     NEVIGRGHFGCVYHGTLLDNDGKKIHCAVKSLNRITDIGEVSQFLTEGIIMKDFSHPNVL
    mutant_tumor   NEVIGRGHFGCVYHGTLLDNDGKKIHCAVKSLNRITDIGEVSQFLTEVIIMKDFSHPNVL
                   *********************************************** ************ 
                1081        .         .         .         .         .         1140 

                1141        .         .         .         .         .         1200 
    wt_healthy     SLLGICLRSEGSPLVVLPYMKHGDLRNFIRNETHNPTVKDLIGFGLQVAKGMKYLASKKF
    mutant_tumor   SLLGICLRSEGSPLVVLPYMKHGDLRNFERNETHNPTVKDLIGFGLQVAKGMKYLASKKF
                   **************************** ******************************* 
                1141        .         .         .         .         .         1200 

                1201        .         .         .         .         .         1260 
    wt_healthy     VHRDLAARNCMLDEKFTVKVADFGLARDMYDKEYYSVHNKTGAKLPVKWMALESLQTQKF
    mutant_tumor   VHRDLAARNCMLDEKFTVKVADFGLARDMRDKEYYSVHNKTGAKLPVKWMALESLQTQKF
                   ***************************** ****************************** 
                1201        .         .         .         .         .         1260 

                1261        .         .         .         .         .         1320 
    wt_healthy     TTKSDVWSFGVLLWELMTRGAPPYPDVNTFDITVYLLQGRRLLQPEYCPDPLYEVMLKCW
    mutant_tumor   YTKSDVWSFGVLLWELMTRGAPPYPDVNTFDITVYLLQGRRLLQPEYCPDPLYEVMLKCW
                    *********************************************************** 
                1261        .         .         .         .         .         1320 

                1321        .         .         .         .         .         1380 
    wt_healthy     HPKAEMRPSFSELVSRISAIFSTFIGEHYVHVNATYVNVKCVAPYPSLLSSEDNADDEVD
    mutant_tumor   HPKAEMRPSFSELVSRISAIFSTFIGEHYVHVNATYVNVKCVAPYPSLLSSEDNADDEVD
                   ************************************************************ 
                1321        .         .         .         .         .         1380 

                1381        1390 
    wt_healthy     TRPASFWETS
    mutant_tumor   TRPASFWETS
                   ********** 
                1381        1390 

    Call:
      read.fasta(file = "A18488012_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 1390 position columns (1390 non-gap, 0 gap) 

    + attr: id, ali, call

``` r
blast<-blast.pdb(fast)
```

    Warning in blast.pdb(fast): Multiple sequences detected - using only the first
    sequence in input object

     Searching ... please wait (updates every 5 seconds) RID = V6F4YWBH016 
     .................................
     Reporting 2245 hits

We could score residue conservations…

``` r
mutation.sites <- which(conserv(fast) <1)


paste(fast$ali[1,mutation.sites],
      mutation.sites,
      fast$ali[2,mutation.sites])
```

    [1] "G 1128 V" "I 1169 E" "Y 1230 R" "T 1261 Y"
