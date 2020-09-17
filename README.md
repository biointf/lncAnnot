# lncAnnot

lncRNA annotation pipeline for GeneChip microarrays

In line with previous considerations, we developed a procedure that allowed to re-annotate CDF files according to the principle 'OPOG' (one probe, one gene).

## CDF from University of Michigan 'BrainArray' repository

The annotations can be dowloaded at http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp from University of Michigan. To date, the repository from the University of Michigan represents a comprehensive and updated source of annotations for Genechip data. For our purpose, the user should download the last release of the package based on GENCODET annotation   containing **protein coding** and **non-coding** transcripts. The compressed archive available at the Michigan website (column "Zip of CDF, Seq, Map, Desc", label "Z") includes (i) mapping, (ii) description, (iii) probes sequence and (iv) CDF files.

**Here, exemplarly, we described the procedure for the last version (v24) of the last-generation Human ClariomD arrays. However, it can be run with any other arrays upon substitution of any string referred to ClariomD with the appropriate arrays.**

1. Go to http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
2. Download GENECODET package, e.g. http://mbni.org/customcdf/24.0.0/gencodet.download/ClariomDHuman_Hs_GENCODET_24.0.0.zip
3. Unzip to ClariomDHuman_Hs_GENCODET_24.0.0 folder

## Probe filter strategy and Flat file creation

We applied 2 filters on the probe/probeset downloaded.

1. Probe mapping multiple ENSGENE element will be filtered out
2. Probesets with less than 4 probes will be filtered out

To apply those filters and to produce a Flat file for the next procedure, for convenience we developed a script named brainArray2Flat.sh that can be dowloaded here:
https://github.com/emacgene/lncAnnot

The script can be run from the shell of any unix-based system with the following command:
```
/brainArray2Flat.sh -p ClariomDHuman_Hs_GENCODET_probe_tab -d ClariomDHuman_Hs_GENCODET_desc.txt -n 4 -o ClariomDHuman_Hs_GENCODET_24
```
The syntax represents the exemplar application of the following function:

```
brainArray2Flat.sh -p <FILE.probe_tab> -d <FILE.desc.txt> -o <OUT_FILENAME> [ -n <MIN_PROBES_NUMBER> ] [-h]

  Flag:
        -p/--probe_tab  FILE    Probe tab file from Brain Array CDF package
        -d/--desc       FILE    Description file from Brain Array CDF package
        -o/--out        PREFIX  Output prefix name

  Optional Flag:
        -n/--probe_th   INT     Min number of probes used as threshold to select a probeset [DEFAULT = 4]
                                Es. -n 4 --> Probeset with < 4 probe will be filtered out
```

where the optional parameter -n represent the _minimum_ (and default) number of probes that should be considered to build the transcript expression signal (we have chosen to set 4 ad default number according to previous evidences that recognize this as the minimum set to generate a reliable signal [REF Ferrari, 2007 BMC Bioinformatics]).

## Creation of flat files

The flat2Cdf function from affxparser package for Bioconductor can now be run to create the newly annotated CDF file.
For convenience, a modified version of the flat2Cdf function that add tags with "." instead of "," is included as source file and can be dowloaded at https://github.com/emacgene/lncAnnot
From R environment, the following functions should be run:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("affxparser")
library(affxparser)
source("./flat2Cdf.R") # set the correct path to the downloaded file
flat2Cdf("~/Downloads/ClariomDHuman_Hs_GENCODET_24.0.0/ClariomD_GENECODET.flat", chipType="ClariomD.ENSG", tag="v24", col.class=c("character","integer","integer","character","character","character"), rows=2572, cols= 2680, xynames=c("X","Y"), ucol=5, gcol=6)
```

1. Concerning this last function, please note that it works correctly only if the correct number of rows and columns is indicated in the parameters. To know the effective value, users can easily assess the header of _any_ CEL files generated for the array of interest with the affxparser::readCelHeader() function.
2. 'ucol' and 'gcol' parameters represent, respectively, units and groups to be considered to associate the probe to a probeset (and consequently to the output). Ideally, starting from GENECODET data users could choose to genereate ENST, i.e. transcript-based, or ENSG, i.e. gene-based annotations.
This is a major aspect of lncRNA annotation, since it is a prerogative of the user to choose at which level the lncRNA transcripts would be investigated.

--> describe difference (advatages and disadvantages of transcript/gene -level analysis)

## Creation of CDF files

Later, a conventional procedure can be applied to transform the .CDF file into usable R package:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("makecdfenv")
library(makecdfenv)
pkgpath=("Set_the_Path_to_CDF_file")
make.cdf.package("ClariomDHuman_Hs_GENCODET.cdf", compress = FALSE, species="Homo_sapiens", unlink=TRUE, cdf.path = pkgpath, package.path = pkgpath)
```
Finally, the library can be installed and loaded into R (optionally, from the Unix shell with R CMD command):
```
# R CMD build --force clariomdhumanhsgencodetcdf
# R CMD INSTALL clariomdhumanhsgencodetcdf_1.64.0.tar.gz
library(clariomdhumanhsgencodetcdf)
```

At this point, users could virtually choose custom procedure to generate the expression data. However, here we prefer to indicate the usage of the conventional affy() package and RMA normalization method, considering this last one as the gold standard for microarray data.  
