closeAllConnections()
rm(list = ls())

print("This version is corrected for breaks in the codon file")

setwd("~/PhyloP_PhastCons/chip_exo")

#library(Biostrings)
#library(devtools)
#library(tools)
#library(rBLAST)
library(stringr)
#library(plyr)

AAcodonfilelist <- read.table("list_of_AA_codon_files")
#AAcodonfilelist <- read.table("temp")
names(AAcodonfilelist) <- c("AAcodonfile")

for (a in 1:length(AAcodonfilelist$AAcodonfile)) {
  print(AAcodonfilelist$AAcodonfile[a])
  #findZFout file open
  inpuputname <-  str_split(AAcodonfilelist$AAcodonfile[a], "_AA")
  uniprotname <-  str_split(AAcodonfilelist$AAcodonfile[a], "_")
  uniprotid <- uniprotname[[1]][2]
  genename <-uniprotname[[1]][3]
  
  #--------------------------------------------------------
  ##read findznf
  dfFZ <-
    read.csv(
      paste0(inpuputname[[1]][1], ".pfamZfs.table.txt", collapse = ""),
      sep = "\t",
      header = T,
      row.names = NULL
    )
  # get the start and end of ZNF helix AA pos coordinates
  dfFZ$start    <-
    str_locate(as.character(dfFZ$ZF.sequence),
               as.character(dfFZ$Helix.sequence))[, 1]
  dfFZ$end      <-
    str_locate(as.character(dfFZ$ZF.sequence),
               as.character(dfFZ$Helix.sequence))[, 2]
  dfFZ$Helixstart <- dfFZ$ZF.start.position + dfFZ$start - 1
  dfFZ$Helixend   <- dfFZ$ZF.start.position + dfFZ$end - 1
  dfFZ$start <- NULL
  dfFZ$end <- NULL
  
  system("rm -f temp.dat")
  system("rm -f temp.bedGraph")
  for (w in 1:length(dfFZ$Helixstart)) {
    system(
      paste(
        "transvar panno --gencode -i '",
        genename,
        ":p.",
        df$Helixstart[w],
        "_",
        df$Helixend[w],
        "' --uniprot > temp.dat;",
        sep = ""
      )
    )
    ##read transvar output
    tempfile <-
      read.csv(
        as.character("temp.dat"),
        sep = "\t",
        header = T,
        row.names = NULL
      )
    dfFZ$coordinates[w] <-
      as.character(tempfile$coordinates.gDNA.cDNA.protein.[1])
    dfFZ$strand[w] <- as.character(tempfile$strand[1])
  }
  
  #--------------------------------------------------------
  inpuputname <-  str_split(AAcodonfilelist$AAcodonfile[a], "_AA")
  dfAA <-
    read.csv(
      as.character(AAcodonfilelist$AAcodonfile[a]),
      sep = ",",
      header = T,
      row.names = NULL
    )
  #--------------------------------------------------------
  
  y = 0
  for (y in 1:length(dfAA$Chrom)) {
    #codon1
    system("rm -f temp.bedGraph")
    system(
      paste(
        "~/soft/kentUtils/bin/bigWigToBedGraph ../data_phyloP46way/",
        dfAA$Chrom[1],
        ".phyloP46way.bw temp.bedGraph -chrom=",
        as.character(dfAA$Chrom[1]),
        " -start=",
        as.character(as.numeric(dfAA$Codon1[y]) - 1),
        " -end=",
        as.character(as.numeric(dfAA$Codon1[y])),
        sep = ""
      )
    )
    bedfile1 <-
      read.csv(
        as.character("temp.bedGraph"),
        sep = "\t",
        header = F,
        row.names = NULL
      )
    
    #codon2
    system("rm -f temp.bedGraph")
    system(
      paste(
        "~/soft/kentUtils/bin/bigWigToBedGraph ../data_phyloP46way/",
        dfAA$Chrom[1],
        ".phyloP46way.bw temp.bedGraph -chrom=",
        as.character(dfAA$Chrom[1]),
        " -start=",
        as.character(as.numeric(dfAA$Codon2[y]) - 1),
        " -end=",
        as.character(as.numeric(dfAA$Codon2[y])),
        sep = ""
      )
    )
    bedfile2 <-
      read.csv(
        as.character("temp.bedGraph"),
        sep = "\t",
        header = F,
        row.names = NULL
      )
    #codon3
    system("rm -f temp.bedGraph")
    system(
      paste(
        "~/soft/kentUtils/bin/bigWigToBedGraph ../data_phyloP46way/",
        dfAA$Chrom[1],
        ".phyloP46way.bw temp.bedGraph -chrom=",
        as.character(dfAA$Chrom[1]),
        " -start=",
        as.character(as.numeric(dfAA$Codon3[y]) - 1),
        " -end=",
        as.character(as.numeric(dfAA$Codon3[y])),
        sep = ""
      )
    )
    bedfile3 <-
      read.csv(
        as.character("temp.bedGraph"),
        sep = "\t",
        header = F,
        row.names = NULL
      )
    
    dfAA$strand[y] <- dfFZ$strand[1]
    dfAA$phyloP[y] <-
      paste(bedfile1$V4, bedfile2$V4, bedfile3$V4, sep = "|")
    
  }
  
  write.csv(dfAA,
            file = paste0(AAcodonfilelist$AAcodonfile[a],
                          "_phyloP.csv",
                          sep = ""))
  
  ##NEWCSV1_Helix_contact_AA
  for (k in 1:nrow(dfFZ)) {
    dfFZ$phyloP[k] <-
      paste0(dfAA$phyloP[dfFZ$Helixstart[k]:dfFZ$Helixend[k]], collapse = "|")
    
  }
  
  write.csv(
    dfFZ,
    file = paste0(
      cbind(unlist(inpuputname))[1],
      "_ZNFinder_output_withA_Allele_Frequency_phyloP.csv",
      sep = ""
    )
  )
  
}
