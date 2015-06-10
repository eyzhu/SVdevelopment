setwd("/wdata/ezh/forestSV_development/test_generate_features_1kGenome")

source("../functions.R")
source("http://bioconductor.org/biocLite.R")
biocLite("scanBamHeader")

myFiles = list.files(path = "./data", "bam$", recursive = T)

for (myF in myFiles) {
      
      basename = gsub("\\..*", "", myF)
      allFeatures(myF,basename)
}