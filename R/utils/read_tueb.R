# code for reading in data pairs is borrowed from 
# Marx, A. and Vreeken, J. Telling Cause from Effect using MDL-based Local and Global Regression.
# In ICDM, 2017

uv <- c(1:46, 48:51,56:69,72:104,106,108)
ref <- read.table("data/tuebingen_benchmark/README_polished_may18.tab", sep="\t", header=F, stringsAsFactors = F)
meta <- read.table("data/tuebingen_benchmark/pairmeta.txt")
ref.uv <- ref[uv, ]

readI = function(i, r = ref){
  f <- paste(c("data/tuebingen_benchmark/", r$V1[i], ".txt"), collapse="")
  t <- read.table(f, sep=" ", header=F, stringsAsFactors = F)
  return(t)
}


