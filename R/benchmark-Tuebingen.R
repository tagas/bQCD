rm(list = ls())

source("R/baselines.R")
source("R/bqcd.R")
source("R/utils/read_tueb.R")

cnt <- 0
pair <- rep("", length(uv))
eps <- rep(0, length(uv))
cds <- rep("--", length(uv))
time <- rep(0, length(uv))

for(i in 1:length(uv)){
  t1 <- Sys.time()
  t <- readI(uv[i])[,1:2]
  res <- QCD_wrap(t[,1], t[,2], "bQCD_qnn")
  #res <- lingamWrapper(cbind(t[,1], t[,2]))
  #res$cd <- ifelse(res$cd==1, "->", "<-")
  t2 <- Sys.time()
  elapsed <- as.numeric(difftime(t2,t1), units="secs")
  pair[i] <- uv[i]
  eps[i] <- res$eps
  cds[i] <- res$cd
  time[i] <- elapsed
  print(paste0("pair = ", i, ", correct = ", ref.uv$V6[i] == cds[i] , 
               ", accuracy = ", mean(ref.uv$V6[1:i] == cds[1:i])))
}

corr <- rep(0,length(uv))
corr[ref.uv$V6 == cds] <- 1

print("Correct causal directions recovered:")
print(sum(corr))

print("Weighted accuracy:")
print(sum(corr*meta$V6[uv])/sum(meta$V6[uv]))

resQCCD = data.frame(Correct = corr, Eps = eps, Cds = cds, T = time)
write.table(resQCCD, file = "./results/bQCD_qnn_m7.tab", row.names = F, quote = F, sep="\t")

