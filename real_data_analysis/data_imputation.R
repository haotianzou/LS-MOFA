library(tidyverse)
library(Amelia)

load('dataset/new/multi-omics.RData')
V = data_all$V
N = data_all$N
Y_all = data_all$Y
metadata_all = data_all$metadata
P_all = c(ncol(Y_all[[1]]), ncol(Y_all[[2]]))  ## 749, 104

seed = 1
set.seed(2024*seed)

##
t.argvals = seq(0, 1, by = 0.01)
lt = length(t.argvals)
uID_all = 1:N

## Impute the missing value: 1 missing value in the second omics data on the 103-th feature
M = 5
long = data.frame(ID = metadata_all[[2]]$ID, time = metadata_all[[2]]$time, Y = Y_all[[2]][, 103])
sum(is.na(long$Y))

long.impute <- amelia(x = long, m = M, idvars = 'ID', ts = 'time', splinetime = 6)
long.impute.dataset <- long.impute$imputations

long.i <- long
J = 1
for (j in 1:J){
  tmp.y <- rep(0, nrow(long))
  for (m in 1:M){
    tmp.y <- tmp.y + long.impute.dataset[[m]]$Y
  }
  long.i[, j+2] <- tmp.y/M
}

Y_all[[2]][, 103] = long.i$Y


data_all = list(Y = Y_all, metadata = metadata_all, long = data_all$long, 
                max.time = data_all$max.time, V = V, N = N)
save(list = 'data_all', file = 'dataset/new/multi-omics_imputed.RData')
