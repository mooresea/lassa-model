# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load library
library(extraDistr)
library(tidyverse)

# load data from iceberg paper
iceberg = read_csv('data/iceberg_data.csv')
iceberg$type = sapply(1:nrow(iceberg),function(ii)
  paste(names(iceberg)[which(!is.na(iceberg[ii,]))],collapse='-'))

# function for calculating Dirichlet multinomial likelihood of data
LL.iceberg = function(par){
  d = subset(iceberg,type=='A-MS')[,c('A','MS')]
  LL = sum(ddirmnom(d,rowSums(d),c(par[1],par[2]),log=T))
  
  d = subset(iceberg,type=='F-MS')[,c('MS','F')]
  LL = LL + sum(ddirmnom(d,rowSums(d),c(par[2],par[3]),log=T))

  return(LL)  
}

# MLE of alpha parameters of Dirichlet distribution
prior.alpha.AMSF = exp(optim(rep(0,3),function(par)-LL.iceberg(exp(par)))$par)

# save output to file
save(prior.alpha.AMSF, file='results/prior_iceberg.RData')
