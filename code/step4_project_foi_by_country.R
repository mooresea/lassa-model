# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load libraries
library(doParallel)
library(tidyverse)

##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100

# load estimated FOI
load(paste0('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData'))
ss$level=admin
##Save locs with foi data for later use
write_csv(ss,file=paste0("results/estimated_foi_adm",admin,'_revrate',as.character(rev_rate*100),".csv"))

# read in demographic data for all adm1s
if(admin==2){
  r = read.csv('data/adm_2_pop_upd.csv')
  r$uid = r$GID_2
}else{
  r = read.csv('data/adm_1_pop_upd.csv')
  r$uid = paste(r$ISO,'.',r$SP_ID_1,'_1',sep='')
}
r = subset(r, r$YEAR>2009)
# r2014 = subset(r, r$YEAR == 2014)
# for(yy in 2015:2023){
#   r = rbind(r, r2014)
# }

# read in adm1s of interest to Lassa and limit to those
if(admin==2){
  c = read.csv('data/lassa_adm2_all_covariates_upd.csv')
  r = subset(r, uid %in% unique(c$GID_2))
}else{
  c = read.csv('data/lassa_adm1_all_covariates_upd.csv')
  r = subset(r, uid %in% unique(c$GID_1))
}


# indices of columns containing population by age data
ind.pop = which(substr(names(r),0,4)=='POP_')

# load projected infections
load(paste0('results/projected_infections_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData'))

# matrix to store FOI values corresponding to projected infections
rr.foi = matrix(NA,nrow(rr.inf),ncol(rr.inf))
row.names(rr.foi) = row.names(rr.inf)

# initiate cluster for parallel computing
cl = makeCluster(4)
registerDoParallel(cl)

# find FOI corresponding to each number of infections
for(ii in 1:nrow(rr)){
  which.ii = which(r$uid==rr$uid[ii])
  print(ii)
  if(rr$uid[ii] %in% ss$uid){
    #Already estimated FOI from serology - use for FOI projections
    ssind=which(ss$uid==rr$uid[ii])
    rr.foi[ii,]=rgamma(ncol(rr.foi),shape=ss$foi.shape[ssind],scale=ss$foi.scale[ssind])
  }else{
    #Project FOI from infections
    rr.foi[ii,] = foreach(jj=1:ncol(rr.inf),.combine=c) %dopar% {
      optimize(f = function(foi){
        abs(sum(as.matrix(((1-(foi/(foi+rev_rate))*(1-exp(-(foi+rev_rate)*(0:99))))*(foi/(foi+rev_rate))*(1-exp(-(foi+rev_rate))))  %*%t(as.matrix(r[which.ii,ind.pop])))) - rr.inf[ii,jj])},
        #abs(sum(as.matrix((exp(-foi*(0:99))*(1-exp(-foi)))%*%t(as.matrix(r[which.ii,ind.pop])))) - rr.inf[ii,jj])},
        interval = c(0,1))$minimum    
    }
  }
}

# stop cluster when finished with parallel computing
stopCluster(cl)

# save projected force of infection and info about each adm1
save(rr,rr.foi,file=paste0('results/projected_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_updated.RData'))

