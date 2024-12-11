# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load package for Dirichlet multinomial probability
library(extraDistr)


##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100

# read in demographic data for all adm1s
if(admin==2){
  r = read.csv('data/adm_2_pop_upd.csv')
  r$uid = r$GID_2
}else{
  r = read.csv('data/adm_1_pop_upd.csv')
  r$uid = paste(r$ISO,'.',r$SP_ID_1,'_1',sep='')
}

r = subset(r, r$YEAR %in% 2010:2023)
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

# read in parameters describing proportions of infection outcomes
load(paste0('results/proportion_by_type_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData'))

# set up data frame with one row for each adm
rr = data.frame(
  uid = sort(unique(r$uid)))

# determine the maximum number of infections that could have possibly
# occurred in each adm1. reasoning for this is that in the limit as
# force of infection goes to infinity, everyone born prior to 1980
# would have already been infected and all children would get infected
# during their first year of life. thus, sum newborns over all years.
rr$max.inf = sapply(rr$uid, function(uu){
  r.tmp = subset(r,uid==uu)
  foi = 10
  floor(sum(as.matrix((exp(-foi*(0:99))*(1-exp(-foi)))%*%t(as.matrix(r.tmp[,ind.pop])))))})

# read in data about outbreaks and aggregate over time
ob = read.csv('data/case_reports_Lassa_CASES.csv')
if(admin==2){
  ob$uid = ob$ADM2_code  
}else{
  ob$uid = ob$ADM1_code
}

ob$year = as.numeric(substr(ob$YEAR,0,4))
ob = subset(ob, ob$year >= 2010)
ob = ob[-which(ob$Human.Human..1...Zoonotic..0...Unknown..NA.==1),]
ob = data.frame(
  uid = sort(unique(ob$uid)),
  cases = aggregate(ob$CASES,list(ob$uid),sum,na.rm=T)[,2],
  deaths = aggregate(ob$DEATHS,list(ob$uid),sum,na.rm=T)[,2])
ob = ob[-which(ob$uid==''),]

# add outbreak data to the data set for all adm1s
rr$cases = 0
rr$deaths = 0
for(ii in 1:nrow(ob)){
  which.rr = which(as.character(rr$uid) == as.character(ob$uid[ii]))
  rr[which.rr,c('cases','deaths')] = ob[ii,c('cases','deaths')]
}

# project the number of infections conditional on cases and deaths and proportion unobserved
rr$ADM0=substring(rr$uid,1,3)
infections = list()
for(ii in 1:nrow(rr)){
  print(ii)
  if(rr$ADM0[ii] %in% row.names(prop.FCU)){
    ii.prop.FCU=prop.FCU[which(row.names(prop.FCU)==rr$ADM0[ii]),] 
  }else{
    ii.prop.FCU=prop.FCU[which(row.names(prop.FCU)=="AVG"),]
  }
  x.in = cbind(rr$deaths[ii],rr$cases[ii],0:(rr$max.inf[ii]-rr$cases[ii]-rr$deaths[ii]))
  size.in = (rr$cases[ii]+rr$deaths[ii]):rr$max.inf[ii]
  probs = ddirmnom(x.in,size.in,ii.prop.FCU,log=F)
  probs = probs / sum(probs)
  infections[[ii]] = probs
}

# draw random samples of the number of infections in each adm1
rr.inf = matrix(NA,nrow(rr),1e3)
for(ii in 1:nrow(rr.inf)){
  rr.inf[ii,] = rr$cases[ii] + rr$deaths[ii] - 1 +
    sample(length(infections[[ii]]),ncol(rr.inf),replace=T,infections[[ii]])
}
row.names(rr.inf) = rr$uid

# save sampled infection numbers and info about each adm1
save(rr,rr.inf,file=paste0('results/projected_infections_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData'))
# save summarized case/death observation data for later use
write_csv(ob,file=paste0('results/observed_cases_deaths_adm',admin,".csv"))