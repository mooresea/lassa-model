# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load libraries
library(BayesianTools)
library(doParallel)
library(foreach)
library(mc2d)

##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100

# load posterior estimates of force of infection
load(paste('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData',sep=''))
##REmove COG from study
ss=ss[which(substring(ss$uid,1,3)!="COG"),]

# read in data on population by age for each year 2010+
if(admin==2){
  d = read.csv('data/adm_2_pop_upd.csv')
  #d$uid = paste(d$ISO,'.',d$SP_ID_1,'_1',sep='')
  d = subset(d, d$YEAR >= 2010) #1980)

}else{
  d = read.csv('data/adm_1_pop_upd.csv')
  d$uid = paste(d$ISO,'.',d$SP_ID_1,'_1',sep='')
  d = subset(d, d$YEAR >= 2010) #1980)
  # d2014 = subset(d, d$YEAR == 2014)
  # for(yy in 2015:2019){
  #   d = rbind(d, d2014)
  # }
}


# indices of columns containing population by age data
ind.pop = which(substr(names(d),0,4)=='POP_')

# set up data frame to store priors for Dirichlet parameters for proportions of infection states
ob = read.csv('data/case_reports_Lassa_CASES.csv')
if(admin==2){
  ob$uid = ob$ADM2_code 
}else{
  ob$uid = ob$ADM1_code 
}
ob$year = as.numeric(substr(ob$YEAR,0,4))
ob = subset(ob, ob$year >= 2010) #1980)
ob = ob[which(ob$uid %in% ss$uid),]

# aggregate outbreak data across 1980+ at adm level
ob = data.frame(
  uid = sort(unique(ob$uid)),
  cases = aggregate(ob$CASES,list(ob$uid),sum,na.rm=T)[,2],
  deaths = aggregate(ob$DEATHS,list(ob$uid),sum,na.rm=T)[,2])
ob$ADM0=substring(ob$uid,1,3)

# add outbreak data to serological dataframe
ss$ob.cases = ss$ob.deaths = 0
for(ii in 1:nrow(ss)){
  if(as.character(ss$uid[ii]) %in% as.character(ob$uid)){
    ss[ii,c('ob.cases','ob.deaths')] = ob[which(as.character(ob$uid)==as.character(ss$uid[ii])),c('cases','deaths')]
  }
}

# matrix of unvaccinated person years by year (row) and age (column)
unvax = as.matrix(sapply(ss$uid,function(uu)colSums(d[which(d$uid%in%uu),ind.pop])),100,nrow(ob))
unvax[is.na(unvax)] = 0
ss$unvax = ceiling(colSums(unvax))

# load priors for probabilities of asymptomatics, cases, and deaths
load('results/prior_iceberg.RData')

##
## TRY TO RUN THIS FOR EACH COUNTRY SEP (will it blow up on countries with no cases reported?)
##
# functions to transform and inverse transform parameters
transform.par = function(par.in){
  c(log(par.in[c(ind.Ui,ind.F)] / (1 - par.in[c(ind.Ui,ind.F)])),
    log(par.in[-c(ind.Ui,ind.F)]))
}
inv.transform.par = function(par.in){
  c(exp(par.in[c(ind.Ui,ind.F)]) / (1 + exp(par.in[c(ind.Ui,ind.F)])),
    exp(par.in[-c(ind.Ui,ind.F)]))
}

# function for log likelihood
logLik = function(par){
  par = inv.transform.par(par)
  
  #IAR.by.age = 1 - exp(-par[ind.FOIi] %*% matrix(0:99,1,100))
  IAR.by.age = (par[ind.FOIi]/(par[ind.FOIi]+rev_rate))*(1 - exp(-(rev_rate+par[ind.FOIi]) %*% matrix(0:99,1,100)))
  
  yes.prev = rowSums(IAR.by.age %*% unvax) / sum(unvax)
  #no.prev.no.curr = (1 - yes.prev) * exp(-par[ind.FOIi])  
  no.prev.no.curr = 1-((par[ind.FOIi]/(par[ind.FOIi]+rev_rate))+( yes.prev-par[ind.FOIi]/(par[ind.FOIi]+rev_rate)) * exp(-(rev_rate+par[ind.FOIi])))
  no.prev.yes.curr = 1 - yes.prev - no.prev.no.curr
  
  sum(dmultinomial(
    cbind(ss$ob.deaths, ss$ob.cases, ss$unvax-ss$ob.deaths-ss$ob.cases),
    ss$unvax,
    cbind(no.prev.yes.curr * (1 - par[ind.Ui]) * par[ind.F],
          no.prev.yes.curr * (1 - par[ind.Ui]) * (1 - par[ind.F]),
          yes.prev + no.prev.no.curr + no.prev.yes.curr * par[ind.Ui]),
    log=T)) +
    sum(dbeta(par[ind.Ui], par[ind.Ualpha], par[ind.Ubeta], log=T))  
}

# function for log prior
logPrior = function(par){
  par = inv.transform.par(par)
  
  dbeta(par[ind.F], prior.alpha.AMSF[3], prior.alpha.AMSF[2], log=T) +
    sum(dgamma(par[ind.FOIi], shape=ss$foi.shape[1:length(ind.FOIi)], scale=ss$foi.scale[1:length(ind.FOIi)]), log=T)
}

##By country
ss$ADM0=substring(ss$uid,1,3)
ss.countries=unique(ss$ADM0)
ss.all=ss
ss.all=ss.all[order(ss.all$uid),]
mcmc=list()
prop.FCU=matrix(NA,nrow=length(ss.countries)+1,ncol=3)
posterior.foi=matrix(NA,nrow=15003,ncol=nrow(ss.all))
for(cc in 1:length(ss.countries)){
  print(cc)
    cc.ind=which(ss.all$ADM0==ss.countries[cc])
    ss=ss.all[cc.ind,]
    # parameters - Ui, F, Ualpha, Ubeta, FOIi
    par = c(
      1 - (ss$ob.cases+ss$ob.deaths) * (sum(prior.alpha.AMSF)/sum(prior.alpha.AMSF[2:3])) / ((1-exp(-ss$foi.shape*ss$foi.scale))*ss$unvax),
      prior.alpha.AMSF[3] / sum(prior.alpha.AMSF[2:3]),
      99,
      1,
      ss$foi.shape*ss$foi.scale)
    ind.Ui = 1:nrow(ss)
    ind.F = nrow(ss) + 1
    ind.Ualpha = nrow(ss) + 2
    ind.Ubeta = nrow(ss) + 3
    ind.FOIi = nrow(ss) + 3 + (1:nrow(ss))
    par.lower = c(
      rep(0.1, length(ind.Ui)),
      0.01,
      10,
      1,
      qgamma(0.001,shape=ss$foi.shape,scale=ss$foi.scale))
    par.upper = c(
      rep(1-1e-6, length(ind.Ui)),
      0.99,
      500,
      50,
      qgamma(0.999,shape=ss$foi.shape,scale=ss$foi.scale))
    
    # perform MCMC to estimate parameters
    bayesianSetup = createBayesianSetup(
      likelihood = function(par){logLik(par) + logPrior(par)},
      lower = transform.par(par.lower), upper = transform.par(par.upper))
    settings = list(iterations = 3e5)
    mcmc[[cc]] = runMCMC(bayesianSetup, settings = settings)
    gelmanDiagnostics(mcmc[[cc]],thin=1e1,start=1e5-5e4,end=NULL)
    post = getSample(mcmc[[cc]],thin=1e1,start=1e5-5e4,end=NULL)
    
    # make draws from the approximate posterior distribution
    for(ii in 1:nrow(post)){
      post[ii,] = inv.transform.par(post[ii,])
    }
    
    # use the approximate posterior to generate a distribution of proportions of deaths, cases, and unreported
    propns = matrix(0,nrow(post),3)
    propns[,3] = rbeta(nrow(propns), post[,ind.Ualpha], post[,ind.Ubeta])
    propns[,1] = (1 - propns[,3]) * post[,ind.F]
    propns[,2] = 1 - rowSums(propns)
    propns = propns[which(apply(propns,1,max)<1),]
    if(cc==1){
      pp=propns
    }else{
      pp=rbind(pp,propns)
    }
    
    # get maximum-likelihood estimate of Dirichlet parameters for deaths, cases, and unreported
    dir.opt = optim(
      par=rep(0,3),
      function(par)-sum(log(ddirichlet(propns,exp(par)))))
    prop.FCU[cc,] = exp(dir.opt$par)
    
    # save posterior estimates of force of infection for each site
    posterior.foi[,cc.ind] = post[,ind.FOIi]
}
colnames(posterior.foi) = ss.all$uid
row.names(prop.FCU)=c(ss.countries,"AVG")

##Get average prop.fcu
ss=ss.all
# parameters - Ui, F, Ualpha, Ubeta, FOIi
par = c(
  1 - (ss$ob.cases+ss$ob.deaths) * (sum(prior.alpha.AMSF)/sum(prior.alpha.AMSF[2:3])) / ((1-exp(-ss$foi.shape*ss$foi.scale))*ss$unvax),
  prior.alpha.AMSF[3] / sum(prior.alpha.AMSF[2:3]),
  99,
  1,
  ss$foi.shape*ss$foi.scale)
ind.Ui = 1:nrow(ss)
ind.F = nrow(ss) + 1
ind.Ualpha = nrow(ss) + 2
ind.Ubeta = nrow(ss) + 3
ind.FOIi = nrow(ss) + 3 + (1:nrow(ss))
par.lower = c(
  rep(0.1, length(ind.Ui)),
  0.01,
  10,
  1,
  qgamma(0.001,shape=ss$foi.shape,scale=ss$foi.scale))
par.upper = c(
  rep(1-1e-6, length(ind.Ui)),
  0.99,
  500,
  50,
  qgamma(0.999,shape=ss$foi.shape,scale=ss$foi.scale))

# perform MCMC to estimate parameters
bayesianSetup = createBayesianSetup(
  likelihood = function(par){logLik(par) + logPrior(par)},
  lower = transform.par(par.lower), upper = transform.par(par.upper))
settings = list(iterations = 3e5)
mcmc[[cc]] = runMCMC(bayesianSetup, settings = settings)
gelmanDiagnostics(mcmc[[cc]],thin=1e1,start=1e5-5e4,end=NULL)
post = getSample(mcmc[[cc]],thin=1e1,start=1e5-5e4,end=NULL)

# make draws from the approximate posterior distribution
for(ii in 1:nrow(post)){
  post[ii,] = inv.transform.par(post[ii,])
}

# use the approximate posterior to generate a distribution of proportions of deaths, cases, and unreported
propns = matrix(0,nrow(post),3)
propns[,3] = rbeta(nrow(propns), post[,ind.Ualpha], post[,ind.Ubeta])
propns[,1] = (1 - propns[,3]) * post[,ind.F]
propns[,2] = 1 - rowSums(propns)
propns = propns[which(apply(propns,1,max)<1),]

# get maximum-likelihood estimate of Dirichlet parameters for deaths, cases, and unreported

dir.opt = optim(
  par=rep(0,3),
  function(par)-sum(log(ddirichlet(propns,exp(par)))))
prop.FCU[nrow(prop.FCU),] = exp(dir.opt$par)

# save parameters describing proportions of infection outcomes
save(mcmc,posterior.foi,prop.FCU,post,file=paste0('results/proportion_by_type_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData'))
