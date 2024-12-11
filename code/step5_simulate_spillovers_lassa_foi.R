# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

# load libraries
library(mc2d)
library(dplyr)


##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100
raw.ind=as.numeric(args[3])

# load force of infection random draws by adm1 (rr.foi)
if(raw.ind==1){
  load(paste0('results/projected_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_updated.RData'))  
}else{
  load(paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_11_outpredict.RData',sep=''))
  load(paste0('results/foi_wtd_adm',admin,'_revrate',as.character(rev_rate*100),'_outpredict.RData'))
  rr.foi=10^foi.all.wtd #[,,1]
  row.names(rr.foi)=c$uid
}


# read in demographic data for all adm1s
# read in demographic data for all adm1s
if(admin==2){
  r = read.csv('data/adm_2_pop_upd.csv')
  r$uid = r$GID_2
}else{
  r = read.csv('data/adm_1_pop_upd.csv')
  r$uid = paste(r$ISO,'.',r$SP_ID_1,'_1',sep='')
}
r = subset(r, r$YEAR==2023)
r = subset(r, r$uid %in% row.names(rr.foi))
r = r[match(row.names(rr.foi),r$uid),]
ind.pop = which(substr(names(r),0,4)=='POP_')

# load proportion of deaths, cases, unobserved (prop.FCU), iceberg
#load(paste0('results/proportion_by_type_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData'))

##Asymp_rate
##Mccormick et al data
n.ill=c(4,1,3,1)
n.pos=c(12,10,15,11)
asymp_mean=sum(n.ill)/sum(n.pos)
asymp_post=rbinom(1000,sum(n.pos),sum(n.ill)/sum(n.pos))/sum(n.pos)

##Relative risk of infection and disease for seropositive individuals - use mean values for now
seropos.risk.inf=0.538
seropos.risk.lf=0.362
serorev.risk.inf=seropos.risk.inf

# set up a matrix to store simulated spillover
# dimension 1 = adm1
# dimension 2 = age
# dimension 2 = FOI replicate
seroprev.age.foi = array(NA,c(nrow(r),100,ncol(rr.foi)))
inf.age.foi=seroprev.age.foi #Track fraction who have ever been infected
spillover.infections = seroprev.age.foi
spillover.infections.inf.risk = seroprev.age.foi
spillover.infections.lf.risk = seroprev.age.foi
spillover.infections.reinf = seroprev.age.foi #Lower risk of disease among reinfected
spillover.infections.inf.risk.reinf = seroprev.age.foi
spillover.infections.lf.risk.reinf = seroprev.age.foi

# set random number seed so that the random draws
# from the Dirichlet distribution are reproducible
set.seed(1234)

# loop across each FOI replicate
for(ii in 1:ncol(rr.foi)){

  # specify the force of infection for each adm1
  foi = rr.foi[,ii]
  
  # specify the age range
  age = 0:99
  
  # matrix of population distribution across adm1s and ages
  pop.by.adm1.by.age = r[,ind.pop]
  
  # probability that a person of a given age in a given
  # adm1 has not yet been infected with Lassa virus
  #pr.susc.by.adm1.by.age.old = 1-exp(-as.matrix(foi)%*%t(as.matrix(age)))
  pr.pos.by.adm1.by.age = (foi/(foi+rev_rate))*(1-exp(-as.matrix(foi+rev_rate)%*%t(as.matrix(age))))
  seroprev.age.foi[,,ii]=pr.pos.by.adm1.by.age
  pr.susc.by.adm1.by.age = 1 - pr.pos.by.adm1.by.age
  
  #probability that a person has ever been infected
  inf.age.foi[,,ii]=1-exp(-as.matrix(foi)%*%t(as.matrix(age)))
  pr.seronaive.by.adm1.by.age = 1 - inf.age.foi[,,ii]
  #fraction who are seronegative but seroreverted
  pr.serorev.by.adm1.by.age=inf.age.foi[,,ii]-seroprev.age.foi[,,ii]
    
  # probability that a person of any age in a given adm1
  # is infected with Lassa virus in a given year
  pr.inf.1.yr = diag(1 - exp(-foi))
  
  # probability that a person of a given age is seronegative
  # and gets infected in a given adm1 in a given year
  pr.inf.by.adm1.by.age.seroneg = pr.inf.1.yr %*% pr.susc.by.adm1.by.age
  
  # probability that a person of a given age is seroreverted
  # and gets infected in a given adm1 in a given year
  pr.inf.by.adm1.by.age.serorev = (pr.inf.1.yr %*% pr.serorev.by.adm1.by.age)*serorev.risk.inf
  # probability that a person of a given age is seronaive
  # and gets infected in a given adm1 in a given year
  pr.inf.by.adm1.by.age.seronaive = pr.inf.1.yr %*% pr.seronaive.by.adm1.by.age
  pr.inf.by.adm1.by.age.serorev.total=pr.inf.by.adm1.by.age.serorev+pr.inf.by.adm1.by.age.seronaive
  
  # probability that a person of a given age is seropositive
  # and gets infected in a given adm1 in a given year
  pr.inf.by.adm1.by.age.seropos.inf.risk = (pr.inf.1.yr %*% pr.pos.by.adm1.by.age)*seropos.risk.inf
  pr.inf.by.adm1.by.age.seropos.lf.risk = (pr.inf.1.yr %*% pr.pos.by.adm1.by.age)*seropos.risk.lf
  
  pr.inf.by.adm1.by.age.no.risk=pr.inf.by.adm1.by.age.seroneg
  pr.inf.by.adm1.by.age.inf.risk=pr.inf.by.adm1.by.age.seroneg+pr.inf.by.adm1.by.age.seropos.inf.risk
  pr.inf.by.adm1.by.age.lf.risk=pr.inf.by.adm1.by.age.seroneg+pr.inf.by.adm1.by.age.seropos.lf.risk
  pr.inf.by.adm1.by.age.inf.risk.serorev.total=pr.inf.by.adm1.by.age.seropos.inf.risk+pr.inf.by.adm1.by.age.serorev.total
  pr.inf.by.adm1.by.age.lf.risk.serorev.total=pr.inf.by.adm1.by.age.seropos.lf.risk+pr.inf.by.adm1.by.age.serorev.total
  
  print(ii)
  # calculate the expected number of infections in each adm1
  for(jj in 1:nrow(pop.by.adm1.by.age)){
    for(kk in 1:ncol(pop.by.adm1.by.age)){
      spillover.infections[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.no.risk[jj,kk])
      spillover.infections.inf.risk[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.inf.risk[jj,kk])
      spillover.infections.lf.risk[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.lf.risk[jj,kk])
      spillover.infections.reinf[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.serorev.total[jj,kk])
      spillover.infections.inf.risk.reinf[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.inf.risk.serorev.total[jj,kk])
      spillover.infections.lf.risk.reinf[jj,kk,ii] = rbinom(1, ceiling(pop.by.adm1.by.age[jj,kk]), pr.inf.by.adm1.by.age.lf.risk.serorev.total[jj,kk])
    }
  }
  # inf = rowSums(inf)
  # 
  # # simulate how many of those infections are deaths, cases,
  # # or unobserved in each adm1
  # FCU.by.adm1 = sapply(inf,function(inf.adm1)rmultinom(1,inf.adm1,prop.FCU.ii))
  # spillover.infections[,ii] = t(FCU.by.adm1[1,] + FCU.by.adm1[2,])
}

# # store simulated spillovers for each adm1
# #adm1s <- select(r, ISO, SP_ID_1)
# adm1s <- subset(r, select=c(ISO, SP_ID_1))
# adm1s <- mutate(adm1s, index=row_number())
# spillover.infections <- as.data.frame(spillover.infections)
# spillover.infections <- mutate(spillover.infections, index=row_number())
# spillover.infections <- merge(spillover.infections, adm1s, by="index")
# #spillover.infections <- select(spillover.infections, -c("index"))
# spillover.infections <- spillover.infections[,-grep("index", colnames(spillover.infections))]
# spillover.infections <- spillover.infections[,c(1001, 1002, 1:1000)]
# write.csv(spillover.infections, file="results/simulated_spillovers_recent.csv")

# save projected force of infection and info about each adm1
save(r,rr.foi,seroprev.age.foi,spillover.infections,spillover.infections.inf.risk,
     spillover.infections.lf.risk,spillover.infections.inf.risk.reinf,
     spillover.infections.lf.risk.reinf,spillover.infections.reinf,
     file=paste0('results/simulated_spillovers_adm',admin,'_revrate',as.character(rev_rate*100),'_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData'))
