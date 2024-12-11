# setwd("~/Repos/cepi_lassa/")
rm(list=ls())

##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100

# load serological data
s = read.csv('data/case_reports_Lassa_IGG.csv')
s$YEAR = as.numeric(substr(s$YEAR,0,4))
s$AGE_LOWER = s$AGE..lower.
s$AGE_UPPER = s$AGE..upper.
s$AGE_LOWER[is.na(s$AGE_LOWER)] = 0
s$AGE_UPPER[is.na(s$AGE_UPPER)] = 99
s$SAMPLE_SIZE = s$SAMPLE.SIZE

# subset to studies performed during the timeframe of interest
s = s[which(as.numeric(substr(s$YEAR,0,3)) >= 198),]
##Subset to studies on general population
s = s[which(s$OCCUPATION==""|is.na(s$OCCUPATION)),]
# subset to adm1 level data
if(admin==2){
  s = subset(s, ADM2_code != '')  
  s$uid = s$ADM2_code
}else{
  s = subset(s, ADM1_code != '')  
  s$uid = s$ADM1_code
}


# load population by age data
a.tmp = read.csv('data/adm_1_pop_upd.csv')
a.tmp$ADM1_code = paste(a.tmp$ISO,'.',a.tmp$SP_ID_1,'_1',sep='')
a = a.tmp
for(ii in 1:nrow(s)){
  a[ii,] = a.tmp[which(paste(a.tmp$ADM1_code,a.tmp$YEAR) == paste(s$ADM1_code,pmin(2014,s$YEAR))[ii]),]
}
a = a[1:nrow(s),]

# indices of columns containing population by age data
ind.pop = which(substr(names(a),0,4)=='POP_')

# unique id by study
ss = data.frame(
  uid = sort(unique(s$uid)))
ss$POSITIVE = aggregate(s$POSITIVE,list(s$uid),sum)[which(sort(unique(s$uid))%in%ss$uid),2]
ss$SAMPLE_SIZE = aggregate(s$SAMPLE_SIZE,list(s$uid),sum)[which(sort(unique(s$uid))%in%ss$uid),2]

# calculate parameters of beta distribution for posterior
# probability of being seropositive in each study
ss$alpha = 1 + ss$POSITIVE
ss$beta = 1 + ss$SAMPLE_SIZE - ss$POSITIVE

# specify a range of force of infection to work with
foi = 10 ^ seq(-6, 1, length.out=200)

# determine the probability of seropositivity associated with each
# force of infection in each study
##UPDATED code to include potential for seroreversion
p.by.foi = sapply(1:length(foi),function(jj)sapply(1:nrow(s),function(ii)weighted.mean((foi[jj]/(foi[jj]+rev_rate))*(1-exp(-(rev_rate+foi[jj])*(s$AGE_LOWER[ii]:s$AGE_UPPER[ii]))),a[ii,ind.pop[s$AGE_LOWER[ii]:s$AGE_UPPER[ii]+1]],na.rm=T)))

p.by.foi.agg = matrix(NA,nrow(ss),ncol(p.by.foi))
for(ii in 1:nrow(p.by.foi.agg)){
  p.by.foi.agg[ii,] = sapply(1:length(foi),function(jj)weighted.mean(p.by.foi[which(s$uid==ss$uid[ii]),jj],s$SAMPLE_SIZE[which(s$uid==ss$uid[ii])]))
}

# fit gamma distributions describing uncertainty in force of infection
# that corresponds to beta distributions describing uncertainty in
# probability of being seropostive in each study
shape.scale = matrix(0,nrow(p.by.foi.agg),2)
for(ii in 1:nrow(p.by.foi.agg)){
  shape.scale[ii,] = exp(optim(par=c(0,0),fn=function(par)sum((pgamma(foi,shape=exp(par[1]),scale=exp(par[2])) - pbeta(p.by.foi.agg[ii,],ss$alpha[ii],ss$beta[ii]))^2))$par)
}
ss$foi.shape = shape.scale[,1]
ss$foi.scale = shape.scale[,2]

ss$foi.median = qgamma(.5,shape=ss$foi.shape,scale=ss$foi.scale)
ss$foi.lo = qgamma(.025,shape=ss$foi.shape,scale=ss$foi.scale)
ss$foi.hi = qgamma(.975,shape=ss$foi.shape,scale=ss$foi.scale)
# save output to file
save(list=ls(),file=paste0('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData'))
