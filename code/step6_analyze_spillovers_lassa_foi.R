library(ggplot2)
library(tidyverse)

admin=2
rev_rate=0.06
raw.ind=1
# save projected force of infection and info about each adm1
load(file=paste0('results/simulated_spillovers_adm',admin,'_revrate',as.character(rev_rate*100),'_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_risk.RData'))

ind.pop=grep("POP",names(r))
r$total_pop=apply(r[,ind.pop],1,sum)
r$pop_18plus=apply(r[,ind.pop[19:100]],1,sum)
r$pop_u12=apply(r[,ind.pop[1:12]],1,sum)
r$pop_12_17=apply(r[,ind.pop[13:18]],1,sum)

if(admin==1){
  adm.names=subset(r, select=c(COUNTRY, ISO, ADMIN_1,uid,total_pop,pop_18plus,pop_u12,pop_12_17))
  adm.names=adm.names %>% rename(Country=COUNTRY,NAME_1=ADMIN_1)
}else{
  adm.names=subset(r, select=c(Country, GID_0, NAME_1,NAME_2,uid,total_pop,pop_18plus,pop_u12,pop_12_17))
}
##order to match rr.foi
adm.names=adm.names[match(row.names(rr.foi),as.character(adm.names$uid)),]
adm.names <- mutate(adm.names, index=row_number())

infs_norisk_mat=matrix(NA,nrow=dim(spillover.infections)[1],ncol = dim(spillover.infections)[3])
infs_infrisk_mat=infs_lfrisk_mat=infs_norisk_mat
infs_norisk_rev_mat=infs_infrisk_rev_mat=infs_lfrisk_rev_mat=infs_norisk_mat
infs_norisk_mat_18=infs_norisk_mat
infs_infrisk_mat_18=infs_lfrisk_mat_18=infs_norisk_mat_18
infs_norisk_rev_mat_18=infs_infrisk_rev_mat_18=infs_lfrisk_rev_mat_18=infs_norisk_mat_18
infs_norisk_mat_17=infs_norisk_mat
infs_norisk_mat_11=infs_norisk_mat
infs_infrisk_mat_17=infs_lfrisk_mat_17=infs_norisk_mat_17
infs_norisk_rev_mat_17=infs_infrisk_rev_mat_17=infs_lfrisk_rev_mat_17=infs_norisk_mat_17
infs_infrisk_mat_11=infs_lfrisk_mat_11=infs_norisk_mat_11
infs_norisk_rev_mat_11=infs_infrisk_rev_mat_11=infs_lfrisk_rev_mat_11=infs_norisk_mat_11
for(i in 1:dim(spillover.infections)[1]){
  if (i%%50==0){print(i)}
  infs_norisk_mat[i,]=apply(spillover.infections[i,,],2,sum)
  infs_infrisk_mat[i,]=apply(spillover.infections.inf.risk[i,,],2,sum)
  infs_lfrisk_mat[i,]=apply(spillover.infections.lf.risk[i,,],2,sum)
  infs_norisk_rev_mat[i,]=apply(spillover.infections.reinf[i,,],2,sum)
  infs_infrisk_rev_mat[i,]=apply(spillover.infections.inf.risk.reinf[i,,],2,sum)
  infs_lfrisk_rev_mat[i,]=apply(spillover.infections.lf.risk.reinf[i,,],2,sum)
  ##18+
  infs_norisk_mat_18[i,]=apply(spillover.infections[i,19:100,],2,sum)
  infs_infrisk_mat_18[i,]=apply(spillover.infections.inf.risk[i,19:100,],2,sum)
  infs_lfrisk_mat_18[i,]=apply(spillover.infections.lf.risk[i,19:100,],2,sum)
  infs_norisk_rev_mat_18[i,]=apply(spillover.infections.reinf[i,19:100,],2,sum)
  infs_infrisk_rev_mat_18[i,]=apply(spillover.infections.inf.risk.reinf[i,19:100,],2,sum)
  infs_lfrisk_rev_mat_18[i,]=apply(spillover.infections.lf.risk.reinf[i,19:100,],2,sum)
  ##12-17 yr old
  infs_norisk_mat_17[i,]=apply(spillover.infections[i,13:18,],2,sum)
  infs_infrisk_mat_17[i,]=apply(spillover.infections.inf.risk[i,13:18,],2,sum)
  infs_lfrisk_mat_17[i,]=apply(spillover.infections.lf.risk[i,13:18,],2,sum)
  infs_norisk_rev_mat_17[i,]=apply(spillover.infections.reinf[i,13:18,],2,sum)
  infs_infrisk_rev_mat_17[i,]=apply(spillover.infections.inf.risk.reinf[i,13:18,],2,sum)
  infs_lfrisk_rev_mat_17[i,]=apply(spillover.infections.lf.risk.reinf[i,13:18,],2,sum)
  ##Under 12
  infs_norisk_mat_11[i,]=apply(spillover.infections[i,1:12,],2,sum)
  infs_infrisk_mat_11[i,]=apply(spillover.infections.inf.risk[i,1:12,],2,sum)
  infs_lfrisk_mat_11[i,]=apply(spillover.infections.lf.risk[i,1:12,],2,sum)
  infs_norisk_rev_mat_11[i,]=apply(spillover.infections.reinf[i,1:12,],2,sum)
  infs_infrisk_rev_mat_11[i,]=apply(spillover.infections.inf.risk.reinf[i,1:12,],2,sum)
  infs_lfrisk_rev_mat_11[i,]=apply(spillover.infections.lf.risk.reinf[i,1:12,],2,sum)
}

infs_norisk_dat <- as.data.frame(infs_norisk_mat)
infs_norisk_dat <- mutate(infs_norisk_dat, index=row_number())
infs_norisk_dat <- merge(adm.names,infs_norisk_dat, by="index")

infs_infrisk_dat <- as.data.frame(infs_infrisk_mat)
infs_infrisk_dat <- mutate(infs_infrisk_dat, index=row_number())
infs_infrisk_dat <- merge(adm.names,infs_infrisk_dat, by="index")

infs_lfrisk_dat <- as.data.frame(infs_lfrisk_mat)
infs_lfrisk_dat <- mutate(infs_lfrisk_dat, index=row_number())
infs_lfrisk_dat <- merge(adm.names,infs_lfrisk_dat, by="index")

infs_norisk_rev_dat <- as.data.frame(infs_norisk_rev_mat)
infs_norisk_rev_dat <- mutate(infs_norisk_rev_dat, index=row_number())
infs_norisk_rev_dat <- merge(adm.names,infs_norisk_rev_dat, by="index")
infs_lfrisk_rev_dat <- as.data.frame(infs_lfrisk_rev_mat)
infs_lfrisk_rev_dat <- mutate(infs_lfrisk_rev_dat, index=row_number())
infs_lfrisk_rev_dat <- merge(adm.names,infs_lfrisk_rev_dat, by="index")
infs_infrisk_rev_dat <- as.data.frame(infs_infrisk_rev_mat)
infs_infrisk_rev_dat <- mutate(infs_infrisk_rev_dat, index=row_number())
infs_infrisk_rev_dat <- merge(adm.names,infs_infrisk_rev_dat, by="index")

#18+
infs_norisk_dat_18 <- as.data.frame(infs_norisk_mat_18)
infs_norisk_dat_18 <- mutate(infs_norisk_dat_18, index=row_number())
infs_norisk_dat_18 <- merge(adm.names,infs_norisk_dat_18, by="index")

infs_infrisk_dat_18 <- as.data.frame(infs_infrisk_mat_18)
infs_infrisk_dat_18 <- mutate(infs_infrisk_dat_18, index=row_number())
infs_infrisk_dat_18 <- merge(adm.names,infs_infrisk_dat_18, by="index")

infs_lfrisk_dat_18 <- as.data.frame(infs_lfrisk_mat_18)
infs_lfrisk_dat_18 <- mutate(infs_lfrisk_dat_18, index=row_number())
infs_lfrisk_dat_18 <- merge(adm.names,infs_lfrisk_dat_18, by="index")

infs_norisk_rev_dat_18 <- as.data.frame(infs_norisk_rev_mat_18)
infs_norisk_rev_dat_18 <- mutate(infs_norisk_rev_dat_18, index=row_number())
infs_norisk_rev_dat_18 <- merge(adm.names,infs_norisk_rev_dat_18, by="index")
infs_lfrisk_rev_dat_18 <- as.data.frame(infs_lfrisk_rev_mat_18)
infs_lfrisk_rev_dat_18 <- mutate(infs_lfrisk_rev_dat_18, index=row_number())
infs_lfrisk_rev_dat_18 <- merge(adm.names,infs_lfrisk_rev_dat_18, by="index")
infs_infrisk_rev_dat_18 <- as.data.frame(infs_infrisk_rev_mat_18)
infs_infrisk_rev_dat_18 <- mutate(infs_infrisk_rev_dat_18, index=row_number())
infs_infrisk_rev_dat_18 <- merge(adm.names,infs_infrisk_rev_dat_18, by="index")

#12-17
infs_norisk_dat_17 <- as.data.frame(infs_norisk_mat_17)
infs_norisk_dat_17 <- mutate(infs_norisk_dat_17, index=row_number())
infs_norisk_dat_17 <- merge(adm.names,infs_norisk_dat_17, by="index")

infs_infrisk_dat_17 <- as.data.frame(infs_infrisk_mat_17)
infs_infrisk_dat_17 <- mutate(infs_infrisk_dat_17, index=row_number())
infs_infrisk_dat_17 <- merge(adm.names,infs_infrisk_dat_17, by="index")

infs_lfrisk_dat_17 <- as.data.frame(infs_lfrisk_mat_17)
infs_lfrisk_dat_17 <- mutate(infs_lfrisk_dat_17, index=row_number())
infs_lfrisk_dat_17 <- merge(adm.names,infs_lfrisk_dat_17, by="index")

infs_norisk_rev_dat_17 <- as.data.frame(infs_norisk_rev_mat_17)
infs_norisk_rev_dat_17 <- mutate(infs_norisk_rev_dat_17, index=row_number())
infs_norisk_rev_dat_17 <- merge(adm.names,infs_norisk_rev_dat_17, by="index")
infs_lfrisk_rev_dat_17 <- as.data.frame(infs_lfrisk_rev_mat_17)
infs_lfrisk_rev_dat_17 <- mutate(infs_lfrisk_rev_dat_17, index=row_number())
infs_lfrisk_rev_dat_17 <- merge(adm.names,infs_lfrisk_rev_dat_17, by="index")
infs_infrisk_rev_dat_17 <- as.data.frame(infs_infrisk_rev_mat_17)
infs_infrisk_rev_dat_17 <- mutate(infs_infrisk_rev_dat_17, index=row_number())
infs_infrisk_rev_dat_17 <- merge(adm.names,infs_infrisk_rev_dat_17, by="index")

#<12
infs_norisk_dat_11 <- as.data.frame(infs_norisk_mat_11)
infs_norisk_dat_11 <- mutate(infs_norisk_dat_11, index=row_number())
infs_norisk_dat_11 <- merge(adm.names,infs_norisk_dat_11, by="index")

infs_infrisk_dat_11 <- as.data.frame(infs_infrisk_mat_11)
infs_infrisk_dat_11 <- mutate(infs_infrisk_dat_11, index=row_number())
infs_infrisk_dat_11 <- merge(adm.names,infs_infrisk_dat_11, by="index")

infs_lfrisk_dat_11 <- as.data.frame(infs_lfrisk_mat_11)
infs_lfrisk_dat_11 <- mutate(infs_lfrisk_dat_11, index=row_number())
infs_lfrisk_dat_11 <- merge(adm.names,infs_lfrisk_dat_11, by="index")

infs_norisk_rev_dat_11 <- as.data.frame(infs_norisk_rev_mat_11)
infs_norisk_rev_dat_11 <- mutate(infs_norisk_rev_dat_11, index=row_number())
infs_norisk_rev_dat_11 <- merge(adm.names,infs_norisk_rev_dat_11, by="index")
infs_lfrisk_rev_dat_11 <- as.data.frame(infs_lfrisk_rev_mat_11)
infs_lfrisk_rev_dat_11 <- mutate(infs_lfrisk_rev_dat_11, index=row_number())
infs_lfrisk_rev_dat_11 <- merge(adm.names,infs_lfrisk_rev_dat_11, by="index")
infs_infrisk_rev_dat_11 <- as.data.frame(infs_infrisk_rev_mat_11)
infs_infrisk_rev_dat_11 <- mutate(infs_infrisk_rev_dat_11, index=row_number())
infs_infrisk_rev_dat_11 <- merge(adm.names,infs_infrisk_rev_dat_11, by="index")

##INF columns
inf.ind=grep("V",colnames(infs_norisk_dat))
infs_norisk_dat$infections_mn=apply(infs_norisk_dat[,inf.ind],1,mean)
infs_norisk_dat$infections_md=apply(infs_norisk_dat[,inf.ind],1,median)
infs_norisk_dat$infections_lo=apply(infs_norisk_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_dat$infections_hi=apply(infs_norisk_dat[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_dat$infections_mn=apply(infs_infrisk_dat[,inf.ind],1,mean)
infs_infrisk_dat$infections_md=apply(infs_infrisk_dat[,inf.ind],1,median)
infs_infrisk_dat$infections_lo=apply(infs_infrisk_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_dat$infections_hi=apply(infs_infrisk_dat[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_dat$infections_mn=apply(infs_lfrisk_dat[,inf.ind],1,mean)
infs_lfrisk_dat$infections_md=apply(infs_lfrisk_dat[,inf.ind],1,median)
infs_lfrisk_dat$infections_lo=apply(infs_lfrisk_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_dat$infections_hi=apply(infs_lfrisk_dat[,inf.ind],1,function(x) quantile(x,.975))

infs_norisk_rev_dat$infections_mn=apply(infs_norisk_rev_dat[,inf.ind],1,mean)
infs_norisk_rev_dat$infections_md=apply(infs_norisk_rev_dat[,inf.ind],1,median)
infs_norisk_rev_dat$infections_lo=apply(infs_norisk_rev_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_rev_dat$infections_hi=apply(infs_norisk_rev_dat[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_rev_dat$infections_mn=apply(infs_infrisk_rev_dat[,inf.ind],1,mean)
infs_infrisk_rev_dat$infections_md=apply(infs_infrisk_rev_dat[,inf.ind],1,median)
infs_infrisk_rev_dat$infections_lo=apply(infs_infrisk_rev_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_rev_dat$infections_hi=apply(infs_infrisk_rev_dat[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_rev_dat$infections_mn=apply(infs_lfrisk_rev_dat[,inf.ind],1,mean)
infs_lfrisk_rev_dat$infections_md=apply(infs_lfrisk_rev_dat[,inf.ind],1,median)
infs_lfrisk_rev_dat$infections_lo=apply(infs_lfrisk_rev_dat[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_rev_dat$infections_hi=apply(infs_lfrisk_rev_dat[,inf.ind],1,function(x) quantile(x,.975))


infs_norisk_dat$inf_ar_mn=infs_norisk_dat$infections_mn/infs_norisk_dat$total_pop
infs_norisk_dat$inf_ar_md=infs_norisk_dat$infections_md/infs_norisk_dat$total_pop
infs_norisk_dat$inf_ar_lo=infs_norisk_dat$infections_lo/infs_norisk_dat$total_pop
infs_norisk_dat$inf_ar_hi=infs_norisk_dat$infections_hi/infs_norisk_dat$total_pop
infs_infrisk_dat$inf_ar_mn=infs_infrisk_dat$infections_mn/infs_infrisk_dat$total_pop
infs_infrisk_dat$inf_ar_md=infs_infrisk_dat$infections_md/infs_infrisk_dat$total_pop
infs_infrisk_dat$inf_ar_lo=infs_infrisk_dat$infections_lo/infs_infrisk_dat$total_pop
infs_infrisk_dat$inf_ar_hi=infs_infrisk_dat$infections_hi/infs_infrisk_dat$total_pop
infs_lfrisk_dat$inf_ar_mn=infs_lfrisk_dat$infections_mn/infs_lfrisk_dat$total_pop
infs_lfrisk_dat$inf_ar_md=infs_lfrisk_dat$infections_md/infs_lfrisk_dat$total_pop
infs_lfrisk_dat$inf_ar_lo=infs_lfrisk_dat$infections_lo/infs_lfrisk_dat$total_pop
infs_lfrisk_dat$inf_ar_hi=infs_lfrisk_dat$infections_hi/infs_lfrisk_dat$total_pop

infs_norisk_rev_dat$inf_ar_mn=infs_norisk_rev_dat$infections_mn/infs_norisk_rev_dat$total_pop
infs_norisk_rev_dat$inf_ar_md=infs_norisk_rev_dat$infections_md/infs_norisk_rev_dat$total_pop
infs_norisk_rev_dat$inf_ar_lo=infs_norisk_rev_dat$infections_lo/infs_norisk_rev_dat$total_pop
infs_norisk_rev_dat$inf_ar_hi=infs_norisk_rev_dat$infections_hi/infs_norisk_rev_dat$total_pop
infs_infrisk_rev_dat$inf_ar_mn=infs_infrisk_rev_dat$infections_mn/infs_infrisk_rev_dat$total_pop
infs_infrisk_rev_dat$inf_ar_md=infs_infrisk_rev_dat$infections_md/infs_infrisk_rev_dat$total_pop
infs_infrisk_rev_dat$inf_ar_lo=infs_infrisk_rev_dat$infections_lo/infs_infrisk_rev_dat$total_pop
infs_infrisk_rev_dat$inf_ar_hi=infs_infrisk_rev_dat$infections_hi/infs_infrisk_rev_dat$total_pop
infs_lfrisk_rev_dat$inf_ar_mn=infs_lfrisk_rev_dat$infections_mn/infs_lfrisk_rev_dat$total_pop
infs_lfrisk_rev_dat$inf_ar_md=infs_lfrisk_rev_dat$infections_md/infs_lfrisk_rev_dat$total_pop
infs_lfrisk_rev_dat$inf_ar_lo=infs_lfrisk_rev_dat$infections_lo/infs_lfrisk_rev_dat$total_pop
infs_lfrisk_rev_dat$inf_ar_hi=infs_lfrisk_rev_dat$infections_hi/infs_lfrisk_rev_dat$total_pop

infs_norisk_dat$uidF=factor(infs_norisk_dat$uid,levels=infs_norisk_dat$uid[order(infs_norisk_dat$inf_ar_md)])
infs_lfrisk_dat$uidF=factor(infs_lfrisk_dat$uid,levels=infs_norisk_dat$uid[order(infs_lfrisk_dat$inf_ar_md)])
infs_infrisk_dat$uidF=factor(infs_infrisk_dat$uid,levels=infs_norisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)])

infs_norisk_rev_dat$uidF=factor(infs_norisk_rev_dat$uid,levels=infs_norisk_rev_dat$uid[order(infs_norisk_rev_dat$inf_ar_md)])
infs_lfrisk_rev_dat$uidF=factor(infs_lfrisk_rev_dat$uid,levels=infs_norisk_rev_dat$uid[order(infs_lfrisk_rev_dat$inf_ar_md)])
infs_infrisk_rev_dat$uidF=factor(infs_infrisk_rev_dat$uid,levels=infs_norisk_rev_dat$uid[order(infs_infrisk_rev_dat$inf_ar_md)])

#18+
infs_norisk_dat_18$infections_mn=apply(infs_norisk_dat_18[,inf.ind],1,mean)
infs_norisk_dat_18$infections_md=apply(infs_norisk_dat_18[,inf.ind],1,median)
infs_norisk_dat_18$infections_lo=apply(infs_norisk_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_dat_18$infections_hi=apply(infs_norisk_dat_18[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_dat_18$infections_mn=apply(infs_infrisk_dat_18[,inf.ind],1,mean)
infs_infrisk_dat_18$infections_md=apply(infs_infrisk_dat_18[,inf.ind],1,median)
infs_infrisk_dat_18$infections_lo=apply(infs_infrisk_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_dat_18$infections_hi=apply(infs_infrisk_dat_18[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_dat_18$infections_mn=apply(infs_lfrisk_dat_18[,inf.ind],1,mean)
infs_lfrisk_dat_18$infections_md=apply(infs_lfrisk_dat_18[,inf.ind],1,median)
infs_lfrisk_dat_18$infections_lo=apply(infs_lfrisk_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_dat_18$infections_hi=apply(infs_lfrisk_dat_18[,inf.ind],1,function(x) quantile(x,.975))

infs_norisk_rev_dat_18$infections_mn=apply(infs_norisk_rev_dat_18[,inf.ind],1,mean)
infs_norisk_rev_dat_18$infections_md=apply(infs_norisk_rev_dat_18[,inf.ind],1,median)
infs_norisk_rev_dat_18$infections_lo=apply(infs_norisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_rev_dat_18$infections_hi=apply(infs_norisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_rev_dat_18$infections_mn=apply(infs_infrisk_rev_dat_18[,inf.ind],1,mean)
infs_infrisk_rev_dat_18$infections_md=apply(infs_infrisk_rev_dat_18[,inf.ind],1,median)
infs_infrisk_rev_dat_18$infections_lo=apply(infs_infrisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_rev_dat_18$infections_hi=apply(infs_infrisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_rev_dat_18$infections_mn=apply(infs_lfrisk_rev_dat_18[,inf.ind],1,mean)
infs_lfrisk_rev_dat_18$infections_md=apply(infs_lfrisk_rev_dat_18[,inf.ind],1,median)
infs_lfrisk_rev_dat_18$infections_lo=apply(infs_lfrisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_rev_dat_18$infections_hi=apply(infs_lfrisk_rev_dat_18[,inf.ind],1,function(x) quantile(x,.975))


infs_norisk_dat_18$inf_ar_mn=infs_norisk_dat_18$infections_mn/infs_norisk_dat_18$pop_18plus
infs_norisk_dat_18$inf_ar_md=infs_norisk_dat_18$infections_md/infs_norisk_dat_18$pop_18plus
infs_norisk_dat_18$inf_ar_lo=infs_norisk_dat_18$infections_lo/infs_norisk_dat_18$pop_18plus
infs_norisk_dat_18$inf_ar_hi=infs_norisk_dat_18$infections_hi/infs_norisk_dat_18$pop_18plus
infs_infrisk_dat_18$inf_ar_mn=infs_infrisk_dat_18$infections_mn/infs_infrisk_dat_18$pop_18plus
infs_infrisk_dat_18$inf_ar_md=infs_infrisk_dat_18$infections_md/infs_infrisk_dat_18$pop_18plus
infs_infrisk_dat_18$inf_ar_lo=infs_infrisk_dat_18$infections_lo/infs_infrisk_dat_18$pop_18plus
infs_infrisk_dat_18$inf_ar_hi=infs_infrisk_dat_18$infections_hi/infs_infrisk_dat_18$pop_18plus
infs_lfrisk_dat_18$inf_ar_mn=infs_lfrisk_dat_18$infections_mn/infs_lfrisk_dat_18$pop_18plus
infs_lfrisk_dat_18$inf_ar_md=infs_lfrisk_dat_18$infections_md/infs_lfrisk_dat_18$pop_18plus
infs_lfrisk_dat_18$inf_ar_lo=infs_lfrisk_dat_18$infections_lo/infs_lfrisk_dat_18$pop_18plus
infs_lfrisk_dat_18$inf_ar_hi=infs_lfrisk_dat_18$infections_hi/infs_lfrisk_dat_18$pop_18plus

infs_norisk_rev_dat_18$inf_ar_mn=infs_norisk_rev_dat_18$infections_mn/infs_norisk_rev_dat_18$pop_18plus
infs_norisk_rev_dat_18$inf_ar_md=infs_norisk_rev_dat_18$infections_md/infs_norisk_rev_dat_18$pop_18plus
infs_norisk_rev_dat_18$inf_ar_lo=infs_norisk_rev_dat_18$infections_lo/infs_norisk_rev_dat_18$pop_18plus
infs_norisk_rev_dat_18$inf_ar_hi=infs_norisk_rev_dat_18$infections_hi/infs_norisk_rev_dat_18$pop_18plus
infs_infrisk_rev_dat_18$inf_ar_mn=infs_infrisk_rev_dat_18$infections_mn/infs_infrisk_rev_dat_18$pop_18plus
infs_infrisk_rev_dat_18$inf_ar_md=infs_infrisk_rev_dat_18$infections_md/infs_infrisk_rev_dat_18$pop_18plus
infs_infrisk_rev_dat_18$inf_ar_lo=infs_infrisk_rev_dat_18$infections_lo/infs_infrisk_rev_dat_18$pop_18plus
infs_infrisk_rev_dat_18$inf_ar_hi=infs_infrisk_rev_dat_18$infections_hi/infs_infrisk_rev_dat_18$pop_18plus
infs_lfrisk_rev_dat_18$inf_ar_mn=infs_lfrisk_rev_dat_18$infections_mn/infs_lfrisk_rev_dat_18$pop_18plus
infs_lfrisk_rev_dat_18$inf_ar_md=infs_lfrisk_rev_dat_18$infections_md/infs_lfrisk_rev_dat_18$pop_18plus
infs_lfrisk_rev_dat_18$inf_ar_lo=infs_lfrisk_rev_dat_18$infections_lo/infs_lfrisk_rev_dat_18$pop_18plus
infs_lfrisk_rev_dat_18$inf_ar_hi=infs_lfrisk_rev_dat_18$infections_hi/infs_lfrisk_rev_dat_18$pop_18plus

infs_norisk_dat_18$uidF=factor(infs_norisk_dat_18$uid,levels=infs_norisk_dat_18$uid[order(infs_norisk_dat_18$inf_ar_md)])
infs_lfrisk_dat_18$uidF=factor(infs_lfrisk_dat_18$uid,levels=infs_norisk_dat_18$uid[order(infs_lfrisk_dat_18$inf_ar_md)])
infs_infrisk_dat_18$uidF=factor(infs_infrisk_dat_18$uid,levels=infs_norisk_dat_18$uid[order(infs_infrisk_dat_18$inf_ar_md)])

infs_norisk_rev_dat_18$uidF=factor(infs_norisk_rev_dat_18$uid,levels=infs_norisk_rev_dat_18$uid[order(infs_norisk_rev_dat_18$inf_ar_md)])
infs_lfrisk_rev_dat_18$uidF=factor(infs_lfrisk_rev_dat_18$uid,levels=infs_norisk_rev_dat_18$uid[order(infs_lfrisk_rev_dat_18$inf_ar_md)])
infs_infrisk_rev_dat_18$uidF=factor(infs_infrisk_rev_dat_18$uid,levels=infs_norisk_rev_dat_18$uid[order(infs_infrisk_rev_dat_18$inf_ar_md)])

#12-17
infs_norisk_dat_17$infections_mn=apply(infs_norisk_dat_17[,inf.ind],1,mean)
infs_norisk_dat_17$infections_md=apply(infs_norisk_dat_17[,inf.ind],1,median)
infs_norisk_dat_17$infections_lo=apply(infs_norisk_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_dat_17$infections_hi=apply(infs_norisk_dat_17[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_dat_17$infections_mn=apply(infs_infrisk_dat_17[,inf.ind],1,mean)
infs_infrisk_dat_17$infections_md=apply(infs_infrisk_dat_17[,inf.ind],1,median)
infs_infrisk_dat_17$infections_lo=apply(infs_infrisk_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_dat_17$infections_hi=apply(infs_infrisk_dat_17[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_dat_17$infections_mn=apply(infs_lfrisk_dat_17[,inf.ind],1,mean)
infs_lfrisk_dat_17$infections_md=apply(infs_lfrisk_dat_17[,inf.ind],1,median)
infs_lfrisk_dat_17$infections_lo=apply(infs_lfrisk_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_dat_17$infections_hi=apply(infs_lfrisk_dat_17[,inf.ind],1,function(x) quantile(x,.975))

infs_norisk_rev_dat_17$infections_mn=apply(infs_norisk_rev_dat_17[,inf.ind],1,mean)
infs_norisk_rev_dat_17$infections_md=apply(infs_norisk_rev_dat_17[,inf.ind],1,median)
infs_norisk_rev_dat_17$infections_lo=apply(infs_norisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_rev_dat_17$infections_hi=apply(infs_norisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_rev_dat_17$infections_mn=apply(infs_infrisk_rev_dat_17[,inf.ind],1,mean)
infs_infrisk_rev_dat_17$infections_md=apply(infs_infrisk_rev_dat_17[,inf.ind],1,median)
infs_infrisk_rev_dat_17$infections_lo=apply(infs_infrisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_rev_dat_17$infections_hi=apply(infs_infrisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_rev_dat_17$infections_mn=apply(infs_lfrisk_rev_dat_17[,inf.ind],1,mean)
infs_lfrisk_rev_dat_17$infections_md=apply(infs_lfrisk_rev_dat_17[,inf.ind],1,median)
infs_lfrisk_rev_dat_17$infections_lo=apply(infs_lfrisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_rev_dat_17$infections_hi=apply(infs_lfrisk_rev_dat_17[,inf.ind],1,function(x) quantile(x,.975))


infs_norisk_dat_17$inf_ar_mn=infs_norisk_dat_17$infections_mn/infs_norisk_dat_17$pop_12_17
infs_norisk_dat_17$inf_ar_md=infs_norisk_dat_17$infections_md/infs_norisk_dat_17$pop_12_17
infs_norisk_dat_17$inf_ar_lo=infs_norisk_dat_17$infections_lo/infs_norisk_dat_17$pop_12_17
infs_norisk_dat_17$inf_ar_hi=infs_norisk_dat_17$infections_hi/infs_norisk_dat_17$pop_12_17
infs_infrisk_dat_17$inf_ar_mn=infs_infrisk_dat_17$infections_mn/infs_infrisk_dat_17$pop_12_17
infs_infrisk_dat_17$inf_ar_md=infs_infrisk_dat_17$infections_md/infs_infrisk_dat_17$pop_12_17
infs_infrisk_dat_17$inf_ar_lo=infs_infrisk_dat_17$infections_lo/infs_infrisk_dat_17$pop_12_17
infs_infrisk_dat_17$inf_ar_hi=infs_infrisk_dat_17$infections_hi/infs_infrisk_dat_17$pop_12_17
infs_lfrisk_dat_17$inf_ar_mn=infs_lfrisk_dat_17$infections_mn/infs_lfrisk_dat_17$pop_12_17
infs_lfrisk_dat_17$inf_ar_md=infs_lfrisk_dat_17$infections_md/infs_lfrisk_dat_17$pop_12_17
infs_lfrisk_dat_17$inf_ar_lo=infs_lfrisk_dat_17$infections_lo/infs_lfrisk_dat_17$pop_12_17
infs_lfrisk_dat_17$inf_ar_hi=infs_lfrisk_dat_17$infections_hi/infs_lfrisk_dat_17$pop_12_17

infs_norisk_rev_dat_17$inf_ar_mn=infs_norisk_rev_dat_17$infections_mn/infs_norisk_rev_dat_17$pop_12_17
infs_norisk_rev_dat_17$inf_ar_md=infs_norisk_rev_dat_17$infections_md/infs_norisk_rev_dat_17$pop_12_17
infs_norisk_rev_dat_17$inf_ar_lo=infs_norisk_rev_dat_17$infections_lo/infs_norisk_rev_dat_17$pop_12_17
infs_norisk_rev_dat_17$inf_ar_hi=infs_norisk_rev_dat_17$infections_hi/infs_norisk_rev_dat_17$pop_12_17
infs_infrisk_rev_dat_17$inf_ar_mn=infs_infrisk_rev_dat_17$infections_mn/infs_infrisk_rev_dat_17$pop_12_17
infs_infrisk_rev_dat_17$inf_ar_md=infs_infrisk_rev_dat_17$infections_md/infs_infrisk_rev_dat_17$pop_12_17
infs_infrisk_rev_dat_17$inf_ar_lo=infs_infrisk_rev_dat_17$infections_lo/infs_infrisk_rev_dat_17$pop_12_17
infs_infrisk_rev_dat_17$inf_ar_hi=infs_infrisk_rev_dat_17$infections_hi/infs_infrisk_rev_dat_17$pop_12_17
infs_lfrisk_rev_dat_17$inf_ar_mn=infs_lfrisk_rev_dat_17$infections_mn/infs_lfrisk_rev_dat_17$pop_12_17
infs_lfrisk_rev_dat_17$inf_ar_md=infs_lfrisk_rev_dat_17$infections_md/infs_lfrisk_rev_dat_17$pop_12_17
infs_lfrisk_rev_dat_17$inf_ar_lo=infs_lfrisk_rev_dat_17$infections_lo/infs_lfrisk_rev_dat_17$pop_12_17
infs_lfrisk_rev_dat_17$inf_ar_hi=infs_lfrisk_rev_dat_17$infections_hi/infs_lfrisk_rev_dat_17$pop_12_17

infs_norisk_dat_17$uidF=factor(infs_norisk_dat_17$uid,levels=infs_norisk_dat_17$uid[order(infs_norisk_dat_17$inf_ar_md)])
infs_lfrisk_dat_17$uidF=factor(infs_lfrisk_dat_17$uid,levels=infs_norisk_dat_17$uid[order(infs_lfrisk_dat_17$inf_ar_md)])
infs_infrisk_dat_17$uidF=factor(infs_infrisk_dat_17$uid,levels=infs_norisk_dat_17$uid[order(infs_infrisk_dat_17$inf_ar_md)])

infs_norisk_rev_dat_17$uidF=factor(infs_norisk_rev_dat_17$uid,levels=infs_norisk_rev_dat_17$uid[order(infs_norisk_rev_dat_17$inf_ar_md)])
infs_lfrisk_rev_dat_17$uidF=factor(infs_lfrisk_rev_dat_17$uid,levels=infs_norisk_rev_dat_17$uid[order(infs_lfrisk_rev_dat_17$inf_ar_md)])
infs_infrisk_rev_dat_17$uidF=factor(infs_infrisk_rev_dat_17$uid,levels=infs_norisk_rev_dat_17$uid[order(infs_infrisk_rev_dat_17$inf_ar_md)])


#<12
infs_norisk_dat_11$infections_mn=apply(infs_norisk_dat_11[,inf.ind],1,mean)
infs_norisk_dat_11$infections_md=apply(infs_norisk_dat_11[,inf.ind],1,median)
infs_norisk_dat_11$infections_lo=apply(infs_norisk_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_dat_11$infections_hi=apply(infs_norisk_dat_11[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_dat_11$infections_mn=apply(infs_infrisk_dat_11[,inf.ind],1,mean)
infs_infrisk_dat_11$infections_md=apply(infs_infrisk_dat_11[,inf.ind],1,median)
infs_infrisk_dat_11$infections_lo=apply(infs_infrisk_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_dat_11$infections_hi=apply(infs_infrisk_dat_11[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_dat_11$infections_mn=apply(infs_lfrisk_dat_11[,inf.ind],1,mean)
infs_lfrisk_dat_11$infections_md=apply(infs_lfrisk_dat_11[,inf.ind],1,median)
infs_lfrisk_dat_11$infections_lo=apply(infs_lfrisk_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_dat_11$infections_hi=apply(infs_lfrisk_dat_11[,inf.ind],1,function(x) quantile(x,.975))

infs_norisk_rev_dat_11$infections_mn=apply(infs_norisk_rev_dat_11[,inf.ind],1,mean)
infs_norisk_rev_dat_11$infections_md=apply(infs_norisk_rev_dat_11[,inf.ind],1,median)
infs_norisk_rev_dat_11$infections_lo=apply(infs_norisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_norisk_rev_dat_11$infections_hi=apply(infs_norisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.975))
infs_infrisk_rev_dat_11$infections_mn=apply(infs_infrisk_rev_dat_11[,inf.ind],1,mean)
infs_infrisk_rev_dat_11$infections_md=apply(infs_infrisk_rev_dat_11[,inf.ind],1,median)
infs_infrisk_rev_dat_11$infections_lo=apply(infs_infrisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_infrisk_rev_dat_11$infections_hi=apply(infs_infrisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.975))
infs_lfrisk_rev_dat_11$infections_mn=apply(infs_lfrisk_rev_dat_11[,inf.ind],1,mean)
infs_lfrisk_rev_dat_11$infections_md=apply(infs_lfrisk_rev_dat_11[,inf.ind],1,median)
infs_lfrisk_rev_dat_11$infections_lo=apply(infs_lfrisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.025))
infs_lfrisk_rev_dat_11$infections_hi=apply(infs_lfrisk_rev_dat_11[,inf.ind],1,function(x) quantile(x,.975))


infs_norisk_dat_11$inf_ar_mn=infs_norisk_dat_11$infections_mn/infs_norisk_dat_11$pop_u12
infs_norisk_dat_11$inf_ar_md=infs_norisk_dat_11$infections_md/infs_norisk_dat_11$pop_u12
infs_norisk_dat_11$inf_ar_lo=infs_norisk_dat_11$infections_lo/infs_norisk_dat_11$pop_u12
infs_norisk_dat_11$inf_ar_hi=infs_norisk_dat_11$infections_hi/infs_norisk_dat_11$pop_u12
infs_infrisk_dat_11$inf_ar_mn=infs_infrisk_dat_11$infections_mn/infs_infrisk_dat_11$pop_u12
infs_infrisk_dat_11$inf_ar_md=infs_infrisk_dat_11$infections_md/infs_infrisk_dat_11$pop_u12
infs_infrisk_dat_11$inf_ar_lo=infs_infrisk_dat_11$infections_lo/infs_infrisk_dat_11$pop_u12
infs_infrisk_dat_11$inf_ar_hi=infs_infrisk_dat_11$infections_hi/infs_infrisk_dat_11$pop_u12
infs_lfrisk_dat_11$inf_ar_mn=infs_lfrisk_dat_11$infections_mn/infs_lfrisk_dat_11$pop_u12
infs_lfrisk_dat_11$inf_ar_md=infs_lfrisk_dat_11$infections_md/infs_lfrisk_dat_11$pop_u12
infs_lfrisk_dat_11$inf_ar_lo=infs_lfrisk_dat_11$infections_lo/infs_lfrisk_dat_11$pop_u12
infs_lfrisk_dat_11$inf_ar_hi=infs_lfrisk_dat_11$infections_hi/infs_lfrisk_dat_11$pop_u12

infs_norisk_rev_dat_11$inf_ar_mn=infs_norisk_rev_dat_11$infections_mn/infs_norisk_rev_dat_11$pop_u12
infs_norisk_rev_dat_11$inf_ar_md=infs_norisk_rev_dat_11$infections_md/infs_norisk_rev_dat_11$pop_u12
infs_norisk_rev_dat_11$inf_ar_lo=infs_norisk_rev_dat_11$infections_lo/infs_norisk_rev_dat_11$pop_u12
infs_norisk_rev_dat_11$inf_ar_hi=infs_norisk_rev_dat_11$infections_hi/infs_norisk_rev_dat_11$pop_u12
infs_infrisk_rev_dat_11$inf_ar_mn=infs_infrisk_rev_dat_11$infections_mn/infs_infrisk_rev_dat_11$pop_u12
infs_infrisk_rev_dat_11$inf_ar_md=infs_infrisk_rev_dat_11$infections_md/infs_infrisk_rev_dat_11$pop_u12
infs_infrisk_rev_dat_11$inf_ar_lo=infs_infrisk_rev_dat_11$infections_lo/infs_infrisk_rev_dat_11$pop_u12
infs_infrisk_rev_dat_11$inf_ar_hi=infs_infrisk_rev_dat_11$infections_hi/infs_infrisk_rev_dat_11$pop_u12
infs_lfrisk_rev_dat_11$inf_ar_mn=infs_lfrisk_rev_dat_11$infections_mn/infs_lfrisk_rev_dat_11$pop_u12
infs_lfrisk_rev_dat_11$inf_ar_md=infs_lfrisk_rev_dat_11$infections_md/infs_lfrisk_rev_dat_11$pop_u12
infs_lfrisk_rev_dat_11$inf_ar_lo=infs_lfrisk_rev_dat_11$infections_lo/infs_lfrisk_rev_dat_11$pop_u12
infs_lfrisk_rev_dat_11$inf_ar_hi=infs_lfrisk_rev_dat_11$infections_hi/infs_lfrisk_rev_dat_11$pop_u12

infs_norisk_dat_11$uidF=factor(infs_norisk_dat_11$uid,levels=infs_norisk_dat_11$uid[order(infs_norisk_dat_11$inf_ar_md)])
infs_lfrisk_dat_11$uidF=factor(infs_lfrisk_dat_11$uid,levels=infs_norisk_dat_11$uid[order(infs_lfrisk_dat_11$inf_ar_md)])
infs_infrisk_dat_11$uidF=factor(infs_infrisk_dat_11$uid,levels=infs_norisk_dat_11$uid[order(infs_infrisk_dat_11$inf_ar_md)])

infs_norisk_rev_dat_11$uidF=factor(infs_norisk_rev_dat_11$uid,levels=infs_norisk_rev_dat_11$uid[order(infs_norisk_rev_dat_11$inf_ar_md)])
infs_lfrisk_rev_dat_11$uidF=factor(infs_lfrisk_rev_dat_11$uid,levels=infs_norisk_rev_dat_11$uid[order(infs_lfrisk_rev_dat_11$inf_ar_md)])
infs_infrisk_rev_dat_11$uidF=factor(infs_infrisk_rev_dat_11$uid,levels=infs_norisk_rev_dat_11$uid[order(infs_infrisk_rev_dat_11$inf_ar_md)])

###
##Generate tables of LF incidence rates for different sero scenarios
### Multiply by 200 to get rate per 1,000 (assuming 20% have symptomatic inf)
###

##Re-organize different results so they can be combined
if(admin==2){
  infs_infrisk_d = infs_infrisk_dat %>% 
    mutate(RevHi_PosInf_m=inf_ar_md*200,
           RevHi_PosInf_l=inf_ar_lo*200,
           RevHi_PosInf_h=inf_ar_hi*200) %>%
    select(Country,NAME_1,NAME_2,uid,
           RevHi_PosInf_m,RevHi_PosInf_l,RevHi_PosInf_h)
}else{
  infs_infrisk_d = infs_infrisk_dat %>% 
    mutate(RevHi_PosInf_m=inf_ar_md*200,
           RevHi_PosInf_l=inf_ar_lo*200,
           RevHi_PosInf_h=inf_ar_hi*200) %>%
    select(Country,NAME_1,uid,
           RevHi_PosInf_m,RevHi_PosInf_l,RevHi_PosInf_h)  
}


infs_norisk_d = infs_norisk_dat %>% 
  mutate(RevHi_PosNo_m=inf_ar_md*200,
         RevHi_PosNo_l=inf_ar_lo*200,
         RevHi_PosNo_h=inf_ar_hi*200) %>%
  select(Country,uid,
         RevHi_PosNo_m,RevHi_PosNo_l,RevHi_PosNo_h)

infs_lfrisk_d = infs_lfrisk_dat %>% 
  mutate(RevHi_PosLF_m=inf_ar_md*200,
         RevHi_PosLF_l=inf_ar_lo*200,
         RevHi_PosLF_h=inf_ar_hi*200) %>%
  select(Country,uid,
         RevHi_PosLF_m,RevHi_PosLF_l,RevHi_PosLF_h)

infs_infrisk_rev_d = infs_infrisk_rev_dat %>% 
  mutate(RevLo_PosInf_m=inf_ar_md*200,
         RevLo_PosInf_l=inf_ar_lo*200,
         RevLo_PosInf_h=inf_ar_hi*200) %>%
  select(Country,uid,
         RevLo_PosInf_m,RevLo_PosInf_l,RevLo_PosInf_h)

infs_norisk_rev_d = infs_norisk_rev_dat %>% 
  mutate(RevLo_PosNo_m=inf_ar_md*200,
         RevLo_PosNo_l=inf_ar_lo*200,
         RevLo_PosNo_h=inf_ar_hi*200) %>%
  select(Country,uid,
         RevLo_PosNo_m,RevLo_PosNo_l,RevLo_PosNo_h)

infs_lfrisk_rev_d = infs_lfrisk_rev_dat %>% 
  mutate(RevLo_PosLF_m=inf_ar_md*200,
         RevLo_PosLF_l=inf_ar_lo*200,
         RevLo_PosLF_h=inf_ar_hi*200) %>%
  select(Country,uid,
         RevLo_PosLF_m,RevLo_PosLF_l,RevLo_PosLF_h)

incid_dat=infs_norisk_d %>% 
  left_join(infs_lfrisk_d,by=c("Country","uid")) %>%
  left_join(infs_infrisk_d,by=c("Country","uid")) %>%
  left_join(infs_norisk_rev_d,by=c("Country","uid")) %>%
  left_join(infs_lfrisk_rev_d,by=c("Country","uid")) %>%
  left_join(infs_infrisk_rev_d,by=c("Country","uid"))
  
incid_fn=paste0("results/LF_incidence_rates_adm",admin,'_revrate',
                as.character(rev_rate*100),
                ifelse(raw.ind==1,'_raw','_modeled_outpredict'),'.csv')
write_csv(incid_dat,file=incid_fn)

##
##Plot of Age-specific Incid Rates 
##
## Use SeroRev-partial, SeroPos-LFrisk as default?
top10_uid=infs_lfrisk_rev_dat$uid[rev(order(infs_lfrisk_rev_dat$inf_ar_md))][1:10]
age11_df=infs_lfrisk_rev_dat_11 %>% filter(uid %in% top10_uid)
age17_df=infs_lfrisk_rev_dat_17 %>% filter(uid %in% top10_uid)
age18_df=infs_lfrisk_rev_dat_18 %>% filter(uid %in% top10_uid)
ageAll_df=infs_lfrisk_rev_dat %>% filter(uid %in% top10_uid)

age11_l = age11_df %>% 
  pivot_longer(cols=starts_with("V"),names_to = "sim",values_to="infections") %>%
  mutate(incid=round(infections*200/pop_u12,1),age="<12 yrs") #%>% 
  #select(uid,age,sim,incid)
age17_l = age17_df %>% 
  pivot_longer(cols=starts_with("V"),names_to = "sim",values_to="infections") %>%
  mutate(incid=round(infections*200/pop_12_17,1),age="12-17 yrs") #%>% 
  #select(uid,age,sim,incid)
age18_l = age18_df %>% 
  pivot_longer(cols=starts_with("V"),names_to = "sim",values_to="infections") %>%
  mutate(incid=round(infections*200/pop_18plus,1),age="18+ yrs") #%>% 
  #select(uid,age,sim,incid)
ageAll_l = ageAll_df %>% 
  pivot_longer(cols=starts_with("V"),names_to = "sim",values_to="infections") %>%
  mutate(incid=round(infections*200/total_pop,1),age="Overall")

ages_all = bind_rows(ageAll_l,age11_l,age17_l,age18_l)
ages_all$name_pretty=paste0(ages_all$Country," -\n ",ages_all$NAME_2)

age_plot<-ggplot(ages_all,aes(x=name_pretty,y=incid,fill=as.factor(age)))+
  geom_boxplot(alpha=0.3)+
  theme_bw()+xlab("")+
  ylab("Annual Lassa Fever\n incidence rate (per 1,000)")+
  guides(fill=guide_legend(title="Age"))+
  theme(axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        axis.text=element_text(size=11),
        axis.text.x=element_text(angle=40,vjust=1,hjust=1))

ggsave(age_plot,file="plots/age_incid_serorev6_revrisk_LFpos.pdf",
       width=8,height=5,units = "in")

##Plot of all scenarios for 2 highest regions?
####Which is default scenario

##Histogram of interannual variability
##Plot of distribution


##Plot of distribution
pdf(paste('plots/top_LF_rates_adm',admin,'_revrate',as.character(rev_rate*100),"_",ifelse(raw.ind==1,'raw','model_outpredict'),'_interannual_variation.pdf',sep=''),width=8,height=7.5,pointsize=14)
infs_lfrisk_rev_dat$NAME_1[grep("Nz",infs_lfrisk_rev_dat$NAME_1)]="Nzerekore"
distj=tail(order(infs_lfrisk_rev_dat$inf_ar_md),9)
par(mfrow=c(3,3))
for(jj in distj){
  hist(200*as.numeric(infs_lfrisk_rev_dat[jj,inf.ind])/infs_lfrisk_rev_dat$total_pop[jj],breaks=40,
       ylab="",xlab="",#"Annual LF incidence rate (per 1000)",
       xlim=c(0,50),
       main=paste(infs_lfrisk_rev_dat$Country[jj],infs_lfrisk_rev_dat$NAME_2[jj],sep=" - "),
       cex.main=1.5,cex.axis=1.4)
  abline(v=200*infs_lfrisk_rev_dat$inf_ar_md[jj],lty=2,lwd=1.5,col="red")
}
dev.off()







###
### OLD PLOTS
###
#plot distributions of force of infection
jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_modeled_outpredict.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(8,4,4,2))
plot(-100,-100,type='l',xlim=c(1346,nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0, 15),
      #max(infs_infrisk_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nModeled FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

axis(1,at=1346:nrow(rr.foi),labels=infs_norisk_dat$uidF[order(infs_norisk_dat$inf_ar_md)][1346:nrow(rr.foi)],las=2)
#axis(1,at=1346:nrow(rr.foi),labels=infs_infrisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)][1346:nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_norisk_dat$inf_ar_lo[order(infs_norisk_dat$inf_ar_md)]*200,
    rev(infs_norisk_dat$inf_ar_hi[order(infs_norisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.76,0.26,0.2))
lines(infs_norisk_dat$inf_ar_md[order(infs_norisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_lfrisk_dat$inf_ar_lo[order(infs_lfrisk_dat$inf_ar_md)]*200,
    rev(infs_lfrisk_dat$inf_ar_hi[order(infs_lfrisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.26,0.26,0.2))
lines(infs_lfrisk_dat$inf_ar_md[order(infs_lfrisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seropos - no risk","Seropos - infection risk","Seropos - LF risk"),lwd=2,col=c(rgb(0.6,0.76,0.26),rgb(0.96,0.76,0.26),rgb(0.9,0.26,0.26)),bty="n")

dev.off()
#write.csv(spillover.infections, file="results/simulated_spillovers_recent.


# plot distributions of force of infection
jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(6,4,4,2))
plot(-100,-100,type='l',xlim=c((nrow(rr.foi)-29),nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,
           max(infs_infrisk_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=(nrow(rr.foi)-29):nrow(rr.foi),labels=infs_infrisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)][(nrow(rr.foi)-29):nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_norisk_dat$inf_ar_lo[order(infs_norisk_dat$inf_ar_md)]*200,
    rev(infs_norisk_dat$inf_ar_hi[order(infs_norisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.76,0.26,0.2))
lines(infs_norisk_dat$inf_ar_md[order(infs_norisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_lfrisk_dat$inf_ar_lo[order(infs_lfrisk_dat$inf_ar_md)]*200,
    rev(infs_lfrisk_dat$inf_ar_hi[order(infs_lfrisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.26,0.26,0.2))
lines(infs_lfrisk_dat$inf_ar_md[order(infs_lfrisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seropos - no risk","Seropos - LF risk","Seropos - infection risk"),lwd=2,col=c(rgb(0.6,0.76,0.26),rgb(0.6,0.26,0.26),rgb(0.96,0.76,0.26)),bty="n")

dev.off()
#write.csv(spillover.infections, file="results/simulated_spillovers_recent.csv")

jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_raw_serorev_risk.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(6,4,4,2))
plot(-100,-100,type='l',xlim=c((nrow(rr.foi)-29),nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,
            max(infs_infrisk_rev_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=(nrow(rr.foi)-29):nrow(rr.foi),labels=infs_infrisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)][(nrow(rr.foi)-29):nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_rev_dat$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_rev_dat$inf_ar_hi[order(infs_infrisk_rev_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(infs_infrisk_rev_dat$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_norisk_rev_dat$inf_ar_lo[order(infs_norisk_rev_dat$inf_ar_md)]*200,
    rev(infs_norisk_rev_dat$inf_ar_hi[order(infs_norisk_rev_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.76,0.26,0.2))
lines(infs_norisk_rev_dat$inf_ar_md[order(infs_norisk_rev_dat$inf_ar_md)]*200,col=rgb(0.6,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_lfrisk_rev_dat$inf_ar_lo[order(infs_lfrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_lfrisk_rev_dat$inf_ar_hi[order(infs_lfrisk_rev_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.26,0.26,0.2))
lines(infs_lfrisk_rev_dat$inf_ar_md[order(infs_lfrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seropos - no risk","Seropos - infection risk","Seropos - LF risk"),lwd=2,col=c(rgb(0.6,0.76,0.26),rgb(0.6,0.26,0.26),rgb(0.96,0.76,0.26)),bty="n")

dev.off()


jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_',ifelse(raw.ind==1,'raw','modeled_outpredict'),'_serorev_compare.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(6,4,4,2))
plot(-100,-100,type='l',xlim=c((nrow(rr.foi)-29),nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,
            max(infs_infrisk_rev_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=(nrow(rr.foi)-29):nrow(rr.foi),labels=infs_infrisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)][(nrow(rr.foi)-29):nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_rev_dat$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_rev_dat$inf_ar_hi[order(infs_infrisk_rev_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.176,0.826,0.5))
lines(infs_infrisk_rev_dat$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.96,0.176,0.826,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_rev_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.16,0.176,0.926,0.2))
lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.16,0.176,0.926,1))

# polygon(
#   c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
#     infs_lfrisk_rev_dat$inf_ar_lo[order(infs_lfrisk_rev_dat$inf_ar_md)]*200,
#     rev(infs_lfrisk_rev_dat$inf_ar_hi[order(infs_lfrisk_rev_dat$inf_ar_md)]*200)),
#   border=NA,col=rgb(0.6,0.26,0.26,0.2))
# lines(infs_lfrisk_rev_dat$inf_ar_md[order(infs_lfrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seroreverted - full risk","Seroreverted - reduced risk"),lwd=2,col=c(rgb(0.16,0.176,0.926),rgb(0.96,0.176,0.826)),bty="n")

dev.off()

##Plot of distribution
jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),"_",ifelse(raw.ind==1,'raw','model_outpredict'),'_interannual_variation.jpeg',sep=''),width=10,height=6.5,units='in',res=300,pointsize=14)
infs_infrisk_rev_dat$NAME_1[grep("Nz",infs_infrisk_rev_dat$NAME_1)]="Nzerekore"
distj=tail(order(infs_infrisk_rev_dat$inf_ar_md),6)
par(mfrow=c(2,3))
for(jj in distj){
  hist(as.numeric(infs_infrisk_rev_dat[jj,inf.ind])/infs_infrisk_rev_dat$total_pop[jj],breaks=40,xlab="Annual attack rate",xlim=c(0,0.22),
       main=paste(infs_infrisk_rev_dat$Country[jj],infs_infrisk_rev_dat$NAME_1[jj],infs_infrisk_rev_dat$NAME_2[jj],sep=" - "),cex=1.2)
  abline(v=infs_infrisk_rev_dat$inf_ar_md[jj],lty=2,lwd=1.5)
  
}
dev.off()

##Plot of distribution
jpeg(paste('plots/top_LF_rates_adm',admin,'_revrate',as.character(rev_rate*100),"_",ifelse(raw.ind==1,'raw','model_outpredict'),'_interannual_variation.jpeg',sep=''),width=8,height=7.5,units='in',res=300,pointsize=14)
infs_infrisk_rev_dat$NAME_1[grep("Nz",infs_infrisk_rev_dat$NAME_1)]="Nzerekore"
distj=tail(order(infs_infrisk_rev_dat$inf_ar_md),9)
par(mfrow=c(3,3))
for(jj in distj){
  hist(200*as.numeric(infs_infrisk_rev_dat[jj,inf.ind])/infs_infrisk_rev_dat$total_pop[jj],breaks=40,xlab="Annual LF incidence rate (per 1000)",
       xlim=c(0,50),
       main=paste(infs_infrisk_rev_dat$Country[jj],infs_infrisk_rev_dat$NAME_1[jj],infs_infrisk_rev_dat$NAME_2[jj],sep=" - "),cex=1.2)
  abline(v=200*infs_infrisk_rev_dat$inf_ar_md[jj],lty=2,lwd=1.5,col="red")
  
}
dev.off()


##
## ADMIN 1 PLOTS
##
# plot distributions of force of infection
jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_modeled_outpredict.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(8,4,4,2))
plot(-100,-100,type='l',xlim=c(145,nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0, 20),
     #max(infs_infrisk_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm1s sorted by median Attack Rate - Top 20\nModeled FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=145:nrow(rr.foi),labels=paste(infs_infrisk_dat$ISO,infs_infrisk_dat$ADMIN_1,sep=" - ")[order(infs_infrisk_dat$inf_ar_md)][145:nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_norisk_dat$inf_ar_lo[order(infs_norisk_dat$inf_ar_md)]*200,
    rev(infs_norisk_dat$inf_ar_hi[order(infs_norisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.76,0.26,0.2))
lines(infs_norisk_dat$inf_ar_md[order(infs_norisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_lfrisk_dat$inf_ar_lo[order(infs_lfrisk_dat$inf_ar_md)]*200,
    rev(infs_lfrisk_dat$inf_ar_hi[order(infs_lfrisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.26,0.26,0.2))
lines(infs_lfrisk_dat$inf_ar_md[order(infs_lfrisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seropos - no risk","Seropos - infection risk","Seropos - LF risk"),lwd=2,col=c(rgb(0.6,0.76,0.26),rgb(0.96,0.76,0.26),rgb(0.9,0.26,0.26)),bty="n")

dev.off()
#write.csv(spillover.infections, file="results/simulated_spillovers_recent.


# plot distributions of force of infection
jpeg(paste('plots/top_attack_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(8,4,4,2))
plot(-100,-100,type='l',xlim=c(187,nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,15),
     # max(infs_infrisk_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm1s sorted by median Attack Rate - Top 20\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=187:nrow(rr.foi),labels=paste(infs_infrisk_dat$ISO,infs_infrisk_dat$ADMIN_1,sep=" - ")[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_norisk_dat$inf_ar_lo[order(infs_norisk_dat$inf_ar_md)]*200,
    rev(infs_norisk_dat$inf_ar_hi[order(infs_norisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.76,0.26,0.2))
lines(infs_norisk_dat$inf_ar_md[order(infs_norisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.76,0.26,1))

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_lfrisk_dat$inf_ar_lo[order(infs_lfrisk_dat$inf_ar_md)]*200,
    rev(infs_lfrisk_dat$inf_ar_hi[order(infs_lfrisk_dat$inf_ar_md)]*200)),
  border=NA,col=rgb(0.6,0.26,0.26,0.2))
lines(infs_lfrisk_dat$inf_ar_md[order(infs_lfrisk_dat$inf_ar_md)]*200,col=rgb(0.6,0.26,0.26,1))
abline(h=10,lty=2)
legend("topleft",c("Seropos - no risk","Seropos - infection risk","Seropos - LF risk"),lwd=2,col=c(rgb(0.6,0.76,0.26),rgb(0.6,0.26,0.26),rgb(0.96,0.76,0.26)),bty="n")

dev.off()
#write.csv(spillover.infections, file="results/simulated_spillovers_recent.csv")


# plot age-specific LF of force of infection
jpeg(paste('plots/top_age-specific_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(6,4,4,2))
plot(-100,-100,type='l',xlim=c((nrow(rr.foi)-29),nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,
            max(infs_infrisk_dat$inf_ar_hi)*200+1),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=(nrow(rr.foi)-29):nrow(rr.foi),labels=infs_infrisk_dat$uid[order(infs_infrisk_dat$inf_ar_md)][(nrow(rr.foi)-29):nrow(rr.foi)],las=2)

# polygon(
#   c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
#     infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
#     rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
#   border=NA,col=rgb(0.96,0.76,0.26,0.5))
# lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1),lwd=2)



polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat_11$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat_11$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.16,0.16,0.906,0.2))
lines(infs_lfrisk_dat_11$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.16,0.16,0.906,1),lwd=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat_17$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat_17$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.76,0.26,0.26,0.2))
lines(infs_infrisk_dat_17$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.76,0.26,0.26,1),lwd=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_dat_18$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
    rev(infs_infrisk_dat_18$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.16,0.76,0.06,0.2))
lines(infs_infrisk_dat_18$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.16,0.76,0.06,1),lwd=2)

abline(h=10,lty=2)
legend("topleft",c("Age 18+","12-17 years","<12 years"),lwd=2,col=c(rgb(0.16,0.76,0.06,1),rgb(0.76,0.26,0.26,1),rgb(0.16,0.16,0.906)),bty="n")

dev.off()


# plot age-specific LF of force of infection
jpeg(paste('plots/top_age-specific_rates_adm',admin,'_revrate',as.character(rev_rate*100),'_raw_serorev_risk.jpeg',sep=''),width=8,height=6.5,units='in',res=300)
par(mar=c(6,4,4,2))
plot(-100,-100,type='l',xlim=c((nrow(rr.foi)-29),nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(0,
            max(infs_infrisk_rev_dat$inf_ar_hi)*200+2),
     xlab='',ylab=' Annual Lassa Fever Attack Rate (per 1000)',xaxs='i',yaxs='i')
mtext(paste0('Adm2s sorted by median Attack Rate - Top 30\nProjected FOI (seroreversion = ',rev_rate*100,'%/yr)'),3,line=0.75)

#axis(1,at=187:nrow(rr.foi),labels=infs_infrisk_dat$uidF[order(infs_infrisk_dat$inf_ar_md)][187:nrow(rr.foi)],las=2)
axis(1,at=(nrow(rr.foi)-29):nrow(rr.foi),labels=infs_infrisk_rev_dat$uid[order(infs_infrisk_dat$inf_ar_md)][(nrow(rr.foi)-29):nrow(rr.foi)],las=2)

# polygon(
#   c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
#     infs_infrisk_dat$inf_ar_lo[order(infs_infrisk_dat$inf_ar_md)]*200,
#     rev(infs_infrisk_dat$inf_ar_hi[order(infs_infrisk_dat$inf_ar_md)])*200),
#   border=NA,col=rgb(0.96,0.76,0.26,0.5))
# lines(infs_infrisk_dat$inf_ar_md[order(infs_infrisk_dat$inf_ar_md)]*200,col=rgb(0.96,0.76,0.26,1),lwd=2)



polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_rev_dat_11$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_rev_dat_11$inf_ar_hi[order(infs_infrisk_rev_rev_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.16,0.16,0.906,0.2))
lines(infs_infrisk_rev_dat_11$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.16,0.16,0.906,1),lwd=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_rev_dat_17$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_rev_dat_17$inf_ar_hi[order(infs_infrisk_rev_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.76,0.26,0.26,0.2))
lines(infs_infrisk_rev_dat_17$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.76,0.26,0.26,1),lwd=2)

polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    infs_infrisk_rev_dat_18$inf_ar_lo[order(infs_infrisk_rev_dat$inf_ar_md)]*200,
    rev(infs_infrisk_rev_dat_18$inf_ar_hi[order(infs_infrisk_rev_dat$inf_ar_md)])*200),
  border=NA,col=rgb(0.16,0.76,0.06,0.2))
lines(infs_infrisk_rev_dat_18$inf_ar_md[order(infs_infrisk_rev_dat$inf_ar_md)]*200,col=rgb(0.16,0.76,0.06,1),lwd=2)

abline(h=10,lty=2)
legend("topleft",c("Age 18+","12-17 years","<12 years"),lwd=2,col=c(rgb(0.16,0.76,0.06,1),rgb(0.76,0.26,0.26,1),rgb(0.16,0.16,0.906)),bty="n")

dev.off()