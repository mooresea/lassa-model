# specify which scenario is being run

##Which level to conduct analysis
args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100

which.crossval = as.numeric(args[3])
# scenarios = expand.grid(crossval=1:11,sero=1:8)
# which.sero = scenarios$sero[which.scenario]
# which.crossval = scenarios$crossval[which.scenario]

# load library
library(mgcv)
library(randomForest)
library(gbm)
library(rgdal)
library(spdep)
library(tidyverse)

# allow for multiple threads for GAMs
ctrl = gam.control(nthreads=10)

##Load observation and FOI data - limit model fitting to those with data (Exclude projections only)
obs_dat=read_csv(file=paste0('results/observed_cases_deaths_adm',admin,'.csv'))
foi_dat=read_csv(file=paste0('results/estimated_foi_adm',admin,'.csv'))
obs_uids=unique(c(obs_dat$uid,foi_dat$uid))

# load analysis of serological data
load(paste('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData',sep=''))

# read in force of infection for sites with serological studies
load(paste('results/proportion_by_type_adm',admin,'_revrate',as.character(rev_rate*100),'_country_upd.RData',sep=''))

# read in projected force of infection for all sites
load(paste('results/projected_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_updated2.RData',sep=''))


# load shape file and subset to locations we are modeling
if(admin==2){
  c = read.csv('data/lassa_adm2_all_covariates_upd.csv')
  c$uid = c$GID_2
  c = subset(c, uid %in% rr$uid)
  g = readOGR('../gadm36_lassa/gadm36_2_lassa.shp')
  g@data$uid = g@data$GID_2
}else{
  c = read.csv('data/lassa_adm1_all_covariates_upd.csv')
  c$uid = c$GID_1
  c = subset(c, uid %in% rr$uid)
  g = readOGR('../gadm36_lassa/gadm36_1_lassa.shp')
  g@data$uid = g@data$GID_1  
}

g = g[which(g@data$uid %in% c$uid),]

# # remove adm1s that cause trouble matching data sets
# g = g[-which(g@data$uid == 'GNB_4'),]
# c = c[-which(c$uid == 'GNB_4'),]
# rr.foi = rr.foi[-setdiff(1:nrow(rr.foi),which(row.names(rr.foi) %in% c$uid)),]

# get indicies of variables that have monthly values
ind.prec = which(substr(names(c),0,8)=='prec_pop')
ind.temp = which(substr(names(c),0,9)=='tmean_pop')
ind.ndvi = which(substr(names(c),0,8)=='ndvi_pop')
c = c[,-c(ind.prec,ind.temp,ind.ndvi)]
ind.ndvi = which(substr(names(c),0,4)=='ndvi')
ind.prec = which(substr(names(c),0,4)=='prec')
ind.temp = which(substr(names(c),0,5)=='tmean')

# find principle components that explain 95% of variation
pca.ndvi = prcomp(c[,ind.ndvi],center=T,scale=T)
# PC1-PC2 account for 95% of variation
pca.prec = prcomp(c[,ind.prec],center=T,scale=T)
# PC1-PC3 account for 95% of variation
pca.temp = prcomp(c[,ind.temp],center=T,scale=T)
# PC1-PC3 account for 95% of variation

# replace monthly variables with dominant principal components
c = c[,-c(ind.ndvi,ind.prec,ind.temp)]
c$ndvi1 = pca.ndvi$x[,1]
c$ndvi2 = pca.ndvi$x[,2]
c$prec1 = pca.prec$x[,1]
c$prec2 = pca.prec$x[,2]
c$prec3 = pca.prec$x[,3]
#c$prec4 = pca.prec$x[,4]
c$temp1 = pca.temp$x[,1]
c$temp2 = pca.temp$x[,2]
c$temp3 = pca.temp$x[,3]

# load in IHME's health access and quality index
tmp = read.csv('data/IHME_GBD_2016_HAQ_INDEX_1990_2016_SCALED_CAUSE_VALUES_Y2018M05D23.CSV')
tmp = subset(tmp,ihme_loc_id %in% unique(c$GID_0))
tmp = subset(tmp,indicator_name == 'Healthcare Access and Quality Index')
##GET MEAN VALUE FROM 1990-2016
aggregate(tmp$val,by=list(ihme_loc_id=tmp$ihme_loc_id),FUN=mean)
tmp = aggregate(tmp$val,by=list(ihme_loc_id=tmp$ihme_loc_id),FUN=mean)
names(tmp) = c('GID_0','HAQI')
c = merge(c,tmp)

# get indices of all variables to be used in FOI prediction
ind.ndvi = which(substr(names(c),0,4)=='ndvi')
ind.prec = which(substr(names(c),0,4)=='prec')
ind.temp = which(substr(names(c),0,4)=='temp')
ind.elev = which(names(c)=='elev')
ind.lon = which(names(c)=='lon')
ind.lat = which(names(c)=='lat')
ind.trav = which(names(c)=='travel_time_unweighted')
ind.forest = which(names(c)=='forest_loss')
#ind.primate_occur = which(names(c)=='primate_occur')
#ind.NHP_count = which(names(c)=='NHP_count')
ind.frontier_pct = which(names(c)=='frontier_pct')
ind.trop_pct = which(names(c)=='trop_pct')
ind.haqi = which(names(c)=='HAQI')
ind.mast_o=which(names(c)=="mastomys_occur")
ind.mast_i=which(names(c)=="mastomys_inf")
ind.house=which(names(c)=="imp_housing")
ind.crop=which(names(c)=="crop_pct")
ind.meat=which(names(c)=="bush_meat")
ind.pov=which(names(c)=="poverty_pct")
ind.malaria=which(names(c)=="malaria_pfpr")


# center and scale predictor variables
c[,ind.elev] = scale(c[,ind.elev],center=T,scale=T)
c[,ind.lon] = scale(c[,ind.lon],center=T,scale=T)
c[,ind.lat] = scale(c[,ind.lat],center=T,scale=T)
c[,ind.trav] = scale(c[,ind.trav],center=T,scale=T)
c[,ind.forest] = scale(c[,ind.forest],center=T,scale=T)
#c[,ind.primate_occur] = scale(c[,ind.primate_occur],center=T,scale=T)
#c[,ind.NHP_count] = scale(c[,ind.NHP_count],center=T,scale=T)
c[,ind.frontier_pct] = scale(c[,ind.frontier_pct],center=T,scale=T)
c[,ind.trop_pct] = scale(c[,ind.trop_pct],center=T,scale=T)
c[,ind.mast_o] = scale(c[,ind.mast_o],center=T,scale=T)
c[,ind.mast_i] = scale(c[,ind.mast_i],center=T,scale=T)
c[,ind.crop] = scale(c[,ind.crop],center=T,scale=T)
c[,ind.meat] = scale(c[,ind.meat],center=T,scale=T)
c[,ind.house] = scale(c[,ind.house],center=T,scale=T)
c[,ind.pov] = scale(c[,ind.pov],center=T,scale=T)
c[,ind.malaria] = scale(c[,ind.malaria],center=T,scale=T)

###
## Remove locations that we don't want to model or missing data
## CMR, CAF, TCD, SEN, GMB? parts of NER?
(c$uid[is.na(c$mastomys_occur)])
c=c[c$GID_0!="TCD",]
c=c[!is.na(c$imp_housing),]
c=c[!is.na(c$mastomys_inf),]

#Only needed for GMRF models
# # figure out neighbors between polygons
# c$uid = as.factor(c$uid)
# #g_fit=g[g@data$uid %in% c_fit$uid,]
# g=g[g@data$uid %in% c$uid,]
# nb = poly2nb(g, row.names = g@data$uid)
# names(nb) = g@data$uid

##
## Separate out locations to estimate and locations to predict
##
c_fit=c[which(c$uid %in% obs_uids),]
c_pred=c[which(!c$uid %in% obs_uids),]


rr.foi=rr.foi[which(rr$uid %in% c_fit$uid),]
rr=rr[which(rr$uid %in% c_fit$uid),]
rr.foi=rr.foi[match(c_fit$uid,rownames(rr.foi)),]

# merge force of infection projections into covariate data set
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('
    c_fit$foi_',ii,' = log(rr.foi[,ii],10)',sep='')))
}
c_fit$m_foi=log(apply(rr.foi,1,median),10)
# partition data for cross-validation
sets = list()
#countries = names(sort(table(c_fit$GID_0)))[c(c(15),c(14),c(13),c(12),c(11),c(1,10),c(2,9),c(3,8),c(4,7),c(5,6))]
#reps_ii=sort(c(1:10,6:10))
s_ids=sample(1:nrow(c_fit),size=nrow(c_fit),replace=F)
nrm=round(nrow(c_fit)/10,0)
ii.cnt=1
for(ii in 1:10){
  ii_ids=sort(s_ids[ii.cnt:(ii.cnt+nrm-1)])
  sets[[ii]] = ii_ids #which(!1:nrow(c_fit) %in% ii_ids) #which(c$GID_0 %in% countries[reps_ii==ii])
  ii.cnt=ii.cnt+nrm
}

# define data sets for training and testing
if(which.crossval %in% 1:10){
  c.test = c_fit[sets[[which.crossval]],]
  c_fit = c_fit[-sets[[which.crossval]],]
} else {
  c.test = c_fit
}

##
##Get explanatory power of variables for mean FOI estimate
##
lm.fit=lm('m_foi~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + HAQI + mastomys_occur + mastomys_inf + crop_pct + bush_meat + imp_housing + malaria_pfpr + 1', data=c_fit)

##Too many parameters to have interactions between all of them
# fit linear models with interactions
if(admin==1){
    lm2.fit = lm('m_foi ~ ndvi1 + ndvi2 + ndvi1:temp1 + ndvi1:temp2 + ndvi1:temp3 + ndvi2:temp1 + ndvi2:temp2 + ndvi2:temp3 + prec1  + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + temp1:lon + temp2:lon + temp3:lon + temp1:lat + temp2:lat + temp3:lat + mastomys_occur + mastomys_inf + imp_housing + malaria_pfpr + HAQI + crop_pct + bush_meat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + 1', data=c_fit)
}else{
   lm2.fit = lm('m_foi ~ (ndvi1 + ndvi2 + prec1 + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + mastomys_occur + mastomys_inf + imp_housing + malaria_pfpr) ^ 2 + HAQI + crop_pct + bush_meat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + 1', data=c_fit)
}


rf.fit = randomForest(m_foi ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + HAQI+ mastomys_occur + mastomys_inf + crop_pct + bush_meat + imp_housing + malaria_pfpr, data=c_fit)
brt.fit = gbm(m_foi ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + HAQI+ mastomys_occur + mastomys_inf + crop_pct + bush_meat + imp_housing + malaria_pfpr, data=c_fit, distribution="gaussian")

importance(rf.fit)

write.csv(as.data.frame(importance(rf.fit)),file=paste0("results/RF_cov_importance_adm",admin,".csv"))
write.csv(as.data.frame(relative.influence(brt.fit)),file=paste0("results/BRT_cov_influence_adm",admin,".csv"))
write.csv(as.data.frame(summary(lm2.fit)$coefficients),file=paste0("results/LM_interact_cov_significance_adm",admin,".csv"))






# fit boosted regression trees
brt.pred = matrix(NA,nrow(c.test),ncol(rr.foi))
brt.pred.out = matrix(NA,nrow(c_pred),ncol(rr.foi))
for(ii in 1:ncol(rr.foi)){
  eval(parse(text=paste('tmp = gbm(foi_',ii,' ~ ndvi1 + ndvi2 + prec1 + prec2 + prec3 + temp1 + temp2 + temp3 + elev + lon + lat + travel_time_unweighted + forest_loss + frontier_pct + trop_pct + HAQI+ mastomys_occur + mastomys_inf + crop_pct + bush_meat + imp_housing + malaria_pfpr, data=c_fit, distribution="gaussian")',sep='')))
  brt.pred.out[,ii] = predict(tmp,newdata=c_pred,n.trees=tmp$n.trees)
  
  if(ii %% 10==0){
    print(ii)
    print("BRT")
  }
}

# save outputs to file
save(
  null.pred,lm.pred,lm2.pred,
  mrkLO.pred,mrkHI.pred,
  mrkLOcovs.pred,mrkHIcovs.pred,
  rf.pred,brt.pred,
  c,c.test,sets,which.crossval,
  c_fit,c_pred,
  null.pred.out,lm.pred.out,lm2.pred.out,
  mrkLO.pred.out,mrkHI.pred.out,
  mrkLOcovs.pred.out,mrkHIcovs.pred.out,
  rf.pred.out,brt.pred.out,
  file=paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_',which.crossval,'_outpredict.RData',sep=''))
