# load libraries
library(extraDistr)

args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100
which.scenario = 11

# allocate storage for optimal model weights
model.wts = matrix(0,9,9+1)

### NOTE: This and other triple-commented lines should be run when making
### the figure showing the likelihood associated with each model (Fig. 6).
# # run model marginal log likelihood calculations
# load('../output/model_wts_regression.RData')
# LL.models = matrix(NA,8,9)

# find optimal model weights for each set of serological data assumptions
#for(which.sero in 1:8){
which.sero=1

  # read in force of infection for sites with serological studies
  load(paste('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData',sep=''))
  
  # read in force of infection for sites with serological studies
  load(paste('results/projected_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_updated2.RData',sep=''))
  
  # read in covariates data set
  if(admin==2){
    c = read.csv('data/lassa_adm2_all_covariates_upd.csv')
    c$uid = c$GID_2
  }else{
    c = read.csv('data/lassa_adm1_all_covariates_upd.csv')
    c$uid = c$GID_1    
  }

  c.all = subset(c, uid %in% rr$uid)
  
  load(paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_11_outpredict.RData',sep=''))
  c.all = c_fit #subset(c.all, uid %in% c_fit$uid)
  
  # c.all
  # # merge force of infection projections into covariate data set
  # for(ii in 1:ncol(rr.foi)){
  #   eval(parse(text=paste('
  #     c.all$foi_',ii,' = log(rr.foi[c.all$uid,ii],10)',sep='')))
  # }
  
  # allocate variables to store the mean and standard deviation of regression predictions
  pred.mean = pred.sd = data.frame(SPID = c.all$uid)
  eval(parse(text=paste('pred.mean$null.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$lm.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$lm2.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrkLO.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrkHI.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrkLOcovs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$mrkHIcovs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
    eval(parse(text=paste('pred.mean$rf.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.mean$brt.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$null.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$lm.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$lm2.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrkLO.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrkHI.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrkLOcovs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$mrkHIcovs.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$rf.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  eval(parse(text=paste('pred.sd$brt.pred_',which.sero,' = rep(NA,nrow(c.all))',sep='')))
  test_uids=c()
  # loop across the ten country partitions
  c.test.all=list()
  for(which.partition in 1:10){
  
    # load the predictions of the 10% of adm1s withheld based on a fit to the other 90%
    load(paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_11_outpredict.RData',sep=''))
    
    file.to.load = paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_',which.partition,'_outpredict.RData',sep='')
    #if(file.to.load %in% system('ls results/model_foi_[0-9]*.RData',intern=T)){
      load(file.to.load)
      print(all(c.all$uids %in% c(c_fit$uid,c.test$uid)))
    c.test.all[[which.partition]]=c.test
      test_uids=c(test_uids,as.character(c.test$uid))
      # calculate the mean log10 FOI predicted by each model for each adm1 from the 10%      
      pred.mean[which(c.all$uid %in% c.test$uid),paste('null.pred_',which.sero,sep='')] =
        apply(null.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('lm.pred_',which.sero,sep='')] =
        apply(lm.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('lm2.pred_',which.sero,sep='')] =
        apply(lm2.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('mrkLO.pred_',which.sero,sep='')] =
        apply(mrkLO.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('mrkHI.pred_',which.sero,sep='')] =
        apply(mrkHI.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('mrkLOcovs.pred_',which.sero,sep='')] =
        apply(mrkLOcovs.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('mrkHIcovs.pred_',which.sero,sep='')] =
        apply(mrkHIcovs.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('rf.pred_',which.sero,sep='')] =
        apply(rf.pred,1,mean)
      pred.mean[which(c.all$uid %in% c.test$uid),paste('brt.pred_',which.sero,sep='')] =
        apply(brt.pred,1,mean)
      
      # calculate the s.d. of log10 FOI predicted by each model for each adm1 from the 10%
      pred.sd[which(c.all$uid %in% c.test$uid),paste('null.pred_',which.sero,sep='')] =
        apply(null.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('lm.pred_',which.sero,sep='')] =
        apply(lm.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('lm2.pred_',which.sero,sep='')] =
        apply(lm2.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('mrkLO.pred_',which.sero,sep='')] =
        apply(mrkLO.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('mrkHI.pred_',which.sero,sep='')] =
        apply(mrkHI.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('mrkLOcovs.pred_',which.sero,sep='')] =
        apply(mrkLOcovs.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('mrkHIcovs.pred_',which.sero,sep='')] =
        apply(mrkHIcovs.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('rf.pred_',which.sero,sep='')] =
        apply(rf.pred,1,sd)
      pred.sd[which(c.all$uid %in% c.test$uid),paste('brt.pred_',which.sero,sep='')] =
        apply(brt.pred,1,sd)
    #}
    print(which.partition)
  }

  ##Might need to figure out why certain UIDs were not part of testing when sets includes all vals
  ##Because didn't set.seed - so sample is different each time (would it be anyway?)
   # # remove any admin units where there are NA predictions
  c.all = c.all[-which(is.na(rowSums(pred.mean[,-1]))),]
  pred.mean = pred.mean[-which(is.na(rowSums(pred.mean[,-1]))),]
  pred.sd = pred.sd[-which(is.na(rowSums(pred.sd[,-1]))),]
 
  ##Plot of model validation fits
  rr.foi.pred=rr.foi[rownames(rr.foi) %in% pred.mean$SPID,]
  rr.foi.pred=rr.foi.pred[match(pred.mean$SPID,rownames(rr.foi.pred)),]
  rr.foi.pred=log(rr.foi.pred,10)
  
  # make scatter plots of predicted and directly estimated force of infection
  jpeg(paste('plots/foi_scatter_adm',admin,'_revrate',as.character(rev_rate*100),'_outpredict_validation.jpeg',sep=''),
       width=10,height=7.5,units='in',res=300)
  
  layout(matrix(1:9,3,3,byrow=T))
  par(oma=rep(0,4),mar=c(5,6,1,2))
  
  
  range.x = range(
    pred.mean[,"null.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"null.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"null.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Intercept-only model (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"lm.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"lm.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"lm.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Linear model (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"lm2.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"lm2.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"lm2.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Linear model with interactions (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"mrkLO.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"mrkLO.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"mrkLO.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Markov random field, 10x10 (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"mrkHI.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"mrkHI.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"mrkHI.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Markov random field, 20x20 (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"mrkLOcovs.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"mrkLOcovs.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"mrkLOcovs.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Markov random field, 10x10 + linear (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"mrkHIcovs.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"mrkHIcovs.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"mrkHIcovs.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Markov random field, 20x20 + linear (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"rf.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"rf.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"rf.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Random forest (R^2 = ',R2,')'),3,cex=0.7)
  
  range.x = range(
    pred.mean[,"brt.pred_1"])
  range.y = range(c(
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.025)),
    apply(rr.foi.pred,1,function(ii)quantile(ii,0.975))))
  plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
       xlab=expression('Regression prediction of '*log[10]*' FOI'),
       ylab=expression('Direct estimate of '*log[10]*' FOI '))
  abline(0,1)
  
  points(pred.mean[,"brt.pred_1"],apply(rr.foi.pred,1,median),
         pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
  #mtext('Linear model',3,cex=0.7)
  #mtext('B',3,cex=0.7,at=range.x[1])
  R2 = summary(lm(apply(rr.foi.pred,1,median)~pred.mean[,"brt.pred_1"]))$r.squared
  R2 = round(1000 * R2) / 1000
  #eval(parse(text=paste("
  #text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
  mtext(paste0('Boosted regression trees (R^2 = ',R2,')'),3,cex=0.7)
  
  dev.off()
   
  # negative marginal log likelihood function pooling all predictions of withheld data
  margLik = function(par){
    wts = head(par,-1)
    par = tail(par,1)
    pred.mean.wtd = as.matrix(pred.mean[,-1])%*%matrix(wts,9,1)
    pred.sd.wtd = sqrt(as.matrix(pred.sd[,-1]^2)%*%matrix(wts^2,9,1)) + par
    LL = log(rowMeans(apply(c.all[,tail((1:ncol(c.all)),1000)],2,function(x)
      dnorm(x,pred.mean.wtd,pred.sd.wtd))))
    -sum(ifelse(is.finite(LL),LL,-750))
  }
  
  # find MLE of model weights
  opt = constrOptim(
    theta = c(rep(1/9,9),0.1), f = margLik, method = 'Nelder-Mead',
    ui = rbind(
      c(rep(1,9),0),
      c(rep(-1,9),0),
      diag(10),
      cbind(-diag(9),rep(0,9))),
    ci = c(0.999,-1.001,rep(0,10),rep(-1,9)))
  model.wts[which.sero,] = opt$par

  # save output
  save(model.wts,pred.mean,pred.sd,file=paste0('results/model_wts_regression_adm',admin,'_revrate',as.character(rev_rate*100),'_outpredict.RData'))
  
  ### NOTE: This and other triple-commented lines should be run when making
  ### the figure showing the likelihood associated with each model (Fig. 6).
  # # calculate marginal log likelihood for each model plus ensemble
  # for(which.model in 1:8){
  #   LL.models[which.sero,which.model] = -margLik(c(diag(8)[which.model,],0))
  # }
  # LL.models[which.sero,9] = -margLik(model.wts[which.sero,])




### NOTE: This and other triple-commented lines should be run when making
### the figure showing the likelihood associated with each model (Fig. 6).
# # plot ensemble model weights and model performance in cross-validation
# jpeg('../figures/model_weights.jpeg',width=6.5,height=5,units='in',res=300)
# layout(1:3,heights=c(0.4,0.4,0.2))
# par(oma=rep(0,4),mar=c(1.75,5,1.5,0.25))
# 
# barplot(
#   -LL.models,col=rgb(0,0,1,seq(0,1,length.out=8)),las=1,beside=T,
#   names.arg=c('Int.','Lin.','Lin.+','MRF10','MRF20','MRF10+','RF','BRT','Ens.'))
# mtext('A',3,at=1,line=0,cex=0.8)
# mtext('Negative marginal log likelihood',2,line=3.75,cex=0.8)
# mtext('Regression model',1,line=2.5,cex=0.8)
# 
# barplot(
#   model.wts[,-9],col=rgb(0,0,1,seq(0,1,length.out=8)),las=1,beside=T,ylim=c(0,0.55),
#   names.arg=c('Int.','Lin.','Lin.+','MRF10','MRF20','MRF10+','RF','BRT'))
# mtext('B',3,at=1,line=-1,cex=0.8)
# mtext('Ensemble weight',2,line=3.75,cex=0.8)
# mtext('Regression model',1,line=2.5,cex=0.8)
# 
# plot.new()
# legend('bottom',legend=1:8,fill=rgb(0,0,1,seq(0,1,length.out=8)),bty='n',horiz=T)
# mtext('Serology scenario',1,cex=0.8)
# 
# dev.off()
