# specify which scenario is being run
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
admin=as.numeric(args[1])
rev_rate=as.numeric(args[2])/100
which.scenario = 11

# load force of infection predictions
load(paste('results/foi_from_sero_adm',admin,'_revrate',as.character(rev_rate*100),'.RData',sep=''))
load(paste('results/projected_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_updated2.RData',sep=''))
load(paste('results/model_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_11_outpredict.RData',sep=''))

# read in covariates data set
if(admin==2){
  c.orig = read.csv('data/lassa_adm2_all_covariates_upd.csv')
  c.orig$uid = c.orig$GID_2
}else{
  c.orig = read.csv('data/lassa_adm1_all_covariates_upd.csv')
  c.orig$uid = c.orig$GID_1
}

c.orig = subset(c.orig, uid %in% rr$uid)

# force of infection that is being modeled
#rr.foi = log(rr.foi[-setdiff(1:nrow(rr.foi),which(row.names(rr.foi) %in% c$uid)),],10)
rr.foi = log(rr.foi,10)

##Remove rr and rr.foi that weren't modeled
rr.foi=rr.foi[which(rr$uid %in% c$uid),]
rr=rr[which(rr$uid %in% c$uid),]
rr.foi.all=rr.foi
rr.all=rr
rr.foi=rr.foi[which(rr$uid %in% c_fit$uid),]
rr=rr[which(rr$uid %in% c_fit$uid),]

##Join model fits and predictions
null.pred.all=rbind(null.pred,null.pred.out)
lm.pred.all=rbind(lm.pred,lm.pred.out)
lm2.pred.all=rbind(lm2.pred,lm2.pred.out)
mrkLO.pred.all=rbind(mrkLO.pred,mrkLO.pred.out)
mrkHI.pred.all=rbind(mrkHI.pred,mrkHI.pred.out)
mrkLOcovs.pred.all=rbind(mrkLOcovs.pred,mrkLOcovs.pred.out)
mrkHIcovs.pred.all=rbind(mrkHIcovs.pred,mrkHIcovs.pred.out)
rf.pred.all=rbind(rf.pred,rf.pred.out)
brt.pred.all=rbind(brt.pred,brt.pred.out)

reord.matF<-function(pred.mat,rrf){
  rownames(pred.mat)=c(c_fit$uid)
  pred.mat=pred.mat[match(rrf$uid,rownames(pred.mat)),]
  return(pred.mat)
}
reord.matF2<-function(pred.mat,rrf){
  rownames(pred.mat)=c(c_fit$uid,c_pred$uid)
  pred.mat=pred.mat[match(rrf$uid,rownames(pred.mat)),]
  return(pred.mat)
}

null.pred=reord.matF(null.pred,rr)
lm.pred=reord.matF(lm.pred,rr)
lm2.pred=reord.matF(lm2.pred,rr)
mrkLO.pred=reord.matF(mrkLO.pred,rr)
mrkHI.pred=reord.matF(mrkHI.pred,rr)
mrkLOcovs.pred=reord.matF(mrkLOcovs.pred,rr)
mrkHIcovs.pred=reord.matF(mrkHIcovs.pred,rr)
rf.pred=reord.matF(rf.pred,rr)
brt.pred=reord.matF(brt.pred,rr)

null.pred.all=reord.matF2(null.pred.all,rr.all)
lm.pred.all=reord.matF2(lm.pred.all,rr.all)
lm2.pred.all=reord.matF2(lm2.pred.all,rr.all)
mrkLO.pred.all=reord.matF2(mrkLO.pred.all,rr.all)
mrkHI.pred.all=reord.matF2(mrkHI.pred.all,rr.all)
mrkLOcovs.pred.all=reord.matF2(mrkLOcovs.pred.all,rr.all)
mrkHIcovs.pred.all=reord.matF2(mrkHIcovs.pred.all,rr.all)
rf.pred.all=reord.matF2(rf.pred.all,rr.all)
brt.pred.all=reord.matF2(brt.pred.all,rr.all)

# make scatter plots of predicted and directly estimated force of infection
jpeg(paste('plots/foi_scatter_adm',admin,'_revrate',as.character(rev_rate*100),'_outpredict.jpeg',sep=''),
     width=10,height=7.5,units='in',res=300)

layout(matrix(1:9,3,3,byrow=T))
par(oma=rep(0,4),mar=c(5,6,1,2))

range.x = range(c(
  apply(null.pred,1,function(ii)quantile(ii,0.025)),
  apply(null.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(null.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(null.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(null.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(null.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(null.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
mtext('Intercept only model',3,cex=0.7)
mtext('A',3,cex=0.7,at=range.x[1])

range.x = range(c(
  apply(lm.pred,1,function(ii)quantile(ii,0.025)),
  apply(lm.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(lm.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(lm.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(lm.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(lm.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(lm.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Linear model',3,cex=0.7)
#mtext('B',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(lm.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#text(x=-5.6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Linear model (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(lm2.pred,1,function(ii)quantile(ii,0.025)),
  apply(lm2.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(lm2.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(lm2.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(lm2.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(lm2.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(lm2.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Linear model with interactions',3,cex=0.7)
#mtext('C',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(lm2.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-6.8,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Linear model with interactions (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(mrkLO.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrkLO.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrkLO.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrkLO.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrkLO.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrkLO.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrkLO.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Markov random field, 10x10',3,cex=0.7)
#mtext('D',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrkLO.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#   text(x=-5.4,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Markov random field, 10x10 (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(mrkHI.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrkHI.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrkHI.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrkHI.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrkHI.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrkHI.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrkHI.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Markov random field, 20x20',3,cex=0.7)
#mtext('E',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrkHI.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#   text(x=-5.3,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Markov random field, 20x20 (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(mrkHIcovs.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrkHIcovs.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrkHIcovs.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrkHIcovs.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrkHIcovs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrkHIcovs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrkHIcovs.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)

#mtext('F',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrkHIcovs.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-5.5,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Markov random field, 20x20 + linear (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(mrkLOcovs.pred,1,function(ii)quantile(ii,0.025)),
  apply(mrkLOcovs.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(mrkLOcovs.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(mrkLOcovs.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(mrkLOcovs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(mrkLOcovs.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(mrkLOcovs.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Markov random field, 10x10 + linear',3,cex=0.7)
#mtext('F',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(mrkLOcovs.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-5.5,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Markov random field, 10x10 + linear (R^2 = ',R2,')'),3,cex=0.7)

range.x = range(c(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(rf.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(rf.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)

#mtext('G',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(rf.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Random forest (R^2 = ',R2,")"),3,cex=0.7)

range.x = range(c(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(brt.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(brt.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Boosted regression trees',3,cex=0.7)
#mtext('H',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(brt.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-5.7,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Boosted regression trees (R^2 = ',R2,")"),3,cex=0.7)
dev.off()

######
##Plot of just RF and BRT (2 best fit)
#####
pdf(paste('plots/foi_scatter_adm',admin,'_revrate',as.character(rev_rate*100),'_outpredict_bestfit.pdf',sep=''),
          width=8,height=4,pointsize=14)

par(mfrow=c(1,2))
range.x = range(c(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(rf.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(rf.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rf.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(rf.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)

#mtext('G',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(rf.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-6,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Random forest (R^2 = ',R2,")"),3,cex=0.7)

range.x = range(c(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,function(ii)quantile(ii,0.975))))
range.y = range(c(
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,function(ii)quantile(ii,0.975))))
plot(-100,-100,xlim=range.x,ylim=range.y,las=1,
     xlab=expression('Regression prediction of '*log[10]*' FOI'),
     ylab=expression('Direct estimate of '*log[10]*' FOI '))
abline(0,1)
segments(
  apply(brt.pred,1,function(ii)quantile(ii,0.025)),
  apply(rr.foi,1,median),
  apply(brt.pred,1,function(ii)quantile(ii,0.975)),
  apply(rr.foi,1,median),
  col=rgb(0,0,0,0.1))
segments(
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.025)),
  apply(brt.pred,1,median),
  apply(rr.foi,1,function(ii)quantile(ii,0.975)),
  col=rgb(0,0,0,0.1))
points(apply(brt.pred,1,median),apply(rr.foi,1,median),
       pch=21,col=rgb(0,0,0,0.75),bg=rgb(0,1,0,1),cex=0.7)
#mtext('Boosted regression trees',3,cex=0.7)
#mtext('H',3,cex=0.7,at=range.x[1])
R2 = summary(lm(apply(rr.foi,1,median)~apply(brt.pred,1,median)))$r.squared
R2 = round(1000 * R2) / 1000
#eval(parse(text=paste("
#  text(x=-5.7,y=0,expression(R^2*' = '*",R2,"))",sep='')))
mtext(paste0('Boosted regression trees (R^2 = ',R2,")"),3,cex=0.7)

dev.off()












# range of force of infection values to be plotted
range.foi = range(c(
  c(apply(rr.foi,1,median)),
  c(apply(lm.pred,1,median)),
  c(apply(lm2.pred,1,median)),
  c(apply(brt.pred,1,median)),
  c(apply(rf.pred,1,median)),
  c(apply(mrkLO.pred,1,median)),
  c(apply(mrkHI.pred,1,median)),
  c(apply(mrkHIcovs.pred,1,median)),
  c(apply(mrkLOcovs.pred,1,median))))
if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
  if(round(range.foi[1])==round(range.foi[2])){
    range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
  }else{
    range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
  }
}

# load shape file and subset to locations we are modeling
library(rgdal)
library(spdep)
library(sf)
if(admin==2){
  g = read_sf('../gadm36_lassa/gadm36_2_lassa.shp')
  g$uid = g$GID_2
}else{
  g = read_sf('../gadm36_lassa/gadm36_1_lassa.shp')
  g$uid = g$GID_1  
}
g0=read_sf('../gadm36_lassa/gadm36_0_lassa_reduced.shp')
g = g[which(g$uid %in% c$uid),]
#g = g[which(g$GID_0!='CMR'),]

g$foi[order(g$uid)]=apply(rr.foi.all,1,median)[order(row.names(rr.foi.all))]
g$foi.lm[order(g$uid)]=apply(lm.pred.all,1,median)[order(c$uid)]
g$foi.lm2[order(g$uid)]=apply(lm2.pred.all,1,median)[order(c$uid)]
g$foi.lo[order(g$uid)]=apply(mrkLO.pred.all,1,median)[order(c$uid)]
g$foi.locovs[order(g$uid)]=apply(mrkLOcovs.pred.all,1,median)[order(c$uid)]
g$foi.hi[order(g$uid)]=apply(mrkHI.pred.all,1,median)[order(c$uid)]
g$foi.hicovs[order(g$uid)]=apply(mrkHIcovs.pred.all,1,median)[order(c$uid)]
g$foi.rf[order(g$uid)]=apply(rf.pred.all,1,median)[order(c$uid)]
g$foi.brt[order(g$uid)]=apply(brt.pred.all,1,median)[order(c$uid)]

# plot locations with positive cases or deaths
jpeg(paste0('plots/map_foi_data_adm',admin,'_revrate',as.character(rev_rate*100),'_upd.jpeg'),width=6.5,height=6,units='in',res=300)

g$foi_data="No"
g$foi_data[which(g$uid %in% ss$uid)]="Yes"

ggplot(g)+geom_sf(aes(fill=foi_data))+theme_bw()+
  #scale_fill_gradient(low="white",high="blue")+
  scale_fill_manual(values=c("white","blue"),breaks=c("No","Yes"))+
  labs(fill="Serology data\navailable")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

dev.off()

# plot locations with positive cases or deaths
jpeg(paste0('plots/map_poscasesdeaths_adm',admin,'_revrate',as.character(rev_rate*100),'_upd.jpeg'),width=6.5,height=6,units='in',res=300)

g$candd=ifelse((rr$cases+rr$deaths)>0,"Yes","No")[
  match(as.character(g$uid),
        as.character(rr$uid))]
g$candd=ifelse(is.na(g$candd),"No",g$candd)
ggplot(g)+geom_sf(aes(fill=candd))+theme_bw()+
  scale_fill_manual(values=c("white","orange"),breaks=c("No","Yes"))+labs(fill="Cases or deaths\nreported")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=6.5,height=6,units='in',res=300)


ggplot(g) + geom_sf(aes(fill=foi),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="red")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

dev.off()



# plot map of standard deviation of force of infection
jpeg(paste('plots/map_foi_sd_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

g$foi.sd=apply(rr.foi,1,sd)[
  match(as.character(g$uid),
        row.names(rr.foi))]

ggplot(g)+geom_sf(aes(fill=foi.sd))+theme_bw()+
  scale_fill_gradient(low="white",high="orange")+labs(fill="Std. Deviation\nFOI")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

dev.off()



# plot distributions of force of infection
jpeg(paste('plots/foi_dist_adm',admin,'_revrate',as.character(rev_rate*100),'_raw.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

par(oma=c(0,0,0,0),mar=c(2.5,3.25,2,0.5))
plot(-100,-100,type='l',xlim=c(-1,nrow(rr.foi)),las=1,xaxt='n',
     ylim=c(min(apply(rr.foi,1,function(x)quantile(x,0.005))),
            max(apply(rr.foi,1,function(x)quantile(x,0.975)))),
     xlab='',ylab='',xaxs='i',yaxs='i')
mtext('Adm1s sorted by median FOI',1,line=0.75)
mtext(expression(log[10]*' FOI'),2,line=2)
polygon(
  c(1:nrow(rr.foi),rev(1:nrow(rr.foi))),c(
    apply(rr.foi,1,function(x)quantile(x,0.025))[order(apply(rr.foi,1,median))],
    rev(apply(rr.foi,1,function(x)quantile(x,0.975))[order(apply(rr.foi,1,median))])),
  border=NA,col=rgb(0.96,0.76,0.26,0.5))
lines(apply(rr.foi,1,median)[order(apply(rr.foi,1,median))],col=rgb(0.96,0.76,0.26,1))

dev.off()

10 ^ quantile(apply(rr.foi,1,median),c(0.005,0.025,0.05,0.5,0.95,0.975,0.995))
# 0%         2.5%           5%          50%          95%        97.5%         100% 
# 4.651126e-07 2.626911e-06 3.872381e-06 3.853664e-05 2.030817e-03 9.069367e-03 1.121432e+00

myvar = function(x){mean((x-mean(x))^2)}
myvar(rowMeans(rr.foi)) / myvar(as.numeric(rr.foi))
#0.4228 # 0.6373473


# 
# # plot maps of force of infection
# jpeg(paste('plots/map_foi_sd_adm',admin,'_revrate',as.character(rev_rate*100),'_null.jpeg',sep=''),width=6.5,height=6,units='in',res=300)
# 
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(null.pred,1,sd)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=1+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,0.5))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(1000 * labs) / 1000
# axis(
#   1,
#   at=1.5,
#   labels=labs[1],cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Standard deviation of '*log[10]*' FOI'),1,cex=1,line=1.5)
# 
# dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_lm.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.lm),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(lm.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_lm2.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.lm2),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(lm2.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_brt.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.brt),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(brt.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_rf.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.rf),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(rf.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_mrkLO.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.lo),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(mrkLO.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_mrkHI.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.hi),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)

# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(mrkHI.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()



# plot maps of force of infection
jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_mrkLOcovs.jpeg',sep=''),width=6.5,height=6,units='in',res=300)

ggplot(g) + geom_sf(aes(fill=foi.locovs),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
  geom_sf(data=g0,fill=NA,color="black",lwd=1)
# par(oma=c(3,0,0,0),mar=c(0.1,0.5,0,0.5))
# layout(matrix(1:2,2,1),heights=c(0.95,0.05))
# 
# fig.vec = apply(mrkLOcovs.pred,1,median)[
#   match(as.character(g$uid),c$uid)]
# fig.vec[is.na(fig.vec)] = min(fig.vec,na.rm=T)
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(quantile(fig.vec,seq(0.2,1,0.2))>fig.vec[ii])))
# col.ind = unname(sapply(1:length(fig.vec),function(ii)
#   which.max(seq(min(fig.vec),max(fig.vec),length.out=9)[-1]>=fig.vec[ii])))
# 
# range.foi = range(col.ind)
# if((floor(range.foi[2])-ceiling(range.foi[1])) <= 0){
#   if(round(range.foi[1])==round(range.foi[2])){
#     range.foi = c(round(range.foi[1])-1,round(range.foi[2])+1)
#   }else{
#     range.foi = c(floor(range.foi[1]),ceiling(range.foi[2]))
#   }
# }
# 
# col = col.ind
# col = (col - range.foi[1]) / diff(range.foi)
# col = ifelse(col<1,col,1-1e-8)
# col = ifelse(col>0,col,1e-8)
# plot(g,col=rgb(0,1,0,col))
# 
# p = sort(unique(col.ind))
# plot(-100,-100,xlim=range.foi+c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',xaxs='i',yaxs='i')
# for(ii in 1:(length(p))){
#   polygon(c(p[ii],p[ii]+1,p[ii]+1,p[ii]),c(-100,-100,100,100),border=NA,col=rgb(0,1,0,(p[ii]-min(p))/diff(range(p))))
# }
# colvals = ceiling(range.foi[1]):floor(range.foi[2])
# labs = seq(min(fig.vec),max(fig.vec),length.out=max(p)+1)
# labs = round(10 * labs) / 10
# axis(
#   1,
#   at=1:(max(p)+1),
#   labels=labs,cex.axis=0.8,mgp=c(3,0.5,0))
# mtext(expression('Median '*log[10]*' FOI'),1,cex=1,line=1.5)

dev.off()

jpeg(paste('plots/map_foi_adm',admin,'_revrate',as.character(rev_rate*100),'_mrkHIcovs.jpeg',sep=''),width=6.5,height=6,units='in',res=300)


ggplot(g) + geom_sf(aes(fill=foi.hicovs),col=NA) + theme_bw()+
  scale_fill_gradient(low="white",high="blue")+labs(fill="FOI (log10)")+
geom_sf(data=g0,fill=NA,color="black",lwd=1)

dev.off()